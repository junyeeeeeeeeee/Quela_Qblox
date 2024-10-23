"""
The RO-freq is located in bare cavity now, use different RO attenuation and RO-amp to fit the cavity freq.\n
Because the RO-freq is the bare cavity, the 'window_shift' should be 0 when the conditions are fior a bare cavity.\n
On the other hand, the 'window_shift' will be a rough dispersive shift for the dressed cavity conditions. 
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.m2_CavitySpec import Cavity_spec, multiplexing_CS_ana
from Modularize.support import Data_manager, QDmanager, cds
from Modularize.support.UserFriend import *
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down
from numpy import linspace

Quality_values = ["Qi_dia_corr", "Qc_dia_corr", "Ql"]
Quality_errors = ["Qi_dia_corr_err", "absQc_err", "Ql_err"]

def preciseCavity_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,run:bool=True, ro_amps:dict={})->dict:
    
    if run:
        cs_ds = Cavity_spec(QD_agent,meas_ctrl,ro_elements,run=True,ro_amps=ro_amps,n_avg=100)
        CS_results = multiplexing_CS_ana(QD_agent, list(ro_elements.keys()), cs_ds)
    else:
        _ = Cavity_spec(QD_agent,meas_ctrl,ro_elements,run=False,ro_amps=ro_amps)
        CS_results = {}
    
    return CS_results

def fillin_PDans(QD_agent:QDmanager,ans:dict):
    """
    Fill in the power dep answer to the quantum_device.\n
    format:\n
    `ans = {"q0":{"dressF_Hz":,"dressP":,"bareF_Hz":},...}`
    """
    for q in ans:
        qubit = QD_agent.quantum_device.get_element(q)
        if ans[q]["dressP"] != "": qubit.measure.pulse_amp(ans[q]["dressP"]) 
        if ans[q]["dressF_Hz"] != "": qubit.clock_freqs.readout(ans[q]["dressF_Hz"])
        if ans[q]["bareF_Hz"] != "": QD_agent.Notewriter.save_bareFreq_for(target_q=q,bare_freq=ans[q]["bareF_Hz"])
        if ans[q]["ro_atte"] != "": QD_agent.Notewriter.save_DigiAtte_For(atte_dB=ans[q]["ro_atte"],target_q=q,mode='ro')

    QD_agent.refresh_log("PD answer stored!")
    QD_agent.QD_keeper()

def BDC_waiter(QD_agent:QDmanager, state:str, qubits:dict, ro_atte:dict, ro_span_Hz:float, fpts:int):
    amps = {}
    ro_elements = {}
    PD_ans = {}
    original_rof = {}
    for q in qubits:
        if state in list(qubits[q].keys()):
            amps[q] = qubits[q][state]['ro_amp']
            
            rof = QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            original_rof[q] = rof
            ro_elements[q] = linspace(rof+qubits[q][state]["window_shift"]-ro_span_Hz, rof+qubits[q][state]["window_shift"]+ro_span_Hz, fpts)
            
        PD_ans[q] = {"dressF_Hz":"","dressP":"","bareF_Hz":"","ro_atte":""}
    init_system_atte(QD_agent.quantum_device,list(qubits.keys()),ro_out_att=ro_atte[state])
    return ro_elements, amps, PD_ans, original_rof


if __name__ == "__main__":

    """ Fill in """
    execution:bool = True 
    sweetSpot:bool = 0     # If true, only support one one qubit
    chip_info_restore:bool = 0
    DRandIP = {"dr":"drke","last_ip":"242"}
    ro_element = {
        "q2":{  "bare" :{"ro_amp":0.1,"window_shift":0e6},
                "dress":{"ro_amp":0.3,"window_shift":0e6}},
        "q3":{  "bare" :{"ro_amp":0.1,"window_shift":0e6},
                "dress":{"ro_amp":0.3,"window_shift":0e6}},
    }
    ro_attes = {"dress":50, "bare":20} # all ro_elements shared

    """ Optional paras"""
    half_ro_freq_window_Hz = 20e6
    freq_data_points = 400


    
   
    CS_results = {}
    for state in ["bare", "dress"]: # bare is always the first !!
        eyeson_print(f"{state} cavities: ")
        """ Preparations """ 
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
        ro_elements, amps, PD_ans, original_rof = BDC_waiter(QD_agent, state, ro_element, ro_attes, half_ro_freq_window_Hz, freq_data_points)
        # Create or Load chip information
        chip_info = cds.Chip_file(QD_agent=QD_agent)


        """ Running """
        if sweetSpot:
            Fctrl[list(ro_element.keys())[0]](QD_agent.Fluxmanager.get_sweetBiasFor(target_q=list(ro_element.keys())[0]))
        
        CS_results[state] = preciseCavity_executor(QD_agent,meas_ctrl,ro_elements,run=execution,ro_amps=amps)
        if sweetSpot:
            Fctrl[list(ro_element.keys())[0]](0)
        cluster.reset()
        for q in ro_elements:
            QD_agent.quantum_device.get_element(q).clock_freqs.readout(original_rof[q])
        
        for qubit in CS_results[state]:
            print(f"{qubit} @ {state} cavity: ")
            if state == "bare":
                PD_ans[qubit]["bareF_Hz"] = CS_results[state][qubit]['fr']
            else:
                PD_ans[qubit]["dressF_Hz"] = CS_results[state][qubit]['fr']
                PD_ans[qubit]["dressP"] = ro_element[qubit][state]["ro_amp"]
                PD_ans[qubit]["ro_atte"] = ro_attes[state]
            for Qua_idx, Qua in enumerate(Quality_values):
                print(f"{Qua[:2]} = {round(float(CS_results[state][qubit][Qua])/1000,2)} 土 {round(float(CS_results[state][qubit][Quality_errors[Qua_idx]])/1000,2)} k")
            


        """ Storing """
        fillin_PDans(QD_agent, PD_ans)

        if chip_info_restore:
            chip_info.update_QD(CS_results)
            if sweetSpot:
                chip_info.update_BDCavityFit_sweet(CS_results)
            else:
                chip_info.update_BDCavityFit(CS_results)


        """ Close """
        print('Cavity quality fit done!')
        shut_down(cluster,Fctrl)




    # # If you want to analyze a cavity nc by ResonatorSpectroscopyAnalysis 
    # # re-analyze a nc
    # from numpy import pi
    # from xarray import open_dataset
    # from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
    # import quantify_core.data.handling as dh
    # from quantify_core.data.handling import set_datadir
    
    # set_datadir('path_to_datadir')
   
    # rs_ds = open_dataset("Modularize/Meas_raw/2024_4_29/DR1q0_CavitySpectro_H20M41S2.nc")
    # x = ResonatorSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
    # print((x.quantities_of_interest["fr"].nominal_value*2*pi/x.quantities_of_interest["Ql"].nominal_value)/(2*pi*1e6))