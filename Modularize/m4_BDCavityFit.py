import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.m2_CavitySpec import Cavity_spec
from Modularize.support import Data_manager, QDmanager, cds
from Modularize.support.UserFriend import *
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down

def preciseCavity_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_amp:float,specific_qubits:str,ro_span_Hz:float=10e6,run:bool=True,f_shifter:float=0,fpts=200):
    # QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(5.92e9)
    rof = {str(specific_qubits):QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout()+f_shifter}
    
    if run:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=True,points=fpts)[specific_qubits]
    else:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=False)[specific_qubits]
    
    QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(rof[specific_qubits]-f_shifter)

    return qb_CSresults

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


if __name__ == "__main__":

    """ Fill in """
    execution:bool = True
    sweetSpot:bool = 0
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        # "q0":{  "bare" :{"ro_amp":1,"ro_atte":30,"window_shift":0},
                # "dress":{"ro_amp":0.01,"ro_atte":30,"window_shift":3}},
        # "q1":{  "bare" :{"ro_amp":1,"ro_atte":30,"window_shift":0},
                #  "dress":{"ro_amp":0.01,"ro_atte":30,"window_shift":0}},
         "q2":{  "bare" :{"ro_amp":1,"ro_atte":30,"window_shift":1},
                  "dress":{"ro_amp":0.03,"ro_atte":30,"window_shift":4}}
        #  "q3":{  "bare" :{"ro_amp":1,"ro_atte":30,"window_shift":0},
                #  "dress":{"ro_amp":0.01,"ro_atte":30,"window_shift":0}},
        #  "q4":{  "bare" :{"ro_amp":1,"ro_atte":30,"window_shift":0},
                #  "dress":{"ro_amp":0.01,"ro_atte":30,"window_shift":0}},

    }
    # 1 = Store
    # 0 = not store
    chip_info_restore = 1

    """ Preparations """ 
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    # Create or Load chip information
    chip_info = cds.Chip_file(QD_agent=QD_agent)


    """ Running """
    CS_results = {}
    PDans = {}
    for qubit in ro_elements:
        CS_results[qubit] = {}
        PDans[qubit] = {"dressF_Hz":0,"dressP":0,"bareF_Hz":0,"ro_atte":0}
        for state in ro_elements[qubit]:
            if ro_elements[qubit][state]["ro_atte"] != '':
                QD_agent.Notewriter.save_DigiAtte_For(ro_elements[qubit][state]["ro_atte"],qubit,'ro')
            if sweetSpot:
                Fctrl[qubit](QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qubit))
            else:
                Fctrl[qubit](0)
            init_system_atte(QD_agent.quantum_device,[qubit],ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
            CS_results[qubit][state] = preciseCavity_executor(QD_agent=QD_agent,meas_ctrl=meas_ctrl,specific_qubits=qubit,ro_amp=ro_elements[qubit][state]["ro_amp"],run = execution, f_shifter=ro_elements[qubit][state]["window_shift"],ro_span_Hz=3e6)
            Fctrl[qubit](0)
            cluster.reset()
            if state == "bare":
                PDans[qubit]["bareF_Hz"] = CS_results[qubit][state].quantities_of_interest['fr'].nominal_value
            else:
                PDans[qubit]["dressF_Hz"] = CS_results[qubit][state].quantities_of_interest['fr'].nominal_value
                PDans[qubit]["dressP"] = ro_elements[qubit][state]["ro_amp"]
                PDans[qubit]["ro_atte"] = QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro')


    """ Storing (future) """
    fillin_PDans(QD_agent, PDans)

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