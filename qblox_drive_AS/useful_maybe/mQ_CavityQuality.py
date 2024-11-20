""" 600 fpoints, 100 avg ~ 0.6 min, Multiplexing"""
""" Do this after m2 """

import json, os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from numpy import arange, pi, mean, array, linspace
from qblox_drive_AS.SOP.CavitySpec import Cavity_spec, multiplexing_CS_ana
from qblox_drive_AS.Configs.ClusterAddress_rec import ip_register
from qblox_instruments import Cluster
from qblox_drive_AS.support.QDmanager import Data_manager, QDmanager
from qblox_drive_AS.support.UserFriend import *
from quantify_core.measurement.control import MeasurementControl
from xarray import Dataset
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, set_atte_for, shut_down
from qblox_drive_AS.aux_measurement.CW_generator import CW_executor, turn_off_sequencer
from quantify_scheduler.helpers.collections import find_port_clock_path


def qualityFiter(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_amp:dict,ro_span_Hz:dict,run:bool=True,fpts=600,n_avg:int=100, data_folder:str="")->tuple[Dataset, dict]:
    rof_origin = {}
    ro_elements = {}
    for qb in list(ro_amp.keys()):
        rof = QD_agent.quantum_device.get_element(qb).clock_freqs.readout()
        rof_origin[qb] = rof
        ro_elements[qb] = linspace(rof-ro_span_Hz[qb], rof+ro_span_Hz[qb], fpts)

    if run:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_elements,run=True,n_avg=n_avg,particular_folder=data_folder,ro_amps=ro_amp)
    else:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_elements,run=False,ro_amps=ro_amp)

    for qb in list(ro_amp.keys()):
        QD_agent.quantum_device.get_element(qb).clock_freqs.readout(rof_origin[qb])
    
    return qb_CSresults, ro_elements

def get_quality_for(QD_agent:QDmanager, CS_result:Dataset, ro_elements:dict, pic_save:bool=False)->dict:
    fit_results = multiplexing_CS_ana(QD_agent, CS_result, ro_elements, save_pic=pic_save)
    return fit_results

def generate_data_folder(target_q:str, additional_name:str="")->str:
    DaTar = Data_manager()
    DaTar.build_folder_today(DaTar.raw_data_dir)
    Quality_folder_path = os.path.join(DaTar.raw_folder,f"{target_q}_CavityQuality_{additional_name}_{DaTar.get_time_now()}")
    if not os.path.isdir(Quality_folder_path):
        os.mkdir(Quality_folder_path)
    eyeson_print(f"Check folder path: {Quality_folder_path}")
    return Quality_folder_path


def QD_RO_init_qualityCase(QD_agent:QDmanager, target_q:str, ro_amp:float=0.3):
    ro_pulse_length:float = 100e-6
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.reset.duration(ro_pulse_length+300e-6)
    qubit.measure.acq_delay(280e-9)
    qubit.measure.pulse_amp(ro_amp) # readout amp set here
    qubit.measure.pulse_duration(ro_pulse_length)
    qubit.measure.integration_time(ro_pulse_length-280e-9-4e-9)

def get_LO_freq(QD_agent:QDmanager,q:str)->float:
    qubit= QD_agent.quantum_device.get_element(q)
    clock=qubit.name + ".ro"
    port='q:res'
    hw_config = QD_agent.quantum_device.hardware_config()
    
    output_path = find_port_clock_path(
        hw_config, port=port, clock= clock)
    
    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    
    return float(hw_config[cluster_key][module_key][output_key]["lo_freq"])

def cavQuality_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,qrm_slot_idx:int,RT_real_atte:int,ro_elements:dict,atte_settings:dict,execution:bool=True,ro_amp=0.3,window_ratio:int=6):
    this_qubit_exp_folder = generate_data_folder('Multiplex',additional_name=f"RTatte{RT_real_atte}dB")
    atte_to_apply = arange(atte_settings["atte_start"], atte_settings['atte_end']+atte_settings["step"], atte_settings["step"])
    SA_dBm = {"q0":5.3,"q1":5.6,"q2":3.7}
    ro_amps = {}
    ro_span_Hzs = {}
    avg_n = 100
    if RT_real_atte >= 100:
        avg_n *= 10


    # # Connect SA get power in dBm
    mark_input("Connect your RO to the SA now! press ENTER to start...")
    for q_idx, q in enumerate(list(ro_elements.keys())):
        if ro_elements[q]["assigned_rof_Hz"] != 0:
            QD_agent.quantum_device.get_element(q).clock_freqs.readout(ro_elements[q]["assigned_rof_Hz"])
        desired_rof = QD_agent.quantum_device.get_element(q).clock_freqs.readout()
        # cluster, connected_module = CW_executor(cluster, qrm_slot_idx, port_idx=0, ro_atte=0, amp=ro_amp, RF_freq=desired_rof, LO_freq=get_LO_freq(QD_agent,q))
        # SA_dBm[q] = float(mark_input(f"input the power (dBm) which is at freq = {round(desired_rof*1e-9,3)} GHz for {q}: "))
        # turn_off_sequencer(cluster, connected_module)

        # if q_idx != len(list(ro_elements.keys())) -1:
        #     conti = mark_input("press ENTER to continue measure the next qubit...")
        # else:
        #     conti = mark_input("Once the cable connected back to DR, press ENTER to continue measure the next qubit...")

        QD_RO_init_qualityCase(QD_agent, q, ro_amp)
        ro_amps[q] = ro_amp
        ro_span_Hzs[q]=7e6
        set_atte_for(QD_agent.quantum_device,min(atte_to_apply),'ro',[q]) # only for the first scan to align the cavity window

    # first multiplexing scan for window align
    eyeson_print("Window align... ")
    positional_result, elements  = qualityFiter(QD_agent,meas_ctrl,ro_amp=ro_amps,ro_span_Hz=ro_span_Hzs,run=execution,n_avg=avg_n)
    align_results = get_quality_for(QD_agent, positional_result, elements)
    for q in ro_elements.keys():
        peak_width = round((align_results[q]["fr"]*2*pi/align_results[q]["Ql"])/(2*pi*1e6),2)*1e6
        ro_span_Hzs[q] = (window_ratio/2)*peak_width

    # start changing atte 
    freqs = {}
    if execution:
        for ro_atte in atte_to_apply:
            highlight_print(f"Now atte = {ro_atte} dB")
            set_atte_for(QD_agent.quantum_device,ro_atte,'ro',list(ro_amps.keys()))
            CS_results, ro_info = qualityFiter(QD_agent,meas_ctrl,ro_amps,run = execution,ro_span_Hz=ro_span_Hzs,data_folder=this_qubit_exp_folder,fpts=window_ratio*100,n_avg=avg_n)
            cluster.reset()
        
        for q in ro_info:
            freqs[q] = array(ro_info[q]).tolist()
        
        return this_qubit_exp_folder, SA_dBm, freqs, list(array(atte_to_apply).astype(float))
    else:
        print(f"Cavity width = {peak_width*1e-6} MHz")
        raise RuntimeError(f"No worries, This execution is set as False ~ ")



if __name__ == "__main__":
    # we fix the readout amp = 0.3, change the ro_atte from 0 to 60 dB
    """ Fill in """
    execution = True
    RT_real_atte = 120
    qrm_slot_idx = 18
    DRandIP = {"dr":"dr4","last_ip":"81"}
    atte_settings:dict = {"atte_start":0, "atte_end":60, "step":10} # atte should be multiples of 2
    ro_elements = {
        "q0":{"assigned_rof_Hz":0}, # if assigned_rof in Hz was given, the following exp will use this RO frequency. 
        "q1":{"assigned_rof_Hz":0},
        "q2":{"assigned_rof_Hz":0},
    }
    

    """ Don't change it if it's unnecessary """
    ro_amp = 0.3
    window_over_peakwidth = 6


    """ Preparations """ 
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')


    """ Running """
    this_qubit_exp_folder, SA_dBm, ro_elements, applied_atte = cavQuality_executor(QD_agent,cluster,meas_ctrl,qrm_slot_idx,RT_real_atte,ro_elements,atte_settings,execution,ro_amp,window_ratio=window_over_peakwidth)

    additional_info = {"RT_atte_dB":float(RT_real_atte),"ro_elements":ro_elements, "SA_dBm":SA_dBm, "applied_atte":applied_atte}

    """ Storing """
    with open(os.path.join(this_qubit_exp_folder,f"Additional_info.json"), "w") as json_file:
        json.dump(additional_info, json_file)


    """ Close """
    print('Cavity quality fit done!')
    shut_down(cluster,Fctrl)