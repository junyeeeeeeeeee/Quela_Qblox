""" 600 fpoints, 100 avg ~ 0.6 min, Multiplexing"""
""" Do this after m2 """

import json, os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from numpy import arange, pi, mean, array, linspace
from Modularize.m2_CavitySpec import Cavity_spec, multiplexing_CS_ana
from Modularize.support.Experiment_setup import ip_register
from qblox_instruments import Cluster
from Modularize.support import Data_manager, QDmanager
from Modularize.support.UserFriend import *
from quantify_core.measurement.control import MeasurementControl
from xarray import Dataset
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, set_atte_for, shut_down
from Modularize.aux_measurement.CW_generator import CW_executor, turn_off_sequencer
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
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.reset.duration(300e-6)
    qubit.measure.acq_delay(280e-9)
    qubit.measure.pulse_amp(ro_amp) # readout amp set here
    qubit.measure.pulse_duration(100e-6)
    qubit.measure.integration_time(100e-6-280e-9-4e-9)

def get_LO_freq(QD_agent:QDmanager,q:str)->float:
    qubit= QD_agent.quantum_device.get_element(q)
    clock=qubit.name + ".ro"
    port=qubit.ports.readout()
    hw_config = QD_agent.quantum_device.hardware_config()
    
    output_path = find_port_clock_path(
        hw_config, port=port, clock= clock)
    
    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    
    return float(hw_config[cluster_key][module_key][output_key]["lo_freq"])

def cavQuality_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,qrm_slot_idx:int,RT_real_atte:int,ro_elements:dict,atte_settings:dict,execution:bool=True,ro_amp=0.3,window_ratio:int=6):
    this_qubit_exp_folder = generate_data_folder('Multiplex',additional_name=f"RTatte{RT_real_atte}dB")
    atte_to_apply = arange(atte_settings["atte_start"], atte_settings['atte_end']+atte_settings["step"], atte_settings["step"])
    SA_dBm = {}
    ro_amps = {}
    ro_span_Hzs = {}
    # # Connect SA get power in dBm
    mark_input("Connect your RO to the SA now! press ENTER to start...")
    for q in ro_elements.keys():
        if ro_elements[q]["assigned_rof_Hz"] != 0:
            QD_agent.quantum_device.get_element(q).clock_freqs.readout(ro_elements[q]["assigned_rof_Hz"])
        desired_rof = QD_agent.quantum_device.get_element(q).clock_freqs.readout()
        cluster, connected_module = CW_executor(cluster, qrm_slot_idx, port_idx=0, ro_atte=0, amp=ro_amp, RF_freq=desired_rof, LO_freq=get_LO_freq(QD_agent,q))
        SA_dBm[q] = float(mark_input(f"input the power (dBm) which is at freq = {round(desired_rof*1e-9,3)} GHz for {q}: "))
        turn_off_sequencer(cluster, connected_module)
        conti = mark_input("Once you have connected the RO cable bact to DR, press ENTER to continue...")
        QD_RO_init_qualityCase(QD_agent, q, ro_amp)
        ro_amps[q] = ro_amp
    
        set_atte_for(QD_agent.quantum_device,mean(atte_to_apply),'ro',[q]) # only for the first scan to align the cavity window

    # first multiplexing scan for window align
    positional_result, elements  = qualityFiter(QD_agent,meas_ctrl,ro_amp=ro_amps,ro_span_Hz=10e6,run=execution)
    align_results = get_quality_for(QD_agent, positional_result, elements)
    for q in ro_elements.keys():
        peak_width = round((align_results[q]["fr"]*2*pi/positional_result[q]["Ql"])/(2*pi*1e6),2)*1e6
        ro_span_Hzs[q] = (window_ratio/2)*peak_width

    # start changing atte 
    rec:dict = {}
    if execution:
        for ro_atte in atte_to_apply:
            highlight_print(f"Now atte = {ro_atte} dB")
            set_atte_for(QD_agent.quantum_device,ro_atte,'ro',list(ro_amps.keys()))
            CS_results, ro_info = qualityFiter(QD_agent,meas_ctrl,ro_amps,run = execution,ro_span_Hz=ro_span_Hzs,data_folder=this_qubit_exp_folder,fpts=window_ratio*100)
            cluster.reset()
            qua = get_quality_for(QD_agent, CS_results, ro_info)
            dBm_output = {}
            for q in qua:
                dBm_output[q] = SA_dBm[q]-RT_real_atte-ro_atte

            rec[ro_atte] = {"output_dBm":dBm_output}|qua
        
        return this_qubit_exp_folder, rec
    else:
        print(f"Cavity width = {peak_width*1e-6} MHz")
        raise RuntimeError(f"No worries, This execution is set as False ~ ")



if __name__ == "__main__":
    # we fix the readout amp = 0.3, change the ro_atte from 0 to 60 dB
    """ Fill in """
    execution = True
    RT_real_atte = 0
    qrm_slot_idx = 6
    DRandIP = {"dr":"dr1","last_ip":"11"}
    atte_settings:dict = {"atte_start":0, "atte_end":60, "step":4} # atte should be multiples of 2
    ro_elements = {
        "q2":{"assigned_rof_Hz":0} # if assigned_rof in Hz was given, the following exp will use this RO frequency. 
    }
    

    """ Don't change it if it's unnecessary """
    ro_amp = 0.3
    window_over_peakwidth = 6


    """ Preparations """ 
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')


    """ Running """
    this_qubit_exp_folder, qubit_rec = cavQuality_executor(QD_agent,cluster,meas_ctrl,qrm_slot_idx,RT_real_atte,ro_elements,atte_settings,execution,ro_amp,window_ratio=window_over_peakwidth)


    """ Storing """
    with open(os.path.join(this_qubit_exp_folder,f"Quality_results_RT{RT_real_atte}db.json"), "w") as json_file:
        json.dump(qubit_rec, json_file)


    """ Close """
    print('Cavity quality fit done!')
    shut_down(cluster,Fctrl)