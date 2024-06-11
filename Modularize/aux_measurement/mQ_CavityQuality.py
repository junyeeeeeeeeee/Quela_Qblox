""" 600 fpoints, 100 avg ~ 0.6 min"""

import json, os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from numpy import arange, pi, mean, array
from Modularize.m2_CavitySpec import Cavity_spec
from Modularize.support.Experiment_setup import ip_register
from qblox_instruments import Cluster
from Modularize.support import Data_manager, QDmanager
from Modularize.support.UserFriend import *
from quantify_core.measurement.control import MeasurementControl
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, set_atte_for, shut_down
from Modularize.aux_measurement.CW_generator import CW_executor, turn_off_sequencer
from quantify_scheduler.helpers.collections import find_port_clock_path


def qualityFiter(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_amp:float,specific_qubits:str,ro_span_Hz:float=10e6,run:bool=True,fpts=600,n_avg:int=100, data_folder:str=""):
    rof = {str(specific_qubits):QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout()}
    
    if run:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=True,points=fpts,n_avg=n_avg,particular_folder=data_folder)[specific_qubits]
    else:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=False)[specific_qubits]
    QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(rof[specific_qubits])
    return qb_CSresults

def get_quality_for(CS_result:ResonatorSpectroscopyAnalysis, target_q:str, mark_label:str="")->dict:
    qi = round(CS_result.quantities_of_interest['Qi'].nominal_value*1e-3,2) # in k
    qi_sd = round(CS_result.quantities_of_interest['Qi'].std_dev*1e-3,2)
    ql = round(CS_result.quantities_of_interest['Ql'].nominal_value*1e-3,2)
    ql_sd = round(CS_result.quantities_of_interest['Ql'].std_dev*1e-3,2)
    qc = round(CS_result.quantities_of_interest['Qc'].nominal_value*1e-3,2)
    qc_sd = round(CS_result.quantities_of_interest['Qc'].std_dev*1e-3,2)
    eyeson_print(f"{target_q}{mark_label}: Qi= {qi}k, Qc= {qc}k, Ql= {ql}k")
    return {"QI":qi*1e3,"QI_sd":qi_sd*1e3,"QC":qc*1e3,"QC_sd":qc_sd*1e3,"QL":ql*1e3,"QL_sd":ql_sd*1e3}

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
    qubit.measure.acq_delay(100e-9)
    qubit.measure.pulse_amp(ro_amp) # readout amp set here
    qubit.measure.pulse_duration(100e-6)
    qubit.measure.integration_time(100e-6-100e-9-4e-9)

def get_LO_freq(QD_agent:QDmanager,q:str)->float:
    qubit= QD_agent.quantum_device.get_element(q)
    clock=qubit.name + ".ro"
    port=qubit.ports.readout()
    hw_config = QD_agent.quantum_device.hardware_config()
    
    output_path = find_port_clock_path(
        hw_config, port=port, clock= clock)
    
    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    
    return float(hw_config[cluster_key][module_key][output_key]["lo_freq"])

def cavQuality_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,target_q:str,qrm_slot_idx:int,RT_real_atte:int,ro_elements:dict,execution:bool=True,ro_amp=0.3,window_ratio:int=6):
    qubit_rec = {target_q:{}}
    this_qubit_exp_folder = generate_data_folder(target_q,additional_name=f"RTatte{RT_real_atte}dB")
    if ro_elements[target_q]["assigned_rof_Hz"] == 0:
        desired_rof = QD_agent.quantum_device.get_element(target_q).clock_freqs.readout()
    else:
        desired_rof = ro_elements[target_q]["assigned_rof_Hz"]
    # # Connect SA get power in dBm
    mark_input("Connect your RO to the SA now! press ENTER to start...")
    cluster, connected_module = CW_executor(cluster, qrm_slot_idx, port_idx=0, ro_atte=0, amp=ro_amp, RF_freq=QD_agent.quantum_device.get_element(target_q).clock_freqs.readout(), LO_freq=get_LO_freq(QD_agent,target_q))
    dBm_bySA = float(mark_input(f"input the power (dBm) which is at freq = {round(desired_rof*1e-9,3)} GHz: "))
    turn_off_sequencer(cluster, connected_module)
    conti = mark_input("Once you have connected the RO cable bact to DR, press ENTER to continue...")
    QD_RO_init_qualityCase(QD_agent, target_q, ro_amp)
    # # align window
    set_atte_for(QD_agent.quantum_device,mean(ro_elements[target_q]["ro_atte_dB"]),'ro',[target_q])
    positional_result = qualityFiter(QD_agent,meas_ctrl,ro_amp=ro_amp,specific_qubits=target_q,ro_span_Hz=10e6,run=execution)
    peak_width = round((positional_result.quantities_of_interest["fr"].nominal_value*2*pi/positional_result.quantities_of_interest["Ql"].nominal_value)/(2*pi*1e6),2)*1e6

    dBm_output = list(dBm_bySA - ro_elements[target_q]["ro_atte_dB"] - RT_real_atte)
    qualities = {"dBm":dBm_output,"Qi":[],"Qi_sd":[],"Qc":[],"Qc_sd":[],"Ql":[],"Ql_sd":[]}
    if execution:
        for ro_atte in ro_elements[target_q]["ro_atte_dB"]:
            highlight_print(f"Now atte = {ro_atte} dB")
            set_atte_for(QD_agent.quantum_device,ro_atte,'ro',[target_q])
            CS_results = qualityFiter(QD_agent=QD_agent,meas_ctrl=meas_ctrl,specific_qubits=target_q,ro_amp=ro_amp,run = execution,ro_span_Hz=(window_ratio/2)*peak_width,data_folder=this_qubit_exp_folder,fpts=window_ratio*100)
            cluster.reset()
            qua = get_quality_for(CS_results,target_q,f", ro_atte_{ro_atte}dB")
            qualities["Qi"].append(qua["QI"])
            qualities["Qc"].append(qua["QC"])
            qualities["Ql"].append(qua["QL"])
            qualities["Qi_sd"].append(qua["QI_sd"])
            qualities["Qc_sd"].append(qua["QC_sd"])
            qualities["Ql_sd"].append(qua["QL_sd"])
        qubit_rec[target_q] = qualities
        return this_qubit_exp_folder, qubit_rec
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
    ro_elements = {
        "q2":{"ro_atte_dB":arange(start=0,stop=60+4,step=4),"assigned_rof_Hz":0} # if assigned_rof in Hz was given, the following exp will use this RO frequency. 
    }
    

    """ Don't change it if it's unnecessary """
    ro_amp = 0.3
    window_over_peakwidth = 6


    """ Preparations """ 
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')


    """ Running """
    for qubit in ro_elements:
        this_qubit_exp_folder, qubit_rec = cavQuality_executor(QD_agent,cluster,meas_ctrl,qubit,qrm_slot_idx,RT_real_atte,ro_elements,execution,ro_amp,window_ratio=window_over_peakwidth)

        """ Storing """
        with open(os.path.join(this_qubit_exp_folder,"Quality_results.json"), "w") as json_file:
            json.dump(qubit_rec, json_file)


    """ Close """
    print('Cavity quality fit done!')
    shut_down(cluster,Fctrl)