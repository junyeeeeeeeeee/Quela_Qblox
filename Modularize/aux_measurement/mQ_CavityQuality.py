import json, os
from numpy import arange, pi, array
from Modularize.m2_CavitySpec import Cavity_spec
from Modularize.support.Experiment_setup import ip_register
from Modularize.support import Data_manager, QDmanager
from Modularize.support.UserFriend import *
from quantify_core.measurement.control import MeasurementControl
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, set_atte_for, shut_down
from Modularize.aux_measurement.CW_generator import CW_executor, turn_off_sequencer
from quantify_scheduler.helpers.collections import find_port_clock_path


def qualityFiter(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_amp:float,specific_qubits:str,ro_span_Hz:float=10e6,run:bool=True,fpts=600,n_avg:int=10, data_folder:str=""):
    rof = {str(specific_qubits):QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout()}
    
    if run:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=True,points=fpts,n_avg=n_avg,particular_folder=data_folder)[specific_qubits]
    else:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=False)[specific_qubits]
    QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(rof[specific_qubits])
    return qb_CSresults

def show_quality_for(CS_result:ResonatorSpectroscopyAnalysis, target_q:str)->dict:
    qi = round(CS_result.quantities_of_interest['Qi'].nominal_value*1e-3,2) # in k
    ql = round(CS_result.quantities_of_interest['Ql'].nominal_value*1e-3,2)
    qc = round(CS_result.quantities_of_interest['Qc'].nominal_value*1e-3,2)

    eyeson_print(f"{target_q}: Qi= {qi}k, Qc= {qc}k, Ql= {ql}k")
    return {"QI":qi,"QC":qc,"QL":ql}

def generate_data_folder(target_q:str, additional_name:str="")->str:
    DaTar = Data_manager()
    DaTar.build_folder_today(DaTar.raw_data_dir)
    Quality_folder_path = os.path.join(DaTar.raw_folder,f"{target_q}_CavityQuality_{additional_name}_{DaTar.get_time_now()}")
    if not os.path.isdir(Quality_folder_path):
        os.mkdir(Quality_folder_path)
    eyeson_print(f"Check folder path: {Quality_folder_path}")
    return Quality_folder_path


def QD_RO_init_qualityCase(QD_agent:QDmanager, target_q:str):
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.reset.duration(300e-6)
    qubit.measure.acq_delay(500e-9)
    qubit.measure.pulse_amp(0.3) # readout amp set here
    qubit.measure.pulse_duration(60e-6)
    qubit.measure.integration_time(60e-6-500e-9-4e-9)

def get_LO_freq(QD_agent:QDmanager,q:str)->float:
    qubit= QD_agent.quantum_device.get_element(q)
    clock=qubit.name + ".ro"
    port=qubit.ports.readout()
    hw_config = QD_agent.quantum_device.hardware_config()
    
    output_path = find_port_clock_path(
        hw_config, port=port, clock= clock)
    
    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    
    return float(hw_config[cluster_key][module_key][output_key]["lo_freq"])



if __name__ == "__main__":
    # we fix the readout amp = 0.3, change the ro_atte from 0 to 60 dB
    """ Fill in """
    execution = True
    RT_real_atte = 0
    qrm_slot_idx = 6
    DRandIP = {"dr":"dr1","last_ip":"11"}
    ro_elements = {
        "q2":{"ro_atte_dB":arange(start=0,stop=10+2,step=2)} # plz connect to SA and get the signal power "dBm_bySA" by the amp=0.3 and with the desired rof
    }
    # 

    """ Preparations """ 
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')

    """ Running """
    qubit_rec = {}
    for qubit in ro_elements:
        this_qubit_exp_folder = generate_data_folder(qubit,additional_name=f"RTatte{RT_real_atte}dB")
        desired_rof = QD_agent.quantum_device.get_element(qubit).clock_freqs.readout()
        mark_input("Connect your RO to the SA now! press ENTER to start...")
        cluster, connected_module = CW_executor(cluster, qrm_slot_idx, port_idx=0, ro_atte=0, amp=0.3, RF_freq=QD_agent.quantum_device.get_element(qubit).clock_freqs.readout(), LO_freq=get_LO_freq(QD_agent,qubit))
        dBm_bySA = float(mark_input(f"input the power (dBm) which is at freq = {round(desired_rof*1e-9,3)} GHz: "))
        turn_off_sequencer(cluster, connected_module)
        conti = mark_input("Once you have connected the RO cable bact to DR, press ENTER to continue...")
        QD_RO_init_qualityCase(QD_agent, qubit)
        dBm_output = list(dBm_bySA - ro_elements[qubit]["ro_atte_dB"] - RT_real_atte)
        qualities = {"dBm":dBm_output,"Qi":[],"Qc":[],"Ql":[]}
        for ro_atte in ro_elements[qubit]["ro_atte_dB"]:
            set_atte_for(QD_agent.quantum_device,ro_atte,'ro',[qubit])
            positional_result = qualityFiter(QD_agent,meas_ctrl,ro_amp=0.3,specific_qubits=qubit,ro_span_Hz=10e6,run=True)
            peak_width = round((positional_result.quantities_of_interest["fr"].nominal_value*2*pi/positional_result.quantities_of_interest["Ql"].nominal_value)/(2*pi*1e6),2)*1e6
            CS_results = qualityFiter(QD_agent=QD_agent,meas_ctrl=meas_ctrl,specific_qubits=qubit,ro_amp=0.3,run = execution,ro_span_Hz=3*peak_width,data_folder=this_qubit_exp_folder)
            cluster.reset()
            qua = show_quality_for(CS_results,qubit)
            qualities["Qi"].append(qua["QI"])
            qualities["Qc"].append(qua["QC"])
            qualities["Ql"].append(qua["QL"])
        qubit_rec[qubit] = qualities

        """ Storing """
        with open(os.path.join(this_qubit_exp_folder,"Quality_results.json"), "w") as json_file:
            json.dump(qubit_rec, json_file)

    """ Close """
    print('Cavity quality fit done!')
    shut_down(cluster,Fctrl)