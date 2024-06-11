import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, cds
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Fit_analysis_plot
from Modularize.m12_T2 import ramsey_executor
from pandas import Series

def xyf_calibrator(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubit:str,max_evoT:float,execution:bool=True,plot_cali_result:bool=False):
    Trustable = False
    avg_N = 500
    _, _, _, detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,specific_qubit,artificial_detune=0e6,freeDura=max_evoT,histo_counts=1,run=execution,plot=False,avg_n=avg_N)
    if execution:
        to_cali_detune = [detune[specific_qubit], -detune[specific_qubit]]
        detune_after = []
        for arti_detune in to_cali_detune:
            slightly_print(f"Ramsey with detuning = {round(arti_detune*1e-6,2)} MHz")
            ramsey_results, _, _, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,specific_qubit,artificial_detune=arti_detune,freeDura=max_evoT,histo_counts=1,run=execution,plot=False,avg_n=avg_N)
            detune_after.append(average_actual_detune[specific_qubit])


        min_detune_after = min(detune_after)
        min_idx = Series(detune_after).idxmin()
        if min_detune_after < detune[specific_qubit]:
            highlight_print(f"{specific_qubit} XYF calibration done, now detuning = {round(min_detune_after*1e-6,2)} MHz")
            original_xyf = QD_agent.quantum_device.get_element(specific_qubit).clock_freqs.f01()
            QD_agent.quantum_device.get_element(specific_qubit).clock_freqs.f01(original_xyf+to_cali_detune[min_idx])
            highlight_print(f"{specific_qubit} correct XYF = {round(QD_agent.quantum_device.get_element(specific_qubit).clock_freqs.f01()*1e-9,3)} GHz")
            if plot_cali_result:
                Fit_analysis_plot(ramsey_results[specific_qubit],P_rescale=False,Dis=None)
            Trustable = True
        else :
            highlight_print('There is no need for detuning!')
    return Trustable
        

if __name__ == "__main__":
    
    """ Fill in """
    execution = 1
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        "q1":{"evoT":45e-6}
    }
    couplers = ['c0', 'c1']
    # 1 = Store
    # 0 = not store
    chip_info_restore = 1

    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    """ Running """
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        freeTime = ro_elements[qubit]["evoT"]

        Trustable = xyf_calibrator(QD_agent,cluster,meas_ctrl,Fctrl,qubit,freeTime,execution)
        
        """ Storing """
        if execution:
            if Trustable:
                QD_agent.QD_keeper()
                if chip_info_restore:
                    chip_info.update_xyfCali(qb=qubit)
        
    """ Close """
    shut_down(cluster,Fctrl,Cctrl)