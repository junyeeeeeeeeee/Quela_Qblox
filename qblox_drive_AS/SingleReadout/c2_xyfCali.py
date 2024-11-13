import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_instruments import Cluster
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import Fit_analysis_plot
from qblox_drive_AS.SingleReadout.m12_T2 import ramsey_executor
from pandas import Series

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = {
        "q0":{"evoT":5e-6}
    }
    couplers = []


    """ Iteration (Do NOT touch!)"""
    abs_detuning = 0
    step = 0 # should be 0, 1, 2 and then break 

    for qubit in ro_elements:
        while True:
            """ Preparation """
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
            

            """ Running """
            Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            if step == 0:
                _, _, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=abs_detuning,freeDura=ro_elements[qubit]["evoT"],run=execution,avg_n=500)
                abs_detuning += average_actual_detune[qubit]
            else:
                trying_detune = ((-1)**(step))*abs_detuning
                _, _, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=trying_detune,freeDura=ro_elements[qubit]["evoT"],run=execution,avg_n=500)
            eyeson_print(f"Detuning = {round(abs_detuning*1e-6,2)} MHz")
            if step > 0 and average_actual_detune[qubit] < abs_detuning:
                highlight_print(f"XYF calibration done successfully! detuning = {round(average_actual_detune[qubit]*1e-6,2)} MHz")
                """ Storing """
                if execution:
                    original_xyf = QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()
                    QD_agent.quantum_device.get_element(qubit).clock_freqs.f01(original_xyf+trying_detune)
                    QD_agent.QD_keeper()
                """ Close """
                shut_down(cluster,Fctrl)
                break
            elif step == 2 :
                warning_print("Didn't find a good XYF !")
                """ Close """
                shut_down(cluster,Fctrl)
                break
            else:
                """ Close """
                shut_down(cluster,Fctrl)

            step += 1
            

        
        
            
                
           