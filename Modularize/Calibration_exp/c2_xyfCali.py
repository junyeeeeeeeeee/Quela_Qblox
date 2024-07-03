import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_instruments import Cluster
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Fit_analysis_plot
from Modularize.m12_T2 import ramsey_executor
from pandas import Series

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        "q0":{"evoT":50e-6},
        # "q1":{"evoT":30e-6},
    }
    couplers = ['c0','c1','c2','c3']


    """ Iteration (Do NOT touch!)"""
    abs_detuning = 0
    step = 0 # should be 0, 1, 2 and then break 

    for qubit in ro_elements:
        while True:
            """ Preparation """
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            

            """ Running """
            Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            if step == 0:
                _, _, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=abs_detuning,freeDura=ro_elements[qubit]["evoT"],run=execution,avg_n=800)
                abs_detuning += average_actual_detune[qubit]
            else:
                trying_detune = ((-1)**(step))*abs_detuning
                _, _, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=trying_detune,freeDura=ro_elements[qubit]["evoT"],run=execution,avg_n=800)
            
            if step > 0 and average_actual_detune[qubit] < abs_detuning:
                highlight_print(f"XYF calibration done successfully! detuning = {round(average_actual_detune[qubit]*1e-6,4)} MHz")
                """ Storing """
                if execution:
                    original_xyf = QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()
                    QD_agent.quantum_device.get_element(qubit).clock_freqs.f01(original_xyf+trying_detune)
                    QD_agent.QD_keeper()
                """ Close """
                shut_down(cluster,Fctrl,Cctrl)
                break
            elif step == 2 :
                warning_print("Didn't find a good XYF !")
                """ Close """
                shut_down(cluster,Fctrl,Cctrl)
                break
            else:
                """ Close """
                shut_down(cluster,Fctrl,Cctrl)

            step += 1
            

        
        
            
                
           