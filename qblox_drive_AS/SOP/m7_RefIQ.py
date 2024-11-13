import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from numpy import NaN, array
from utils.tutorial_utils import show_args
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_Qubit_SS_sche, Single_shot_ref_fit_analysis, pulse_preview, Single_shot_fit_plot

def Single_shot_ref_spec(QD_agent:QDmanager,ro_elements:dict,shots:int=1000,run:bool=True,Experi_info:dict={}):
    print("Single shot start")
    sche_func = multi_Qubit_SS_sche   
    analysis_result = {}

    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        
        qubit_info.measure.pulse_duration(2e-6)
        qubit_info.measure.integration_time(1.5e-6)
        qubit_info.reset.duration(250e-6)
        qubit_info.measure.pulse_amp(ro_elements[q]["ro_amp_factor"]*float(qubit_info.measure.pulse_amp()))
        if qubit_info.rxy.amp180() is NaN:
            qubit_info.rxy.amp180(0)
        if qubit_info.rxy.duration() is NaN:
            qubit_info.rxy.duration(0)

    sched_kwargs = dict(   
        ini_state='g',
        pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements,'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
    )
    exp_kwargs= dict(shots=shots,)
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(shots)
        iq_tuples = gettable.get() # tuple?
        for q_idx, q in enumerate(ro_elements):
            IQ_array = array([iq_tuples[2*q_idx],iq_tuples[2*q_idx+1]])
            analysis_result[q] = Single_shot_ref_fit_analysis(IQ_array)
        
        # Data_manager().save_raw_data(QD_agent=QD_agent,ds=ss_ds,qb=q,exp_type='SS')
        

        
    else:
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        
        
    return analysis_result

def refIQ_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict,ro_elements:dict,run:bool=True,shots_num:int=50000):
    
    if run:
        for qu in ro_elements:
            if ro_elements[qu]["specified_bias"] == 0:
                Fctrl[qu](float(QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=qu)))
            else:
                Fctrl[qu](ro_elements[qu]["specified_bias"])
        analysis_result = Single_shot_ref_spec(QD_agent,ro_elements,shots=shots_num)
        for qu in ro_elements:
            Fctrl[qu](0.0)
        cluster.reset()
        try :
            for qu in ro_elements:
                I_ref, Q_ref= analysis_result[qu]['fit_pack'][0],analysis_result[qu]['fit_pack'][1]
                QD_agent.memo_refIQ({str(qu):[I_ref,Q_ref]})
                Single_shot_fit_plot(analysis_result[qu],qu)
        except:
            shut_down(cluster,Fctrl)
            raise ValueError ("Analysis goes wrong!")
    else:
        analysis_result = Single_shot_ref_spec(QD_agent,ro_elements,shots=shots_num,run=False)
    
    
if __name__ == "__main__":
    
    """ Fill in """
    execution = True
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        'q0':{"ro_amp_factor":1,"specified_bias":0},
        # 'q1':{"ro_amp_factor":1,"specified_bias":0}
    } # If 'specified_bias' is not zero, it will bias with it while the exp running.
                
    couplers = []
    shots:int = 50000

    
    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    check_mark = False

    """ Running """
    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
        if ro_elements[qubit]["ro_amp_factor"] != 1: check_mark = True
    
    refIQ_executor(QD_agent,cluster,Fctrl,ro_elements,run=execution,shots_num=shots)
    
    if check_mark:
        keep = mark_input(f"Keep these RO amp for all the qubits ?[y/n]")
    else:
        keep = 'y'


    """ Storing """
    if execution:
        if keep.lower() in ["y", "yes"]:
            QD_agent.refresh_log("After IQ ref checking!")
            QD_agent.QD_keeper()


    """ Close """
    print('IQ ref checking done!')
    shut_down(cluster,Fctrl)

