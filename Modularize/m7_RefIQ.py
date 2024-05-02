import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.Pulse_schedule_library import Qubit_SS_sche, Single_shot_ref_fit_analysis, pulse_preview, Single_shot_fit_plot

def Single_shot_ref_spec(QD_agent:QDmanager,shots:int=1000,run:bool=True,q:str='q1',Experi_info:dict={},want_state:str='g',ro_amp_scaling:float=1):
    print("Single shot start")
    sche_func = Qubit_SS_sche   
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    qubit_info.measure.pulse_amp(ro_amp_scaling*float(qubit_info.measure.pulse_amp()))

    # qubit_info.clock_freqs.readout(5.7225e9)
    if want_state == 'g':
        XYL = 0
    else:
        XYL = qubit_info.rxy.amp180()
    sched_kwargs = dict(   
        q=q,
        ini_state=want_state,
        pi_amp={str(q):XYL},
        pi_dura={str(q):4e-9},
        R_amp={str(q):qubit_info.measure.pulse_amp()}, #
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
    )
    exp_kwargs= dict(shots=shots,
                     )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(shots)
        ss_ds= gettable.get() # tuple?
        # Data_manager().save_raw_data(QD_agent=QD_agent,ds=ss_ds,qb=q,exp_type='SS')
        analysis_result[q] = Single_shot_ref_fit_analysis(ss_ds)
        
        show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        
        show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result

def refIQ_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict,specific_qubits:str,run:bool=True,ro_amp_adj:float=1,shots_num:int=7000):

    if run:
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_tuneawayBiasFor(target_q=specific_qubits)))
        analysis_result = Single_shot_ref_spec(QD_agent,q=specific_qubits,want_state='g',shots=shots_num,ro_amp_scaling=ro_amp_adj)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
        try :
            I_ref, Q_ref= analysis_result[specific_qubits]['fit_pack'][0],analysis_result[specific_qubits]['fit_pack'][1]
            QD_agent.memo_refIQ({str(specific_qubits):[I_ref,Q_ref]})
            Single_shot_fit_plot(analysis_result[specific_qubits])
        except:
            shut_down(cluster,Fctrl)
            raise ValueError ("Analysis goes wrong!")
        
    else:
        analysis_result = Single_shot_ref_spec(QD_agent,q=specific_qubits,want_state='g',shots=shots_num,run=False)

    


if __name__ == "__main__":
    
    """ Fill in """
    execution = True
    DRandIP = {"dr":"dr1","last_ip":"11"}
    ro_elements = {'q0':{"ro_amp_factor":1}}


    """ Preparations """
    QD_path = "Modularize/QD_backup/2024_4_29/DR1#11_SumInfo-44G.pkl"#find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    if ro_elements == 'all':
        ro_elements = list(Fctrl.keys())


    """ Running """
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
        refIQ_executor(QD_agent,cluster,Fctrl,specific_qubits=qubit,run=execution,ro_amp_adj=ro_elements[qubit]["ro_amp_factor"])
        
        if ro_elements[qubit]["ro_amp_factor"] !=1:
            keep = input(f"Keep this RO amp for {qubit}?[y/n]")
        else:
            keep = 'y'


    """ Storing """
    if execution:
        if keep:
            QD_agent.refresh_log("After IQ ref checking!")
            QD_agent.QD_keeper(special_path=QD_path)


    """ Close """
    print('IQ ref checking done!')
    shut_down(cluster,Fctrl)

