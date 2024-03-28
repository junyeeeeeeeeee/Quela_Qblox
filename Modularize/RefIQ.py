
from utils.tutorial_utils import show_args
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.Pulse_schedule_library import Qubit_SS_sche, Single_shot_ref_fit_analysis, pulse_preview, Single_shot_fit_plot

def Single_shot_ref_spec(QD_agent:QDmanager,shots:int=1000,run:bool=True,q:str='q1',Experi_info:dict={},want_state:str='g'):
    print("Single shot start")
    sche_func = Qubit_SS_sche   
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    # qubit_info.clock_freqs.readout(5.7225e9)
    if want_state == 'g':
        XYL = 0
    else:
        XYL = qubit_info.rxy.amp180()
    sched_kwargs = dict(   
        q=q,
        ini_state=want_state,
        pi_amp={str(q):XYL},
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

def refIQ_executor(QD_agent:QDmanager,Fctrl:dict,specific_qubits:str,run:bool=True):

    if run:
        init_system_atte(QD_agent.quantum_device,list([specific_qubits]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'ro'))
        
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_sweetBiasFor(target_q=specific_qubits)))
        analysis_result = Single_shot_ref_spec(QD_agent,q=specific_qubits,want_state='g',shots=10000)
        Fctrl[specific_qubits](0.0)
        try :
            I_ref, Q_ref= analysis_result[specific_qubits]['fit_pack'][0],analysis_result[specific_qubits]['fit_pack'][1]
            QD_agent.memo_refIQ({str(specific_qubits):[I_ref,Q_ref]})
            Single_shot_fit_plot(analysis_result[specific_qubits])
        except:
            shut_down(cluster,Fctrl)
            raise ValueError ("Analysis goes wrong!")
        
    else:
        analysis_result = Single_shot_ref_spec(QD_agent,q=specific_qubits,want_state='g',shots=10000,run=False)

    


if __name__ == "__main__":
    
    """ Fill in """
    execution = True
    QD_path = 'Modularize/QD_backup/2024_3_28/DR2#171_SumInfo.pkl'
    ro_elements = ["q4"]


    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')


    """ Running """
    for qubit in ro_elements:
        refIQ_executor(QD_agent,Fctrl,specific_qubits=qubit,run=execution)


    """ Storing """
    if execution:
        QD_agent.refresh_log("After IQ ref checking!")
        QD_agent.QD_keeper()


    """ Close """
    print('IQ ref checking done!')
    shut_down(cluster,Fctrl)

