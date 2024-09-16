import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Qubit_SS_sche, Single_shot_ref_fit_analysis, pulse_preview, Single_shot_fit_plot

def Single_shot_ref_spec(QD_agent:QDmanager,shots:int=1000,run:bool=True,q:str='q1',Experi_info:dict={},want_state:str='g',ro_amp_scaling:float=1, thermal_pop_mode:bool=True):
    print("Single shot start")
    sche_func = Qubit_SS_sche   
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    if thermal_pop_mode:
        qubit_info.measure.pulse_duration(2e-6)
        qubit_info.measure.integration_time(1.5e-6)
        qubit_info.reset.duration(250e-6)
    else:
        qubit_info.measure.pulse_duration(100e-6)
        qubit_info.measure.integration_time(100e-6-1e-6)
        qubit_info.reset.duration(250e-6)
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

def refIQ_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict,specific_qubits:str,run:bool=True,ro_amp_adj:float=1,shots_num:int=50000,want_see_p01:bool=True,specify_bias=0):

    if run:
        if specify_bias == 0:
            Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=specific_qubits)))
        else:
            Fctrl[specific_qubits](specify_bias)
        analysis_result = Single_shot_ref_spec(QD_agent,q=specific_qubits,want_state='g',shots=shots_num,ro_amp_scaling=ro_amp_adj,thermal_pop_mode=want_see_p01)
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
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = {'q0':{"ro_amp_factor":1.2},}
                
    couplers = []


    for qubit in ro_elements:
        """ Preparations """
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
        

        """ Running """
        Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
        refIQ_executor(QD_agent,cluster,Fctrl,specific_qubits=qubit,run=execution,ro_amp_adj=ro_elements[qubit]["ro_amp_factor"])
        
        if ro_elements[qubit]["ro_amp_factor"] !=1:
            keep = mark_input(f"Keep this RO amp for {qubit}?[y/n]")
        else:
            keep = 'y'


        """ Storing """
        if execution:
            if keep.lower() in ["y", "yes"]:
                QD_agent.refresh_log("After IQ ref checking!")
                QD_agent.QD_keeper()


        """ Close """
        print('IQ ref checking done!')
        shut_down(cluster,Fctrl,Cctrl)

