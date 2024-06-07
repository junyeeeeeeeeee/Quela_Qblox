import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import matplotlib.pyplot as plt
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from Modularize.support.UserFriend import *
from qcodes.parameters import ManualParameter
from Modularize.m7_RefIQ import refIQ_executor
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from numpy import linspace, array, where, max, ndarray, sqrt, arctan2
from Modularize.support import QDmanager, Data_manager, init_meas, shut_down, init_system_atte,coupler_zctrl
from Modularize.support.Pulse_schedule_library import ROF_Cali_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array


def rofCali(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_span_Hz:float=3e6,IF:int=150e6,n_avg:int=500,f_points:int=101,run:bool=True,q='q1',Experi_info:dict={},data_folder:str=''):
    
    analysis_result = {}
    qubit = QD_agent.quantum_device.get_element(q)

    ro_f_origin= qubit.clock_freqs.readout()
    LO= ro_f_origin+IF+ro_span_Hz
    from numpy import NaN
    qubit.clock_freqs.readout(NaN)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='readout',LO_frequency=LO)
    
    option_rof = ManualParameter(name="freq", unit="Hz", label="Frequency")
    option_rof.batched = True
    ro_f_samples = linspace(ro_f_origin-ro_span_Hz,ro_f_origin+ro_span_Hz,f_points)
    
    sche_func= ROF_Cali_sche

    def state_dep_sched(ini_state:str):
        sched_kwargs = dict(
            q=q,
            ro_freq=option_rof,
            ini_state=ini_state,
            pi_amp={str(q):qubit.rxy.amp180()},
            pi_dura={str(q):qubit.rxy.duration()},
            R_amp={str(q):qubit.measure.pulse_amp()},
            R_duration={str(q):qubit.measure.pulse_duration()},
            R_integration={str(q):qubit.measure.integration_time()},
            R_inte_delay=qubit.measure.acq_delay(),
            )
        exp_kwargs= dict(sweep_ROF=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]])
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func,
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
            )
            
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            meas_ctrl.gettables(gettable)
            meas_ctrl.settables(option_rof)
            meas_ctrl.setpoints(ro_f_samples)
            
            
            rfs_ds = meas_ctrl.run("Rof-Calibrate")
            # Save the raw data into netCDF
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=rfs_ds,qb=q,label=ini_state,exp_type='RofCali',specific_dataFolder=data_folder)
            
            I,Q= dataset_to_array(dataset=rfs_ds,dims=1)
            
            show_args(exp_kwargs, title="RofCali_kwargs: Meas.qubit="+q)

        else:
            n_s = 2
            sweep_para= array(ro_f_samples[:n_s])
            sched_kwargs['ro_freq']= sweep_para.reshape(sweep_para.shape or (1,))
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
            

            show_args(exp_kwargs, title="RofCali_kwargs: Meas.qubit="+q)
            if Experi_info != {}:
                show_args(Experi_info(q))
            I, Q = [], []
        
        return array(I), array(Q)
    
    slightly_print("Running |1>")
    I_e, Q_e = array(state_dep_sched('e'))
    slightly_print("Running |0>")
    I_g, Q_g = array(state_dep_sched('g'))
    I_diff = I_e-I_g
    Q_diff = Q_e-Q_g
    dis_diff = sqrt((I_diff)**2+(Q_diff)**2)
    if run:
        great_rof = anal_rof_cali(I_e,Q_e,I_g,Q_g,dis_diff,ro_f_samples,ro_f_origin)
    else:
        great_rof = 0
    qubit.clock_freqs.readout(ro_f_origin)
    return great_rof


def anal_rof_cali(I_e:ndarray,Q_e:ndarray,I_g:ndarray,Q_g:ndarray,dis_diff:ndarray,ro_f_samples:ndarray,original_rof=float):
    max_dif_idx = where(dis_diff==max(dis_diff))[0][0]
    optimized_rof = ro_f_samples[max_dif_idx]
    fig, ax = plt.subplots(3,1)
    amp_e = sqrt(I_e**2+Q_e**2)
    amp_g = sqrt(I_g**2+Q_g**2)
    ax[0].plot(ro_f_samples,amp_e,label='|1>')
    ax[0].plot(ro_f_samples,amp_g,label='|0>')
    ax[0].set_xlabel('ROF (Hz)')
    ax[0].set_ylabel("Amplitude (V)")
    ax[0].legend()

    pha_e = arctan2(Q_e,I_e)
    pha_g = arctan2(Q_g,I_g)
    ax[1].plot(ro_f_samples,pha_e,label='|1>')
    ax[1].plot(ro_f_samples,pha_g,label='|0>')
    ax[1].set_xlabel('ROF (Hz)')
    ax[1].set_ylabel("phase")
    ax[1].legend()

    ax[2].plot(ro_f_samples,dis_diff,label='diff')
    ax[2].vlines(x=array([optimized_rof]),ymin=min(dis_diff),ymax=max(dis_diff),colors='black',linestyles='--',label='optimal')
    ax[2].vlines(x=array([original_rof]),ymin=min(dis_diff),ymax=max(dis_diff),colors='#DCDCDC',linestyles='--',label='original')
    ax[2].set_xlabel('ROF (Hz)')
    ax[2].set_ylabel("diff (V)")
    ax[2].legend()
    plt.show()

    return optimized_rof
    

def rofCali_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,execution:bool=True,ro_f_span:float=2e6,fpts:int=100):
    if execution:
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        optimal_rof = rofCali(QD_agent,meas_ctrl,ro_span_Hz=ro_f_span,q=specific_qubits,f_points=fpts,run=execution)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
    else:
        optimal_rof = rofCali(QD_agent,meas_ctrl,ro_span_Hz=ro_f_span,q=specific_qubits,f_points=fpts,run=execution)

    return optimal_rof



if __name__ == '__main__':

    """ Fill in """
    execute = True
    DRandIP = {"dr":"dr1","last_ip":"11"}
    ro_elements = {'q0':{"span_Hz":8e6}}
    couplers = ['c0','c1']

    """ Preparation """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    


    """ Running """
    keep = False
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        ro_span = ro_elements[qubit]["span_Hz"]

        optimal_rof = rofCali_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,execution=execute,ro_f_span=ro_span)
        if execute:
            if mark_input(f"Update the optimal ROF for {qubit}?[y/n]").lower() in ['y', 'yes']:
                QD_agent.quantum_device.get_element(qubit).clock_freqs.readout(optimal_rof)
                refIQ_executor(QD_agent,cluster,Fctrl,qubit)
                keep = True

        """ Storing """ 
        if execute:
            if keep:
                QD_agent.QD_keeper()
                keep = False 

    """ Close """    
    shut_down(cluster,Fctrl,Cctrl)

