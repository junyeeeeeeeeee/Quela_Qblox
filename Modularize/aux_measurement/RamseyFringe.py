import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
import matplotlib.pyplot as plt
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from numpy import std, arange, array, average, mean, ndarray, linspace
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot, Fit_T2_cali_analysis_plot, T1_fit_analysis
from Modularize.analysis.RabiChevAna import plot_fringe

def RamseyFringe(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,detuning:int,freq_pts:int,IF:int=250e6,n_avg:int=1000,points:int=101,run:bool=True,q='q1', ref_IQ:list=[0,0],Experi_info:dict={},exp_idx:int=0,data_folder:str=''):
    
    qubit_info = QD_agent.quantum_device.get_element(q)

    
    # qubit_info.reset.duration(qubit_info.reset.duration()*2)
    print("Integration time ",qubit_info.measure.integration_time()*1e6, "µs")
    print("Reset time ", qubit_info.reset.duration()*1e6, "µs")
    # Manually change f01
    # f01 = qubit.clock_freqs.f01()
    # qubit.clock_freqs.f01(f01-2.47e6)
    
    xyf = ManualParameter(name="xyf", unit="Hz", label="XY Frequency")
    xyf.batched = False
    f01_samples = linspace(qubit_info.clock_freqs.f01()-detuning,qubit_info.clock_freqs.f01()+detuning,freq_pts)
    
    LO= f01_samples[-1]+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    
    gap = (freeduration)*1e9 // points + (((freeduration)*1e9 // points) %4) # multiple by 8 ns
    
    
    samples = arange(0,freeduration,gap*1e-9)
    samples = modify_time_point(samples, 4e-9)
    
    sche_func= Ramsey_sche
    sched_kwargs = dict(
        q=q,
        pi_amp={str(q):qubit_info.rxy.amp180()},
        New_fxy=xyf,
        freeduration=Para_free_Du,
        pi_dura=qubit_info.rxy.duration(),
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        echo_pi_num=0
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]])
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
        meas_ctrl.settables([Para_free_Du,xyf])
        meas_ctrl.setpoints_grid((samples,f01_samples))
        
        
        ramsey_ds = meas_ctrl.run('Ramsey')

        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=ramsey_ds,label=exp_idx,qb=q,exp_type='fringe',specific_dataFolder=data_folder,get_data_loc=True)
        
        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        # n_s = 2
        sweep_para= array([samples[0],samples[-1]])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        nc_path = ''
    return nc_path


def modify_time_point(ary:ndarray,factor1:int, factor2:int=0):
    x = []
    for i in ary:
        if i % factor1 == 0:
            ii = i

        else:
            multiples = i // factor1
            ii = factor1*(multiples+1)
        
        if factor2 != 0:
            if ii % factor2 == 0: 
                    if ii not in x :
                        x.append(ii)
                    else:
                        pass
            else:
                multiples = ii // factor2
                multiples_of_factor = factor2*(multiples+1)

                if multiples_of_factor % factor1 == 0:
                    if multiples_of_factor not in x :
                        x.append(multiples_of_factor)
                    else:
                        pass
        else:
            x.append(ii)


    return array(x)


def fringe_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,detune:float,freq_pts:int,freeDura:float=30e-6,ith:int=1,run:bool=True,specific_folder:str='',pts:int=100, avg_n:int=800, IF:float=250e6):
    if run:
        qubit_info = QD_agent.quantum_device.get_element(specific_qubits)
        ori_reset = qubit_info.reset.duration()
        qubit_info.reset.duration(qubit_info.reset.duration()+freeDura)
    
        slightly_print(f"The {ith}-th T2:")
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        
        nc_path = RamseyFringe(QD_agent,meas_ctrl,detuning=detune,freq_pts=freq_pts,freeduration=freeDura,n_avg=avg_n,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=pts,run=True,exp_idx=ith,data_folder=specific_folder,IF=IF)
        Fctrl[specific_qubits](0.0)
        
        cluster.reset()
        qubit_info.reset.duration(ori_reset)
        
    else:
        nc_path = RamseyFringe(QD_agent,meas_ctrl,detuning=detune,freq_pts=freq_pts,freeduration=freeDura,n_avg=1000,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=False)


    return nc_path


if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = {
        "q0":{"detune":2e6,"evoT":10e-6},#-0.174e6
    }
    couplers = []

    """ Optional paras """
    time_data_points = 100
    freq_pts:int = 20
    avg_n = 1000
    xy_IF = 250e6
    adj_freq = 0e6


    """ Iteration """
    ncs = {}
    for q_idx, qubit in enumerate(list(ro_elements.keys())):
        
        start_time = time.time()
        """ Preparations """
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
        QD_agent.quantum_device.get_element(qubit).clock_freqs.f01(QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()+adj_freq)
        
        """ Running """
        Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        slightly_print(f"Ramsey with detuning = {round(ro_elements[qubit]['detune']*1e-6,2)} MHz")
        ncs[qubit] = fringe_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,detune=ro_elements[qubit]["detune"],freq_pts=freq_pts,freeDura=ro_elements[qubit]["evoT"],ith=0,run=execution,pts=time_data_points,avg_n=avg_n,IF=xy_IF)
        
        """ Storing """
        if q_idx == len(ro_elements)-1:
            for q in ncs:
                plot_fringe(QD_agent,q,ncs[q])

        """ Close """
        print('Ramsey Fringe done!')
        shut_down(cluster,Fctrl,Cctrl)
        end_time = time.time()
        slightly_print(f"time cost: {round(end_time-start_time,1)} secs")
            
        
        
            
            
        


    