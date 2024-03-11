
from numpy import linspace, array
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Pulse_schedule_library import T1_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T1_fit_analysis

def T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float=80e-6,IF:int=150e6,n_avg:int=300,points:int=200,run:bool=True,q='q1',times:int=1, Experi_info:dict={},ref_IQ:list=[0,0]):
    analysis_result= {}
    T1_us = {}
    analysis_result[q]= []
    T1_us[q] = []
    sche_func=T1_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    LO= qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    samples = linspace(0,freeduration,points)
    
    sched_kwargs = dict(
        q=q,
        pi_amp={str(q):qubit_info.rxy.amp180()},
        freeduration=Para_free_Du,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                     )
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
        meas_ctrl.settables(Para_free_Du)
        meas_ctrl.setpoints(samples)
        for i in range(times):
            T1_ds = meas_ctrl.run('T1')
            # Save the raw data into netCDF
            Data_manager.save_raw_data(QD_agent,T1_ds,histo_label=i,qb=q,exp_type='T1')
            I,Q= dataset_to_array(dataset=T1_ds,dims=1)
            data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
            data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=8e-6)
            analysis_result[q].append(data_fit)
            T1_us[q].append(data_fit.attrs['T1_fit']*1e6)
             
        show_args(exp_kwargs, title="T1_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="T1_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
  
    return analysis_result, T1_us




if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Pulse_schedule_library import Fit_analysis_plot
    from numpy import mean
    import os
    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_5/DR1#170_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    
    # Set system attenuation
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=20)
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50) 

    execute = True
    error_log = []
    for qb in Fctrl:
        print(f"{qb} are under the measurement ...")
        T1_results, T1_hist = T1(QD_agent,meas_ctrl,q=qb,times=10,ref_IQ=QD_agent.refIQ[qb])
        if T1_results == {}:
            error_log.append(qb)
        else:
            if execute:
                linecut=0
                qubit = QD_agent.quantum_device.get_element(qb)
                Fit_analysis_plot(T1_results[qb][linecut],P_rescale=False,Dis=None)
                # set the histo save path
                Data_manager.save_histo_pic(QD_agent,T1_hist,qb,mode="t1")
                print(f'{qb} average_T1= {mean(array(T1_hist[qb]))} us')
    
    print('T1 done!')
    shut_down(cluster,Fctrl)