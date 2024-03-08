
from numpy import linspace, array, average
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, save_raw_data
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis


def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,arti_detune:int=0,IF:int=130e6,n_avg:int=300,points:int=201,run:bool=True,q='q1',times=1, ref_IQ:list=[0,0],Experi_info:dict={}):
    analysis_result= {}
    T2_us = {}
    analysis_result[q]= []
    T2_us[q] = []
    Real_detune= {}
    Real_detune[q]= []
    qubit = QD_agent.quantum_device.get_element(q)
    New_fxy= qubit.clock_freqs.f01()+arti_detune
    LO= New_fxy+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    samples = linspace(0,freeduration,points)
    sche_func= Ramsey_sche
    sched_kwargs = dict(
        q=q,
        pi_amp={str(q):qubit.rxy.amp180()},
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        R_amp={str(q):qubit.measure.pulse_amp()},
        R_duration={str(q):qubit.measure.pulse_duration()},
        R_integration={str(q):qubit.measure.integration_time()},
        R_inte_delay=qubit.measure.acq_delay(),
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                     f_xy='%E' %sched_kwargs['New_fxy'],
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
            ramsey_ds = meas_ctrl.run('Ramsey')

            # Save the raw data into netCDF
            save_raw_data(QD_agent,ramsey_ds,histo_label=i,qb=q,exp_type='T2')

            I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
            data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])
            data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=8e-6)
            analysis_result[q].append(data_fit)
            T2_us[q].append(data_fit.attrs['T2_fit']*1e6)
            Real_detune[q].append(data_fit.attrs['f']-arti_detune)
        
        average_fit_detune= average(array(Real_detune[q]))
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    return analysis_result, T2_us, average_fit_detune


if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, save_histo_pic, shut_down, reset_offset
    from Pulse_schedule_library import Fit_analysis_plot
    from numpy import mean

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
        linecut = 0
        Ramsey_results,T2_hist, average_actual_detune= Ramsey(QD_agent,meas_ctrl,freeduration=30e-6,n_avg=500,q=qb,times=10,ref_IQ=QD_agent.refIQ[qb])
        Fit_analysis_plot(Ramsey_results[qb][linecut],P_rescale=False,Dis=None)
        # set the histo save path
        save_histo_pic(QD_agent,T2_hist,qb,mode="t2")
        print('average_actual_detune=', average_actual_detune*1e-6,'MHz')
        print('average_T2=', mean(array(T2_hist[qb])),'us')

    
    print('T2 done!')
    shut_down(cluster,Fctrl)