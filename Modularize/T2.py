def Ramsey(quantum_device:QuantumDevice,ro_elements:list,pi_amp:dict,detuning:float,freeduration:float,n_avg:int=300,points:int=201,run:bool=True,q='q1',times=1):
    analysis_result= {}
    T2_us = {}
    analysis_result[q]= []
    T2_us[q] = []
    Real_detune= {}
    Real_detune[q]= []
    New_fxy= f01[q]+detuning
    LO= New_fxy+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    samples = np.linspace(0,freeduration,points)
    sche_func= Ramsey_sche
    sched_kwargs = dict(
        q=q,
        pi_amp=pi_amp,
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                     f_xy='%E' %sched_kwargs['New_fxy'],
                     )
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
        )
        
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Para_free_Du)
        meas_ctrl.setpoints(samples)
        
        for i in range(times):
            ramsey_ds = meas_ctrl.run('Ramsey')
            ramsey_ds
            I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
            data= IQ_data_dis(I,Q,ref_I=I_ref,ref_Q=Q_ref)
            data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=8e-6)
            analysis_result[q].append(data_fit)
            T2_us[q].append(data_fit.attrs['T2_fit']*1e6)
            Real_detune[q].append(data_fit.attrs['f']-detuning)
        
        average_fit_detune= np.average(np.array(Real_detune[q]))
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    else:
        sweep_para= np.array(samples[:n_s])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    return analysis_result, T2_us, average_fit_detune