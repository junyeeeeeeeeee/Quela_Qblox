
from Pulse_schedule_library import *

n_s=2 # number of schedule in pulse_preview
#%%

def Cavity_spec(quantum_device:QuantumDevice,ro_bare_guess:dict,ro_span_Hz:int=5e6,n_avg:int=1000,points:int=200,run:bool=True,q:str='q1'):

    sche_func = One_tone_sche
        
    analysis_result = {}
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        powerDep=False,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     )
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(ro_f_samples)
        
        
        
        rs_ds = meas_ctrl.run("One-tone")
        rs_ds
        analysis_result[q] = ResonatorSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
        print(f"{q} Cavity:")
        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    else:
        sweep_para= np.array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    return analysis_result


execute = True

if execute:
    Cavityspec_results = Cavity_spec(quantum_device,
                               ro_bare_guess=ro_bare,
                               ro_span_Hz=2e6,
                               n_avg=500,
                               points=51,
                               run=True,
                               q='q1')
#%%
plot_data =True

if plot_data:
    for q in Cavityspec_results:
        print('Qi=',Cavityspec_results[q].quantities_of_interest['Qi'])
        print('Qc=',Cavityspec_results[q].quantities_of_interest['Qc'])
        print('Ql=',Cavityspec_results[q].quantities_of_interest['Ql'])
        Cavityspec_results[q].display_figs_mpl()
        Amp_phase_plot(quantum_device, Cavityspec_results[q],title="One_tone_spec")
    

#%%

def Cavity_powerDep_spec(quantum_device:QuantumDevice,ro_elements:list,ro_bare_guess:dict,ro_span_Hz:int=5e6,ro_p_min:float=0.1,ro_p_max:float=0.7,n_avg:int=1000,f_points:int=200,p_points:int=200,run:bool=True,q:str='q1'):

    sche_func = One_tone_sche
        
    analysis_result = {}
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,f_points)
    ro_p_samples = linspace(ro_p_min,ro_p_max,p_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
    ro_pulse_amp.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp=ro_pulse_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        powerDep=True,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     Power=['start '+'%E' %ro_p_samples[0],'end '+'%E' %ro_p_samples[-1]])
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((ro_f_samples,ro_p_samples))
        
        
        
        rp_ds = meas_ctrl.run("One-tone-powerDep")
        rp_ds
        analysis_result[q] = Basic2DAnalysis(tuid=rp_ds.attrs["tuid"], dataset=rp_ds).run()
        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    else:
        sweep_para1= np.array(ro_f_samples[:n_s])
        sweep_para2= np.array(ro_p_samples[:2])
        spec_sched_kwargs['frequencies']= sweep_para1.reshape(sweep_para1.shape or (1,))
        spec_sched_kwargs['R_amp']= {q:sweep_para2.reshape(sweep_para2.shape or (1,))[0]}
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    return analysis_result


execute = True

if execute:
    Cavity_power_spec_results = Cavity_powerDep_spec(quantum_device,
                               ro_element,
                               ro_bare_guess=ro_bare,
                               ro_span_Hz=5e6,
                               ro_p_min=0.1,
                               ro_p_max=0.7,
                               n_avg=500,
                               f_points=101,
                               p_points=21,
                               run=True,
                               q='q1')
#%%    
plot_data =True
update= True  #update readout amplitude
desired_R_amp=0.1

if plot_data:
    for q in Cavity_power_spec_results:
        Cavity_power_spec_results[q].display_figs_mpl()
        
if update:
    R_amp[q]= desired_R_amp
            
#%%

def Cavity_FluxDep_spec(quantum_device:QuantumDevice,ro_elements:list,ro_bare_guess:dict,ro_span_Hz:int=5e6,flux_span:float=0.3,n_avg:int=1000,f_points:int=200,flux_points:int=200,run:bool=True,q:str='q1'):

    sche_func = One_tone_sche
        
    analysis_result = {}
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,f_points)
    flux_samples = linspace(-flux_span,flux_span,flux_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        powerDep=False,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     Flux=['start '+'%E' %flux_samples[0],'end '+'%E' %flux_samples[-1]])
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,flux_settable_map[q]])
        meas_ctrl.setpoints_grid((ro_f_samples,flux_samples))
        
        
        
        rfs_ds = meas_ctrl.run("One-tone-Flux")
        rfs_ds
        analysis_result[q] = ResonatorFluxSpectroscopyAnalysis(tuid=rfs_ds.attrs["tuid"], dataset=rfs_ds).run()
        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    else:
        sweep_para= np.array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    return analysis_result


execute = True

if execute:
    FD_results = Cavity_FluxDep_spec(quantum_device,
                               ro_element,
                               ro_bare_guess=ro_bare,
                               ro_span_Hz=2e6,
                               flux_span=0.2, #center at zero bias
                               n_avg=300,
                               f_points=51,
                               flux_points=21,
                               run=True,
                               q='q1')
    
#%%    
plot_data =True
bias_at_sweet_spot=True 

if plot_data:
    for q in FD_results:
        Amp_phase_plot(quantum_device, FD_results[q],title="One_tone_FluxDep_spec")
        FD_results[q].display_figs_mpl()
        if bias_at_sweet_spot:
            sweet_spot_offset= FD_results[q].quantities_of_interest["offset_0"].nominal_value
            sweet_spot_R_F= FD_results[q].quantities_of_interest["freq_0"]
            qubit = quantum_device.get_element(q)
            qubit.clock_freqs.readout(sweet_spot_R_F)
            flux_settable_map[q](sweet_spot_offset)
            bias[q]=sweet_spot_offset
            R_F[q]= sweet_spot_R_F
        show_args(Experi_info(q),title='Update parameters')   
         
#%%

def Single_shot_ref_check(quantum_device:QuantumDevice,shots:int=1000,run:bool=True,q:str='q1'):
    
    sche_func = Qubit_SS_sche   
    analysis_result = {}

    sched_kwargs = dict(   
        q=q,
        ini_state='g',
        pi_amp=pi_amp,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay
    )
    exp_kwargs= dict(shots=shots,
                     )
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(shots)
        ss_ds= gettable.get()
        
        analysis_result[q] = Single_shot_ref_fit_analysis(ss_ds)
        
        show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    else:
        pulse_preview(quantum_device,sche_func,sched_kwargs)
        
        show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    return analysis_result


execute = True

if execute: 
    SS_ref_result= Single_shot_ref_check(quantum_device,
                shots=40000,
                run=True,
                q='q1')



#%%
plot_data =True
update= True
if plot_data:
   for q in SS_ref_result:
       Single_shot_fit_plot(SS_ref_result[q])
       if update:
           I_ref, Q_ref= SS_ref_result[q]['fit_pack'][0],SS_ref_result[q]['fit_pack'][1]
           show_args(Experi_info(q),title='Update parameters')



#%%

def Two_tone_spec(quantum_device:QuantumDevice,ro_elements:list,f01_guess:dict,xyf_span_Hz:int=200e6,xyamp:float=0.02,n_avg:int=1000,points:int=200,run:bool=True,q:str='q1'):
    
    sche_func = Two_tone_sche   
    analysis_result = {}
    f01_high = f01_guess[q]+100*1e6
    f01_samples = linspace(f01_high-xyf_span_Hz,f01_high,points)
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=f01_high)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True

    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        spec_amp=xyamp,
        spec_Du=50*1e-6,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(f01_samples)
        
        qs_ds = meas_ctrl.run("Two-tone")
        qs_ds
        I,Q= dataset_to_array(dataset=qs_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=I_ref,ref_Q=Q_ref)
        analysis_result[q] = QS_fit_analysis(data,f=f01_samples)
        
        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    else:
        sweep_para= np.array(f01_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    return analysis_result


execute = True

if execute: #f01_guess[q] is the predicted f01 for q, f_LO = f01_guess[q]+100*1e6
    QS_results = Two_tone_spec(quantum_device,
                               ro_element,
                               f01_guess=f01_guess,
                               xyf_span_Hz=200e6,
                               xyamp=0.01,
                               n_avg=1000,
                               points=101,
                               run=True,
                               q='q1')
    
#%%
plot_data =True
update= False  #update f01 here   
if plot_data:
    for q in QS_results:
        Fit_analysis_plot(QS_results[q],P_rescale=True,Dis=0.003048161)
        if update:
            Revised_f01= QS_results[q].attrs['f01_fit']
            qubit = quantum_device.get_element(q)
            qubit.clock_freqs.f01(Revised_f01)
            f01[q]= Revised_f01
        show_args(Experi_info(q),title='Update parameters')      
        
#%%

def Zgate_two_tone_spec(quantum_device:QuantumDevice,ro_elements:list,Z_amp_start:float,Z_amp_end:float,xyf_highest:float,xyf_span_Hz:float=200e6,xyamp:float=0.02,n_avg:int=1000,Z_points:int=200,f_points:int=200,run:bool=True,q:str='q1'):
    
    sche_func = Z_gate_two_tone_sche
        
    analysis_result = {}
    data={}
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=xyf_highest)
    f01_samples = linspace(xyf_highest-xyf_span_Hz,xyf_highest,f_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples = linspace(Z_amp_start,Z_amp_end,Z_points)
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        Z_amp=Z_bias,
        spec_amp=xyamp,
        spec_Du=50*1e-6,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,Z_bias])
        meas_ctrl.setpoints_grid((f01_samples,Z_samples))
        qs_ds = meas_ctrl.run("Zgate-two-tone")
        qs_ds
        
        analysis_result[q] = QubitFluxSpectroscopyAnalysis(tuid=qs_ds.attrs["tuid"], dataset=qs_ds).run()
        
        
        show_args(exp_kwargs, title="Zgate_two_tone_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    else:
        sweep_para1= np.array(f01_samples[:n_s])
        sweep_para2= np.array(Z_samples[:2])
        spec_sched_kwargs['frequencies']= sweep_para1.reshape(sweep_para1.shape or (1,))
        spec_sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        
        
        show_args(exp_kwargs, title="Zgate_two_tone_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    return analysis_result


execute = True

if execute:
    QS_Z_results = Zgate_two_tone_spec(quantum_device,
                               ro_element,
                               Z_amp_start=-0.05,
                               Z_amp_end=0.05,
                               xyf_highest=4.2e9,
                               xyf_span_Hz=300e6,
                               xyamp=0.01,
                               n_avg=500,
                               Z_points=6,
                               f_points=51,
                               run=True,
                               q='q1')
    
    
 #%%
plot_data =True
update= False   # Update the flux sweetspot
if plot_data:
    for q in QS_Z_results:
        Amp_phase_plot(quantum_device, QS_Z_results[q],title='Zgate_two_tone_spec')
        QS_Z_results[q].display_figs_mpl()
        if update: 
            flux_settable(qfs_analysis.quantities_of_interest["offset_0"].nominal_value)
        show_args(Experi_info(q),title='Update parameters') 



        