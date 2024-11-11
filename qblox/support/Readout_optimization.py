from Pulse_schedule_library import *

n_s=2 # number of schedule in pulse_preview

IF= 200*1e6 #IF frequency used in time-domain
#%% Readout update panel

update= True
desired_R_amp=0.12
desired_R_duration= 2* 1e-6
desired_R_integration=2* 1e-6
desired_R_inte_delay=0#192*1e-9
q='q1'

if update:
    R_amp[q]= desired_R_amp
    R_integration[q]= desired_R_integration
    R_duration[q]= desired_R_duration
    R_inte_delay= desired_R_inte_delay
    show_args(Experi_info(q),title='Update parameters')

#%%

def Qubit_state_avg_timetrace(quantum_device:QuantumDevice,trace_recordlength:float,n_avg:int=1000,run:bool=True,q:str='q1'):
    
    sche_func = Trace_sche 
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    data = {}
    result = {}

    def state_dep_sched(ini_state):
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp=pi_amp,
            R_amp=R_amp,
            R_duration=R_duration,
            R_integration=R_integration,
            R_inte_delay=R_inte_delay,
            trace_recordlength=trace_recordlength,
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
            t_ds= gettable.get()
            
            data[ini_state] = t_ds
            
            #show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
    
        else:
            pulse_preview(quantum_device,sche_func,sched_kwargs)
            
            #show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
            
    data['trace_recordlength']= trace_recordlength   
    state_dep_sched('g')
    state_dep_sched('e')
    result[q]=  data  
    Current_Readout_IF= 5.95*1e9-R_F[q]
    print('Current_Readout_IF=',Current_Readout_IF/1e6,'MHz')
    return result


execute = True

if execute: 
    Avgtimetrace_result= Qubit_state_avg_timetrace(quantum_device,
                                                   trace_recordlength=8*1e-6,
                                                   n_avg=100,
                                                   run=True,
                                                   q='q1')

#%%
# Trace acquisition default gives the IF-band data  
# Program will perform digital down-coversion and filtering to baseband 
plot_data =True
if plot_data:              
   for q in Avgtimetrace_result:   
       Baseband_data= Qubit_state_Avgtimetrace_plot(Avgtimetrace_result[q],fc=100*1e6,Digital_downconvert=True, IF=5.95*1e9-R_F[q])
R_inte_delay= 192*1e-9

#%%
T1= 7.5e-6 #np.mean(T1_hist)


def Qubit_state_single_shot(quantum_device:QuantumDevice,shots:int=1000,run:bool=True,q:str='q1'):
    
    sche_func = Qubit_SS_sche  
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    data = {}
    analysis_result = {}
    exp_kwargs= dict(shots=shots,
                     )
    
    def state_dep_sched(ini_state):
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp=pi_amp,
            R_amp=R_amp,
            R_duration=R_duration,
            R_integration=R_integration,
            R_inte_delay=R_inte_delay,
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
            
            data[ini_state] = ss_ds
            
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
    
        else:
            pulse_preview(quantum_device,sche_func,sched_kwargs)
            
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
            
    tau= R_integration[q]        
    state_dep_sched('g')
    state_dep_sched('e')
        
    analysis_result[q]= Qubit_state_single_shot_fit_analysis(data,T1=T1,tau=tau)
        
    return analysis_result


execute = True

if execute: 
    SS_result= Qubit_state_single_shot(quantum_device,
                shots=20000,
                run=True,
                q='q1')



#%%
plot_data =True
if plot_data:               #y_scale='log' , 'linear'
   for q in SS_result:     #Plot_type='g_only', 'e_only','both'  
       Qubit_state_single_shot_plot(SS_result[q],Plot_type='g_only',y_scale='log')
       show_args(SS_result[q]['error_pack'],title='Error budget')

#%%

def Readout_F_calibration(quantum_device:QuantumDevice,ro_bare_guess:dict,ro_span_Hz:int=5e6,n_avg:int=1000,points:int=200,run:bool=True,q:str='q1'):
    

    sche_func = Qubit_state_heterodyne_spec_sched_nco
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)    
    analysis_result = {}
    data = {} 
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     )
    def state_dep_sched(ini_state):
        spec_sched_kwargs = dict(   
            frequencies=freq,
            q= q,
            ini_state= ini_state,
            pi_amp=pi_amp,
            R_amp=R_amp,
            R_duration=R_duration,
            R_integration=R_integration,
            R_inte_delay=R_inte_delay
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
            rs_ds = meas_ctrl.run("Readout_F_cal")
            rs_ds
            mag,pha= dataset_to_array(dataset=rs_ds,dims=1)
            data[ini_state] = tuple([mag,pha])
            show_args(exp_kwargs, title="Readout_F_cal_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
            
        else:
            sweep_para= np.array(ro_f_samples[:n_s])
            spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
            pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
            show_args(exp_kwargs, title="Readout_F_cal_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
            
    state_dep_sched('g')   
    state_dep_sched('e')  
    analysis_result[q]= Readout_F_opt_analysis(data,ro_f_samples)
    f_g= analysis_result[q]['f_g']
    f_e= analysis_result[q]['f_e']
    f_eff_bare= (f_g+f_e)/2
    print('')
    print('f_g=',f_g,' GHz')  
    print('f_e=',f_e,' GHz')
    print('f_eff_bare=',f_eff_bare,' GHz')
    

    
    return analysis_result, f_eff_bare



execute = True

if execute:
    Readout_F_cal_results,f_eff_bare = Readout_F_calibration(quantum_device,
                               ro_bare_guess=ro_bare,
                               ro_span_Hz=3e6,
                               n_avg=500,
                               points=201,
                               run=True,
                               q='q1')
    
#%%    
plot_data =True
update= True   

if plot_data:
    for q in Readout_F_cal_results:
        Readout_F_opt_Plot(quantum_device,Readout_F_cal_results[q])
        if update:
            qubit = quantum_device.get_element(q)
            qubit.clock_freqs.readout(f_eff_bare)
            R_F[q]= f_eff_bare
            show_args(Experi_info(q),title='Update parameters')       

#%%
T1= 7.5e-6 #np.mean(T1_hist)


def Readout_amp_opt(quantum_device:QuantumDevice,R_amp_low:float,R_amp_high:float,points:int=2,shots:int=1000,run:bool=True,q:str='q1'):
    
    sche_func = Qubit_amp_SS_sche
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO) 
    analysis_result = {}
    analysis_result[q]=[]
    amp_samples = linspace(R_amp_low,R_amp_high,points)

    exp_kwargs= dict(sweep_amp=['start '+'%E' %amp_samples[0],'end '+'%E' %amp_samples[-1]],
                     )
    
    def state_dep_sched(ini_state,amp,data):
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp=pi_amp,
            R_amp=amp,
            R_duration=R_duration,
            R_integration=R_integration,
            R_inte_delay=R_inte_delay,
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
            R_amp_ds= gettable.get()
            
            data[ini_state]=R_amp_ds
            
            return data
    
        else:
            sweep_para= np.array(amp_samples[:n_s])
            spec_sched_kwargs['R_amp']= sweep_para.reshape(sweep_para.shape or (1,))

            return spec_sched_kwargs
    
    tau= R_integration[q] 
    if run:
        for i in range(points):
            data = {}
            state_dep_sched('g',amp_samples[i],data)
            state_dep_sched('e',amp_samples[i],data)
            analysis_result[q].append(Qubit_state_single_shot_fit_analysis(data,T1=T1,tau=tau))
    else:        
        spec_sched_kwargs= state_dep_sched('g')
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)    
        
    show_args(exp_kwargs, title="Single_shot_amp_opt_kwargs: Meas.qubit="+q)
    show_args(Experi_info(q))    
        
    return  analysis_result, amp_samples


execute = True

if execute: 
    R_amp_opt_results, para= Readout_amp_opt(quantum_device,
                                      R_amp_low=0.05,
                                      R_amp_high=0.15,
                                      points=6,
                                      shots=20000,
                                      run=True,
                                      q='q1')



#%%
plot_data =True
linecut=5
print('R_amp=',amp_samples[linecut])

if plot_data:               #y_scale='log' , 'linear'
   for q in R_amp_opt_results:     #Plot_type='g_only', 'e_only','both'  
       Qubit_state_single_shot_plot(R_amp_opt_results[q][linecut],Plot_type='both',y_scale='linear')
       show_args(R_amp_opt_results[q][linecut]['error_pack'],title='Error budget')



#%%
T1= 7.5e-6 #np.mean(T1_hist)


def Readout_integration_opt(quantum_device:QuantumDevice, trace_recordlength=8*1e-6,shots:int=1000,run:bool=True,q:str='q1'):
    
    sche_func = Trace_sche   
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    data = {}
    data_trace = []
    analysis_result = {}
    exp_kwargs= dict(shots=shots,
                     )
    
    def state_dep_sched(ini_state):
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp=pi_amp,
            R_amp=R_amp,
            R_duration=R_duration,
            R_integration=R_integration,
            R_inte_delay=R_inte_delay,
            trace_recordlength=trace_recordlength,
        )
        
        if run:
            gettable = ScheduleGettable(
                quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
            )
            for i in range(shots):
                quantum_device.cfg_sched_repetitions(1)
                ss_ds= gettable.get()
                data_trace.append(ss_ds)
                
            data[ini_state] = data_trace
            
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
    
        else:
            pulse_preview(quantum_device,sche_func,sched_kwargs)
            
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q))
            
    tau= R_integration[q]        
    state_dep_sched('g')
    #state_dep_sched('e')
        
    #analysis_result[q]= Qubit_state_single_shot_fit_analysis(data,T1=T1,tau=tau)
        
    return data #analysis_result


execute = True

if execute: 
    SSTrace_results= Readout_integration_opt(quantum_device,
                                       trace_recordlength=8*1e-6,
                                       shots=100,
                                       run=True,
                                       q='q1')



#%%
plot_data =True
if plot_data:               #y_scale='log' , 'linear'
   for q in SSTrace_results:     #Plot_type='g_only', 'e_only','both'  
       Qubit_state_single_shot_plot(SSTrace_results[q],Plot_type='both',y_scale='linear')
       show_args(SSTrace_results[q]['error_pack'],title='Error budget')
       
       
       
       