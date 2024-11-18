from Modularize.support.Pulse_schedule_library import *

n_s=2 # number of schedule in pulse_preview

IF= 200*1e6 #IF frequency used in time-domain


def Rabi(quantum_device:QuantumDevice,ro_elements:list,XY_amp:float, XY_duration:float,n_avg:int=300,points:int=200,run:bool=True,XY_theta:str='X_theta',Rabi_type:str='TimeRabi',q:str='q1'):
    analysis_result = {}
    sche_func= Rabi_sche
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    if Rabi_type=='TimeRabi':
       Para_XY_amp= XY_amp
       Sweep_para=Para_XY_Du = ManualParameter(name="XY_Duration", unit="s", label="Time")
       str_Rabi= 'XY_duration'
       Sweep_para.batched = True
       samples = np.linspace(0, XY_duration,points)
       exp_kwargs= dict(sweep_duration=[Rabi_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                        Amp='%E' %XY_amp,
                        )
    elif Rabi_type=='PowerRabi':
        Sweep_para= Para_XY_amp= ManualParameter(name="XY_amp", unit="V", label="Voltage")
        str_Rabi= 'XY_amp'
        Para_XY_amp.batched = True
        Para_XY_Du = XY_duration
        samples = np.linspace(0,XY_amp,points) 
        exp_kwargs= dict(sweep_amp=[Rabi_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                         Duration='%E' %XY_duration,
                         )
    else: raise KeyError ('Typing error: Rabi_type')
    
    sched_kwargs = dict(
        q=q,
        XY_amp=Para_XY_amp,
        XY_duration=Para_XY_Du,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        XY_theta='X_theta',
        Rabi_type=Rabi_type,
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
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(samples)
    
       
        rabi_ds = meas_ctrl.run(Rabi_type)
        rabi_ds
        I,Q= dataset_to_array(dataset=rabi_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=I_ref,ref_Q=Q_ref)
        data_fit= Rabi_fit_analysis(data=data,samples=samples,Rabi_type=Rabi_type)
        analysis_result[q]= data_fit
        
        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    else:
       sweep_para= np.array(samples[:n_s])
       sched_kwargs[str_Rabi]= sweep_para.reshape(sweep_para.shape or (1,))
       pulse_preview(quantum_device,sche_func,sched_kwargs)
       

       show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
       show_args(Experi_info(q))
    return analysis_result
    
executor = True

if executor:
    Rabi_results = Rabi(quantum_device,
                        ro_elements=ro_element,
                        XY_amp=0.5, 
                        XY_duration=20*1e-9,
                        n_avg=1000,
                        points=201,
                        run=True,
                        XY_theta='X_theta',
                        Rabi_type='PowerRabi',
                        q='q1',) # Rabi_type='TimeRabi' or 'PowerRabi'

#%%   pi pulse duration is 20ns
plot_data =True
update= True 
  
if plot_data:       
    for q in Rabi_results:
        Fit_analysis_plot(Rabi_results[q],P_rescale=False,Dis=None)
        if update:
           pi_amp[q]= Rabi_results[q].attrs['pi_2']
        show_args(Experi_info(q),title='Update parameters')
        
#%%

def Zgate_Rabi(quantum_device:QuantumDevice,ro_elements:list,XY_amp:float, XY_duration:float, Z_amp_min:float,Z_amp_max:float,n_avg:int=300,Rabi_points:int=200,Z_points:int=201,run:bool=True,XY_theta:str='X_theta',Rabi_type:str='TimeRabi',q:str='q1'):
    analysis_result = {}
    sche_func= Zgate_Rabi_sche
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples = linspace(Z_amp_min,Z_amp_max,Z_points)
    
    if Rabi_type=='TimeRabi':
       Para_XY_amp= XY_amp
       Sweep_para=Para_XY_Du = ManualParameter(name="XY_Duration", unit="s", label="Time")
       str_Rabi= 'XY_duration'
       Sweep_para.batched = True
       samples = np.linspace(0, XY_duration,Rabi_points)
       exp_kwargs= dict(sweep_duration=[Rabi_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                        Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                        Amp='%E' %XY_amp,
                        )
    elif Rabi_type=='PowerRabi':
        Sweep_para= Para_XY_amp= ManualParameter(name="XY_amp", unit="V", label="Voltage")
        str_Rabi= 'XY_amp'
        Para_XY_amp.batched = True
        Para_XY_Du = XY_duration
        samples = np.linspace(0,XY_amp,Rabi_points) 
        exp_kwargs= dict(sweep_amp=[Rabi_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                         Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                         Duration='%E' %XY_duration,
                         )
    else: raise KeyError ('Typing error: Rabi_type')  
    
    sched_kwargs = dict(
        q=q,
        XY_amp=Para_XY_amp,
        XY_duration=Para_XY_Du,
        Z_amp=Z_bias,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        XY_theta='X_theta',
        Rabi_type=Rabi_type,
        )
    
    
    if run:
        gettable = ScheduleGettable(
        quantum_device,
        schedule_function=sche_func,
        schedule_kwargs=sched_kwargs,
        real_imag=False,
        batched=True,
        )
    
   
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(samples)
    
       
        rabi_ds = meas_ctrl.run('Zgate_'+Rabi_type)
        rabi_ds
        analysis_result[q]= Basic2DAnalysis(tuid=rabi_ds.attrs["tuid"], dataset=rabi_ds).run()
        show_args(exp_kwargs, title="Zgate_Rabi_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    else:
       sweep_para1= np.array(samples[:n_s])
       sweep_para2= np.array(Z_samples[:2])
       sched_kwargs[str_Rabi]= sweep_para1.reshape(sweep_para1.shape or (1,))
       sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
       pulse_preview(quantum_device,sche_func,sched_kwargs)
       show_args(exp_kwargs, title="Zgate_Rabi_kwargs: Meas.qubit="+q)
       show_args(Experi_info(q))
    return analysis_result
    
executor = True

if executor:
    Zgate_Rabi_results = Zgate_Rabi(quantum_device,
                        ro_elements=ro_element,
                        XY_amp=0.5, 
                        XY_duration=200*1e-9,
                        Z_amp_min=0,
                        Z_amp_max=0.3,
                        n_avg=1000,
                        Rabi_points=11,
                        Z_points=21,
                        run=False,
                        XY_theta='X_theta',
                        Rabi_type='PowerRabi',
                        q='q1',) # Rabi_type='TimeRabi' or 'PowerRabi'

   
plot_data =False

  
if plot_data:       
    for q in Zgate_Rabi_results:
        Amp_phase_plot(quantum_device, Zgate_Rabi_results[q],title='Zgate_Rabi')
        Zgate_Rabi_results[q].display_figs_mpl()

           


#%%
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
    
executor = True

if executor:
    Ramsey_results,T2_hist, average_actual_detune= Ramsey(quantum_device,
                            ro_elements=ro_element,
                            pi_amp=pi_amp, 
                            detuning=0.5e6,
                            freeduration=20e-6,
                            n_avg=500,
                            points=101,
                            run=True,
                            q='q1',
                            times=1,
                            ) 

qubit = quantum_device.get_element(q)
f_01_old= qubit.clock_freqs.f01()
print('f_01_old=', f_01_old,'GHz')  
  
#%%    
linecut=0
plot_data =True
update= True
  
if plot_data:       
    for q in Ramsey_results:
        Fit_analysis_plot(Ramsey_results[q][linecut],P_rescale=False,Dis=None)
        hist_plot(q,T2_hist ,title=r"$T_{2}\  (\mu$s)")
        print('average_actual_detune=', average_actual_detune*1e-6,'MHz')
        print('average_T2=', np.mean(np.array(T2_hist[q])),'us')
        if update:
           Revised_f01= f_01_old-average_actual_detune
           qubit.clock_freqs.f01(Revised_f01)
           f01[q]= Revised_f01
           print('f_01_revised=', Revised_f01,'GHz')  
        show_args(Experi_info(q),title='Update parameters')
        


        
#%%
def Zgate_Ramsey(quantum_device:QuantumDevice,ro_elements:list,pi_amp:dict,detuning:float,freeduration:float, Z_amp_min:float,Z_amp_max:float,n_avg:int=300,freeDu_points:int=201,Z_points:int=201,run:bool=True,q='q1',times=1):
    analysis_result= {}
    T2_us = {}
    analysis_result[q]= []
    T2_us[q] = []
    Real_detune= {}
    Real_detune[q]=[]
    New_fxy= f01[q]+detuning
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    freeDu_samples = np.linspace(0,freeduration,freeDu_points)
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples = linspace(Z_amp_min,Z_amp_max,Z_points)
    sche_func= Zgate_Ramsey_sche
    sched_kwargs = dict(
        q=q,
        pi_amp=pi_amp,
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        Z_amp=Z_bias,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %freeDu_samples[0],'end '+'%E' %freeDu_samples[-1]],
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
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
        meas_ctrl.settables([Para_free_Du,Z_bias])
        meas_ctrl.setpoints_grid((freeDu_samples,Z_samples))
        
 
        for i in range(times):
            ramsey_ds = meas_ctrl.run('Zgate_Ramsey')
            ramsey_ds
            for j in range(Z_points):
                I,Q= dataset_to_array(dataset=ramsey_ds,dims=2)
                data= IQ_data_dis(I,Q,ref_I=I_ref,ref_Q=Q_ref).transpose()
                data_fit= T2_fit_analysis(data=data[j],freeDu=freeDu_samples,T2_guess=8e-6)
                analysis_result[q].append(data_fit)
                T2_us[q].append(data_fit.attrs['T2_fit']*1e6)
                Real_detune[q].append(data_fit.attrs['f']-detuning)
                
        T2_us['plot_parameters']=[times,Z_samples]
        Real_detune['plot_parameters']=[times,Z_samples]
        show_args(exp_kwargs, title="Zgate_Ramsey_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    else:
        sweep_para1= np.array(freeDu_samples[:n_s])
        sweep_para2= np.array(Z_samples[:2])
        sched_kwargs['freeduration']= sweep_para1.reshape(sweep_para1.shape or (1,))
        sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Zgate_Ramsey_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    return analysis_result, T2_us, Real_detune
    
executor = True

if executor:
    Zgate_Ramsey_results,T2_bp,Real_detune_bp = Zgate_Ramsey(quantum_device,
                            ro_elements=ro_element,
                            pi_amp=pi_amp, 
                            detuning=0.1e6,
                            freeduration=16e-6,
                            Z_amp_min=-0.005,
                            Z_amp_max=0.01,
                            n_avg=100,
                            freeDu_points=501,
                            Z_points=11,
                            run=True,
                            q='q1',
                            times=30) 
    
#%%    
   
plot_data =True
linecut= 10

if plot_data:    
    for q in Zgate_Ramsey_results:
        Fit_analysis_plot(Zgate_Ramsey_results[q][linecut],P_rescale=True,Dis=0.004389268544379189)
        #T2_array= Z_bias_error_bar_plot(q,T2_bp,title=r"$T_{2}\  (\mu$s)") 
        #Ramsey_F_array= Ramsey_F_Z_bias_error_bar_plot(q,Real_detune_bp)
        
#%%


def T1(quantum_device:QuantumDevice,ro_elements:list,pi_amp:dict,freeduration:float,n_avg:int=300,points:int=200,run:bool=True,q='q1',times=1):
    analysis_result= {}
    T1_us = {}
    analysis_result[q]= []
    T1_us[q] = []
    sche_func=T1_sche
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    samples = np.linspace(0,freeduration,points)
    
    sched_kwargs = dict(
        q=q,
        pi_amp=pi_amp,
        freeduration=Para_free_Du,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
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
            T1_ds = meas_ctrl.run('T1')
            T1_ds
            I,Q= dataset_to_array(dataset=T1_ds,dims=1)
            data= IQ_data_dis(I,Q,ref_I=I_ref,ref_Q=Q_ref)
            data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=8e-6)
            analysis_result[q].append(data_fit)
            T1_us[q].append(data_fit.attrs['T1_fit']*1e6)
             
        show_args(exp_kwargs, title="T1_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    else:
        sweep_para= np.array(samples[:n_s])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="T1_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
  
    return analysis_result, T1_us
    



execute = True

if execute:   
    T1_results, T1_hist = T1(quantum_device,
                    ro_elements=ro_element,
                    pi_amp=pi_amp,
                    freeduration=80e-6,
                    n_avg=500,
                    points=401,
                    run=True,
                    q='q1',
                    times=1,
                    )  
#%%
plot_data =True
linecut=0
  
if plot_data:       
    for q in T1_results:
        Fit_analysis_plot(T1_results[q][linecut],P_rescale=True,Dis=0.0030489)
        hist_plot(q,T1_hist ,title=r"$T_{1}\  (\mu$s)")
        print('average_T1=', np.mean(np.array(T1_hist[q])),'us')
#%%

def Zgate_T1(quantum_device:QuantumDevice,ro_elements:list,pi_amp:dict,freeduration:float, Z_amp_min:float,Z_amp_max:float,n_avg:int=300,freeDu_points:int=201,Z_points:int=201,run:bool=True,q='q1',times=1):
    analysis_result= {}
    T1_us = {}
    analysis_result[q]= []
    T1_us[q] = []
    LO= f01[q]+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    freeDu_samples = np.linspace(0,freeduration,freeDu_points)
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples = linspace(Z_amp_min,Z_amp_max,Z_points)
    sche_func= Zgate_T1_sche
    sched_kwargs = dict(
        q=q,
        pi_amp=pi_amp,
        freeduration=Para_free_Du,
        Z_amp=Z_bias,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %freeDu_samples[0],'end '+'%E' %freeDu_samples[-1]],
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
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
        meas_ctrl.settables([Para_free_Du,Z_bias])
        meas_ctrl.setpoints_grid((freeDu_samples,Z_samples))

        
        for i in range(times):
            T1_ds = meas_ctrl.run('Zgate_T1')
            T1_ds
            for j in range(Z_points):
                I,Q= dataset_to_array(dataset=T1_ds,dims=2)
                data= IQ_data_dis(I,Q,ref_I=I_ref,ref_Q=Q_ref).transpose()
                data_fit= T1_fit_analysis(data=data[j],freeDu=freeDu_samples,T1_guess=8e-6)
                analysis_result[q].append(data_fit)
                T1_us[q].append(data_fit.attrs['T1_fit']*1e6)
        T1_us['plot_parameters']=[times,Z_samples]
        show_args(exp_kwargs, title="Zgate_T1_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
    else:
        sweep_para1= np.array(freeDu_samples[:n_s])
        sweep_para2= np.array(Z_samples[:2])
        sched_kwargs['freeduration']= sweep_para1.reshape(sweep_para1.shape or (1,))
        sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Zgate_T1_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q))
        
    return analysis_result, T1_us
    
executor = True

if executor:
    Zgate_T1_results,T1_bp = Zgate_T1(quantum_device,
                            ro_elements=ro_element,
                            pi_amp=pi_amp, 
                            freeduration=80e-6,
                            Z_amp_min=-0.1,
                            Z_amp_max=0.1,
                            n_avg=100,
                            freeDu_points=101,
                            Z_points=21,
                            run=True,
                            q='q1',
                            times=50,
                            ) 
#%%    
plot_data =True
linecut=0

if plot_data:    
    for q in Zgate_T1_results:
        Fit_analysis_plot(Zgate_T1_results[q][linecut],P_rescale=True,Dis=0.0032053422727767833)
        Z_bias_error_bar_plot(q,T1_bp,title=r"$T_{1}\  (\mu$s)")
 

           