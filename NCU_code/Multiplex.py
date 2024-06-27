from Pulse_schedule_library import *

n_s=2 # number of schedule in pulse_preview

IF= 200*1e6 #IF frequency used in time-domain


#%%

def iSwap(quantum_device:QuantumDevice,pi_amp:dict,target_pi:str='q0',target_Z1:str='q0',target_Z2:str='q0',target_Zc:str='qc0',Z1_on_off:bool=True,Z2_on_off:bool=True,Zc_on_off:bool=True,sweep_bias:str='Z',freeduration:float=5*1e-6, Z1_amp_min:float=0,Z1_amp_max:float=0, Z2_amp_min:float=0,Z2_amp_max:float=0, Zc_amp_min:float=0,Zc_amp_max:float=0,n_avg:int=300,freeDu_points:int=201,Z_points:int=201,run:bool=True,meas_q:list=['q0']):

    data={}
    data_save=[]
    LO= f01[target_pi]+IF
    hw_c= set_LO_frequency(quantum_device,q=target_pi,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    freeDu_samples = np.linspace(0,freeduration,freeDu_points)
    
    if sweep_bias=='Zc':
        str_bias= 'Zc_amp'
        Z_samples = linspace(Zc_amp_min,Zc_amp_max,Z_points)
        Z1_bias= Z1_amp_min
        Z2_bias= Z2_amp_min
        Para_Z=Zc_bias = ManualParameter(name="Zc", unit="V", label="Zc bias")
        Para_Z.batched = False

    elif sweep_bias=='Z1':
        str_bias= 'Z1_amp'
        Z_samples = linspace(Z1_amp_min,Z1_amp_max,Z_points)
        Zc_bias= Zc_amp_min
        Z2_bias= Z2_amp_min
        Para_Z=Z1_bias = ManualParameter(name="Z1", unit="V", label="Z1 bias")
        Para_Z.batched = False
    elif sweep_bias=='Z2':
        str_bias= 'Z2_amp'
        Z_samples = linspace(Z2_amp_min,Z2_amp_max,Z_points)
        Zc_bias= Zc_amp_min
        Z1_bias= Z1_amp_min
        Para_Z=Z2_bias = ManualParameter(name="Z2", unit="V", label="Z2 bias")
        Para_Z.batched = False
    else: 
        raise KeyError ('Zc and Z should not be both sweeping or non-sweeping')
        
    Z_on_off= dict(Z1_on_off=Z1_on_off,Z2_on_off=Z2_on_off,Zc_on_off=Zc_on_off)
    
    sche_func= iSwap_sche
    sched_kwargs = dict(
        q=meas_q,
        pi_amp=pi_amp,
        target_pi=target_pi,
        freeduration=Para_free_Du,
        target_Z1=target_Z1,
        target_Z2=target_Z2,
        target_Zc=target_Zc,
        Z1_amp=Z1_bias,
        Z2_amp=Z2_bias,
        Zc_amp=Zc_bias,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        Z_on_off=Z_on_off,
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %freeDu_samples[0],'end '+'%E' %freeDu_samples[-1]],
                     freeDu_points=freeDu_points,
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                     Z_points=Z_points,
                     )
    if run:

    
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(meas_q),
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([Para_free_Du,Para_Z])
        meas_ctrl.setpoints_grid((freeDu_samples,Z_samples))
        iswap_ds = meas_ctrl.run("iSwap")
        iswap_ds
        
        I,Q=  Multi_dataset_to_array(dataset=iswap_ds,dims=2,Q=meas_q)

        
        
        for q in meas_q:
            data[q]= IQ_data_dis(I[q],Q[q],ref_I=I_ref[q],ref_Q=Q_ref[q]).transpose()
            
            show_args(exp_kwargs, title="iSwap_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q,globals()))
            
        data_save= [freeDu_samples,Z_samples,data]
        

    else:
        sweep_para1= np.array(freeDu_samples[:n_s])
        sweep_para2= np.array(Z_samples[:2])
        sched_kwargs['freeduration']= sweep_para1.reshape(sweep_para1.shape or (1,))
        sched_kwargs[str_bias]= sweep_para2.reshape(sweep_para2.shape or (1,))[0]
        pulse_preview(quantum_device,sche_func,sched_kwargs)
        
        
        for i in meas_q:
            show_args(exp_kwargs, title="iSwap_kwargs: Meas.qubit="+i)
            show_args(Experi_info(i,globals()))
        
    return data_save


execute = True

if executor:
    iSwap_results =iSwap(quantum_device,
                        pi_amp=pi_amp, 
                        target_pi='q0',
                        target_Z1='q0',
                        target_Z2='q1',
                        target_Zc='qc0',
                        Z1_on_off=False,
                        Z2_on_off=True,
                        Zc_on_off=False,
                        sweep_bias='Z1',
                        freeduration=300e-6,
                        Z1_amp_min=-0.2,
                        Z1_amp_max=0.2,
                        Z2_amp_min=-0.5,
                        Z2_amp_max=0.5,
                        Zc_amp_min=-0.1,
                        Zc_amp_max=0.1,
                        n_avg=500,
                        freeDu_points=51,
                        Z_points=7,
                        run=True,
                        meas_q=['q0','q1','q2'],
                        ) 
    
#%%
linecut=1

plot_data =True

if plot_data:
    for q in iSwap_results[2]:
        plot_2D(iSwap_results[0]*1e6,iSwap_results[1],iSwap_results[2][q],
                label=[r"$t_{f}\ (\mu$s)"," Z bias (V)"],
                title="iSwap_"+q,
                readout_qubit_info=False,
                P_rescale=False,
                Dis=None,
                color_bound=False,
                bound_value=[14,16],
                plot_linecut=True,
                linecut=linecut,)
        show_args(Experi_info(q,globals()),title='Update parameters') 

    print("\033[34m Current linecut value= %.3f \033[0m" %iSwap_results[1][linecut])    
    
