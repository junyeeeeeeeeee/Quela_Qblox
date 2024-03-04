from numpy import linspace, array
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QuantumDevice, get_time_now
from quantify_core.measurement.control import MeasurementControl
from Pulse_schedule_library import Rabi_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, Rabi_fit_analysis

def Rabi(quantum_device:QuantumDevice,meas_ctrl:MeasurementControl,XY_amp:float, XY_duration:float=20e-9, IF:int=150e6,n_avg:int=300,points:int=200,run:bool=True,XY_theta:str='X_theta',Rabi_type:str='PowerRabi',q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0]):
    analysis_result = {}
    sche_func= Rabi_sche
    qubit_info = quantum_device.get_element(q)
    LO= qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    if Rabi_type=='TimeRabi':
       Para_XY_amp= XY_amp
       Sweep_para=Para_XY_Du = ManualParameter(name="XY_Duration", unit="s", label="Time")
       str_Rabi= 'XY_duration'
       Sweep_para.batched = True
       samples = linspace(0, XY_duration,points)
       exp_kwargs= dict(sweep_duration=[Rabi_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                        Amp='%E' %XY_amp,
                        )
    elif Rabi_type=='PowerRabi':
        Sweep_para= Para_XY_amp= ManualParameter(name="XY_amp", unit="V", label="Voltage")
        str_Rabi= 'XY_amp'
        Para_XY_amp.batched = True
        Para_XY_Du = XY_duration
        samples = linspace(0,XY_amp,points) 
        exp_kwargs= dict(sweep_amp=[Rabi_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                         Duration='%E' %XY_duration,
                         )
    else: raise KeyError ('Typing error: Rabi_type')
    
    sched_kwargs = dict(
        q=q,
        XY_amp=Para_XY_amp,
        XY_duration=Para_XY_Du,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        XY_theta=XY_theta,
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
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        data_fit= Rabi_fit_analysis(data=data,samples=samples,Rabi_type=Rabi_type)
        analysis_result[q]= data_fit
        
        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sched_kwargs[str_Rabi]= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,sched_kwargs)
       

        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))

    return analysis_result
    