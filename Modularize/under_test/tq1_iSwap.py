
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from Modularize.support import uw
from numpy import array, linspace, arange
from Modularize.support import cds
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import Data_manager, QDmanager, compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support import init_meas, init_system_atte, shut_down, QRM_nco_init
from Modularize.support.Pulse_schedule_library import iSwap_sche, set_LO_frequency, pulse_preview, Multi_dataset_to_array, IQ_data_dis


def iSwap(QD_agent:QDmanager, meas_ctrl:MeasurementControl, who_with_pi:str='q0', target_Z1:str='q0', 
          target_Z2:str='q0', target_Zc:str='qc0', Z1_on_off:bool=True, Z2_on_off:bool=True,
          Zc_on_off:bool=True, sweep_bias:str='Zc', freeduration:float=5*1e-6, Z1_amp_min:float=0,
          Z1_amp_max:float=0, Z2_amp_min:float=0, Z2_amp_max:float=0, Zc_amp_min:float=0,
          Zc_amp_max:float=0, n_avg:int=300, freeDu_points:int=201, Z_points:int=201, run:bool=True, meas_q:list=['q0']):
    
    IF = 150e6
    contributors = [target_Z1,target_Z2,target_Zc]

    data={}
    data_save=[]
    LO= QD_agent.quantum_device.get_element(who_with_pi).clock_freqs.f01() + IF
    hw_c= set_LO_frequency(QD_agent.quantum_device,q=who_with_pi,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    freeDu_samples = linspace(0,freeduration,freeDu_points)
    
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
        pi_amp={str(who_with_pi):QD_agent.quantum_device.get_element(who_with_pi).rxy.amp180()},
        pi_dura={str(who_with_pi):QD_agent.quantum_device.get_element(who_with_pi).rxy.duration()},
        target_pi=who_with_pi,
        freeduration=Para_free_Du,
        target_Z1=target_Z1,
        target_Z2=target_Z2,
        target_Zc=target_Zc,
        Z1_amp=Z1_bias,
        Z2_amp=Z2_bias,
        Zc_amp=Zc_bias,
        R_amp=compose_para_for_multiplexing(QD_agent,meas_q,1),
        R_duration=compose_para_for_multiplexing(QD_agent,meas_q,3),
        R_integration=compose_para_for_multiplexing(QD_agent,meas_q,4),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,meas_q,2),
        Z_on_off=Z_on_off,
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %freeDu_samples[0],'end '+'%E' %freeDu_samples[-1]],
                     freeDu_points=freeDu_points,
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                     Z_points=Z_points,
                     )
    if run:

    
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(meas_q),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([Para_free_Du,Para_Z])
        meas_ctrl.setpoints_grid((freeDu_samples,Z_samples))
        iswap_ds = meas_ctrl.run("iSwap")
        Data_manager().save_2Qraw_data(QD_agent,iswap_ds,contributors,exp_type='iswap')
        
        I,Q=  Multi_dataset_to_array(dataset=iswap_ds,dims=2,Q=meas_q)

        
        
        for q in meas_q:
            data[q]= IQ_data_dis(I[q],Q[q],ref_I=I_ref[q],ref_Q=Q_ref[q]).transpose()
            
            show_args(exp_kwargs, title="iSwap_kwargs: Meas.qubit="+q)
            show_args(Experi_info(q,globals()))
            
        data_save= [freeDu_samples,Z_samples,data]
        

    else:
        n_s = 2
        sweep_para1= array(freeDu_samples[:n_s])
        sweep_para2= array(Z_samples[:2])
        sched_kwargs['freeduration']= sweep_para1.reshape(sweep_para1.shape or (1,))
        sched_kwargs[str_bias]= sweep_para2.reshape(sweep_para2.shape or (1,))[0]
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        
        
        for i in meas_q:
            show_args(exp_kwargs, title="iSwap_kwargs: Meas.qubit="+i)
            show_args(Experi_info(i,globals()))
        
    return data_save