import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN, array, arange, minimum, maximum, cos, sin, pi, linspace
from xarray import Dataset
import matplotlib.pyplot as plt
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_Qubit_SS_sche, Single_shot_ref_fit_analysis, pulse_preview

def Single_shot_ref_spec(QD_agent:QDmanager,ro_elements:dict,shots:int=1000,run:bool=True):
    print("Single shot start")
    sche_func = multi_Qubit_SS_sche   

    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        
        eyeson_print(f"Inte_time= {round(qubit_info.measure.integration_time()*1e6,1)} µs")
        eyeson_print(f"Reset_time= {round(qubit_info.reset.duration()*1e6,1)} µs")
       
        qubit_info.measure.pulse_amp(ro_elements[q]*float(qubit_info.measure.pulse_amp()))
        if qubit_info.rxy.amp180() is NaN:
            qubit_info.rxy.amp180(0)
        if qubit_info.rxy.duration() is NaN:
            qubit_info.rxy.duration(0)

    sched_kwargs = dict(   
        ini_state='g',
        waveformer=QD_agent.Waveformer,
        pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements,'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
    )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(shots)
        iq_tuples = gettable.get() # tuple?
        dict_ = {}
        for q_idx, q in enumerate(ro_elements):
            IQ_array = array([iq_tuples[2*q_idx],iq_tuples[2*q_idx+1]])
            dict_[q] = (["mixer","shots"],IQ_array)
        ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"shots":arange(shots)})
        
    else:
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        ds = ""
        
        
    return ds

def Single_shot_fit_plot(results:dict,title_qubit:str=None,save_pic_folder:str=None):
    c_I,c_Q,sig=1000*results['fit_pack'][0],1000*results['fit_pack'][1],1000*results['fit_pack'][2]
    I,Q= results['data'][0],results['data'][1]
    radius = sig*2
    theta = linspace(0, 2 * pi, 720)
    x = c_I + radius * cos(theta)
    y = c_Q + radius * sin(theta)

    fig, ax = plt.subplots(nrows =1,figsize =(8,8),dpi =200)
    ax.scatter(1000*I, 1000*Q, color="blue", alpha=0.5, s=5)       
    ax.scatter(c_I,c_Q,c='k',s=15)
    ax.plot(x, y, linewidth=0.8, linestyle='--', c='red')
    ax.set_xlabel(r"$I\ $(mV)",size ='15')
    ax.set_ylabel(r"$Q\ $(mV)",size ='15')
    ax.set_title(f'{title_qubit if title_qubit is not None else ""} Single shot raw data')
    ax.set_xlim(1000*minimum(min(I),min(Q)),1000*maximum(max(I),max(Q)))
    ax.set_ylim(1000*minimum(min(I),min(Q)),1000*maximum(max(I),max(Q)))
    fig.tight_layout()
    plt.grid()
    if save_pic_folder is not None:
        plt.savefig(os.path.join(save_pic_folder,f"RefIQ_{title_qubit}.png"))
        plt.close()
    else:
        plt.show()

def IQ_ref_ana(ds:Dataset, q:str, save_pic_folder:str=None):
    analysis_result = Single_shot_ref_fit_analysis(ds[q])
    ref_IQ = array([analysis_result['fit_pack'][0],analysis_result['fit_pack'][1]])
    Single_shot_fit_plot(analysis_result,q,save_pic_folder)
    return ref_IQ
    
