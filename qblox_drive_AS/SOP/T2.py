import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qcodes.parameters import ManualParameter
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from quantify_scheduler.gettables import ScheduleGettable
from numpy import arange, array, arange
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing, QDmanager, Data_manager
from qblox_drive_AS.support.Pulse_schedule_library import multi_ramsey_sche, pulse_preview

def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,time_samples:dict, spin_num:dict={},repeat:int=1,n_avg:int=1000,run:bool=True, second_phase:str='x'):
    
    sche_func= multi_ramsey_sche
    for q in time_samples:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        time_data_idx = arange(time_samples[q].shape[0])
        if q not in list(spin_num.keys()):
            spin_num[q] = 0

    repeat_data_idx = arange(repeat)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
    Para_repeat.batched = False
    

    
   
    sched_kwargs = dict(
        freeduration=time_samples,
        pi_amp=compose_para_for_multiplexing(QD_agent,time_samples,'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,time_samples,'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,time_samples,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,time_samples,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,time_samples,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,time_samples,'r2'),
        echo_pi_num=spin_num,
        second_pulse_phase=second_phase
        )

    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(time_samples.keys())),
        )
        
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([Para_free_Du,Para_repeat])
        meas_ctrl.setpoints_grid((time_data_idx,repeat_data_idx))
        
        
        ds = meas_ctrl.run('Ramsey')
        
        dict_ = {}
        for q_idx, q in enumerate(time_samples):
            i_data = array(ds[f'y{2*q_idx}']).reshape(repeat,time_samples[q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(repeat,time_samples[q].shape[0])
            dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
            time_values = list(time_samples[q])*2*repeat
            dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,repeat,time_samples[q].shape[0]))
        
        ramsey_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":repeat_data_idx,"idx":time_data_idx})
        for var in [i for i in ramsey_ds.data_vars if "_" not in i]:
            ramsey_ds[var].attrs["spin_num"] = spin_num[var]
        ramsey_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        ramsey_ds.attrs["execution_time"] = Data_manager().get_time_now()

        
    else:
        preview_para = {}
        for q in time_samples:
            preview_para[q] = array([time_samples[q][0],time_samples[q][-1]])
        sched_kwargs['freeduration']= preview_para
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        ramsey_ds = ""
    return ramsey_ds










        
    
    
        
            
        


    