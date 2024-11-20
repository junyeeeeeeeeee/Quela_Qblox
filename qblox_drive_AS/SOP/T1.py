import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from numpy import array, arange
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_T1_sche, pulse_preview


def T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,time_samples:dict,repeat:int=1,n_avg:int=300,run:bool=True):
    sche_func= multi_T1_sche

    for q in time_samples:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        time_data_idx = arange(time_samples[q].shape[0])
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
    Para_repeat.batched = False
    repeat_data_idx = arange(repeat)

    sched_kwargs = dict(
        freeduration=time_samples,
        pi_amp=compose_para_for_multiplexing(QD_agent,time_samples,'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,time_samples,'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,time_samples,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,time_samples,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,time_samples,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,time_samples,'r2'),
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

        ds = meas_ctrl.run('T1')
        dict_ = {}
        for q_idx, q in enumerate(time_samples):
            i_data = array(ds[f'y{2*q_idx}']).reshape(repeat,time_samples[q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(repeat,time_samples[q].shape[0])
            dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
            time_values = list(time_samples[q])*2*repeat
            dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,repeat,time_samples[q].shape[0]))
        
        T1_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":repeat_data_idx,"idx":time_data_idx})
        T1_ds.attrs["execution_time"] = Data_manager().get_time_now()
        T1_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

    else:
        preview_para = {}
        for q in time_samples:
            preview_para[q] = array([time_samples[q][0],time_samples[q][-1]])
        sched_kwargs['freeduration']= preview_para
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        T1_ds = ''
    
    return T1_ds


