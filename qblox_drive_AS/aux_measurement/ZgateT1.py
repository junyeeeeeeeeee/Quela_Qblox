import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from numpy import array, arange, ndarray
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_Zgate_T1_sche, pulse_preview
from qblox_drive_AS.SOP.FluxQubit import z_pulse_amp_OVER_const_z

def Zgate_T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,time_samples:dict,z_samples:ndarray,n_avg:int=500,run:bool=True,no_pi_pulse:bool=False):
    sche_func= multi_Zgate_T1_sche
    origin_pi_amp = {}
    for q in time_samples:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        origin_pi_amp[q] = qubit_info.rxy.amp180()
        if no_pi_pulse:
            qubit_info.rxy.amp180(0)
        time_data_idx = arange(time_samples[q].shape[0])
    

    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples:ndarray = z_samples
    
    sched_kwargs = dict(
        freeduration=time_samples,
        Z_amp=Z_bias,
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
        meas_ctrl.settables([Para_free_Du,Z_bias])
        meas_ctrl.setpoints_grid((time_data_idx,Z_samples*z_pulse_amp_OVER_const_z))


        ds = meas_ctrl.run('Zgate_T1')
        dict_ = {}

        for q_idx, q in enumerate(time_samples):   
            i_data = array(ds[f'y{2*q_idx}']).reshape(Z_samples.shape[0],time_samples[q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(Z_samples.shape[0],time_samples[q].shape[0])
            dict_[q] = (["mixer","z_voltage","idx"],array([i_data,q_data]))

            time_values = list(time_samples[q])*2*Z_samples.shape[0]
            dict_[f"{q}_time"] = (["mixer","z_voltage","time"],array(time_values).reshape(2,Z_samples.shape[0],array(time_samples[q]).shape[0]))
            
        zT1_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"z_voltage":Z_samples,"time":arange(time_data_idx.shape[0])})

        for q in time_samples:
            zT1_ds.attrs[f"{q}_ref_bias"] = round(QD_agent.Fluxmanager.get_proper_zbiasFor(q),3)

        zT1_ds.attrs["execution_time"] = Data_manager().get_time_now()
        zT1_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        zT1_ds.attrs["prepare_excited"] = 0 if no_pi_pulse else 1
        zT1_ds.attrs["method"] = "Average"
        zT1_ds.attrs["system"] = "qblox"

    else:
        preview_para = {}
        for q in time_samples:
            preview_para[q] = time_samples[q][:2]
        sweep_para2 = array([z_samples[0],z_samples[-1]])
        sched_kwargs['freeduration']= preview_para
        sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        zT1_ds = ''
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

    for q in origin_pi_amp:
        QD_agent.quantum_device.get_element(q).rxy.amp180(origin_pi_amp[q])
    
    return zT1_ds

