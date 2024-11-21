"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from qblox_drive_AS.support.UserFriend import *
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import linspace, array, arange, NaN, ndarray
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing, reset_offset
from qblox_drive_AS.support.Pulse_schedule_library import multi_PowerRabi_sche, multi_TimeRabi_sche, set_LO_frequency, pulse_preview
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from xarray import Dataset

#? The way to merge two dict a and b : c = {**a,**b}

def round_to_nearest_multiple_of_multipleNS(x, specific_multiple:int=None):
    if specific_multiple is None:
        specific_multiple = 4
    return specific_multiple * round(x / specific_multiple)


def PowerRabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,pi_amp:dict,pi_dura:dict,n_avg:int=300,run:bool=True,OSmode:bool=False):
    
    sche_func= multi_PowerRabi_sche

    for q in pi_amp:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
        pi_amp_idxes = arange(0,pi_amp[q].shape[0])
        
    
    Sweep_para = ManualParameter(name="PowerRabi", unit="Volt", label="Amp")
    Sweep_para.batched = True
    
    
    sched_kwargs = dict(
        pi_amp = pi_amp,
        pi_dura=pi_dura,
        R_amp=compose_para_for_multiplexing(QD_agent,pi_amp,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,pi_amp,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,pi_amp,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,pi_amp,'r2'),
        XY_theta='X_theta',
        OS_or_not=OSmode
        )
    
    
    if run:
        gettable = ScheduleGettable(
        QD_agent.quantum_device,
        schedule_function=sche_func,
        schedule_kwargs=sched_kwargs,
        real_imag=True,
        batched=True,
        num_channels=len(list(pi_amp.keys())),
        )
        
   
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(pi_amp_idxes)
    
       
        rabi_ds = meas_ctrl.run("RabiOsci")

        dict_ = {}
        for idx, q in enumerate(list(pi_amp.keys())):
            I = array(rabi_ds[f'y{2*idx}'])
            Q = array(rabi_ds[f'y{2*idx+1}'])
            dict_[q] = (['mixer','pi_amp'],array([I,Q]))
            dict_[f"{q}_piamp"] = (['mixer','pi_amp'],array(2*list(pi_amp[q])).reshape(2,pi_amp[q].shape[0]))
            
        ds = Dataset(dict_, coords={"mixer":array(["I","Q"]),"pi_amp":pi_amp_idxes})
        ds.attrs["execution_time"] = Data_manager().get_time_now()
        ds.attrs["rabi_type"] = "PowerRabi"
        ds.attrs["OS_mode"] = OSmode
        for q in pi_dura:
            ds.attrs[f"{q}_pidura"] = pi_dura[q]

     
    else:

        preview_para_pi_amp = {}
        preview_para_pi_dura = {}
        for q in pi_amp:
            preview_para_pi_amp[q] = array([pi_amp[q][0],pi_amp[q][-1]])
            preview_para_pi_dura[q] = pi_dura[q]
        
        sched_kwargs['pi_amp']= preview_para_pi_amp
        sched_kwargs['pi_dura']= preview_para_pi_dura
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        ds = ""

    return ds
    
def TimeRabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,pi_amp:dict,pi_dura:dict,n_avg:int=300,run:bool=True,OSmode:bool=False):
    
    sche_func= multi_TimeRabi_sche

    for q in pi_dura:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
        pi_dura_idxes = arange(0,pi_dura[q].shape[0])
        
    
    Sweep_para = ManualParameter(name="TimeRabi", unit="Sec", label="time")
    Sweep_para.batched = True
    
    
    sched_kwargs = dict(
        pi_amp = pi_amp,
        pi_dura=pi_dura,
        R_amp=compose_para_for_multiplexing(QD_agent,pi_dura,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,pi_dura,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,pi_dura,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,pi_dura,'r2'),
        XY_theta='X_theta',
        OS_or_not=OSmode
        )
    
    
    if run:
        gettable = ScheduleGettable(
        QD_agent.quantum_device,
        schedule_function=sche_func,
        schedule_kwargs=sched_kwargs,
        real_imag=True,
        batched=True,
        num_channels=len(list(pi_dura.keys())),
        )
        
   
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(pi_dura_idxes)
    
       
        rabi_ds = meas_ctrl.run("RabiOsci")

        dict_ = {}
        for idx, q in enumerate(list(pi_dura.keys())):
            I = array(rabi_ds[f'y{2*idx}'])
            Q = array(rabi_ds[f'y{2*idx+1}'])
            dict_[q] = (['mixer','pi_dura'],array([I,Q]))
            dict_[f"{q}_pidura"] = (['mixer','pi_dura'],array(2*list(pi_dura[q])).reshape(2,pi_dura[q].shape[0]))
        ds = Dataset(dict_, coords={"mixer":array(["I","Q"]),"pi_dura":pi_dura_idxes})
        ds.attrs["execution_time"] = Data_manager().get_time_now()
        ds.attrs["rabi_type"] = "TimeRabi"
        ds.attrs["OS_mode"] = OSmode
        for q in pi_amp:
            ds.attrs[f"{q}_piamp"] = pi_amp[q]

     
    else:

        preview_para_pi_amp = {}
        preview_para_pi_dura = {}
        for q in pi_amp:
            preview_para_pi_dura[q] = array([pi_dura[q][0],pi_dura[q][-1]])
            preview_para_pi_amp[q] = pi_amp[q]
        
        sched_kwargs['pi_amp']= preview_para_pi_amp
        sched_kwargs['pi_dura']= preview_para_pi_dura
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        ds = ""

    return ds

def conditional_update_qubitInfo(QD_agent:QDmanager,fit_results:Dataset,target_q:str):
    if fit_results.attrs['pi_2'] >= min(fit_results.coords['samples']) and fit_results.attrs['pi_2'] <= max(fit_results.coords['samples']) :
        qubit = QD_agent.quantum_device.get_element(target_q)
        match str(fit_results.attrs['Rabi_type']).lower():
            case 'powerrabi':
                qubit.rxy.amp180(fit_results.attrs['pi_2'])
                qubit.rxy.duration(fit_results.attrs["fix_variable"])
            case 'timerabi':
                qubit.rxy.amp180(fit_results.attrs["fix_variable"])
                qubit.rxy.duration(fit_results.attrs['pi_2'])
    else:
        warning_print(f"Results for {target_q} didn't satisfy the update condition !")

