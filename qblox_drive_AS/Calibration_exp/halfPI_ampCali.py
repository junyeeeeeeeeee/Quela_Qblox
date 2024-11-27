"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, copy
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from numpy import array, arange, moveaxis
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_hPI_amp_cali_sche, pulse_preview


def half_pi_amp_cali(QD_agent:QDmanager,meas_ctrl:MeasurementControl, coef_samples:dict, pi_pair_num:list=[3,5],n_avg:int=300,run:bool=True):
    results = {}
    sche_func= multi_hPI_amp_cali_sche
    
    dataset_2_nc = ""
    for q in coef_samples:
        results[q], results[f"{q}_HalfPIcoef"] = [], []
        data_sample_idx = arange(coef_samples[q].shape[0])
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")

    
    Sweep_para = ManualParameter(name="XY_Amp_coef")
    Sweep_para.batched = True
    
    def pi_pair_dep_exe(pi_pair_num:int, a:QDmanager)->dict:
        
        dict_ = {}
        sched_kwargs = dict(
            pi_amp_coefs=coef_samples,
            pi_pair_num=pi_pair_num,
            waveformer=a.Waveformer,
            XY_amp=compose_para_for_multiplexing(QD_agent,coef_samples,'d1'),
            XY_duration=compose_para_for_multiplexing(QD_agent,coef_samples,'d3'),
            R_amp=compose_para_for_multiplexing(QD_agent,coef_samples,'r1'),
            R_duration=compose_para_for_multiplexing(QD_agent,coef_samples,'r3'),
            R_integration=compose_para_for_multiplexing(QD_agent,coef_samples,'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,coef_samples,'r2'),
            )
        
        
        if run:
            gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(coef_samples.keys())),
            )
            
    
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            meas_ctrl.gettables(gettable)
            meas_ctrl.settables(Sweep_para)
            meas_ctrl.setpoints(data_sample_idx)
        
        
            ds = meas_ctrl.run("Pi amp calibration")
            
            for q_idx, q in enumerate(coef_samples):
                I_data, Q_data = array(ds[f"y{2*q_idx}"]).tolist(), array(ds[f"y{2*q_idx+1}"]).tolist()
                dict_[q] = [I_data,Q_data] # shape in (mixer, pi_amp)
                dict_[f"{q}_HalfPIcoef"] = [list(coef_samples[q])]*2
    
        else:
            preview_para = {}
            for q in coef_samples:
                preview_para[q] = array([coef_samples[q][0],coef_samples[q][-1]])
            sched_kwargs['pi_amp_coefs']= preview_para
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

        return dict_
    
    for idx, pi_pair_number in enumerate(pi_pair_num):
        slightly_print(f"{pi_pair_number}*4 halfpi-pulses ...")
        
            
        dataDict = pi_pair_dep_exe(pi_pair_number,QD_agent) # {"q0":[],"q0_PIcoef":[], ...}
        for var in dataDict:
            results[var].append(dataDict[var])
            if idx == len(pi_pair_num) - 1: results[var] = (["mixer", "PiPairNum", "HalfPiCoef"],moveaxis(array(results[var]),0,1)) # shape (pi-pair_num, mixer, rof) -> (mixer, pi-pair_num, rof)
                
    if run:
        dataset_2_nc = Dataset(results,coords={"mixer":array(["I","Q"]),"PiPairNum":array(pi_pair_num),"HalfPiCoef":data_sample_idx})
        dataset_2_nc.attrs["execution_time"] = Data_manager().get_time_now()
        
   
    return dataset_2_nc

  