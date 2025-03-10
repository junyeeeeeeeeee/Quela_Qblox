import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from numpy import ndarray
from xarray import Dataset
from numpy import array, arange, sqrt
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_Z_gate_two_tone_sche, pulse_preview


z_pulse_amp_OVER_const_z = sqrt(2)/2.5


def Zgate_two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,XYFs:dict,Bias_element:list,Bias_samples:ndarray,n_avg:int=1000,run:bool=True):
    print("Zgate 2tone start")
    
    sche_func = multi_Z_gate_two_tone_sche

    original_xyfs = {}
    for q in XYFs:
        freq_datapoint_idx = arange(0,XYFs[q].shape[0])
        qubit_info = QD_agent.quantum_device.get_element(q)  
        original_xyfs[q] = qubit_info.clock_freqs.f01()
        qubit_info.clock_freqs.f01(NaN)
        spec_pulse_amp = QD_agent.Notewriter.get_2tone_piampFor(q)
    
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    Z_bias = ManualParameter(name="Flux", unit="V", label="Z bias")
    Z_bias.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=XYFs,
        bias_qs=Bias_element,
        Z_amp=Z_bias,
        spec_amp=spec_pulse_amp,
        spec_Du=10e-6,
        R_amp=compose_para_for_multiplexing(QD_agent,XYFs,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,XYFs,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,XYFs,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,XYFs,'r2'),
    )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(XYFs.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,Z_bias])
        meas_ctrl.setpoints_grid((freq_datapoint_idx,Bias_samples*z_pulse_amp_OVER_const_z))
        qs_ds = meas_ctrl.run("Zgate-two-tone")
        
        dict_ = {}
        for idx, q in enumerate(XYFs):
            freq_values = 2*Bias_samples.shape[0]*list(XYFs[q])
            i_data = array(qs_ds[f'y{2*idx}']).reshape(Bias_samples.shape[0],array(XYFs[q]).shape[0])
            q_data = array(qs_ds[f'y{2*idx+1}']).reshape(Bias_samples.shape[0],array(XYFs[q]).shape[0])
            dict_[q] = (["mixer","bias","freq"],array([i_data,q_data]))
    
            dict_[f"{q}_freq"] = (["mixer","bias","freq"],array(freq_values).reshape(2,Bias_samples.shape[0],array(XYFs[q]).shape[0]))

        rfs_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":Bias_samples,"freq":freq_datapoint_idx})
        rfs_ds.attrs["execution_time"] = Data_manager().get_time_now()
        rfs_ds.attrs["method"] = "Average"
        rfs_ds.attrs["system"] = "qblox"

    else:
        n_s = 2
        preview_para = {}
        for q in XYFs:
            preview_para[q] = XYFs[q][:n_s]
        sweep_para2 = array(Bias_samples[:2])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        rfs_ds = ""


        
    for q in XYFs:
        qubit_info = QD_agent.quantum_device.get_element(q)   
        qubit_info.clock_freqs.f01(original_xyfs[q])

    return rfs_ds


def update_by_fluxQubit(QD_agent:QDmanager,correct_results:dict,target_q:str):
    """
    correct_results dict in the form: {"xyf":float,"sweet_bias":float}
    """
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.clock_freqs.f01(correct_results["xyf"])
    QD_agent.Fluxmanager.check_offset_and_correctFor(target_q=target_q,new_offset=correct_results["sweet_bias"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=target_q,bias=correct_results["sweet_bias"])


