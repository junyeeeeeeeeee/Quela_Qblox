
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, NaN, arange
from xarray import Dataset
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.QuFluxFit import calc_Gcoef_inFbFqFd, calc_g
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_Two_tone_sche, pulse_preview


def Two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,xyf_guess:dict,xyl_guess:list,n_avg:int=500,run:bool=True,drive_read_overlap:bool=False)->Dataset:
    sche_func = multi_Two_tone_sche   
    
    original_qubit_info = {}
    for q in xyf_guess:
        original_qubit_info[q] = {}
        qubit_info = QD_agent.quantum_device.get_element(q)
        original_qubit_info[q]["ROW"] = qubit_info.measure.pulse_duration()
        original_qubit_info[q]["ITW"] = qubit_info.measure.integration_time()
        original_qubit_info[q]["XYF"] = qubit_info.clock_freqs.f01()
        freq_datapoint_idx = arange(0,xyf_guess[q].shape[0])
            
        if not drive_read_overlap:
            drive_pulse_ref_pt = 'start'
            drive_pulse_length = 100e-6
        else:
            drive_pulse_ref_pt = 'end'
            drive_pulse_length = 100e-6
            qubit_info.measure.pulse_duration(drive_pulse_length)
            qubit_info.measure.integration_time(drive_pulse_length)

        
        qubit_info.reset.duration(280e-9+drive_pulse_length+0.5e-6) 
        qubit_info.clock_freqs.f01(NaN)
        eyeson_print(f"Inte_time= {round(qubit_info.measure.integration_time()*1e6,1)} µs")
        eyeson_print(f"Reset_time= {round(qubit_info.reset.duration()*1e6,1)} µs")

    freq = ManualParameter(name="XYfreq", unit="Hz", label="Frequency")
    freq.batched = True
    xy_amp = ManualParameter(name="XYamp", unit="V", label="Amplitude")
    xy_amp.batched = False

    spec_sched_kwargs = dict(   
        frequencies=xyf_guess,
        spec_amp=xy_amp,
        spec_Du=drive_pulse_length,
        R_amp=compose_para_for_multiplexing(QD_agent,xyf_guess,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,xyf_guess,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,xyf_guess,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,xyf_guess,'r2'),
        ref_pt=drive_pulse_ref_pt
    )

    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(xyf_guess.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables((freq,xy_amp))
        meas_ctrl.setpoints_grid((freq_datapoint_idx,xyl_guess))
        
        ds = meas_ctrl.run("Two-tone")
        
        dict_ = {}
        for idx, q in enumerate(xyf_guess):
            freq_values = 2*len(xyl_guess)*list(xyf_guess[q])
            i_data = array(ds[f'y{2*idx}']).reshape(len(xyl_guess),array(xyf_guess[q]).shape[0])
            q_data = array(ds[f'y{2*idx+1}']).reshape(len(xyl_guess),array(xyf_guess[q]).shape[0])
            dict_[q] = (["mixer","xy_amp","freq"],array([i_data,q_data]))
            dict_[f"{q}_freq"] = (["mixer","xy_amp","freq"],array(freq_values).reshape(2,len(xyl_guess),array(xyf_guess[q]).shape[0]))

        
        ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"xy_amp":array(xyl_guess),"freq":freq_datapoint_idx})
        ds.attrs["execution_time"] = Data_manager().get_time_now()
        
        
    else:
        n_s = 2
        preview_para = {}
        for q in xyf_guess:
            preview_para[q] = xyf_guess[q][:n_s]
        sweep_para2 = array(xyl_guess[0])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['spec_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[0]
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        ds = ""
    
    # restore
    for q in original_qubit_info:
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.measure.pulse_duration(original_qubit_info[q]["ROW"])
        qubit_info.measure.integration_time(original_qubit_info[q]["ITW"])
        qubit_info.clock_freqs.f01(original_qubit_info[q]["XYF"])

    return ds

def update_2toneResults_for(QD_agent:QDmanager,qb:str,QS_results:dict,XYL:float):
    qubit = QD_agent.quantum_device.get_element(qb)
    Revised_f01 = QS_results[qb].attrs['f01_fit']
    fb = float(QD_agent.Notewriter.get_bareFreqFor(target_q=qb))*1e-6
    fd = QD_agent.quantum_device.get_element(qb).clock_freqs.readout()*1e-6
    A = calc_Gcoef_inFbFqFd(fb,Revised_f01*1e-6,fd)
    sweet_g = calc_g(fb,Revised_f01*1e-6,A)
    qubit.clock_freqs.f01(Revised_f01)
    QD_agent.Notewriter.save_2tone_piamp_for(qb,XYL)
    QD_agent.Notewriter.save_CoefInG_for(target_q=qb,A=A)
    QD_agent.Notewriter.save_sweetG_for(target_q=qb,g_Hz=sweet_g*1e6)

def tune_away_setup(QD_agent:QDmanager,Fctrl:dict,bias_setup:dict,ro_qs:list,zero_mode:bool=False):
    for q in ro_qs:
        if zero_mode:
            Fctrl[q](0.0)
        else:
            offset = QD_agent.Fluxmanager.get_proper_zbiasFor(q)
            if q not in list(bias_setup.keys()):
                bias_setup[q] = 0
            want_bias = offset+bias_setup[q]
            if bias_setup[q] != 0:
                rof = QD_agent.Fluxmanager.sin_for_cav(q,array([want_bias]))[0]
                QD_agent.quantum_device.get_element(q).clock_freqs.readout(rof)
                QD_agent.Fluxmanager.save_tuneawayBias_for('manual',q,want_bias)
                warning_print(f"meas under flux = {round(offset,3)}+{round(bias_setup[q],3)} V")
            Fctrl[q](want_bias) 
