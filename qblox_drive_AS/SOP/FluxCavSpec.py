import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, arange, NaN, ndarray
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Pulse_schedule_library import One_tone_multi_sche, pulse_preview
from qblox_drive_AS.SOP.FluxQubit import z_pulse_amp_OVER_const_z
from xarray import Dataset

def FluxCav_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,flux_ctrl:dict,ro_elements:dict,flux_samples:ndarray,n_avg:int=300,run:bool=True)->Dataset:
    sche_func = One_tone_multi_sche
    original_rof = {}
    flux_dura = 0
    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        flux_dura = qubit_info.reset.duration()+qubit_info.measure.integration_time()
        original_rof[q] = qubit_info.clock_freqs.readout()
        qubit_info.clock_freqs.readout(NaN)
        freq_datapoint_idx = arange(ro_elements[q].shape[0])

    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    bias = ManualParameter(name="bias", unit="V", label="Flux voltage")
    bias.batched = False
    
    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
        powerDep=False,
        bias = bias,
        bias_dura = flux_dura
    )

    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables((freq,bias))
        meas_ctrl.setpoints_grid((freq_datapoint_idx,flux_samples*z_pulse_amp_OVER_const_z)) # x0, x1
        
        ds = meas_ctrl.run("One-tone-Flux")
        dict_ = {}
        for idx, q in enumerate(ro_elements):
            freq_values = 2*flux_samples.shape[0]*list(ro_elements[q])
            i_data = array(ds[f'y{2*idx}']).reshape(flux_samples.shape[0],array(ro_elements[q]).shape[0])
            q_data = array(ds[f'y{2*idx+1}']).reshape(flux_samples.shape[0],array(ro_elements[q]).shape[0])
            dict_[q] = (["mixer","bias","freq"],array([i_data,q_data]))
            dict_[f"{q}_freq"] = (["mixer","bias","freq"],array(freq_values).reshape(2,flux_samples.shape[0],array(ro_elements[q]).shape[0]))

        
        rfs_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":flux_samples,"freq":arange(flux_samples.shape[0])})
        rfs_ds.attrs["execution_time"] = Data_manager().get_time_now()
        

        # reset flux bias
        for q in ro_elements:
            flux_ctrl[q](0.0)
        
    else:
        n_s = 2
        preview_para = {}
        for q in ro_elements:
            preview_para[q] = ro_elements[q][:n_s]
        sweep_para2= array(flux_samples[:2])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['bias']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)

        rfs_ds = ''

    for q in ro_elements:
        QD_agent.quantum_device.get_element(q).clock_freqs.readout(original_rof[q])
    
    return rfs_ds

def update_flux_info_in_results_for(QD_agent:QDmanager,qb:str,FD_results:dict):
    qubit = QD_agent.quantum_device.get_element(qb)
    qubit.clock_freqs.readout(FD_results["sweet_freq"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=FD_results["sweet_flux"])
    QD_agent.Fluxmanager.save_period_for(target_q=qb, period=1/FD_results['f'])
    QD_agent.Fluxmanager.save_tuneawayBias_for(target_q=qb,mode='auto')
    QD_agent.Fluxmanager.save_cavFittingParas_for(target_q=qb,
        f=FD_results['f'],
        amp=FD_results['A'],
        phi=FD_results['phi'],
        offset=FD_results['offset']
    )

