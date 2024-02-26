from numpy import array, linspace
from Modularize.support import build_folder_today
from Modularize.path_book import meas_raw_dir
from Modularize.Pulse_schedule_library import One_tone_sche, pulse_preview
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QuantumDevice, get_time_now, PDresults_alignPlot
from quantify_core.measurement.control import MeasurementControl
import os

def PowerDep_spec(quantum_device:QuantumDevice,meas_ctrl:MeasurementControl,ro_span_Hz:int=3e6,ro_p_min:float=0.1,ro_p_max:float=0.7,n_avg:int=100,f_points:int=30,p_points:int=30,run:bool=True,q:str='q1',Experi_info:dict={})->dict:

    sche_func = One_tone_sche
        
    analysis_result = {}
    qubit_info = quantum_device.get_element(q)
    ro_f_center = qubit_info.clock_freqs.readout()
    # avoid frequency conflicts 
    from numpy import NaN
    qubit_info.clock_freqs.readout(NaN)

    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,f_points)
    ro_p_samples = linspace(ro_p_min,ro_p_max,p_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
    ro_pulse_amp.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp=ro_pulse_amp,
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        powerDep=True,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     Power=['start '+'%E' %ro_p_samples[0],'end '+'%E' %ro_p_samples[-1]])
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((ro_f_samples,ro_p_samples))
        
        
        
        rp_ds = meas_ctrl.run("One-tone-powerDep")
        # Save the raw data into netCDF
        exp_timeLabel = get_time_now()
        raw_folder = build_folder_today(meas_raw_dir)
        rp_ds.to_netcdf(os.path.join(raw_folder,f"{q}_PowCavSpec_{exp_timeLabel}.nc"))
        print(f"Raw exp data had been saved into netCDF with the time label '{exp_timeLabel}'")

        analysis_result[q] = Basic2DAnalysis(tuid=rp_ds.attrs["tuid"], dataset=rp_ds).run()
        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s = 2
        sweep_para1= array(ro_f_samples[:n_s])
        sweep_para2= array(ro_p_samples[:2])
        spec_sched_kwargs['frequencies']= sweep_para1.reshape(sweep_para1.shape or (1,))
        spec_sched_kwargs['R_amp']= {q:sweep_para2.reshape(sweep_para2.shape or (1,))[0]}
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Modularize.Experiment_setup import get_FluxController

    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_2_26/SumInfo.pkl'
    QDmanager, cluster, meas_ctrl, ic = init_meas(QuantumDevice_path=QD_path,mode='l')
    Fctrl = get_FluxController(cluster)
    # default the offset in circuit
    reset_offset(Fctrl)
    # Set system attenuation
    init_system_atte(QDmanager.quantum_device,list(Fctrl.keys()))
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    error_log = []
    for qb in Fctrl:
        print(qb)
        PD_results = PowerDep_spec(QDmanager.quantum_device,meas_ctrl,q=qb)
        if PD_results == {}:
            error_log.append(qb)
        else:
            # TODO: Once the analysis for power dependence completed, fill in the answer to the quantum device here.
            pass
    if error_log != []:
        print(f"Power dependence error qubit: {error_log}")
    QDmanager.refresh_log('after PowerDep')
    QDmanager.QD_keeper()
    print('Power dependence done!')
    shut_down(cluster,Fctrl)