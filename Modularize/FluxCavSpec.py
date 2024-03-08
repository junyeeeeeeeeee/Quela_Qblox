from numpy import array, linspace
from Modularize.support import build_folder_today, QDmanager, get_time_now
from Modularize.path_book import meas_raw_dir
from Modularize.Pulse_schedule_library import One_tone_sche, pulse_preview
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
import os

def FluxCav_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,flux_ctrl:dict,ro_span_Hz:int=3e6,flux_span:float=0.3,n_avg:int=500,f_points:int=30,flux_points:int=40,run:bool=True,q:str='q1',Experi_info:dict={}):

    sche_func = One_tone_sche
        
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    ro_f_center = qubit_info.clock_freqs.readout()
    # avoid frequency conflicts 
    from numpy import NaN
    qubit_info.clock_freqs.readout(NaN)

    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,f_points)
    flux_samples = linspace(-flux_span,flux_span,flux_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        powerDep=False,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     Flux=['start '+'%E' %flux_samples[0],'end '+'%E' %flux_samples[-1]])
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,flux_ctrl[q]])
        meas_ctrl.setpoints_grid((ro_f_samples,flux_samples))
        
        
        
        rfs_ds = meas_ctrl.run("One-tone-Flux")
        # Save the raw data into netCDF
        exp_timeLabel = get_time_now()
        raw_folder = build_folder_today(meas_raw_dir)
        dr_loc = QD_agent.Identity.split("#")[0]
        rfs_ds.to_netcdf(os.path.join(raw_folder,f"{dr_loc}{q}_FluxCavSpec_{exp_timeLabel}.nc"))
        print(f"Raw exp data had been saved into netCDF with the time label '{exp_timeLabel}'")

        analysis_result[q] = ResonatorFluxSpectroscopyAnalysis(tuid=rfs_ds.attrs["tuid"], dataset=rfs_ds).run()
        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        # reset flux bias
        flux_ctrl[q](0.0)
        qubit_info.clock_freqs.readout(ro_f_center)
        
    else:
        n_s = 2
        sweep_para= array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Modularize.Experiment_setup import get_FluxController
    from numpy import pi
    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_8/DR1#170_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic = init_meas(QuantumDevice_path=QD_path,cluster_ip='170',mode='l')
    Fctrl = get_FluxController(cluster,ip_label=QD_agent.Identity.split("#")[-1])
    # default the offset in circuit
    reset_offset(Fctrl)
    # Set system attenuation
    # init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=36)
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50) 

    execute = True
    error_log = []
    for qb in ["q4"]:
        print(f"{qb} are under the measurement ...")
        FD_results = FluxCav_spec(QD_agent,meas_ctrl,Fctrl,ro_span_Hz=0.8e6,q=qb,flux_span=0.4,run=execute)
        if FD_results == {}:
            error_log.append(qb)
        else:
            if execute:
                qubit = QD_agent.quantum_device.get_element(qb)
                qubit.clock_freqs.readout(FD_results[qb].quantities_of_interest["freq_0"])
                QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=FD_results[qb].quantities_of_interest["offset_0"].nominal_value)
                QD_agent.Fluxmanager.save_period_for(target_q=qb, period=2*pi/FD_results[qb].quantities_of_interest["frequency"].nominal_value)
                QD_agent.Fluxmanager.save_tuneawayBias_for(target_q=qb,mode='auto')
                QD_agent.Fluxmanager.save_cavFittingParas_for(target_q=qb,
                    f=FD_results[qb].quantities_of_interest["frequency"].nominal_value,
                    amp=FD_results[qb].quantities_of_interest["amplitude"].nominal_value,
                    phi=FD_results[qb].quantities_of_interest["shift"].nominal_value,
                    offset=FD_results[qb].quantities_of_interest["offset"].nominal_value
                )
    if error_log != []:
        print(f"Flux dependence error qubit: {error_log}")
    if execute:
        QD_agent.refresh_log("after FluxDep")
        QD_agent.QD_keeper()
    print('Flux dependence done!')
    shut_down(cluster,Fctrl)
    