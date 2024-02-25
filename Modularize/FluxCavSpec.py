from numpy import array, linspace
from Modularize.support import build_folder_today
from Modularize.path_book import meas_raw_dir
from Modularize.Pulse_schedule_library import One_tone_sche, pulse_preview

from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QuantumDevice, get_time_now
from quantify_core.measurement.control import MeasurementControl
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
import os

def FluxCav_spec(quantum_device:QuantumDevice,meas_ctrl:MeasurementControl,flux_ctrl:dict,ro_span_Hz:int=3e6,flux_span:float=0.3,n_avg:int=500,f_points:int=30,flux_points:int=40,run:bool=True,q:str='q1',Experi_info:dict={}):

    sche_func = One_tone_sche
        
    analysis_result = {}
    qubit_info = quantum_device.get_element(q)
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
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,flux_ctrl[q]])
        meas_ctrl.setpoints_grid((ro_f_samples,flux_samples))
        
        
        
        rfs_ds = meas_ctrl.run("One-tone-Flux")
        # Save the raw data into netCDF
        exp_timeLabel = get_time_now()
        raw_folder = build_folder_today(meas_raw_dir)
        rfs_ds.to_netcdf(os.path.join(raw_folder,f"{q}_FluxCavSpec_{exp_timeLabel}.nc"))
        print(f"Raw exp data had been saved into netCDF with the time label '{exp_timeLabel}'")

        analysis_result[q] = ResonatorFluxSpectroscopyAnalysis(tuid=rfs_ds.attrs["tuid"], dataset=rfs_ds).run()
        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        # reset flux bias
        flux_ctrl[q](0.0)
        
    else:
        n_s = 2
        sweep_para= array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down
    from numpy import pi
    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_2_25/SumInfo.pkl'
    QDmanager, cluster, meas_ctrl, ic = init_meas(QuantumDevice_path=QD_path,mode='l')
    
    Fctrl: callable = {
        "q0":cluster.module2.out0_offset,
        "q1":cluster.module2.out1_offset,
        "q2":cluster.module2.out2_offset,
        "q3":cluster.module2.out3_offset,
        # "q4":cluster.module10.out0_offset
    }
    
    # default the offset in circuit
    for i in Fctrl:
        Fctrl[i](0.0)
    # Set system attenuation
    init_system_atte(QDmanager.quantum_device,list(Fctrl.keys()))
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50) 

    error_log = []
    for qb in ["q0","q1","q3"]:
        print(qb)
        FD_results = FluxCav_spec(QDmanager.quantum_device,meas_ctrl,Fctrl,q=qb)
        if FD_results == {}:
            error_log.append(qb)
        else:
            
            qubit = QDmanager.quantum_device.get_element(qb)
            print(f"@ {qb} fitting ans: {FD_results[qb].quantities_of_interest}")
            qubit.clock_freqs.readout(FD_results[qb].quantities_of_interest["freq_0"])
            QDmanager.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=FD_results[qb].quantities_of_interest["offset_0"].nominal_value)
            QDmanager.Fluxmanager.save_period_for(target_q=qb, period=2*pi/FD_results[qb].quantities_of_interest["frequency"].nominal_value)
            QDmanager.Fluxmanager.save_tuneawayBias_for(target_q=qb,mode='auto')
            QDmanager.Fluxmanager.save_cavFittingParas_for(target_q=qb,
                f=FD_results[qb].quantities_of_interest["frequency"].nominal_value,
                amp=FD_results[qb].quantities_of_interest["amplitude"].nominal_value,
                phi=FD_results[qb].quantities_of_interest["shift"].nominal_value,
                offset=FD_results[qb].quantities_of_interest["offset"].nominal_value
            )

    print(f"Flux dependence error qubit: {error_log}")
    QDmanager.refresh_log("after FluxDep")
    QDmanager.QD_keeper()
    print('Flux dependence done!')
    print("Flux dict: ",QDmanager.Fluxmanager.get_bias_dict())
    shut_down(cluster,Fctrl)
    