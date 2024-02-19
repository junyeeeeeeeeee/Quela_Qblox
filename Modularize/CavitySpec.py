from numpy import array, linspace
from Modularize.support import build_folder_today
from Modularize.path_book import meas_raw_dir
from Modularize.Pulse_schedule_library import One_tone_sche, pulse_preview
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QuantumDevice, get_time_now
from quantify_core.measurement.control import MeasurementControl
import os

# TODO: Need a test on machine @ 02/17
def Cavity_spec(quantum_device:QuantumDevice,meas_ctrl:MeasurementControl,ro_bare_guess:dict,ro_span_Hz:int=5e6,n_avg:int=1000,points:int=200,run:bool=True,q:str='q1',Experi_info:dict={})->dict:
    """
        Do the cavity search by the given QuantumDevice with a given target qubit q. \n
        Please fill up the initial value about measure for qubit in QuantumDevice first, like: amp, duration, integration_time and acqusition_delay! 
    """
    sche_func = One_tone_sche
    qubit_info = quantum_device.get_element(q)
    analysis_result = {}
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,points)
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
                     )
    
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
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(ro_f_samples)
        
        rs_ds = meas_ctrl.run("One-tone")
        analysis_result[q] = ResonatorSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
        # save the xarrry into netCDF
        exp_timeLabel = get_time_now()
        raw_folder = build_folder_today(meas_raw_dir)
        rs_ds.to_netcdf(os.path.join(raw_folder,f"CavitySpectro_{exp_timeLabel}.nc"))
        print(f"Raw exp data had been saved into netCDF with the time label '{exp_timeLabel}'")
        print(f"{q} Cavity:")
        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s=2 
        sweep_para= array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result


if __name__ == "__main__":
    from Modularize.support import QD_reloader, QD_keeper, connect_clusters, configure_measurement_control_loop, shut_down
    from qblox_instruments import Cluster
    from Modularize.path_book import qdevice_backup_dir
    # Reload the QuantumDevice or build up a new one
    QD_path = '/Users/ratiswu/Documents/GitHub/Quela_Qblox/Modularize/QD_backup/2024_2_19/academia_sinica_device_2024-02-19_06-36-36_UTC.json'
    if QD_path != '':
        quantum_device = QD_reloader(QD_path)
        print("Quantum device successfully reloaded!")
    else:
        from Modularize.Experiment_setup import hardware_cfg
        from Modularize.support import create_quantum_device
        quantum_device = create_quantum_device(hardware_cfg, num_qubits=5)
    # Connect to the Qblox cluster
    connect, ip = connect_clusters()
    cluster = Cluster(name = "cluster0", identifier = ip.get(connect.value))
    meas_ctrl, instrument_coordinator = configure_measurement_control_loop(quantum_device, cluster)
    
    # TODO: keep the particular bias into chip spec, and it should be fast loaded
    flux_settable_map: callable = {
        "q1":cluster.module2.out0_offset,
        "q2":cluster.module2.out1_offset,
        "q3":cluster.module2.out2_offset,
        "q4":cluster.module2.out3_offset,
        "q5":cluster.module10.out0_offset
    }
    # default the offset in circuit
    for i in flux_settable_map:
        flux_settable_map[i](0.00)
    
    # guess
    ro_bare=dict(
        q1 = 5.720 * 1e9,
        q2 = 5.8 * 1e9,
        q3 = 5.9 * 1e9,
        q4 = 6 * 1e9,
        q5 = 6.1 * 1e9,
    )
    error_log = []
    for qb in ro_bare:
        print(qb)
        if QD_path == '':
            qubit = quantum_device.get_element(qb)
            qubit.reset.duration(150e-6)
            qubit.measure.acq_delay(0)
            qubit.measure.pulse_amp(0.1)
            qubit.measure.pulse_duration(2e-6)
            qubit.measure.integration_time(2e-6)
        CS_results = Cavity_spec(quantum_device,meas_ctrl,ro_bare,q=qb)
        if CS_results != {}:
            print(f'Cavity @ {CS_results[qb].quantities_of_interest["fr"].nominal_value} Hz')
            quantum_device.get_element(qb).clock_freqs.readout(CS_results[qb].quantities_of_interest["fr"].nominal_value)
        else:
            error_log.append(qb)
    print(f"Cavity Spectroscopy error qubit: {error_log}")
    exp_timeLabel = get_time_now()
    qd_folder = build_folder_today(qdevice_backup_dir)
    QD_keeper(quantum_device,qd_folder)
    print('CavitySpectro done!')
    
    # shut_down(cluster,flux_settable_map)
