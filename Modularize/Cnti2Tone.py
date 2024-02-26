from numpy import array, linspace
from Modularize.Pulse_schedule_library import Two_tone_sche, set_LO_frequency, pulse_preview, IQ_data_dis, QS_fit_analysis, dataset_to_array
from Modularize.support import build_folder_today
from Modularize.path_book import meas_raw_dir
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QuantumDevice, get_time_now
from quantify_core.measurement.control import MeasurementControl
import os


def Two_tone_spec(quantum_device:QuantumDevice,meas_ctrl:MeasurementControl,f01_guess:dict,xyf_span_Hz:int=200e6,xyamp:float=0.02,n_avg:int=1000,points:int=200,run:bool=True,q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0]):
    
    sche_func = Two_tone_sche   
    analysis_result = {}
    qubit_info = quantum_device.get_element(q)
    f01_high = f01_guess[q]+100*1e6
    f01_samples = linspace(f01_high-xyf_span_Hz,f01_high,points)
    set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=f01_high)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True

    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        spec_amp=xyamp,
        spec_Du=50*1e-6,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(f01_samples)
        
        qs_ds = meas_ctrl.run("Two-tone")
        # Save the raw data into netCDF
        exp_timeLabel = get_time_now()
        raw_folder = build_folder_today(meas_raw_dir)
        qs_ds.to_netcdf(os.path.join(raw_folder,f"{q}_TwTone_{exp_timeLabel}.nc"))
        print(f"Raw exp data had been saved into netCDF with the time label '{exp_timeLabel}'")

        I,Q= dataset_to_array(dataset=qs_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1]) # data = dis for plotting
        analysis_result[q] = QS_fit_analysis(data,f=f01_samples)
        
        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s = 2
        sweep_para= array(f01_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))

    return analysis_result







if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Modularize.Pulse_schedule_library import Fit_analysis_plot
    from Modularize.Experiment_setup import get_FluxController
    from numpy import NaN

    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_2_26/SumInfo.pkl'
    QDmanager, cluster, meas_ctrl, ic = init_meas(QuantumDevice_path=QD_path,mode='l')
    Fctrl = get_FluxController(cluster)
    reset_offset(Fctrl)
    # Set system attenuation
    init_system_atte(QDmanager.quantum_device,list(Fctrl.keys()))
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    XYf_guess=dict(
        q0 = 3.9e9,
        q1 = 4.232e9,
        q2 = 3.841e9,
        q3 = 4.022e9,
    )
    # q4 = 2.5738611635902258 * 1e9,
    
    for qb in ["q0"]:
        qubit = QDmanager.quantum_device.get_element(qb)
        qubit.clock_freqs.f01(NaN)
        for i in Fctrl:
            if i != qb:
                tuneaway = QDmanager.Fluxmanager.get_tuneawayBiasFor(i)
                if abs(tuneaway) <= 0.3:
                    Fctrl[i](tuneaway)
                else:
                    raise ValueError(f"tuneaway bias wrong! = {tuneaway}")


        Fctrl[qb](QDmanager.Fluxmanager.get_sweetBiasFor(qb))
        QS_results = Two_tone_spec(QDmanager.quantum_device,meas_ctrl,XYf_guess,q=qb,ref_IQ=QDmanager.refIQ[qb],xyf_span_Hz=100e6)

        reset_offset(Fctrl)
        Fit_analysis_plot(QS_results[qb],P_rescale=False,Dis=0)
        Revised_f01= QS_results[qb].attrs['f01_fit']
        qubit = QDmanager.quantum_device.get_element(qb)
        qubit.clock_freqs.f01(Revised_f01)

    QDmanager.refresh_log("After continuous 2-tone!")
    QDmanager.QD_keeper()
    print('2-tone done!')
    shut_down(cluster,Fctrl)

