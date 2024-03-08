from numpy import array, linspace
from Modularize.Pulse_schedule_library import Two_tone_sche, set_LO_frequency, pulse_preview, IQ_data_dis, QS_fit_analysis, dataset_to_array, twotone_comp_plot
from Modularize.support import build_folder_today
from Modularize.path_book import meas_raw_dir
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QDmanager, get_time_now
from quantify_core.measurement.control import MeasurementControl
from numpy import NaN
import os

def Two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,f01_guess:int=0,xyf_span_Hz:int=400e6,xyamp:float=0.02,n_avg:int=500,points:int=200,run:bool=True,q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0]):
    
    sche_func = Two_tone_sche   
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    # qubit_info.reset.duration(0)
    if f01_guess != 0:
        f01_high = f01_guess+100e6
    else:
        f01_high = qubit_info.clock_freqs.f01()+100e6
    # if xyamp == 0:
    #     xyamp = qubit_info.rxy.amp180(XYL)
    # Avoid warning
    qubit_info.clock_freqs.f01(NaN)

    f01_samples = linspace(f01_high-xyf_span_Hz,f01_high,points)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=f01_high)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True

    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        spec_amp=xyamp,
        spec_Du=100e-6,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        ref_pt="start"
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
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
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        analysis_result[q] = {}

    return analysis_result, f01_high-100e6



if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Modularize.Pulse_schedule_library import Fit_analysis_plot
    from Modularize.Experiment_setup import get_FluxController
    from numpy import arange

    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_8/DR1#170_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic = init_meas(QuantumDevice_path=QD_path,cluster_ip='170',mode='l')
    Fctrl = get_FluxController(cluster,ip_label=QD_agent.Identity.split("#")[-1])
    reset_offset(Fctrl)
    # Set system attenuation
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=30)
    # for i in range(6):
    #     getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
    #     getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    Fctrl["q4"](float(QD_agent.Fluxmanager.get_sweetBiasFor("q4")))

    XYf_guess=dict(
        q4 = [3.1e9]
    )
    # q4 = 2.5738611635902258 * 1e9,
    
    xyamp_dict = dict(
        q4 = linspace(0,0.6,13)
    )
    for qb in XYf_guess:
        print(f"{qb} is under the measurement!")
        for XYF in XYf_guess[qb]:
            ori_data = []
            for XYL in xyamp_dict[qb]:
                print(f"XYL now ={XYL}, XYF={XYF*1e-9-0.2}~{XYF*1e-9+0.1}GHz")
                qubit = QD_agent.quantum_device.get_element(qb)
                # for i in Fctrl:
                #     if i != qb:
                #         tuneaway = QD_agent.Fluxmanager.get_tuneawayBiasFor(i)
                #         if abs(tuneaway) <= 0.3:
                #             Fctrl[i](tuneaway)
                #         else:
                #             raise ValueError(f"tuneaway bias wrong! = {tuneaway}")

                print(f"bias = {QD_agent.Fluxmanager.get_sweetBiasFor(qb)}")
                QS_results, origin_f01 = Two_tone_spec(QD_agent,meas_ctrl,xyamp=XYL,f01_guess=XYF,q=qb,xyf_span_Hz=500e6,points=200,n_avg=1000,run=True,ref_IQ=QD_agent.refIQ[qb]) # 
                # Fit_analysis_plot(QS_results[qb],P_rescale=False,Dis=0)
                if XYL != 0:
                    twotone_comp_plot(QS_results[qb], ori_data, True)
                else:
                    twotone_comp_plot(QS_results[qb], ori_data, False)
                    ori_data = QS_results[qb].data_vars['data']
                Revised_f01= QS_results[qb].attrs['f01_fit']
                qubit = QD_agent.quantum_device.get_element(qb)
                # qubit.clock_freqs.f01(Revised_f01)
                # temporarily save the xy amp for a clear f01 fig
                # qubit.rxy.amp180(XYL)

    QD_agent.refresh_log("After continuous 2-tone!")
    # QD_agent.QD_keeper()
    print('2-tone done!')
    shut_down(cluster,Fctrl)

