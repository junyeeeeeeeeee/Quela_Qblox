from numpy import array, linspace
from Modularize.support.Pulse_schedule_library import Z_gate_two_tone_sche, set_LO_frequency, pulse_preview
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QDmanager, Data_manager
from quantify_core.measurement.control import MeasurementControl
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from numpy import NaN
import os


def Zgate_two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Z_amp_start:float,Z_amp_end:float,IF:int=100e6,xyf:float=0e9,xyf_span_Hz:float=300e6,n_avg:int=1000,RO_z_amp:float=0,Z_points:int=20,f_points:int=30,run:bool=True,q:str='q1',Experi_info={},get_data_path:bool=False):
    print("Zgate 2tone start")
    trustable = True
    sche_func = Z_gate_two_tone_sche
        
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    if xyf == 0:
        xyf_highest = qubit_info.clock_freqs.f01()+IF
    else:
        xyf_highest = xyf + IF
    qubit_info.clock_freqs.f01(NaN)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf_highest)
    f01_samples = linspace(xyf_highest-xyf_span_Hz,xyf_highest,f_points)
    
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    
    # temperature quard
    if Z_amp_end > 0.4:
        Z_amp_end = 0.4
    elif Z_amp_end < -0.4:
        Z_amp_end = -0.4
    else:
        pass

    if Z_amp_start > 0.4:
        Z_amp_start = 0.4
    elif Z_amp_start < -0.4:
        Z_amp_start = -0.4
    else:
        pass 

    Z_samples = linspace(Z_amp_start,Z_amp_end,Z_points)
    print(Z_samples[0], Z_samples[-1])
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        Z_amp=Z_bias,
        spec_amp=qubit_info.rxy.amp180(),
        spec_Du=50*1e-6,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        Z_ro_amp=RO_z_amp,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
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
        meas_ctrl.settables([freq,Z_bias])
        meas_ctrl.setpoints_grid((f01_samples,Z_samples))
        qs_ds = meas_ctrl.run("Zgate-two-tone")
        print(qs_ds)
        # Save the raw data into netCDF
        if get_data_path:
            path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='F2tone',get_data_loc=get_data_path)
        else:
            path = ''
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='F2tone',get_data_loc=get_data_path)

        try:
            analysis_result[q] = QubitFluxSpectroscopyAnalysis(tuid=qs_ds.attrs["tuid"], dataset=qs_ds).run()
        except:
            analysis_result[q] = {}
            print("Qb vs Flux fitting failed! Raw data had been saved.")
            trustable = False
        
        
        show_args(exp_kwargs, title="Zgate_two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s = 2
        sweep_para1= array(f01_samples[:n_s])
        sweep_para2= array(Z_samples[:2])
        spec_sched_kwargs['frequencies']= sweep_para1.reshape(sweep_para1.shape or (1,))
        spec_sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        
        
        show_args(exp_kwargs, title="Zgate_two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        path = ''

    return analysis_result, path, trustable

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from numpy import absolute as abs

    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_21/DR2#171_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    execute = True
    for qb in ["q0"]:
        # for i in Fctrl:
        #     if i != qb:
        #         tuneaway = QDmanager.Fluxmanager.get_tuneawayBiasFor(i)
        #         if abs(tuneaway) <= 0.3:
        #             Fctrl[i](tuneaway)
        #         else:
        #             raise ValueError(f"tuneaway bias wrong! = {tuneaway}")
        init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qb,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qb,'xy'))
        center = QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qb)
        half_period = QD_agent.Fluxmanager.get_PeriodFor(target_q=qb)/6
        window_shifter = 0
        results, _, trustable= Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=center-half_period+window_shifter,Z_points=10,Z_amp_end=center+half_period+window_shifter,q=qb,run=execute)
        reset_offset(Fctrl)
        if trustable:
            qubit = QD_agent.quantum_device.get_element(qb)
            qubit.clock_freqs.f01(results[qb].quantities_of_interest["freq_0"].nominal_value)
            QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=results[qb].quantities_of_interest["offset_0"].nominal_value)
            QD_agent.QD_keeper()
    
    print('Flux qubit done!')
    shut_down(cluster,Fctrl)
    


