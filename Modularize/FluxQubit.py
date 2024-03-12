from numpy import array, linspace
from Modularize.Pulse_schedule_library import Z_gate_two_tone_sche, set_LO_frequency, pulse_preview
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support import QDmanager, Data_manager
from quantify_core.measurement.control import MeasurementControl
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from numpy import NaN
import os


def Zgate_two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Z_amp_start:float,Z_amp_end:float,xyf:float=0e9,xyf_span_Hz:float=300e6,n_avg:int=500,Z_points:int=26,f_points:int=26,run:bool=True,q:str='q1',Experi_info={}):
    
    sche_func = Z_gate_two_tone_sche
        
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    if xyf == 0:
        xyf_highest = qubit_info.clock_freqs.f01()+100e6
    else:
        xyf_highest = xyf + 100e6
    qubit_info.clock_freqs.f01(NaN)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf_highest)
    f01_samples = linspace(xyf_highest-xyf_span_Hz,xyf_highest,f_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    
    # temperature quard
    if Z_amp_end > 0.25:
        Z_amp_end = 0.25
    elif Z_amp_end < -0.25:
        Z_amp_end = -0.25
    else:
        pass

    if Z_amp_start > 0.25:
        Z_amp_start = 0.25
    elif Z_amp_start < -0.25:
        Z_amp_start = -0.25
    else:
        pass 

    Z_samples = linspace(Z_amp_start,Z_amp_end,Z_points)
    
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
        # Save the raw data into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='2tone')
        
        analysis_result[q] = QubitFluxSpectroscopyAnalysis(tuid=qs_ds.attrs["tuid"], dataset=qs_ds).run()
        
        
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
    return analysis_result, xyf_highest-100e6

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from numpy import absolute as abs

    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_11/DR2#171_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    
    # Set system attenuation
    # init_system_atte(QDmanager.quantum_device,list(Fctrl.keys()))
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    execute = True
    QD_agent.quantum_device.get_element("q0").rxy.amp180(0.19)
    for qb in ["q0"]:
        # for i in Fctrl:
        #     if i != qb:
        #         tuneaway = QDmanager.Fluxmanager.get_tuneawayBiasFor(i)
        #         if abs(tuneaway) <= 0.3:
        #             Fctrl[i](tuneaway)
        #         else:
        #             raise ValueError(f"tuneaway bias wrong! = {tuneaway}")

        center = QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qb)
        half_period = QD_agent.Fluxmanager.get_PeriodFor(target_q=qb)/8
        window_shifter = 0
        results, origin_f01 = Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=center-half_period+window_shifter,Z_amp_end=center+half_period+window_shifter,q=qb,run=execute)
        reset_offset(Fctrl)
        if execute:
            qubit = QD_agent.quantum_device.get_element(qb)
            # qubit.clock_freqs.f01(origin_f01)
    print('Flux qubit done!')
    shut_down(cluster,Fctrl)
    


