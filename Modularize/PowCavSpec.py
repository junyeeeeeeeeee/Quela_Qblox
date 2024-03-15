from numpy import array, linspace
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import Data_manager, QDmanager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from Modularize.Pulse_schedule_library import One_tone_sche, pulse_preview

def PowerDep_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_span_Hz:int=3e6,ro_p_min:float=0.1,ro_p_max:float=0.7,n_avg:int=100,f_points:int=60,p_points:int=30,run:bool=True,q:str='q1',Experi_info:dict={})->dict:

    sche_func = One_tone_sche
        
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
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
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((ro_f_samples,ro_p_samples))
        
        
        
        rp_ds = meas_ctrl.run("One-tone-powerDep")
        # Save the raw data into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=rp_ds,qb=q,exp_type='PD')

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
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    
    qubit_info.clock_freqs.readout(ro_f_center)
    return analysis_result

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down

    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_14/DR2#171_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    # Set system attenuation
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=28)
    
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    error_log = []
    for qb in ["q1"]:
        print(qb)
        # Fctrl[qb](QD_agent.Fluxmanager.get_sweetBiasFor(qb))
        PD_results = PowerDep_spec(QD_agent,meas_ctrl,q=qb, ro_span_Hz=2.5e6)
        if PD_results == {}:
            error_log.append(qb)
        else:
            # TODO: Once the analysis for power dependence completed, fill in the answer to the quantum device here.
            pass
        Fctrl[qb](0.0)
    if error_log != []:
        print(f"Power dependence error qubit: {error_log}")
    QD_agent.refresh_log('after PowerDep')
    QD_agent.QD_keeper()
    print('Power dependence done!')
    shut_down(cluster,Fctrl)