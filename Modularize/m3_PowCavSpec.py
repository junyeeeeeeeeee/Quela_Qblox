import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, linspace
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import Data_manager, QDmanager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.Pulse_schedule_library import One_tone_sche, pulse_preview

def PowerDep_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_span_Hz:int=4e6,ro_p_min:float=0.01,ro_p_max:float=0.5,n_avg:int=100,f_points:int=10,p_points:int=20,run:bool=True,q:str='q1',Experi_info:dict={})->dict:

    sche_func = One_tone_sche
        
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    ro_f_start = qubit_info.clock_freqs.readout()-1e6
    # avoid frequency conflicts 
    from numpy import NaN
    qubit_info.clock_freqs.readout(NaN)

    ro_f_samples = linspace(ro_f_start,ro_f_start+2*ro_span_Hz,f_points)
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
    
    qubit_info.clock_freqs.readout(ro_f_start+1e6)
    return analysis_result


def powerCavity_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,ro_span_Hz:float=3e6,max_power:float=0.7,run:bool=True,sweet_spot:bool=False,fpts:int=30,ppts:int=30,avg_n:int=100):

    if run:
        if sweet_spot:
            Fctrl[specific_qubits](QD_agent.Fluxmanager.get_sweetBiasFor(target_q=specific_qubits))
        PD_results = PowerDep_spec(QD_agent,meas_ctrl,q=specific_qubits, ro_span_Hz=ro_span_Hz,ro_p_max=max_power,f_points=fpts,p_points=ppts,n_avg=avg_n)
        Fctrl[specific_qubits](0.0)
        if PD_results == {}:
            print(f"Power dependence error qubit: {specific_qubits}")
     
    else:
        PD_results = PowerDep_spec(QD_agent,meas_ctrl,q=specific_qubits, ro_span_Hz=ro_span_Hz,run=False,ro_p_max=max_power)

    

if __name__ == "__main__":
    
    """ fill in """
    execution = True
    sweetSpot_dispersive = True
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {    # measurement target q from this dict 
        "q2": {"ro_atte":30},
    }

    """ preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')

    """ Running """
    for qubit in ro_elements:
        QD_agent.Notewriter.save_DigiAtte_For(ro_elements[qubit]["ro_atte"],qubit,'ro')
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
        powerCavity_executor(QD_agent,meas_ctrl,Fctrl,specific_qubits=qubit,run=execution,sweet_spot=sweetSpot_dispersive,max_power=0.5,ro_span_Hz=5e6, fpts=60)
        cluster.reset()
        if not execution:
            break
    
    QD_agent.refresh_log('after PowerDep')
    
    """ Storing """
    if execution: 
        QD_agent.QD_keeper()
    
    """ Close """
    print('Power dependence done!')
    shut_down(cluster,Fctrl)

    