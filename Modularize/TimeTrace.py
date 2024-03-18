
from numpy import array, linspace
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Pulse_schedule_library import Trace_sche, set_LO_frequency, pulse_preview


def Qubit_state_avg_timetrace(QD_agent:QDmanager,trace_recordlength:float,IF:float=150e6,n_avg:int=1000,run:bool=True,q:str='q1',Experi_info:dict={}):
    qubit_info = QD_agent.quantum_device.get_element(q)
    sche_func = Trace_sche 
    LO= 3e9+IF#qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    data = {}
    result = {}

    def state_dep_sched(ini_state:str,qubit_info):
        if ini_state.lower() == 'g':
            XYL = {str(q):0}
        elif ini_state.lower() == 'e':
            XYL = {str(q):0} #qubit_info.rxy.amp180(),
        else:
            raise KeyError(f"ini_state='{ini_state}' can not be recognized!")
        print(XYL)
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp= XYL,
            R_amp={str(q):qubit_info.measure.pulse_amp()},
            R_duration={str(q):qubit_info.measure.pulse_duration()},
            R_integration={str(q):qubit_info.measure.integration_time()},
            R_inte_delay=qubit_info.measure.acq_delay(),
            trace_recordlength=trace_recordlength,
        )
        
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
            )
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            t_ds= gettable.get()
            
            #show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            if Experi_info != {}:
                show_args(Experi_info(q))
            
            return t_ds
    
        else:
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
            
            #show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            if Experi_info != {}:
                show_args(Experi_info(q))
            return 0
            
        
    data['trace_recordlength']= trace_recordlength   
    data['g'] = state_dep_sched('g',qubit_info)
    data['e'] = state_dep_sched('e',qubit_info)
    result[q]=  data  
    # Current_Readout_IF= 5.95*1e9-R_F[q]
    # print('Current_Readout_IF=',Current_Readout_IF/1e6,'MHz')
    return result


if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Modularize.QuFluxFit import calc_Gcoef_inFbFqFd, calc_g, calc_fq_g_excluded
    from Pulse_schedule_library import Qubit_state_Avgtimetrace_plot


    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_17/DR2#171_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')

    # Set system attenuation
    # init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=10)
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    Fctrl['q2'](QD_agent.Fluxmanager.get_sweetBiasFor('q2'))
    Avgtimetrace_result= Qubit_state_avg_timetrace(QD_agent,
                                                   trace_recordlength=8*1e-6,
                                                   n_avg=100,
                                                   run=True,
                                                   q='q2')
    Fctrl['q2'](0.0)
    plot_data =True
    if plot_data:              
        for q in Avgtimetrace_result:   
            Baseband_data= Qubit_state_Avgtimetrace_plot(Avgtimetrace_result[q],fc=100*1e6,Digital_downconvert=False, IF=150e6)
