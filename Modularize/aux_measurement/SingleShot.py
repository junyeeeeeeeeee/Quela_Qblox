import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from numpy import array, linspace
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager,init_system_atte
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Pulse_schedule_library import Qubit_SS_sche, set_LO_frequency, pulse_preview, Qubit_state_single_shot_fit_analysis


def Qubit_state_single_shot(QD_agent:QDmanager,shots:int=1000,run:bool=True,q:str='q1',IF:float=150e6,Experi_info:dict={},ro_amp_factor:float=1,T1:float=15e-6,parent_datafolder:str=''):
    qubit_info = QD_agent.quantum_device.get_element(q)
    sche_func = Qubit_SS_sche  
    LO= qubit_info.clock_freqs.f01()+IF
    qubit_info.measure.pulse_amp(ro_amp_factor*qubit_info.measure.pulse_amp())
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    data = {}
    analysis_result = {}
    exp_kwargs= dict(shots=shots,
                     )
    
    def state_dep_sched(ini_state:str):
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp={str(q):qubit_info.rxy.amp180()},
            pi_dura={str(q):qubit_info.rxy.duration()},
            R_amp={str(q):qubit_info.measure.pulse_amp()},
            R_duration={str(q):qubit_info.measure.pulse_duration()},
            R_integration={str(q):qubit_info.measure.integration_time()},
            R_inte_delay=qubit_info.measure.acq_delay(),
        )
        
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
            )
            QD_agent.quantum_device.cfg_sched_repetitions(shots)
            ss_ds= gettable.get()
            
            data[ini_state] = ss_ds
            print(type(ss_ds), ss_ds)
            # Data_manager().save_raw_data(QD_agent=QD_agent,ds=ss_ds,qb=q,label=ini_state.upper(),exp_type='ss',specific_dataFolder=parent_datafolder)
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            if Experi_info != {}:
                show_args(Experi_info(q))
    
        else:
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
            
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            if Experi_info != {}:
                show_args(Experi_info(q))
            
    tau= qubit_info.measure.integration_time()        
    state_dep_sched('g')
    state_dep_sched('e')
        
    analysis_result[q]= Qubit_state_single_shot_fit_analysis(data,T1=T1,tau=tau)
        
    return analysis_result


def SS_executor(QD_agent:QDmanager,Fctrl:dict,target_q:str,shots:int=1e5,execution:bool=True,data_folder='',plot:bool=True,roAmp_modifier:float=1):
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(target_q,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(target_q,'ro'))
    Fctrl[target_q](float(QD_agent.Fluxmanager.get_sweetBiasFor(target_q)))
    SS_result= Qubit_state_single_shot(QD_agent,
                shots=5000,
                run=execution,
                q=target_q,
                parent_datafolder=data_folder,
                ro_amp_factor=roAmp_modifier)
    Fctrl[target_q](0.0)
    if plot:
        Qubit_state_single_shot_plot(SS_result[target_q],Plot_type='both',y_scale='log')



if __name__ == '__main__':
    from Modularize.support import init_meas, shut_down
    from Modularize.support.Pulse_schedule_library import Qubit_state_single_shot_plot

    # Reload the QuantumDevice or build up a new one
    """ Fill in """
    execute = True
    QD_path = 'Modularize/QD_backup/2024_4_24/DR1#11_SumInfo.pkl'
    ro_elements = {'q0':{"roAmp_factor":1.5}}


    """ Preparation """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    


    """ Running """
    for qubit in ro_elements:
        ro_amp_scaling = ro_elements[qubit]["roAmp_factor"]
        SS_executor(QD_agent,Fctrl,qubit,execution=execute,roAmp_modifier=ro_amp_scaling)

    """ Storing """ 
    QD_agent.QD_keeper() 

    """ Close """    
    shut_down(cluster,Fctrl)