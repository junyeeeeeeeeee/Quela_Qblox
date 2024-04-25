import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from xarray import Dataset
from numpy import array, linspace
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager,init_system_atte, init_meas, shut_down
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Pulse_schedule_library import Qubit_state_single_shot_plot
from Modularize.support.Pulse_schedule_library import Qubit_SS_sche, set_LO_frequency, pulse_preview, Qubit_state_single_shot_fit_analysis


def Qubit_state_single_shot(QD_agent:QDmanager,shots:int=1000,run:bool=True,q:str='q1',IF:float=150e6,Experi_info:dict={},ro_amp_factor:float=1,T1:float=15e-6,exp_idx:int=0,parent_datafolder:str=''):
    qubit_info = QD_agent.quantum_device.get_element(q)
    sche_func = Qubit_SS_sche  
    LO= qubit_info.clock_freqs.f01()+IF
    rof = QD_agent.quantum_device.get_element(q).clock_freqs.readout()-0.2e6
    QD_agent.quantum_device.get_element(q).clock_freqs.readout(rof)
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
    SS_dict = {
        "e":{"dims":("I","Q"),"data":array(data['e'])},
        "g":{"dims":("I","Q"),"data":array(data['g'])},
    }
    SS_ds = Dataset.from_dict(SS_dict)
    Data_manager().save_raw_data(QD_agent=QD_agent,ds=SS_ds,qb=q,exp_type='ss',label=exp_idx,specific_dataFolder=parent_datafolder)
    if parent_datafolder =='':
        analysis_result[q] = Qubit_state_single_shot_fit_analysis(data,T1=T1,tau=tau) 
    else:
        analysis_result[q] = []
    return analysis_result


def SS_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict,target_q:str,shots:int=5000,execution:bool=True,data_folder='',plot:bool=True,roAmp_modifier:float=1,exp_label:int=0):
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(target_q,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(target_q,'ro'))
    Fctrl[target_q](float(QD_agent.Fluxmanager.get_tuneawayBiasFor(target_q)))
    SS_result= Qubit_state_single_shot(QD_agent,
                shots=shots,
                run=execution,
                q=target_q,
                parent_datafolder=data_folder,
                ro_amp_factor=roAmp_modifier,
                exp_idx=exp_label)
    Fctrl[target_q](0.0)
    cluster.reset()
    if plot:
        Qubit_state_single_shot_plot(SS_result[target_q],Plot_type='both',y_scale='log')



if __name__ == '__main__':
    

    """ Fill in """
    execute = True
    QD_path = 'Modularize/QD_backup/2024_4_25/DR1#11_SumInfo.pkl'
    ro_elements = {'q0':{"roAmp_factor":1}}
    

    """ Preparation """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    
    # QD_agent.Notewriter.save_DigiAtte_For(QD_agent.Notewriter.get_DigiAtteFor('q0','ro')-2,'q0','ro')

    """ Running """
    for qubit in ro_elements:
        ro_amp_scaling = ro_elements[qubit]["roAmp_factor"]
        
        SS_executor(QD_agent,cluster,Fctrl,qubit,execution=execute,roAmp_modifier=ro_amp_scaling)
        if ro_amp_scaling !=1:
            keep = input(f"Keep this RO amp for {qubit}?[y/n]")
        else:
            keep = 'y'

    """ Storing """ 
    if execute:
        if keep.lower() == 'y':
            QD_agent.QD_keeper() 
            
    """ Close """    
    shut_down(cluster,Fctrl)