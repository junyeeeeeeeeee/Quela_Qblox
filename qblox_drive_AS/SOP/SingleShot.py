import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from numpy import array, moveaxis, arange
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import QDmanager, Data_manager, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_Qubit_SS_sche, pulse_preview


def Qubit_state_single_shot(QD_agent:QDmanager,ro_elements:list,shots:int=1000,run:bool=True):
    sche_func = multi_Qubit_SS_sche 

    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"{q} Integration time: {round(qubit_info.measure.integration_time()*1e6,1)} µs")

    
    folder = []

    def state_dep_sched(ini_state:str):
        slightly_print(f"Shotting for |{ini_state}>")
        sched_kwargs = dict(   
            ini_state=ini_state,
            waveformer = QD_agent.Waveformer,
            pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'d1'),
            pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements,'d3'),
            R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
            R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
            R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
        )
        
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(ro_elements),
            )
            QD_agent.quantum_device.cfg_sched_repetitions(shots)
            ss_da= gettable.get() # DataArray (2*ro_q, shots)
            reshaped_data = list(array(ss_da).reshape(len(ro_elements),2,shots)) # (ro_q, IQ, shots)
            folder.append(reshaped_data) # (state, ro_q, IQ, shots)

        else:
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

            
    state_dep_sched('g')
    state_dep_sched('e')
    folder = moveaxis(array(folder),0,2) # (state,ro_q,IQ,shots) -> (ro_q, IQ, state, shots)
    output_dict = {}
    for q_idx, q_name in enumerate(ro_elements):
        output_dict[q_name] = (["mixer","prepared_state","index"],folder[q_idx])

    SS_ds = Dataset(output_dict, coords= {"mixer":array(["I","Q"]), "prepared_state":array([0,1]),"index":arange(shots)})
    SS_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    SS_ds.attrs["execution_time"] = Data_manager().get_time_now()
    
    return SS_ds






        
        

        
    
