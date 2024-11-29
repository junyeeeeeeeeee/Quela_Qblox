import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from numpy import array, moveaxis, arange
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import QDmanager, Data_manager, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import Gate_Test_SS_sche, pulse_preview
from qblox_drive_AS.support.WaveformCtrl import GateGenesis

def XGateError_single_shot(QD_agent:QDmanager,ro_elements:list, max_gate_num:int,shots:int=1000, untrained:bool=False,run:bool=True):
    sche_func = Gate_Test_SS_sche 
    sample_pts = 25 # fix

    if max_gate_num < 3:
        raise ValueError("Max gates number must be more than 3 gates.")
    odd_numbers = list(range(3, max_gate_num + 1, 2))
    if len(odd_numbers) > sample_pts:
        step = len(odd_numbers) // (sample_pts-1)  # Select step to get 24 values
        selected_values = odd_numbers[::step][:(sample_pts-1)]  # Select 24 evenly spaced values
        selected_values.append(odd_numbers[-1])  # Add the closest number to max_gates_num
    else:
        selected_values = odd_numbers

    pulse_repeats = array(list([0, 1]) + list(selected_values))

    for q in ro_elements:
        pi_dura = QD_agent.quantum_device.get_element(q).rxy.duration()
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.reset.duration(qubit_info.reset.duration()+max(pulse_repeats)*pi_dura)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"{q} Integration time: {round(qubit_info.measure.integration_time()*1e6,1)} µs")

    
    folder = []
    if untrained:
        slightly_print("Use default pi-pulse !")
        wf_packs = GateGenesis(q_num=len(list(QD_agent.quantum_device.elements())),c_num=0)
    else:
        wf_packs = QD_agent.Waveformer


    def repeat_dep_sched(seq_num:int):
        slightly_print(f"Shotting for {seq_num} pi-pulses")
        sched_kwargs = dict(   
            pulse_num=seq_num,
            waveformer = wf_packs,
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

    for i in pulse_repeats:        
        repeat_dep_sched(i)
    
    folder = moveaxis(array(folder),0,2) # (state,ro_q,IQ,shots) -> (ro_q, IQ, state, shots)
    output_dict = {}
    for q_idx, q_name in enumerate(ro_elements):
        output_dict[q_name] = (["mixer","pulse_num","index"],folder[q_idx])

    SS_ds = Dataset(output_dict, coords= {"mixer":array(["I","Q"]), "pulse_num":pulse_repeats,"index":arange(shots)})
    SS_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    SS_ds.attrs["execution_time"] = Data_manager().get_time_now()
    
    return SS_ds




        
        

        
    
