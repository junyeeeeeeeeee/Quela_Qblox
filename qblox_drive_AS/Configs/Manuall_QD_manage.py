""" Set the readout parameters here about (1) integration time, (2) reset time (3) Driving IF """
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support import QDmanager, init_system_atte
from qblox_drive_AS.support.Pulse_schedule_library import set_LO_frequency


def set_reset_time_by_T1(QD_agent:QDmanager,target_q:str)->float:
    T1 = float(QD_agent.Notewriter.get_T1For(target_q))
    if T1 != 0:
        reset_time = 10*T1 
        reset_time = (reset_time*1e9 - (reset_time*1e9)%4)*1e-9
    else:
        reset_time = 250e-6
    
    return reset_time

def set_roLOfreq_atte(QD_agent:QDmanager,target_q:str,LO_Hz:float,ro_atte:int=None):
    """ ## *Warning*: 
        Set the LO or ro atte for those qubits who shares the same readout module with the `target_q`.
    """
    if LO_Hz is not None:
        set_LO_frequency(QD_agent.quantum_device,target_q,'readout',LO_Hz)
    if ro_atte is not None:
        init_system_atte(QD_agent.quantum_device,list([target_q]),ro_out_att=ro_atte)

if __name__ == "__main__":
    #// Reset time can be set in str 'auto', which will use 10*T1 if there is the T1 record in the QD_agent. Otherwise, set it as 250e-6. 
    #! If you don't want update an item like 'reset_time', please fill in with None.
    #! If you fill in `reset_time` with 'auto', we will calculate 10*T1 for you. 
    Revisers:dict = {
        "q0":{
            "integration_time":None, "reset_time":300e-6, "XY_if":-150e6
        },
        "q1":{
            "integration_time":2e-6, "reset_time":'auto', "XY_if":-150e6
        }
    }   

    QD_path = ""
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    for q in Revisers:
        qubit_info = QD_agent.quantum_device.get_element(q)
        for item in Revisers[q]:
            if Revisers[q][item] is not None:
                match str(item).lower():
                    case "integration_time":
                        qubit_info.measure.pulse_duration(Revisers[q][item])
                        qubit_info.measure.integration_time(Revisers[q][item])
                    case "reset_time":
                        qubit_info.reset.duration(Revisers[q][item] if Revisers[q][item] != 0 else set_reset_time_by_T1(QD_agent,q))
                    case "xy_if":
                        QD_agent.Notewriter.save_xyIF_for(q,Revisers[q][item])
                    case _:
                        print(f"Unknown item was given for {q} as name = {item}")

    # if you want changhe the RO LO freq, set LO in Hz as you want.
    # if you want change the RO attenuation, set ro_atte as you want with multiples of 2 and in the range [0,60]
    
    set_roLOfreq_atte(QD_agent, 'q0', LO_Hz=None, ro_atte=None)
    QD_agent.QD_keeper()