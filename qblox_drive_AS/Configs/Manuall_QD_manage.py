""" Set the Quantum device parameters here about (1) integration time, (2) reset time (3) Driving IF """
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support import QDmanager
from qblox_drive_AS.support.UserFriend import *
from quantify_scheduler.helpers.collections import find_port_clock_path
from qblox_drive_AS.support.Pulse_schedule_library import set_LO_frequency
from qblox_drive_AS.support.QDmanager import find_path_by_port


def set_reset_time_by_T1(QD_agent:QDmanager,target_q:str)->float:
    T1 = float(QD_agent.Notewriter.get_T1For(target_q))
    if T1 != 0:
        reset_time = 10*T1 
        reset_time = (reset_time*1e9 - (reset_time*1e9)%4)*1e-9
    else:
        reset_time = 250e-6
    
    return reset_time

def QD_modifier(QD_agent:QDmanager, toModify_dict:dict):
    for q in toModify_dict:
        qubit_info = QD_agent.quantum_device.get_element(q)
        for item in toModify_dict[q]:
            if toModify_dict[q][item] is not None:
                match str(item).lower():
                    case "integration_time":
                        qubit_info.measure.pulse_duration(toModify_dict[q][item])
                        qubit_info.measure.integration_time(toModify_dict[q][item])
                    case "reset_time":
                        qubit_info.reset.duration(toModify_dict[q][item] if toModify_dict[q][item] != 0 else set_reset_time_by_T1(QD_agent,q))
                    case "xy_if":
                        QD_agent.Notewriter.save_xyIF_for(q,toModify_dict[q][item])
                    case _:
                        print(f"Unknown item was given for {q} as name = {item}")

def set_roLOfreq(QD_agent:QDmanager,LO_Hz:float,target_q:str='q0'):
    """ ## *Warning*: 
        Set the LO for those qubits who shares the same readout module with the `target_q`.
    """
    if LO_Hz is not None:
        slightly_print(f"Set {find_port_clock_path(QD_agent.quantum_device.hardware_config(),'q:res',f'{target_q}.ro')[1]} RO-LO = {round(LO_Hz*1e-9),2} GHz")
        set_LO_frequency(QD_agent.quantum_device,target_q,'readout',LO_Hz)

def set_roAtte(QD_agent:QDmanager,ro_atte:int, target_q:str='q0'):
    """ ## *Warning*: 
        * Set the readout attenuation for those qubits who shares the same readout module with the `target_q`.\n
        ## Args:\n
        ro_atte (int): multiple of 2.
    """
    if ro_atte is not None:
        slightly_print(f"Set {find_port_clock_path(QD_agent.quantum_device.hardware_config(),'q:res',f'{target_q}.ro')[1]} RO-attenuation = {int(ro_atte)} dB")
        who_onTheSameBoat:dict = find_path_by_port(QD_agent.quantum_device.hardware_config(),"q:res")
        for q in who_onTheSameBoat:
            if who_onTheSameBoat[q][0] == who_onTheSameBoat[target_q][0]:
                print(f"set {q} ...")
                QD_agent.Notewriter.save_DigiAtte_For(ro_atte,q,'ro')

def set_sweet_g(QD_agent:QDmanager,g_dict:dict=None):
    """ save g is the `g_dict` for q_key, g unit in Hz"""
    if g_dict != {} or g_dict is not None:
        for q in g_dict:
            QD_agent.Notewriter.save_sweetG_for(g_dict[q],q)

def update_coupler_bias(QD_agent:QDmanager,cp_elements:dict):
    """
    Update the idle bias in Fluxmanager for couplers.\n
    --------------------------
    ### Args:\n
    cp_elements = {"c0":0.2}
    """
    if cp_elements != {} or cp_elements is not None:
        for cp in cp_elements:
            slightly_print(f"Set coupler {cp} at {cp_elements[cp]} V.")
            QD_agent.Fluxmanager.save_idleBias_for(cp, cp_elements[cp])

if __name__ == "__main__":
    #// Reset time can be set in str 'auto', which will use 10*T1 if there is the T1 record in the QD_agent. Otherwise, set it as 250e-6. 
    #! If you don't want update an item 'reset_time' for example, please fill in with None.
    #! If you fill in `reset_time` with 'auto', we will calculate 10*T1 for you. 
    
    QD_path = "qblox_drive_AS/QD_backup/20241118/DR2#10_SumInfo.pkl"
    Modifications:dict = {
        "q0":{
            "integration_time":None, "reset_time":300e-6, "XY_if":-150e6
        },
        "q1":{
            "integration_time":2e-6, "reset_time":'auto', "XY_if":-150e6
        }
    }   

    
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    """ re-set Modifications """
    # QD_modifier(QD_agent, Modifications)

    """ Set RO-LO, RO-atte """
    # if you want changhe the RO LO freq, set LO in Hz as you want.
    # if you want change the RO attenuation, set ro_atte as you want with multiples of 2 and in the range [0,60]
    # If you don't want to modify it, set it with None.
    
    set_roLOfreq(QD_agent, LO_Hz=None, target_q='q0') # LO is global in the same QRM-RF module
    set_roAtte(QD_agent, ro_atte=None, target_q='q0') # RO-attenuation is global in the same QRM-RF module

    """ Set coupler bias """
    # Memorize the coupler offset
    # If you don't want to modify it, set it with None or {}
    update_coupler_bias(QD_agent,cp_elements={})  # cp_elements = {"c0":0.1, "c2":0.05, ...}
    
    """ Set grq for sweet spot """
    set_sweet_g(QD_agent,g_dict={"q0":90e6,"q1":90e6})
    
    QD_agent.QD_keeper()