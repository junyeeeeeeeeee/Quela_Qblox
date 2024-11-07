""" Set the readout parameters here about (1) integration time, (2) reset time """
from Modularize.support.QDmanager import QDmanager

def get_manully_integration_time(target_q:str)->float:
    return Readers[target_q]["integration_time"]

def get_manully_reset_time(target_q:str)->float:
    return Readers[target_q]["reset_time"] if str(Readers[target_q]["reset_time"]).lower() != 'auto' else 0

def get_manully_rotate_angle(target_q:str)->float:
    return Readers[target_q]["rotate_on_I_angle_degree"]


def set_reset_time_by_T1(QD_agent:QDmanager,target_q:str)->float:
    T1 = float(QD_agent.Notewriter.get_T1For(target_q))
    if T1 != 0:
        reset_time = 10*T1 
        reset_time = (reset_time*1e9 - (reset_time*1e9)%4)*1e-9
    else:
        reset_time = 250e-6
    
    return reset_time

#// Reset time can be set in str 'auto', which will use 10*T1 if there is the T1 record in the QD_agent. Otherwise, set it as 250e-6. 
Readers:dict = {
    "q0":{
        "integration_time":2e-6, "reset_time":300e-6, "rotate_on_I_angle_degree":0
    },
    "q1":{
        "integration_time":2e-6, "reset_time":250e-6, "rotate_on_I_angle_degree":0
    }
}
