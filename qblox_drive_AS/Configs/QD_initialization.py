""" Thsi script helps you build a new QD_file """

from qblox_drive_AS.support import QDmanager
from qblox_drive_AS.support.QDmanager import hcfg_composer

cluster_IP:str = "192.168.1.11"
dr_name:str = "dr1"
qubit_number_onChip:int = 4
coupler_number_onChip:int = 0
chip_name:str = "8SQAlOS250402_3"
chip_type:str = "5Q4C"


Hcfg = [
    {"name":"q1:mw", "slot":4, "port":1},
    {"name":"q2:mw", "slot":4, "port":0},
    {"name":"q1:res", "slot":6, "port":0},
    {"name":"q2:res", "slot":6, "port":0},
]


QD_agent = QDmanager()
QD_agent.build_new_QD(qubit_number_onChip,coupler_number_onChip,hcfg_composer(Hcfg, dr_name),cluster_IP,dr_name,chip_name,chip_type)

    
QD_agent.QD_keeper()


    


