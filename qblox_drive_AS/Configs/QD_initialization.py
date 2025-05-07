""" Thsi script helps you build a new QD_file """
from qblox_drive_AS.support.QDmanager import QDmanager, hcfg_composer
####### Port Name rules #######
# 1. Driving port: ":mw", like 'q1:mw', 'q2:mw', ...
# 2. Readout port: ":res", like 'q0:res', 'q3:res', ...
# 3. Flux bias port: ":fl", like 'q12:fl', 'q999:fl', ...


cluster_IP:str = "192.168.1.11"
dr_name:str = "dr1"
qubit_number_onChip:int = 4
coupler_number_onChip:int = 0
chip_name:str = "Test_8SQ"
chip_type:str = "5Q4C"


Hcfg = [
    {"name":"q1:mw", "slot":4, "port":0},
    {"name":"q3:mw", "slot":4, "port":1},
    {"name":"q1:res", "slot":6, "port":0},
    {"name":"q2:res", "slot":6, "port":0},
]




""" Do NOT touch !"""
QD_agent = QDmanager()
QD_agent.build_new_QD(qubit_number_onChip,coupler_number_onChip,hcfg_composer(Hcfg, dr_name),cluster_IP,dr_name,chip_name,chip_type)
QD_agent.QD_keeper()


    


