""" Thsi script helps you build a new QD_file """
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support import QDmanager
from qblox_drive_AS.support.QDmanager import hcfg_composer


cluster_IP:str = "192.168.1.11"
dr_name:str = "dr1"
qubit_number_onChip:int = 4
coupler_number_onChip:int = 0
chip_name:str = "8SQAlOS250402_3"
chip_type:str = "5Q4C"
old_QD_path:str = "qblox_drive_AS/QD_backup/20250417/DR1#11_SumInfo.pkl" # set the path in string When you want to update the Hcfg. Otherwise, set it None


Hcfg = {
    "connectivity": {
        "graph": [
            [
                f"cluster{dr_name}.module4.complex_output_1",
                "q1:mw"
            ],
            [
                f"cluster{dr_name}.module4.complex_output_0",
                "q2:mw"
            ],
            [
                f"cluster{dr_name}.module6.complex_output_0",
                "q1:res"
            ],
            [
                f"cluster{dr_name}.module6.complex_output_0",
                "q2:res"
            ],
        ]
    },
}


if old_QD_path is None or str(old_QD_path) == "" :
    QD_agent = QDmanager()
    QD_agent.build_new_QD(qubit_number_onChip,coupler_number_onChip,hcfg_composer(Hcfg),cluster_IP,dr_name,chip_name,chip_type)
else:
    QD_agent = QDmanager(old_QD_path)
    QD_agent.QD_loader(new_Hcfg=hcfg_composer(Hcfg))
    
QD_agent.QD_keeper()


    


