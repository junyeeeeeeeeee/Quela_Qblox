from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import XGateErrorTest
#// test okay.

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr1","last_ip":"11"}
target_qs:list = ["q4"]
un_trained_pulse:bool = True
max_gate_num = 300
shots:int = 10000


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = XGateErrorTest(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(target_qs,shots,max_gate_num,Execution,un_trained_pulse)
EXP.WorkFlow()
EXP.RunAnalysis()