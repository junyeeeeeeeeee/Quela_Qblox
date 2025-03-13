from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import nSingleShot
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
target_qs:list = ["q0","q2"]
shots:int = 10000
histo_counts:int = 1 # use only when the fitting won't go wrong.

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = nSingleShot(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(target_qs,histo_counts,shots,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()