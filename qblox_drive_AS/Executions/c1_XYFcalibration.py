from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import XYFcali
#// test okay


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
target_qs:list = ["q0", "q1"]
AVG:int = 500

""" try change it ONLY when fitting goes wrong """
evoT = 0.5e-6

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = XYFcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(target_qs,evoT,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()
