from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import EnergyRelaxation


''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr1","last_ip":"11"}
max_evo_time:float = 1000e-6
target_qs:list = ["q1","q3"]
time_sampling_func:str = "linspace"
time_ptsORstep:int|float = 100
AVG:int = 400
histo_counts:int = 1

#?? Notes: While you're using one-shot method, time samples (time_pts) * shots (AVG) must less than or equal to 131000. 
use_OneShot:bool = 0

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = EnergyRelaxation(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(max_evo_time,target_qs,time_sampling_func,time_ptsORstep,histo_counts,AVG,Execution,use_OneShot)
EXP.WorkFlow()
EXP.RunAnalysis()