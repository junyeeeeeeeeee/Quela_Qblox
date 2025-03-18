from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import CPMG
#// 0.9.2 okay


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
time_range:dict = {"q0":[0,100e-6],"q1":[0,100e-6]}
pi_num:int = 2
time_sampling_func:str = "linspace"
time_ptsORstep:int|float = 100
AVG:int = 1000
histo_counts:int = 1

#?? Notes: While you're using one-shot method, time samples (time_pts) * shots (AVG) must less than or equal to 131000. 
use_OneShot:bool = True

''' Don't Touch '''
pi_num_dict = {}
for q in time_range:
    pi_num_dict[q] = pi_num
save_dir = Data_manager().build_packs_folder()
EXP = CPMG(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,pi_num_dict,time_sampling_func,time_ptsORstep,histo_counts,AVG,Execution,OSmode=use_OneShot)
EXP.WorkFlow()
EXP.RunAnalysis()
