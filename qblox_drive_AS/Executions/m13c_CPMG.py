from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import CPMG
#T// test okay


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
time_range:dict = {"q0":[0,60e-6],"q1":[0,50e-6]}
pi_num:int = 2
time_sampling_func:str = "linspace"
time_ptsORstep:int|float = 100
AVG:int = 500
histo_counts:int = 3

''' Don't Touch '''
pi_num_dict = {}
for q in time_range:
    pi_num_dict[q] = pi_num
save_dir = Data_manager().build_packs_folder()
EXP = CPMG(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,pi_num_dict,time_sampling_func,time_ptsORstep,histo_counts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()
