from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerRabiOsci
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
pi_amp_range:dict = {"q4":[-0.9,0.9], "q5":[-0.9,0.9]}    # [pi_amp_start, pi_amp_end]
pi_amp_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange
pi_amp_ptsORstep:int|float = 100  # Depends on the sampling func you use, 'linspace' or 'logspace' set pts in int, 'arange' set step in float
AVG:int = 600

#?? Notes: While you're using one-shot method, amp samples (pi_amp_pts) * shots (AVG) must less than or equal to 131000. 
use_OneShot:bool = False

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerRabiOsci(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(pi_amp_range,pi_amp_sampling_func,pi_amp_ptsORstep,AVG,Execution,OSmode=use_OneShot)
EXP.WorkFlow()
EXP.RunAnalysis()