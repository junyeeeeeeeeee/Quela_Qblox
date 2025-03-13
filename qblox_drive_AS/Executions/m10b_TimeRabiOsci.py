""" This exp doesn't update any info for QDmanager """
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import TimeRabiOsci
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
pi_dura_range:dict = {"q4":[0,100e-9], "q5":[0,100e-9]}    # [pi_amp_start, pi_amp_end]
pi_dura_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange
pi_dura_ptsORstep:int|float = 100  # Depends on the sampling func you use, 'linspace' or 'logspace' set pts in int, 'arange' set step in float
AVG:int = 500

#?? Notes: While you're using one-shot method, time samples (pi_dura_pts) * shots (AVG) must less than or equal to 131000.
use_OneShot:bool = False

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = TimeRabiOsci(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(pi_dura_range,pi_dura_sampling_func,pi_dura_ptsORstep,AVG,Execution,OSmode=use_OneShot)
EXP.WorkFlow()
EXP.RunAnalysis()