from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerRabiOsci
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
pi_amp_range:dict = {"q0":[-0.6,0.6], "q1":[-0.5,0.5]}    # [pi_amp_start, pi_amp_end]
pi_amp_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

pi_dura:dict = {"q0":40e-9, "q1":46e-9}

pi_amp_ptsORstep:int|float = 100  # Depends on the sampling func you use, 'linspace' or 'logspace' set pts in int, 'arange' set step in float
AVG:int = 300

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerRabiOsci(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(pi_amp_range,pi_dura,pi_amp_sampling_func,pi_amp_ptsORstep,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()