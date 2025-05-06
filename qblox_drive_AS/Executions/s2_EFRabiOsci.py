from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.SecExcitedChara.EFRabi import PowerRabi_12_Osci
#// 0.9.2 okay
#// Initialize the EF settings in Manuall_QD_manage (QMaster.init_12_settings(switch="ON")) before you start it.

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr1","last_ip":"11"}
pi_amp_range:dict = {"q1":[-0.9, 0.9], "q3":[-0.9, 0.9]}    # [pi_amp_start, pi_amp_end]
pi_amp_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange
pi_amp_ptsORstep:int|float = 100  # Depends on the sampling func you use, 'linspace' or 'logspace' set pts in int, 'arange' set step in float
AVG:int = 1000


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerRabi_12_Osci(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(pi_amp_range,pi_amp_sampling_func,pi_amp_ptsORstep,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()