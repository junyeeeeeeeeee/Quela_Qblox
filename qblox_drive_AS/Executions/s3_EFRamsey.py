from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.SecExcitedChara.EFRamsey import EF_Ramsey, f12_calibration
#// 0.9.2 okay

''' fill in '''
Execution:bool = 1
Calibrate_mode:bool = 0    # Calibrate f12 or not
DRandIP = {"dr":"dr1","last_ip":"11"}
max_evo_time:float = 10e-6
target_qs:list = ["q3"]
time_sampling_func:str = "linspace"
time_ptsORstep:int|float = 100
AVG:int = 1000
histo_counts:int = 1



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
if not Calibrate_mode:
    EXP = EF_Ramsey(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
else:
    EXP = f12_calibration(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)

EXP.SetParameters(max_evo_time, target_qs,time_sampling_func,time_ptsORstep,histo_counts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()
