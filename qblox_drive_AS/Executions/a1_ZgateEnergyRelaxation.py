from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ZgateEnergyRelaxation
#// test okay

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
time_range:dict = {"q0":[0,80e-6],"q1":[0,70e-6]}    # {"q0":[start, end], ...}
time_sampling_func:str = "linspace"                  # 'linspace', 'logspace', 'arange'
time_ptsORstep:int|float = 75                       # pts/step, depends on time_sampling_func you use

z_amp_range:list = [-0.08, 0.08, 10]                # [start, end, pts/step], depends on bias_sampling_func you use
z_sampling_func:str = 'linspace'                     # 'linspace', 'arange'

time_monitor:bool = False                            # True will use while loop, then analyze it in qblox_drive_AS.analysis.raw_data_demolisher with save_dir path and a same QD_path
prepare_excited:bool = True                          # False will not drive qubit

AVG:int = 500


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ZgateEnergyRelaxation(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,time_sampling_func,z_amp_range,prepare_excited,z_sampling_func,time_ptsORstep,time_monitor,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis(time_dep_plot=False)