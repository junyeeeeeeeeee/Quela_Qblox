from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ZgateEnergyRelaxation
#// 0.9.2 okay

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr4","last_ip":"81"}
max_evo_time:float = 40e-6
target_qs:list = ["q0", "q1"]
time_sampling_func:str = "linspace"                  # 'linspace', 'logspace', 'arange'
time_ptsORstep:int|float = 100                     # pts/step, depends on time_sampling_func you use

z_amp_range:list = [-0.08, 0.08, 10]                # [start, end, pts/step], depends on bias_sampling_func you use
z_sampling_func:str = 'linspace'                     # 'linspace', 'arange'
histo_counts:int = None                              # If it's set a number instead of None, repeat that number no matter what value the `time_monitor` is.
time_monitor:bool = False                            # True will use while loop, then analyze it in qblox_drive_AS.analysis.raw_data_demolisher with save_dir path and a same QD_path
prepare_state:int = 0                          # 0 for |g> ; 1 for |e> 

AVG:int = 500


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ZgateEnergyRelaxation(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(max_evo_time, target_qs,time_sampling_func,z_amp_range,prepare_state,z_sampling_func,time_ptsORstep,time_monitor,AVG,Execution)
EXP.WorkFlow(histo_counts)
EXP.RunAnalysis(time_dep_plot=False)