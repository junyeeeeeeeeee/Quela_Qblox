from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerConti2tone
#// Okay v0.9.2

''' fill in '''
Execution:bool = 1
RO_XY_overlap:bool = False
DRandIP = {"dr":"dr1","last_ip":"11"}
freq_range:dict = {"q0":[3.35e9,3.45e9], "q1":[3.35e9,3.45e9]}    # [freq_start, freq_end] use linspace, or [0] system calculate fq for you.
xyl_range:list = [0.6]                                 # driving power [from, end, pts/step]
xyl_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 100
AVG:int = 300

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerConti2tone(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,xyl_range,xyl_sampling_func,freq_pts,AVG,RO_XY_overlap,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()