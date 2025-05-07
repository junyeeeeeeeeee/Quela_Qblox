from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerConti2tone
#// 0.9.2 okaㄗ

''' fill in '''
Execution:bool = 1 # , "q0":[4.5e9, 5e9], 
RO_XY_overlap:bool = False
DRandIP = {"dr":"dr1","last_ip":"11"}
freq_range:dict = {"q3":[3.24e9, 3.28e9]}    # [freq_start, freq_end] use linspace, or [0] system calculate fq for you.
xyl_range:list = [0.6]                                 # driving power [from, end, pts/step]
xyl_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 200
AVG:int = 300

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerConti2tone(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.fit_half_f02 = True
EXP.SetParameters(freq_range,xyl_range,xyl_sampling_func,freq_pts,AVG,RO_XY_overlap,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()