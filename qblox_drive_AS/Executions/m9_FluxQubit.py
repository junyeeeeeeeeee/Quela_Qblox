from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import FluxQubit
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
freq_span_range:dict = {"q0":[-350e6,100e6], "q1":[-300e6,150e6]}    # [freq_span_start, freq_span_end] use linspace, total span should <= 500 MHz
bias_elements:list = ['q0', 'q1']
z_amp_range:list = [-0.08, 0.08, 60]                                 # z-pulse amplitude [from, end, pts/step]
z_amp_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 40
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = FluxQubit(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,bias_elements,z_amp_range,z_amp_sampling_func,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()