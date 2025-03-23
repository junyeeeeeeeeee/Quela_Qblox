from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import FluxCoupler
#// 0.9.2 okay.

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr4","last_ip":"81"}
freq_span_range:dict = {"q0":[-5e6,+5e6]}    # np.linspace(rof+span, rof+span, freq_pts)
bias_elements:list = ["c0"]
flux_range:list = [-0.3, 0.3, 20]                                 # flux [from, end, pts/step]
flux_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 30
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = FluxCoupler(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,bias_elements,flux_range,flux_sampling_func,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()