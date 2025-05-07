from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.SecExcitedChara.GEF_ROfreq_cali import GEF_ROFcali
#// 0.9.2 okay

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
freq_span_range:dict = {"q1":[-1e6,+1e6],"q3":[-1e6,1e6]}
freq_pts = 100
AVG:int = 500



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = GEF_ROFcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,freq_pts,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()
