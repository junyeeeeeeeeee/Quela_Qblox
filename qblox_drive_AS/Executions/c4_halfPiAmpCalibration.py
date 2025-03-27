from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import hPiAcali
#// 0.9.2 okay


''' fill in '''
Execution:bool = 1

DRandIP = {"dr":"dr4","last_ip":"81"}
amp_range:list = [0.9, 1.1]
target_qs:list = ["q0", "q1"]

coef_sampling_func:str = 'linspace'
half_pi_quadruple_num:list = [3,5]
coef_ptsORstep:int|float = 50
AVG:int = 300



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = hPiAcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)

EXP.SetParameters(amp_range,target_qs,coef_sampling_func,coef_ptsORstep,half_pi_quadruple_num,avg_n=AVG,execution=Execution)

EXP.WorkFlow()
EXP.RunAnalysis()