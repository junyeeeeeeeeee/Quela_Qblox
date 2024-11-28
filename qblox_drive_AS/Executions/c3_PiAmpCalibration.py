from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PiAcali
#// test okay


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
pi_power_coef_range:dict = {"q4":[0.95,1.05]}
coef_sampling_func:str = 'linspace'
pi_pair_num:list = [12,15]
coef_ptsORstep:int|float = 100
AVG:int = 500



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PiAcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(pi_power_coef_range,coef_sampling_func,coef_ptsORstep,pi_pair_num,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()