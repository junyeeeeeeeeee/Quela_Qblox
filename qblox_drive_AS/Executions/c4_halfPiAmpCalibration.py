from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import hPiAcali
#// test okay


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
ro_power_coef_range:dict = {"q0":[0.45,0.55], "q1":[0.3,0.6]}    # half pi pulse coef, it should be around 0.5
coef_sampling_func:str = 'linspace'
half_pi_quadruple_num:list = [3,5]
coef_ptsORstep:int|float = 100
AVG:int = 500



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = hPiAcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(ro_power_coef_range,coef_sampling_func,coef_ptsORstep,half_pi_quadruple_num,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()