from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import DragCali
#// 0.9.2 okay


''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr4","last_ip":"81"}
drag_coef_range:list = [-0.5,0.5]
target_qs = ["q2"]
coef_sampling_func:str = 'linspace'
coef_ptsORstep:int|float = 50
AVG:int = 500


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = DragCali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(drag_coef_range,target_qs,coef_sampling_func,coef_ptsORstep,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()