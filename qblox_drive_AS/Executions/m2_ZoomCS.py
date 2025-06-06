from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.Calibration_exp.TofCali import TofCalirator
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// 0.9.2 okay

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr4","last_ip":"81"}
freq_range:dict = {"q0":[6.03e9, 6.07e9],
                   "q1":[6.055e9, 6.095e9],
                   "q2":[5.97e9, 5.99e9]
                   }
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()


