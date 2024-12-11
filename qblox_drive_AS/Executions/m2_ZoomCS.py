from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// Okay v0.9.2

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
freq_range:dict = {"q0":[4.534e9, 4.544e9],
                   "q1":[4.724e9, 4.734e9],}
                #    "q2":[6.125e9, 6.135e9],
                #    "q3":[6.223e9, 6.233e9],
                #    "q4":[6.315e9, 6.325e9],
                #    "q5":[6.414e9, 6.424e9],}
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()

