from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.Calibration_exp.TofCali import TofCalirator
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
freq_range:dict = {"q1":[4.9e9, 4.92e9],
                   "q3":[5.1e9, 5.12e9],
                   }
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()


### Initialization measurement to get `time_of_flight` and `nco_prop_delay`
CAL = TofCalirator(save_dir)
CAL.QD_path = find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"])
CAL.q = list(freq_range.keys())[0]
CAL.WorkFlow()
CAL.RunAnalysis()
