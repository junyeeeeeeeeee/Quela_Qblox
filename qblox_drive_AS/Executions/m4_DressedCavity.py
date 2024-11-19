""" Go qblox_drive_AS.Configs.Manuall_QG_manage.py set your dressed state readout attenuation first """
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import Dressed_CavitySearching
#// test okay.


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
freq_range:dict = {"q0":[5.9925e9, 5.9975e9], "q1":[6.075e9, 6.080e9]}    # np.linspace(rof+span, rof+span, freq_pts)
ro_amp = {"q0":0.2, "q1":0.2}

freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Dressed_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,ro_amp,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()