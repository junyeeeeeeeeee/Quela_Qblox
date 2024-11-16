from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerCavity

#TODO: TO test


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
freq_span_range:dict = {"q0":[-4e6,+5e6], "q1":[-5e6,+5e6]}    # np.linspace(rof+span, rof+span, freq_pts)
ro_amp_range:list = [1e-2, 0.6, 20]                                 # amp [from, end, pts/step]
ro_amp_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerCavity(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,ro_amp_range,ro_amp_sampling_func,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()