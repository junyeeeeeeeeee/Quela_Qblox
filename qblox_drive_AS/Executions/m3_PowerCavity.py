from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerCavity
# // test okay.


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr4","last_ip":"81"}
freq_span_range:dict = {"q0":[-6e6,+6e6], "q1":[-6e6,+6e6],
                        "q2":[-6e6,+6e6], "q3":[-6e6,+6e6],
                        "q4":[-6e6,+6e6], "q5":[-6e6,+6e6],}    # np.linspace(rof+span, rof+span, freq_pts)
ro_amp_range:list = [0, 0.6, 40]                                 # amp [from, end, pts/step]
ro_amp_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 40
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerCavity(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,ro_amp_range,ro_amp_sampling_func,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()