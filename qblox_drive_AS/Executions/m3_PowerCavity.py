from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import PowerCavity


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr1","last_ip":"11"}
freq_span_range:dict = {"q3":[-1e6,+1e6], "q1":[-1e6,+1e6], "q0":[-1e6,+1e6], "q2":[-1e6,+1e6]}    # np.linspace(rof+span, rof+span, freq_pts)
RO_atteuations:list = [0, 60, 2]                               # atte [from, end, step] , end MUST be less than 60, end step MUST be even.  

freq_pts:int = 50
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = PowerCavity(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,RO_atteuations,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()

""" Tip """
""" If fitting for dress freq is acceptible, go set a right attenuation manually and then SKIP m4. """