from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ParitySwitch

#// 0.9.2 okay. Here only run the measurement

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr2","last_ip":"10"}
time_range:dict = {"q0":806e-9,"q1":806e-9}
Shots:int = 10000
histo_counts:int = 1


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ParitySwitch(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,histo_counts,Shots,Execution,OSmode=True)
EXP.WorkFlow()

# Please go qblox_drive_AS.analysis.ParityAna to analyze the measurement results. 