from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ParitySwitch

#// test okay. Here only run the measurement

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr1","last_ip":"11"}
time_range:dict = {"q0":806e-9}
Shots:int = 10000
histo_counts:int = 100


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ParitySwitch(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,histo_counts,Shots,Execution,OSmode=True)
EXP.WorkFlow()

# Please go qblox_drive_AS.analysis.ParityAna to analyze the measurement results. 