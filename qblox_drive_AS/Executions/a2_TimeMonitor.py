from datetime import datetime
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager, UserFriend
from qblox_drive_AS.support.ExpFrames import CPMG, EnergyRelaxation, SingleShot
#TODO: need a test


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"dr2","last_ip":"10"}
# T1 settings
T1_time_range:dict = {"q0":[0,60e-6],"q1":[0,50e-6]} # skip if empty.
# T2 settings
T2_time_range:dict = {"q0":[0,60e-6],"q1":[0,50e-6]} # skip if empty.
echo_pi_num:dict = {"q0":2,"q1":3}
a_little_detune_Hz:float = 0.3e6                     # for all qubit in T2_time_range. Set None if you do SpinEcho or CPMG
# SingleShot settings, skip if 0 shot
Shots:int = 10000     
target_qs:list = ["q0","q1"]


# T1, T2 shared settings
time_sampling_func:str = "linspace"
time_ptsORstep:int|float = 100
AVG:int = 500


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
start_time = datetime.now()
while True:
    if T1_time_range is not None or T1_time_range != {}:
        EXP = EnergyRelaxation(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
        EXP.SetParameters(T1_time_range,time_sampling_func,time_ptsORstep,1,AVG,Execution)
        EXP.WorkFlow()

    if T2_time_range is not None or T2_time_range != {}:
        EXP = CPMG(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
        EXP.SetParameters(T2_time_range,echo_pi_num,time_sampling_func,time_ptsORstep,1,AVG,Execution)
        EXP.WorkFlow(freq_detune_Hz=a_little_detune_Hz)

    if target_qs is not None or len(target_qs) != 0 or Shots != 0:
        EXP = SingleShot(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
        EXP.SetParameters(target_qs,1,Shots,Execution)
        EXP.WorkFlow()

    end_time = datetime.now()
    time_difference = end_time - start_time
    UserFriend.slightly_print(f"It's {round(time_difference.total_seconds()/3600,2)} hrs recorded.")

# Analyze it in qblox_drive_AS.analysis.TimeTraceAna.py with your save_dir and a same QD_path.