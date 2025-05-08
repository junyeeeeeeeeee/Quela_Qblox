from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import QubitMonitor
#// 0.9.2 okay

DRandIP = {"dr":"dr1","last_ip":"11"}

QM = QubitMonitor(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),save_dir=Data_manager().build_packs_folder(),execution=True)


# T1 settings
QM.T1_target_qs = ["q1", "q3"]   # skip T1 exp if empty. 
QM.T1_max_evo_time = 1000e-6

# T2 settings
QM.T2_target_qs = ["q1", "q3"]   # skip T1 exp if empty. 
QM.T2_max_evo_time = 100e-6
QM.echo_pi_num = [0,1]                               # list [0, 1, 2, ....], it will do the T2 according to the pi-number you assigned

# SingleShot settings, skip if 0 shot
QM.OS_shots = 10000     
QM.OS_target_qs = ["q1","q3"]

# T1, T2 shared settings
QM.time_sampling_func = "linspace"
QM.time_ptsORstep = 100
QM.AVG = 500

#? If use one-shot method, QM.time_ptsORstep * QM.AVG < 131000
QM.OSmode = False


""" While looping start """
QM.StartMonitoring()
