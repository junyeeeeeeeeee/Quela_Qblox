from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import QubitMonitor
#// 0.9.2 okay

DRandIP = {"dr":"dr4","last_ip":"81"}

QM = QubitMonitor(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),save_dir=Data_manager().build_packs_folder(),execution=True)


# T1 settings
QM.T1_time_range = {"q0":[0,40e-6],"q1":[0,40e-6],"q2":[0,40e-6]} # skip if empty.
# T2 settings
QM.T2_time_range = {"q0":[0,20e-6],"q1":[0,20e-6],"q2":[0,20e-6]} # skip if empty.
QM.echo_pi_num = [0,1,2]                               # list [0, 1, 2, ....], it will do the T2 according to the pi-number you assigned
# SingleShot settings, skip if 0 shot
QM.OS_shots = 10000     
QM.OS_target_qs = ["q0","q1","q2"]


# T1, T2 shared settings
QM.time_sampling_func = "linspace"
QM.time_ptsORstep = 100
QM.AVG = 1000

#? If use one-shot method, QM.time_ptsORstep * QM.AVG < 131000
QM.OSmode = True

""" While looping start """
QM.StartMonitoring()
