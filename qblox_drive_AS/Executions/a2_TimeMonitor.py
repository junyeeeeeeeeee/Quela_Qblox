from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import QubitMonitor
#// test okay

DRandIP = {"dr":"dr1","last_ip":"11"}

QM = QubitMonitor(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),save_dir=Data_manager().build_packs_folder(),execution=True)


# T1 settings
QM.T1_time_range = {"q4":[0,30e-6],"q5":[0,30e-6]} # skip if empty.
# T2 settings
QM.T2_time_range = {"q4":[0,15e-6],"q5":[0,15e-6]} # skip if empty.
QM.echo_pi_num = [0,1,2,3]                               # list [0, 1, 2, ....], it will do the T2 according to the pi-number you assigned

# SingleShot settings, skip if 0 shot
QM.OS_shots = 10000     
QM.OS_target_qs = ["q4","q5"]


# T1, T2 shared settings

QM.time_sampling_func = "linspace"
QM.time_ptsORstep = 100
QM.AVG = 500

#? If use one-shot method, QM.time_ptsORstep * QM.AVG < 131000
QM.OSmode = True

""" While looping start """
QM.StartMonitoring()
