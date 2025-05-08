from qblox_drive_AS.Calibration_exp.SQRB import SQRB
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
save_dir = Data_manager().build_packs_folder()
EXP = SQRB(data_folder=save_dir)


# Execution, bool
EXP.execution = 0

# QD_path setting, str
EXP.QD_path = find_latest_QD_pkl_for_dr(which_dr="dr1", ip_label="11")

# measurement qubit
EXP.qs = ["q1", "q3"]

# maximum gate number
EXP.max_gate_num = 100  # 300 is maximum, max usage of sequencer memory  
# gate sampliing number
EXP.gate_pts = 75
# how many random circuits
EXP.circuits_num = 2

# AVG times
EXP.avg_n = 500



""" Do NOT touch !"""
EXP.SetParameters()
EXP.WorkFlow()
EXP.RunAnalysis()