from qblox_drive_AS.Calibration_exp.SQRB import SQRB
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// single readout okay

save_dir = Data_manager().build_packs_folder()
EXP = SQRB(data_folder=save_dir)

# Execution, bool
EXP.execution = 1

# QD_path setting, str
EXP.QD_path = find_latest_QD_pkl_for_dr(which_dr="dr4", ip_label="81")

# measurement qubit
EXP.q = "q1"

# maximum gate number
EXP.max_gate_num = 300
# gate sampliing number
EXP.gate_pts = 100
# AVG times
EXP.avg_n = 1000



""" Do NOT touch !"""
EXP.SetParameters()
EXP.WorkFlow()
EXP.RunAnalysis()