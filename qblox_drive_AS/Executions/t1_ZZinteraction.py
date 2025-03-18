from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from numpy import linspace
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.TQoperations.ZZinteraction import ZZinteraction

save_dir = Data_manager().build_packs_folder()
EXP = ZZinteraction(data_folder=save_dir)

# Execution, bool
EXP.execution = False
# QD_path setting, str
EXP.QD_path = find_latest_QD_pkl_for_dr(which_dr="dr2", ip_label="10")
# method oneshot or not, bool
# AVG or Shot number, int
EXP.OSmode = False
EXP.avg_n = 500
# Evolution time, ndarray
EXP.time_samples = linspace(0, 100e-6, 100)
# Bias samples, ndarray
EXP.z_amp_samples = linspace(-0.05, 0.05, 20)
# Readout qubits, list
EXP.read_qs = ["q0"]
# Probe pi-pulse qubits, list
EXP.probe_qs = ["q1"]
# Biased couplers, list
EXP.bias_cs = ['q0']
# histogram counts, int
EXP.histo_count = 1

EXP.SetParameters()
EXP.WorkFlow()
EXP.RunAnalysis()

