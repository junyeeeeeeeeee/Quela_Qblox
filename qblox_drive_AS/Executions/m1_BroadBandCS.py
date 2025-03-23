from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import BroadBand_CavitySearching
#// 0.9.2 okay

""" Fill in """
DRandIP = {"dr":"dr4","last_ip":"81"}
target_q = "q"
freq_sample_linspaced = [5.9e9, 6.2e9, 600]
save_dir = Data_manager().build_packs_folder()


""" Don't touch """
EXP = BroadBand_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_sample_linspaced[0],freq_sample_linspaced[1],freq_sample_linspaced[2], target_q)
EXP.WorkFlow()
EXP.RunAnalysis()

