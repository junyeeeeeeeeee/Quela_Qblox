from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import BroadBand_CavitySearching


""" Fill in """
DRandIP = {"dr":"dr2","last_ip":"10"}
target_qs = ["q0"]
freq_sample_linspaced = [5.8e9, 6.3e9, 1000]
save_dir = Data_manager().build_packs_folder()


""" Don't touch """
EXP = BroadBand_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(target_qs,freq_sample_linspaced[0],freq_sample_linspaced[1],freq_sample_linspaced[2])
EXP.WorkFlow()
EXP.RunAnalysis()

