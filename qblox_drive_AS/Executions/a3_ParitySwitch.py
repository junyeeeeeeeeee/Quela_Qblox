from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ParitySwitch
from qblox_drive_AS.support.QDmanager import QDmanager
#// test okay

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr1","last_ip":"11"}
time_range:dict = {"q0":806e-9}
AVG:int = 10000
histo_counts:int = 100


''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ParitySwitch(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,histo_counts,AVG,Execution,OSmode=True)
# EXP.WorkFlow()
EXP.RunAnalysis(
    new_QD_path="/Users/machenhsun/Documents/GitHub/OneShot_T1/qblox_drive_AS/QD_backup/20250307/DR1#11_SumInfo.pkl",
    new_file_path="/Users/machenhsun/Documents/GitHub/OneShot_T1/qblox_drive_AS/Meas_raw/20250307/H00M58S42/ParitySwitch_20250307005955.nc",  # 確保是 .nc 檔案
    new_pic_save_place="/Users/machenhsun/Documents/GitHub/OneShot_T1/qblox_drive_AS/Meas_raw/20250307/H00M58S42"
)