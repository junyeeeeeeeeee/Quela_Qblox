from qblox_drive_AS.support.ExpFrames import ParitySwitch

QD_path = "qblox_drive_AS/QD_backup/20250220/DR1#11_SumInfo.pkl"
file_path = "/Users/ratiswu/Downloads/ParitySwitch_20250307005955.nc"

EXP = ParitySwitch(QD_path)
EXP.execution = True
EXP.RunAnalysis(new_file_path=file_path)