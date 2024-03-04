# Make sure you have tkinter library

import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import Modularize.chip_data_store as cds

# ============= manually set Parameters =================
QD_path = "D:\HW\量子元件實驗室\Qblox\Quela codes\Quela_Qblox\Modularize\QD_backup\2024_2_27\SumInfo.pkl"
# ======================put in MeasFlow====================

chip_file = cds.Chip_file()
chip_file.update_QD(QD_path)