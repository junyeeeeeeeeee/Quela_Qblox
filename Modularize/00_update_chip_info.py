# Make sure you have tkinter library

import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
print(os.getcwd())

import Modularize.chip_data_store as cds

# ============= manually set Parameters =================
QD_path = ''
# ======================put in MeasFlow====================

chip_file = cds.Chip_file()
chip_file.update_QD(QD_path = QD_path)