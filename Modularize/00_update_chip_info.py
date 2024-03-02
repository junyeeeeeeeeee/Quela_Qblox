# Make sure you have tkinter library

import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import Modularize.chip_data_store as cds

# ============= manually set Parameters =================
chip_name = "5Qtest0_6"
QD_path = ''
# ======================put in MeasFlow====================

chip_file = cds.Chip_file(chip_name)
chip_file.update_QD(QD_path)