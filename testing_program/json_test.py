import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from Modularize import chip_data_store as cds

# ============= manually set Parameters =================
chip_name = "5Qtest6"
chip_type = "5Q"
QD_path = ''
ro_atte = 20
xy_atte = 0

# ======================put in MeasFlow====================

chip_file = cds.Chip_file(chip_name, chip_type, QD_path, ro_atte, xy_atte)
print(chip_file.get_chip_dict())