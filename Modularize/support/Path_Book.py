"""
    Path_book keeps the path indecated to the directories, and the object will be saved in the folder named by date.
"""
import os, datetime
# The directory for measurement raw data
meas_raw_dir = 'Modularize/Meas_raw'
# The directory for qauntum device
qdevice_backup_dir = 'Modularize/QD_backup'

def find_latest_QD_folder(which_dr:str):
    folders = [name for name in os.listdir(qdevice_backup_dir) if os.path.isdir(os.path.join(qdevice_backup_dir,name))]
    date = []
    for name in folders:
        year_str = name.split("_")[0]
        mon_str  = "0"+name.split("_")[1] if len(name.split("_")[1])==1 else name.split("_")[1]
        day_str  = "0"+name.split("_")[2] if len(name.split("_")[2])==1 else name.split("_")[2]
        date.append(datetime.datetime.strptime(day_str+mon_str+year_str, "%d%m%Y").date())


if __name__ == "__main__":
    x = ["24","3","1"]
    print(x[0]+x[1]+x[2])
