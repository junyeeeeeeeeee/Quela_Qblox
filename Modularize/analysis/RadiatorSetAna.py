import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis, T2_fit_analysis
from Modularize.support.Path_Book import meas_raw_dir
from xarray import Dataset, open_dataset # Dataset.to_dict(SS_ds)
from numpy import array
import matplotlib.pyplot as plt


sub_folders = [name for name in os.listdir(meas_raw_dir) if (os.path.isdir(os.path.join(meas_raw_dir,name)) and name[:8]='Radiator')]

