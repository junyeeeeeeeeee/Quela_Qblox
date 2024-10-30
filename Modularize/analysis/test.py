import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import matplotlib.pyplot as plt
import pandas as pd
from numpy import array, sqrt, linspace, diff, arange, arctan2, sin, cos, hstack, vstack, ndarray, pi, column_stack, moveaxis, empty_like
import os 
import time
from  datetime import datetime
from xarray import open_dataset, Dataset, DataArray
import quantify_core.data.handling as dh
from Modularize.support import QDmanager, Data_manager
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.analysis.raw_data_demolisher import Rabi_dataReducer
from Modularize.support.QuFluxFit import plot_QbFlux_multiVersn
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis, Fit_analysis_plot
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from Modularize.support import rotate_onto_Inphase, rotate_data
# ds = open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M03S06/DR2q0_zT1(0)_H13M05S10.nc")
# time_1 = ds.attrs["end_time"]
# dss= open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M39S47/DR2q0_zT1(0)_H13M44S30.nc")
# time_2 = dss.attrs["end_time"]

# total_sec_diff = (datetime.strptime(time_2,"%Y-%m-%d %H:%M:%S")-datetime.strptime(time_1,"%Y-%m-%d %H:%M:%S")).total_seconds()# print(time.strptime(end_time,"%Y-%m-%d %H:%M:%S")-time.strptime(start_time,"%Y-%m-%d %H:%M:%S"))
# print(total_sec_diff)
# file = "Modularize/Meas_raw/20241029/DR2multiQ_SingleShot(0)_H18M40S57.nc"
# QD_file = "Modularize/QD_backup/20241029/DR2#10_SumInfo.pkl"
# QD_agent = QDmanager(QD_file)
# QD_agent.QD_loader()

# ds = open_dataset(file)
# for var in ds.data_vars:
#     ds[var] = ds[var] * 1000
# data = ds['q0']
# gmm2d_fidelity = GMMROFidelity()
# gmm2d_fidelity._import_data(data)
# gmm2d_fidelity._start_analysis()
# g1d_fidelity = gmm2d_fidelity.export_G1DROFidelity()

# transi_freq = None


# _, angle = rotate_onto_Inphase(gmm2d_fidelity.centers[0],gmm2d_fidelity.centers[1])





# z = moveaxis(array(data),0,1) # (IQ, state, shots) -> (state, IQ, shots)
# rotated_data = empty_like(array(data))
# for state_idx, state_data in enumerate(z):
#     rotated_data[state_idx] = rotate_data(state_data,angle)

# da = DataArray(moveaxis(rotated_data,0,1), coords= [("mixer",["I","Q"]), ("prepared_state",[0,1]), ("index",arange(array(data).shape[2]))] )
# gmm2d_fidelity._import_data(da)
# gmm2d_fidelity._start_analysis()
# g1d_fidelity = gmm2d_fidelity.export_G1DROFidelity()
# print(f"!! center:\n{gmm2d_fidelity.centers}")
# plot_readout_fidelity(da, gmm2d_fidelity, g1d_fidelity,transi_freq,None)
# plt.close()

name = "DR2q0_SingleShot(0)_H11M01S48"
print(name.split("_")[1].split("(")[1][:-1])