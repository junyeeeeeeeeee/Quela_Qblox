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
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis,T2_fit_analysis, Fit_analysis_plot
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from Modularize.support import rotate_onto_Inphase, rotate_data
# ds = open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M03S06/DR2q0_zT1(0)_H13M05S10.nc")
# time_1 = ds.attrs["end_time"]
# dss= open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M39S47/DR2q0_zT1(0)_H13M44S30.nc")
# time_2 = dss.attrs["end_time"]

# total_sec_diff = (datetime.strptime(time_2,"%Y-%m-%d %H:%M:%S")-datetime.strptime(time_1,"%Y-%m-%d %H:%M:%S")).total_seconds()# print(time.strptime(end_time,"%Y-%m-%d %H:%M:%S")-time.strptime(start_time,"%Y-%m-%d %H:%M:%S"))
# print(total_sec_diff)
file = "Modularize/Meas_raw/20241030/DR2multiQ_T2(0)_H16M55S45.nc"
QD_file = "Modularize/QD_backup/20241030/DR2#10_SumInfo.pkl"
QD_agent = QDmanager(QD_file)
QD_agent.QD_loader()

ds = open_dataset(file)
for var in ds.data_vars:
    if var.split("_")[-1] != 'x':
        print(ds[var].)
        time_data = array(ds[f"{var}_x"])[0][0]
        reshaped = moveaxis(array(ds[var]),0,1)  # (repeat, IQ, idx)
        T2_fit = []
        for idx, data in enumerate(reshaped):
            echo:bool=False if ds[var].attrs["spin_num"] == 0 else True
            if len(QD_agent.refIQ[var]) == 1:
                data = rotate_data(data,QD_agent.refIQ[var][0])[0]
            else:
                data = sqrt((data[0]-QD_agent.refIQ[var][0])**2+(data[1]-QD_agent.refIQ[var][1])**2)
            ans = T2_fit_analysis(data,time_data)
            
            if reshaped.shape[0] == 1:
                Fit_analysis_plot(ans,P_rescale=False,Dis=None,spin_echo=echo)
            else:
                T2_fit.append(ans.attrs["T2"])


    # print(array(ds.x0).shape)
    # print(ds.y1.data)
