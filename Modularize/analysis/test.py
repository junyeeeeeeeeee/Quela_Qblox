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
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis,T2_fit_analysis, Fit_analysis_plot,T1_fit_analysis
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from Modularize.support import rotate_onto_Inphase, rotate_data
from Modularize.support.QuFluxFit import plot_QbFlux

def plot_FreqBiasIQ(f:ndarray,z:ndarray,I:ndarray,Q:ndarray,refIQ:list=[], ax:plt.Axes=None)->plt.Axes:
    """
        I and Q shape in (z, freq)
    """
    if len(refIQ) != 2:
        refIQ = [0,0]
        
    data = sqrt((I-refIQ[0])**2+(Q-refIQ[1])**2)
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,8))
        ax:plt.Axes
    c = ax.pcolormesh(z, f*1e-9, data.transpose(), cmap='RdBu')
    ax.set_xlabel("Flux Pulse amp (V)", fontsize=20)
    ax.set_ylabel("Frequency (GHz)", fontsize=20)
    fig.colorbar(c, ax=ax, label='Contrast (V)')
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    
    return ax
# ds = open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M03S06/DR2q0_zT1(0)_H13M05S10.nc")
# time_1 = ds.attrs["end_time"]
# dss= open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M39S47/DR2q0_zT1(0)_H13M44S30.nc")
# time_2 = dss.attrs["end_time"]

# total_sec_diff = (datetime.strptime(time_2,"%Y-%m-%d %H:%M:%S")-datetime.strptime(time_1,"%Y-%m-%d %H:%M:%S")).total_seconds()# print(time.strptime(end_time,"%Y-%m-%d %H:%M:%S")-time.strptime(start_time,"%Y-%m-%d %H:%M:%S"))
# print(total_sec_diff)
file = "Modularize/Meas_raw/20241104/DR2multiQ_FluxCavity_H12M37S17.nc"
QD_file = "Modularize/QD_backup/20241102/DR2#10_SumInfo.pkl"
QD_agent = QDmanager(QD_file)
QD_agent.QD_loader()

ds = open_dataset(file)
cntrled_couplers = ds.attrs["cntrl_couplers"].replace("_", " & ")

for var in ds.data_vars:
    if var.split("_")[-1] != 'freq':
        freq = array(ds.data_vars[f"{var}_freq"])[0][0]
        bias = array(ds.coords["bias"])
        I_data = array(ds.data_vars[var])[0]
        Q_data = array(ds.data_vars[var])[1]
        ax = plot_FreqBiasIQ(freq,bias,I_data,Q_data,refIQ=[0,0])
        ax.set_title(f"{var} Readout",fontsize=20)
        ax.set_xlabel(f"couplers {cntrled_couplers} bias amplitude (V)")
        plt.grid()
        plt.tight_layout()
        plt.show()


