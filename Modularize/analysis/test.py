import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import matplotlib.pyplot as plt
import pandas as pd
from numpy import array, sqrt, linspace, diff, arange, arctan2, sin, cos, hstack, vstack, ndarray, pi, column_stack, moveaxis, empty_like, where, mean
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



def rof_cali_anaplot(ds:Dataset):
    for var in ds.data_vars:
        if var.split("_")[-1] != "rof":
            IQ_data = moveaxis(array(ds[var]),1,0) # shape (mixer, state, rof) -> (state, mixer, rof)
            rof = moveaxis(array(ds[f"{var}_rof"]),0,1)[0][0]
            I_diff = IQ_data[1][0]-IQ_data[0][0]
            Q_diff = IQ_data[1][1]-IQ_data[0][1]
            dis_diff = sqrt((I_diff)**2+(Q_diff)**2)
            # optimal_rof = anal_rof_cali(IQ_data[1][0],IQ_data[1][1],IQ_data[0][0],IQ_data[0][1],dis_diff,rof,mean(rof))
              



file = "Modularize/Meas_raw/20241106/DR2multiQ_RofCali(0)_H14M43S18.nc"
QD_file = "Modularize/QD_backup/20241106/DR2#10_SumInfo.pkl"
QD_agent = QDmanager(QD_file)
QD_agent.QD_loader()

ds = open_dataset(file)
rof_cali_anaplot(ds,QD_agent.refIQ)