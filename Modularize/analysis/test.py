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
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis,T2_fit_analysis, Fit_analysis_plot,T1_fit_analysis, IQ_data_dis, cos_fit_analysis
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from Modularize.support import rotate_onto_Inphase, rotate_data
from Modularize.support.QuFluxFit import plot_QbFlux




file = "Modularize/Meas_raw/20241107/DR2multiQ_HalfPiCali(0)_H12M39S40.nc"
QD_file = "Modularize/QD_backup/20241107/DR2#10_SumInfo.pkl"
QD_agent = QDmanager(QD_file)
QD_agent.QD_loader()

refIQ = QD_agent.refIQ
""" for pi-amp """
# ds = open_dataset(file)
# for var in ds.data_vars:
#     if var.split("_")[-1] != "PIcoef":
#         pi_amp_coef =  moveaxis(array(ds[f"{var}_PIcoef"]),1,0)[0][0]
#         pi_pair_num = array(ds.coords["PiPairNum"])
#         data = moveaxis(array(ds[var]),1,0)
#         refined_data_folder = []
#         for PiPairNum_dep_data in data:
#             if len(refIQ[var]) == 2:
#                 refined_data = IQ_data_dis(PiPairNum_dep_data[0],PiPairNum_dep_data[1],refIQ[var][0],refIQ[var][1])
#             else:
#                 refined_data = rotate_data(PiPairNum_dep_data,refIQ[var][0])[0]
#             refined_data_folder.append(cos_fit_analysis(refined_data,pi_amp_coef))
#         fig, ax = plt.subplots()
#         for idx, refined_data in enumerate(refined_data_folder):
#             x = refined_data.coords['freeDu']
#             x_fit = refined_data.coords['para_fit']  
#             ax.plot(x,refined_data.data_vars['data'],'--',label=f"{pi_pair_num[idx]} PI pairs", alpha=0.8, ms=4)
#             ax.plot(x_fit,refined_data.data_vars['fitting'],'-', alpha=1, lw=2)    
            
#         #     plt.plot(pi_amp_coef,refined_data,label=f"{pi_pair_num[idx]} pairs pi-pulses")
#         plt.title(f"{var} PI-pulse amp coef calibration")
#         plt.legend()
#         plt.grid()
#         plt.show()

""" for half pi """
ds = open_dataset(file)
for var in ds.data_vars:
    if var.split("_")[-1] != "HalfPIcoef":
        pi_amp_coef =  moveaxis(array(ds[f"{var}_HalfPIcoef"]),1,0)[0][0]
        pi_pair_num = array(ds.coords["PiPairNum"])
        data = moveaxis(array(ds[var]),1,0)
        refined_data_folder = []
        for PiPairNum_dep_data in data:
            if len(refIQ[var]) == 2:
                refined_data = IQ_data_dis(PiPairNum_dep_data[0],PiPairNum_dep_data[1],refIQ[var][0],refIQ[var][1])
            else:
                refined_data = rotate_data(PiPairNum_dep_data,refIQ[var][0])[0]
            refined_data_folder.append(cos_fit_analysis(refined_data,pi_amp_coef))
        fig, ax = plt.subplots()
        for idx, refined_data in enumerate(refined_data_folder):
            x = refined_data.coords['freeDu']
            x_fit = refined_data.coords['para_fit']  
            ax.plot(x,refined_data.data_vars['data'],'--',label=f"{pi_pair_num[idx]} halfPI quadruples", alpha=0.8, ms=4)
            ax.plot(x_fit,refined_data.data_vars['fitting'],'-', alpha=1, lw=2)    
            
        #     plt.plot(pi_amp_coef,refined_data,label=f"{pi_pair_num[idx]} pairs pi-pulses")
        plt.title(f"{var} PI-pulse amp coef calibration")
        plt.legend()
        plt.grid()
        plt.show()
