import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import matplotlib.pyplot as plt
import pandas as pd
from numpy import array, sqrt, linspace, diff, arange
import os 
import time
from  datetime import datetime
from xarray import open_dataset, Dataset
import quantify_core.data.handling as dh
from Modularize.support import QDmanager, Data_manager
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.analysis.raw_data_demolisher import Rabi_dataReducer
from Modularize.support.QuFluxFit import plot_QbFlux_multiVersn
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis, Fit_analysis_plot

# ds = open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M03S06/DR2q0_zT1(0)_H13M05S10.nc")
# time_1 = ds.attrs["end_time"]
# dss= open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M39S47/DR2q0_zT1(0)_H13M44S30.nc")
# time_2 = dss.attrs["end_time"]

# total_sec_diff = (datetime.strptime(time_2,"%Y-%m-%d %H:%M:%S")-datetime.strptime(time_1,"%Y-%m-%d %H:%M:%S")).total_seconds()# print(time.strptime(end_time,"%Y-%m-%d %H:%M:%S")-time.strptime(start_time,"%Y-%m-%d %H:%M:%S"))
# print(total_sec_diff)
# file = "Modularize/Meas_raw/20241028/DR2multiQ_Rabi_H15M39S17.nc"
# QD_file = "Modularize/QD_backup/20241027/DR2#10_SumInfo.pkl"
# QD_agent = QDmanager(QD_file)
# QD_agent.QD_loader()

# dss = Rabi_dataReducer(file)
# for q in dss:
#     ds = dss[q]
#     x_data = array(ds['x0'])
#     title = ds['x0'].attrs["long_name"]
#     x_axis_name = ds['x0'].attrs["name"]
#     x_axis_unit = ds['x0'].attrs["unit"]

#     i_data = array(ds['y0'])
#     q_data = array(ds['y1'])

#     contrast = sqrt((i_data-QD_agent.refIQ[q][0])**2+(q_data-QD_agent.refIQ[q][1])**2)
#     fit_pack = Rabi_fit_analysis(contrast,x_data,title)
#     Fit_analysis_plot(fit_pack,P_rescale=None,Dis=None,q=q)
    # plt.plot(x_data,contrast)
    # plt.xlabel(f"{x_axis_name} ({x_axis_unit})")
    # plt.title(f"{q}_{title}")
    # plt.ylabel("Contrast (mV)")
    # plt.show()


# for q in dss:
#     QubitFluxSpectroscopyAnalysis(tuid=dss[q].attrs["tuid"], dataset=dss[q]).run()



# for q_idx, q in enumerate(ds.attrs['RO_qs'].split(" ")[1:]):
#         x0, x1 = array(ds[f"x{2*q_idx}"]), array(ds[f"x{2*q_idx+1}"])
#         y0, y1 = array(ds[f"y{2*q_idx}"]), array(ds[f"y{2*q_idx+1}"])
#         x1 = x1.reshape(int(x1.shape[0]/x0.shape[0]),x0.shape[0]).transpose()[0]
#         y0 = y0.reshape(x0.shape[0],x1.shape[0])
#         y1 = y1.reshape(x0.shape[0],x1.shape[0])
#         amp = sqrt(y0**2+y1**2)
#         plt.pcolormesh(x1,x0,amp)
#         plt.show()

class tools_a():
    def __init__(self):
        self.a = 10
    def plus(self, number:int):
        self.a += number
        
class tools_b():
    def __init__(self):
        self.b = 7
    def minus(self, number:int):
        self.b -= number
        


class Worker(tools_a,tools_b):
    def __init__(self):
        tools_a.__init__(self)
        tools_b.__init__(self)
    
    def works(self,number:int):
        self.minus(number)
        self.plus(number)
        print(self.a)
        print(self.b)


if __name__ == "__main__":
    dad = Worker()
    dad.works(5)
