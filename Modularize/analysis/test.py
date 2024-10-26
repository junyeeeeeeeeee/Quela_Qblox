import matplotlib.pyplot as plt
import pandas as pd
from numpy import array, sqrt
import os 
import time
from  datetime import datetime
from xarray import open_dataset



# ds = open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M03S06/DR2q0_zT1(0)_H13M05S10.nc")
# time_1 = ds.attrs["end_time"]
# dss= open_dataset("Modularize/Meas_raw/20241025/ZgateT1_q0_H13M39S47/DR2q0_zT1(0)_H13M44S30.nc")
# time_2 = dss.attrs["end_time"]

# total_sec_diff = (datetime.strptime(time_2,"%Y-%m-%d %H:%M:%S")-datetime.strptime(time_1,"%Y-%m-%d %H:%M:%S")).total_seconds()# print(time.strptime(end_time,"%Y-%m-%d %H:%M:%S")-time.strptime(start_time,"%Y-%m-%d %H:%M:%S"))
# print(total_sec_diff)
file = "Modularize/Meas_raw/20241026/DR2q1_2tone_H11M52S57.nc"

def twotone_ana(nc_path:str, plot:bool=True, refIQ:dict={})->dict:
    fit_packs = {}
    ds = open_dataset(nc_path)
    joint_qs = ds.attrs["RO_qs"].split(" ")[1:] # the first element is " "

    for idx, q in enumerate(joint_qs):
        xyf = array(ds[f"x{2*idx}"])
        xyl = array(ds[f"x{2*idx+1}"],dtype=float)
        
        xyl = xyl.reshape(int(xyl.shape[0]/xyf.shape[0]),xyf.shape[0]).transpose()[0]

        ii = array(ds[f"y{2*idx}"])
        qq = array(ds[f"y{2*idx+1}"])

        if q not in refIQ:
            refIQ[q] = [0,0]
        
        contrast = sqrt((ii-refIQ[q][0])**2+(qq-refIQ[q][1])**2).reshape(xyl.shape[0],xyf.shape[0])
        
        if plot:
            if xyl.shape[0] != 1:
                plt.pcolormesh(xyf,xyl,contrast)
                plt.title(f"{q} power-dep 2tone")
                plt.ylabel("XY power (V)")
                
            else:
                plt.plot(xyf,contrast[0])
                plt.ylabel("Contrast")
                plt.title(f"{q} 2tone with XY power {xyl[0]} V")
            plt.xlabel("XY frequency (Hz)")
            plt.grid()
            plt.show()
        
        if xyl.shape[0] == 1:
            fit_packs[q] = {"xyf_data":xyf,"contrast":contrast[0]}
        
    return fit_packs