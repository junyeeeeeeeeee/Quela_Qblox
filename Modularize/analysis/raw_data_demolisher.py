import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import open_dataset, Dataset
from numpy import ndarray, array, sqrt
import quantify_core.data.handling as dh
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.support.QDmanager import Data_manager
import matplotlib.pyplot as plt

def fluxCav_dataReductor(nc_file_path:str)->dict:
    """
    For flux cavity meas, each qubit will get 2 x-data and 2 y-data, which are 'x0', 'x1', 'y0', 'y1' accordingly.
    """
    ds = open_dataset(nc_file_path)
    ordered_q_labels = ds.attrs['RO_qs'].split(" ")[1:]
    datasets = {}
    for q_idx, q in enumerate(ordered_q_labels):
        x0, x1 = ds[f"x{2*q_idx}"], ds[f"x{2*q_idx+1}"]
        y0, y1 = ds[f"y{2*q_idx}"], ds[f"y{2*q_idx+1}"]


        new_ds = Dataset(
            data_vars = dict(y0=(["dim_0"],y0.data),y1=(["dim_0"],y1.data)),
            coords = dict(x0=(["dim_0"],x0.data),x1=(["dim_0"],x1.data))
        )
        
        new_ds.attrs = ds.attrs
        new_ds.attrs["tuid"] = new_ds.attrs["tuid"].split("-")[0]+"-"+new_ds.attrs["tuid"].split("-")[1]+"-"+f'{(int(new_ds.attrs["tuid"].split("-")[2])+q_idx):03d}'+"-"+new_ds.attrs["tuid"].split("-")[3]
        
        Data_manager().build_tuid_folder(new_ds.attrs["tuid"],f"{ordered_q_labels[q_idx]}FluxCav")
        to_copy_array_attr = [x0, x1, y0, y1]
        for idx, item in enumerate([new_ds.x0, new_ds.x1, new_ds.y0, new_ds.y1]):
            item.attrs = to_copy_array_attr[idx].attrs

        datasets[ordered_q_labels[q_idx]] = new_ds

    return datasets

def fluxQub_dataReductor(nc_file_path:str)->dict:
    ds = open_dataset(nc_file_path)
    ordered_q_labels = ds.attrs['RO_qs'].split(" ")[1:]
    datasets = {}
    for q_idx, q in enumerate(ordered_q_labels):
        x0, x1 = ds[f"x{2*q_idx}"], ds[f"x{2*q_idx+1}"]
        y0, y1 = ds[f"y{2*q_idx}"], ds[f"y{2*q_idx+1}"]

        print(array(x0).shape)
        print(array(x1).shape)

        new_ds = Dataset(
            data_vars = dict(y0=(["dim_0"],y0.data),y1=(["dim_0"],y1.data)),
            coords = dict(x0=(["dim_0"],x0.data),x1=(["dim_0"],x1.data))
        )
        
        new_ds.attrs = ds.attrs
        new_ds.attrs["tuid"] = new_ds.attrs["tuid"].split("-")[0]+"-"+new_ds.attrs["tuid"].split("-")[1]+"-"+f'{(int(new_ds.attrs["tuid"].split("-")[2])+q_idx):03d}'+"-"+new_ds.attrs["tuid"].split("-")[3]
        
        Data_manager().build_tuid_folder(new_ds.attrs["tuid"],f"{ordered_q_labels[q_idx]}FluxQubit")
        to_copy_array_attr = [x0, x1, y0, y1]
        for idx, item in enumerate([new_ds.x0, new_ds.x1, new_ds.y0, new_ds.y1]):
            item.attrs = to_copy_array_attr[idx].attrs

        datasets[ordered_q_labels[q_idx]] = new_ds

    return datasets

def twotone_ana(nc_path:str, plot:bool=True, refIQ:dict={} , fit_func:callable=None)->dict:
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
        fit_f01s = []
        fif_amps = []
        if fit_func is not None and xyl.shape[0] != 1:
            for xyl_idx, an_amp_data in enumerate(contrast):
                if xyl[xyl_idx] != 0:
                    if min(xyf) <= fit_func(an_amp_data,xyf).attrs['f01_fit'] and fit_func(an_amp_data,xyf).attrs['f01_fit'] <= max(xyf):
                        fit_f01s.append(fit_func(an_amp_data,xyf).attrs['f01_fit']) # Hz
                        fif_amps.append(xyl[xyl_idx])

        if plot:
            if xyl.shape[0] != 1:
                plt.pcolormesh(xyf,xyl,contrast)
                if len(fit_f01s) != 0 :
                    plt.scatter(fit_f01s,fif_amps,marker="*",c='red')
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


if __name__ == "__main__":
    dh.set_datadir(meas_raw_dir)
    # file = "Modularize/Meas_raw/20241024/DR2q1_FluxCavity_H14M35S23.nc"
    file = "Modularize/Meas_raw/2024_10_16/DR4q4_FluxCavity_H13M31S11.nc"
    dss = fluxCav_dataReductor(file, ['q0'])
    # for q in dss:
    #     ResonatorFluxSpectroscopyAnalysis(tuid=dss[q].attrs["tuid"], dataset=dss[q]).run(sweetspot_index=0)
    