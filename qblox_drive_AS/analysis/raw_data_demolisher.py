import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import open_dataset, Dataset
from numpy import ndarray, array, sqrt, cos, sin, transpose, deg2rad, real, imag
import quantify_core.data.handling as dh
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from qblox_drive_AS.support.Path_Book import meas_raw_dir
from qblox_drive_AS.support.QDmanager import Data_manager
import matplotlib.pyplot as plt
import re

class MultiplexingDataReducer():
    def __init__(self, nc:str|Dataset):
        if isinstance(nc,str):
            self.ds = open_dataset(nc)
        elif isinstance(nc,Dataset):
            self.ds = nc
        else:
            raise TypeError("You should give a nc_file path or a xr.Dataset")
        
def fluxCoupler_dataReducer(nc_file_path:str)->Dataset:
    dataset = open_dataset(nc_file_path)

    return dataset


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

        new_ds = Dataset(
            data_vars = dict(y0=(["dim_0"],y0.data),y1=(["dim_0"],y1.data)),
            coords = dict(x0=(["dim_0"],x0.data),x1=(["dim_0"],x1.data))
        )
        
        new_ds.attrs = ds.attrs
        new_ds.attrs["target_q"] = q
        new_ds.attrs["tuid"] = new_ds.attrs["tuid"].split("-")[0]+"-"+new_ds.attrs["tuid"].split("-")[1]+"-"+f'{(int(new_ds.attrs["tuid"].split("-")[2])+q_idx):03d}'+"-"+new_ds.attrs["tuid"].split("-")[3]
        
        # Data_manager().build_tuid_folder(new_ds.attrs["tuid"],f"{ordered_q_labels[q_idx]}FluxQubit")
        to_copy_array_attr = [x0, x1, y0, y1]
        for idx, item in enumerate([new_ds.x0, new_ds.x1, new_ds.y0, new_ds.y1]):
            item.attrs = to_copy_array_attr[idx].attrs

        datasets[ordered_q_labels[q_idx]] = new_ds

    return datasets

def Rabi_dataReducer(nc_file_path:str):
    ds = open_dataset(nc_file_path)
    
    joint_q = [x for x in ds.attrs["RO_qs"].split(" ") if x.startswith('q')] 
    dataset = {}
    for q_idx, q in enumerate(joint_q):
        attr = ds.attrs
        x0 = ds[f"x{q_idx}"]
        y0 = ds[f"y{2*q_idx}"]
        y1 = ds[f"y{2*q_idx+1}"]

        new_ds = Dataset(
            data_vars = dict(y0=(["Mixer_I"],y0.data),y1=(["Mixer_Q"],y1.data)),
            coords = dict(x0=(["Varable_1"],x0.data))
        )
        new_ds.attrs = attr
        new_ds.attrs['target_q'] = q
        to_copy_array_attr = [x0, y0, y1]
        for idx, item in enumerate([new_ds.x0, new_ds.y0, new_ds.y1]):
            item.attrs = to_copy_array_attr[idx].attrs
        
        dataset[q] = new_ds
    
    return dataset

def OneShot_dataReducer(nc_file_path:str):
    dataset = open_dataset(nc_file_path)
    for var in dataset.data_vars:
        dataset[var] = dataset[var] * 1000
    
    return dataset

def T2_dataReducer(nc_file_path:str):
    dataset = open_dataset(nc_file_path)

    return dataset

def T1_dataReducer(nc_file_path:str):
    dataset = open_dataset(nc_file_path)

    return dataset

def ZgateT1_dataReducer(raw_data_folder:str)->dict:
    
    datasets = []
    # Iterate directory
    for path in os.listdir(raw_data_folder):
        # check if current path is a file
        if os.path.isfile(os.path.join(raw_data_folder, path)):
            file_path = os.path.join(raw_data_folder,path)
            if file_path.split(".")[-1] == 'nc':
                datasets.append(open_dataset(file_path))
    # make VIP folder for each qubit
    vip_folders:dict = {}
    for q in [var for var in datasets[0].data_vars if var.split("_")[-1] != "time"]:
        vip_folders[q] = os.path.join(raw_data_folder, f"{q}_ZgateT1")
        if not os.path.exists(vip_folders[q]): os.mkdir(vip_folders[q])

    # make a new dataset with a new dimension with the attrs["end_time"] of each dataset.
    summarized_nc_paths = {}
    for var in vip_folders:
        end_times = []
        data = []
        for dataset in datasets:
            end_times.append(dataset.attrs["end_time"])
            data.append(array(dataset[var]).tolist()) # shape in (mixer, bias, evo-time)
            time_data = [list(dataset[f"{var}_time"])]*len(datasets)
            

        dict_ = {var:(["end_time","mixer","z_voltage","time"],array(data)),f"{var}_time":(["end_time","mixer","z_voltage","time"],array(time_data))}
        zT1_ds = Dataset(dict_,coords={"end_time":array(end_times),"mixer":array(["I","Q"]),"z_voltage":array(dataset.coords["z_voltage"]),"time":array(dataset.coords["time"])})
        zT1_ds.attrs["z_offset"] = [float(re.search(rf"{var}_(\d+\.\d+)", datasets[0].attrs["ref_bias"]).group(1))]
        zT1_ds.attrs["prepare_excited"] = datasets[0].attrs["prepare_excited"]
        summarized_nc_paths[var] = os.path.join(vip_folders[var],f"{var}_Summaized_zT1.nc")
        zT1_ds.to_netcdf(summarized_nc_paths[var])
    
    # return dict contains nc_path with the q_name as its key 
    return summarized_nc_paths


def Conti2tone_dataReducer(nc_file_path:str):
    ds = open_dataset(nc_file_path)
    ordered_q_labels = ds.attrs['RO_qs'].split(" ")[1:]
    datasets = {}
    for q_idx, q in enumerate(ordered_q_labels):
        x0, x1 = ds[f"x{2*q_idx}"], ds[f"x{2*q_idx+1}"]
        y0, y1 = ds[f"y{2*q_idx}"], ds[f"y{2*q_idx+1}"]

        new_ds = Dataset(
            data_vars = dict(y0=(["dim_0"],y0.data),y1=(["dim_0"],y1.data)),
            coords = dict(x0=(["dim_0"],array(list(x0.data)*int(len(list(x1.data))/len(list(x0.data))))),x1=(["dim_0"],x1.data))
        )
        
        new_ds.attrs = ds.attrs
        new_ds.attrs['target_q'] = q 
        to_copy_array_attr = [x0, x1, y0, y1]
        for idx, item in enumerate([new_ds.x0, new_ds.x1, new_ds.y0, new_ds.y1]):
            item.attrs = to_copy_array_attr[idx].attrs

        datasets[ordered_q_labels[q_idx]] = new_ds

    return datasets

def rofcali_dataReducer(nc_file_path:str)->Dataset:
    ds = open_dataset(nc_file_path)
    return ds

def piampcali_dataReducer(nc_file_path:str)->Dataset:
    ds = open_dataset(nc_file_path)
    return ds

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
    x = []
    for i in range(10):
        x.append(i)
        print(x[-1])
    