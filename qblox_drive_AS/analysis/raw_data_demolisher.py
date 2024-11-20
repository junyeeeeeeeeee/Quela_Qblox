import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import open_dataset, Dataset
from numpy import array

class MultiplexingDataReducer():
    def __init__(self, nc:str|Dataset):
        if isinstance(nc,str):
            self.ds = open_dataset(nc)
        elif isinstance(nc,Dataset):
            self.ds = nc
        else:
            raise TypeError("You should give a nc_file path or a xr.Dataset")
        








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
        zT1_ds.attrs["z_offset"] = [float(datasets[0].attrs[f"{var}_ref_bias"])]
        zT1_ds.attrs["prepare_excited"] = datasets[0].attrs["prepare_excited"]
        summarized_nc_paths[var] = os.path.join(vip_folders[var],f"{var}_Summaized_zT1.nc")
        zT1_ds.to_netcdf(summarized_nc_paths[var])
    
    # return dict contains nc_path with the q_name as its key 
    return summarized_nc_paths






if __name__ == "__main__":
    x = []
    for i in range(10):
        x.append(i)
        print(x[-1])
    