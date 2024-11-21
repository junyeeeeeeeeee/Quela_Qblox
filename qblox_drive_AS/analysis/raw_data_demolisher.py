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
                print(datasets[-1].attrs)
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
    from qblox_drive_AS.support.QDmanager import QDmanager
    from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
    folder_path = 'qblox_drive_AS/Meas_raw/20241121/H10M16S21'
    QD_path = 'qblox_drive_AS/QD_backup/20241121/DR2#10_SumInfo.pkl'
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    nc_paths = ZgateT1_dataReducer(folder_path)
    for q in nc_paths:
        if QD_agent.rotate_angle[q][0] != 0:
            ref = QD_agent.rotate_angle[q]
        else:
            print(f"{q} rotation angle is 0, use contrast to analyze.")
            ref = QD_agent.refIQ[q]

        ds = open_dataset(nc_paths[q])
        ANA = Multiplex_analyzer("auxA")
        ANA._import_data(ds,var_dimension=2,refIQ=ref)
        ANA._start_analysis(time_sort=False)
        ANA._export_result(nc_paths[q])
        ds.close()

    