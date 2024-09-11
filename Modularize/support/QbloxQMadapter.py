"""
Aim to transform the dataset generated from Qblox into a analyzable format by the analysis prgrom in QM side.
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import xarray as xr
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from Modularize.support.Pulse_schedule_library import dataset_to_array
from Modularize.support.QDmanager import Data_manager
from numpy import moveaxis, array

class DataTransformer:
    """
    Transform the dataset from Qblox into other system analyzable shape.
    """
    def __init__(self,output_data_system:str):
        if output_data_system.lower() in ["qm"]:
            self.output_system__ =  "qm"
        else:
            raise KeyError(f"Now we can only transform the data into QM's shape, chech your given `output_data_system` = {output_data_system}")

    def zgateT1_adapter(self, zgateT1_nc_file_path:str, ref_z_offset:float, ref_IQ:list=[], save_path:str='')->xr.Dataset:
        """
        trnaslate the given raw data form into the shape (IQ, Z, evoTime) and save the dataset if the `save_path` was given.\n
        -------
        ## Return:\n
        The dataset with the coordinates: 'mixer', 'time', 'z_voltage'

        """
        evo_time, fluxes, _, _ = convert_netCDF_2_arrays(zgateT1_nc_file_path)
        ds = xr.open_dataset(zgateT1_nc_file_path)
        i, q = dataset_to_array(dataset=ds,dims=2)
        I, Q = moveaxis(i,0,-1), moveaxis(q,0,-1)
        if self.output_system__ == "qm":
            output_data = {}
            r_name = os.path.split(zgateT1_nc_file_path)[-1].split("_")[0][-2:]
            
            output_data[r_name] = ( ["mixer","z_voltage","time"],array([I.tolist(), Q.tolist()]) )
            dataset = xr.Dataset(
                output_data,
                coords={ "mixer":array(["I","Q"]), "time": evo_time, "z_voltage":fluxes }
            )
            dataset.attrs["z_offset"] = [ref_z_offset]
            dataset.attrs["ref_IQ"] = ref_IQ

        if save_path != "":
            path = os.path.join(save_path,f"To{self.output_system__.upper()}_{os.path.split(zgateT1_nc_file_path)[-1]}")
            self.save_transformed_data(dataset, path)
        
        return dataset
    
    def save_transformed_data(self, dataset:xr.Dataset, save_path:str, label:str="", format:str='nc'):
        if format.lower() == 'nc':
            if label != "":
                save_path = save_path.split(".")[0] + f"({label})" + ".nc"
            dataset.to_netcdf(save_path)
        else:
            raise KeyError(f"Un-supported data format was given = {format}")
        
        return save_path
    


if __name__ == "__main__" :
    # files = ["Modularize/Meas_raw/2024_6_27/DR1SCAq0_zT1(1)_H14M8S35.nc"]
    raw_folder = "Modularize/Meas_raw/z_gate_T1_test/z_gate_T1_pi_False"
    files = [os.path.join(raw_folder,name) for name in os.listdir(raw_folder) if (os.path.isfile(os.path.join(raw_folder,name)) and name.split(".")[-1]=='nc')]
    for file in files:
        QM_folder = os.path.join(os.path.split(file)[0],"ToQM")
        if not os.path.exists( QM_folder ):
            os.mkdir( QM_folder )
        ds = DataTransformer(output_data_system='qm').zgateT1_adapter(file,0,save_path=QM_folder)


