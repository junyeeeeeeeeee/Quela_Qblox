"""
Aim to transform the dataset generated from Qblox into a analyzable format by the analysis prgrom in QM side.
"""

import xarray as xr
import os
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from numpy import moveaxis, array

def zgateT1_Qblox2QM_adapter(zgateT1_nc_file_path:str, ref_z_offset:float, save_path:str='')->xr.Dataset:
    """
    trnaslate the given raw data form into the shape (IQ, Z, evoTime) and save into the QM analyzable shape\n
    -------
    ## Return:\n
    The dataset with the coordinates: 'mixer', 'time', 'z_voltage'

    """
    
    evo_time, fluxes, I, Q = convert_netCDF_2_arrays(zgateT1_nc_file_path)
    output_data = {}
    r_name = os.path.split(zgateT1_nc_file_path)[-1].split("_")[0][-2:]
    print(I.shape)
    print(Q.shape)
    print(fluxes.shape)
    print(evo_time.shape)
    
    output_data[r_name] = ( ["mixer","z_voltage","time"],[[I.reshape(-1).tolist(), Q.reshape(-1).tolist()], fluxes.tolist(), evo_time.tolist()] )
    dataset = xr.Dataset(
        output_data,
        coords={ "mixer":array(["I","Q"]), "time": evo_time, "z_voltage":fluxes }
    )
    dataset.attrs["z_offset"] = list(ref_z_offset)
    if save_path != "":
        name = os.path.split(zgateT1_nc_file_path)[-1]
        path = os.path.join(save_path,f"ToQM_{name}")
        dataset.to_netcdf(path)
    
    return dataset


if __name__ == "__main__" :
    ds = zgateT1_Qblox2QM_adapter("Modularize/Meas_raw/2024_6_27/DR1SCAq0_zT1(1)_H14M8S35.nc",0)
    di = ds.to_dict()

