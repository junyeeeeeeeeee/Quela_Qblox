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
    want = []
    evo_time, fluxes, i, Q = convert_netCDF_2_arrays(zgateT1_nc_file_path)

    for che in [i, Q]:
        want.append(moveaxis(che,0,-1))  # make shape in [flux, free-time]
    
    output_data = {}
    r_name = os.path.split(zgateT1_nc_file_path)[-1].split("_")[0][-2:]
    
    output_data[r_name] = ( ["mixer","z_voltage","time"],array([want[0], want[1], fluxes.tolist(), evo_time.tolist()]) )
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