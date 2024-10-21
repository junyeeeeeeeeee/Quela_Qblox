import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import open_dataset, Dataset
from numpy import ndarray, array
import quantify_core.data.handling as dh
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.support.QDmanager import Data_manager

def fluxCav_dataReductor(nc_file_path:str, ordered_q_labels:ndarray)->Dataset:
    """
    For flux cavity meas, each qubit will get 2 x-data and 2 y-data, which are 'x0', 'x1', 'y0', 'y1' accordingly.
    """
    ds = open_dataset(nc_file_path)
    datasets = {}
    for q_idx in range(array(ordered_q_labels).shape[0]):
        x0, x1 = ds[f"x{2*q_idx}"], ds[f"x{2*q_idx+1}"]
        y0, y1 = ds[f"y{2*q_idx}"], ds[f"y{2*q_idx+1}"]
        new_ds = Dataset(
            data_vars = dict(y0=(["dim_0"],y0.data),y1=(["dim_0"],y1.data)),
            coords = dict(x0=(["dim_0"],x0.data),x1=(["dim_0"],x1.data))
        )

        new_ds.attrs = ds.attrs
        new_ds.attrs["tuid"] = new_ds.attrs["tuid"].split("-")[0]+"-"+new_ds.attrs["tuid"].split("-")[1]+str(int(new_ds.attrs["tuid"].split("-")[2])+q_idx)+"-"+new_ds.attrs["tuid"].split("-")[3]
        Data_manager().build_tuid_folder(new_ds.attrs["tuid"],f"{ordered_q_labels[q_idx]}FluxCav")
        to_copy_array_attr = [x0, x1, y0, y1]
        for idx, item in enumerate([new_ds.x0, new_ds.x1, new_ds.y0, new_ds.y1]):
            item.attrs = to_copy_array_attr[idx].attrs

        datasets[ordered_q_labels[q_idx]] = new_ds

    return datasets


if __name__ == "__main__":
    dh.set_datadir(meas_raw_dir)
    file = "Modularize/Meas_raw/2024_10_16/DR4q2_FluxCavity_H13M46S37.nc"

    dss = fluxCav_dataReductor(file, ['q0'])
    for q in dss:
        ResonatorFluxSpectroscopyAnalysis(tuid=dss[q].attrs["tuid"], dataset=dss[q]).run(sweetspot_index=0)
    