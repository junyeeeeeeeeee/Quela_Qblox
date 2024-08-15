import os, json, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import xarray as xr
import quantify_core.data.handling as dh
from Modularize.support.QDmanager import QDmanager
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from numpy import sqrt, array, moveaxis, ndarray, cos, sin ,deg2rad, real, imag, linspace
from xarray import open_dataset
import matplotlib.pyplot as plt
from Modularize.support.Pulse_schedule_library import dataset_to_array


if __name__ == '__main__':
    import os
    
    QD_agent = QDmanager('Modularize/QD_backup/2024_8_13/DR1#11_SumInfo.pkl')
    QD_agent.QD_loader()
    qs = ['q4']
    for q in qs:
        print(q,":")
        qubit = QD_agent.quantum_device.get_element(q)
        print(f"bare= {QD_agent.Notewriter.get_bareFreqFor(q)*1e-9} GHz")
        print(f"ROF = {qubit.clock_freqs.readout()*1e-9} GHz")
        print(f"XYF = {qubit.clock_freqs.f01()*1e-9} GHz")
        print(f"x = {(qubit.clock_freqs.readout()-QD_agent.Notewriter.get_bareFreqFor(q))*1e-6} MHz")
        print(f"g = {QD_agent.Notewriter.get_sweetGFor(q)*1e-6} MHz")

    
    
    
    
    
    