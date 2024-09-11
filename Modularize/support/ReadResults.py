import os, json, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import xarray as xr
import quantify_core.data.handling as dh
from Modularize.support.QDmanager import QDmanager
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from numpy import sqrt, array, moveaxis, ndarray, cos, sin ,deg2rad, real, imag, linspace
from xarray import open_dataset
import matplotlib.pyplot as plt
from Modularize.support.Pulse_schedule_library import dataset_to_array, IQ_data_dis, T1_fit_analysis, Fit_analysis_plot, T2_fit_analysis


if __name__ == '__main__':
    import os
    
    QD_agent = QDmanager('Modularize/QD_backup/2024_8_29/DR4#81_SumInfo.pkl')
    QD_agent.QD_loader()
    # qs = ['q1']
    # # for q in qs:
    # #     print(q,":")
    # #     qubit = QD_agent.quantum_device.get_element(q)
    # #     print(f"bare= {QD_agent.Notewriter.get_bareFreqFor(q)*1e-9} GHz")
    # #     print(f"ROF = {qubit.clock_freqs.readout()*1e-9} GHz")
    # #     print(f"XYF = {qubit.clock_freqs.f01()*1e-9} GHz")
    # #     print(f"x = {(qubit.clock_freqs.readout()-QD_agent.Notewriter.get_bareFreqFor(q))*1e-6} MHz")
    # #     print(f"g = {QD_agent.Notewriter.get_sweetGFor(q)*1e-6} MHz")

    # file = 'Modularize/Meas_raw/2024_8_29/DR4q1_T2(0)_H12M58S58.nc'
    # nc = open_dataset(file)
    # I,Q= dataset_to_array(dataset=nc,dims=1)
    # data= IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[qs[0]][0],ref_Q=QD_agent.refIQ[qs[0]][-1])
    # data_fit= T2_fit_analysis(data=data,freeDu=array(nc.x0),T2_guess=5e-6)
    # Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)
    
    
    
    
    