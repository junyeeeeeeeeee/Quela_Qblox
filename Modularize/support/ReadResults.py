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
# from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
# from quantify_core.analysis.base_analysis import Basic2DAnalysis



def plot_QbFlux(nc_path:str, target_q:str):
    ref = [0,0] #Qmanager.refIQ[target_q]
    # plot flux-qubit 
    f,z,i,q = convert_netCDF_2_arrays(nc_path)
    amp = array(sqrt((i-array(ref)[0])**2+(q-array(ref)[1])**2)).transpose()
    print(z)
    print(f)
    print(amp)
    fig, ax = plt.subplots()
    ax:plt.Axes
    c = ax.pcolormesh(z, f, amp, cmap='RdBu')
    fig.colorbar(c, ax=ax)
    plt.show()


    


# from quantify_scheduler.helpers.collections import find_port_clock_path
# QD_agent = QDmanager('Modularize/QD_backup/2024_3_19/DR2#171_SumInfo.pkl')
# QD_agent.QD_loader()
# print(QD_agent.Fluxmanager.get_sweetBiasFor('q2'))
# print(QD_agent.Fluxmanager.get_PeriodFor('q2'))
# qd = QD_agent.quantum_device
# hcfg = QD_agent.Hcfg
# qubit = qd.get_element('q4')
# output_path = find_port_clock_path(
#         hcfg, port=qubit.ports.readout(), clock=qubit.name + ".ro"
#     )
# cluster_key, module_key, output_key, _, _ = tuple(output_path)
# readout_module = hcfg[cluster_key][module_key]
# print(readout_module[output_key]["output_att"])


# meas_datadir = 'tempt'
# dh.set_datadir(meas_datadir)
# print(ds.attrs["tuid"])
# ana = Basic2DAnalysis(tuid=ds.attrs["tuid"], dataset=ds).run()
# print(ana.quantities_of_interest)



# from Modularize.support import QDmanager
# from numpy import array
# QD_path = 'Modularize/QD_backup/2024_3_31/DR2#171_SumInfo.pkl'
# QD_agent = QDmanager(QD_path)
# QD_agent.QD_loader()
# for i in ["q0","q1","q2","q3","q4"]:
#     qu = QD_agent.quantum_device.get_element(i)
#     rof = qu.clock_freqs.readout()
#     xyf = qu.clock_freqs.f01()
#     xyl = qu.rxy.amp180()
#     T1 = QD_agent.Notewriter.get_T1For(target_q=i)
#     T2 = QD_agent.Notewriter.get_T2For(target_q=i)
#     print(f"***{i}:")
#     print(f"rof : {rof}")
#     print(f"xyf : {xyf}")
#     print(f"xyl : {xyl}")
#     print(f"T1 : {T1}")
#     print(f"T2 : {T2}")
#     print("===========================\n")
# def aprx_fq(disper_MHz:float,bareF_MHz:float,g_MHz:float=45.0):
#     return bareF_MHz-(g_MHz**2)/disper_MHz


# print("aprx fq=",aprx_fq(x,bare))

if __name__ == '__main__':
    import os
    from numpy import moveaxis
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


    # plot_QbFlux("Modularize/Meas_raw/2024_8_12/DR4q0_Flux2tone_H13M27S29.nc",'q0')
    
    
    
    
    
    