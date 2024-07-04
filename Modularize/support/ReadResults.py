import os
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



def plot_QbFlux(Qmanager:QDmanager, nc_path:str, target_q:str):
    ref = Qmanager.refIQ[target_q]
    # plot flux-qubit 
    f,z,i,q = convert_netCDF_2_arrays(nc_path)
    amp = array(sqrt((i-array(ref)[0])**2+(q-array(ref)[1])**2)).transpose()
    fig, ax = plt.subplots()
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
    # QD_agent = QDmanager('Modularize/QD_backup/2024_6_11/DR1SCA#11_SumInfo.pkl')
    # QD_agent.QD_loader()
    # qubit = QD_agent.quantum_device.get_element("q0")
    # qubit.measure.acq_delay(280e-9)
    # qubit.measure.pulse_duration(1e-6)
    # qubit.measure.integration_time(500e-9)
    # QD_agent.QD_keeper()
    from matplotlib.figure import Figure
    def save_fig(fig:Figure, all_path:str):
        fig.savefig(all_path)
        plt.close()

    file = "Modularize/Meas_raw/2024_7_4/DR3q4_CavitySpectro_H18M39S34.nc"
    ds = open_dataset(file)

    freq = dict(
        q0=5973306404,
        q1=6083525235,
        q2=5920038579,
        q3=6099702007,
        q4=6010586222
    )
    ro_elements = {}
    for qb in list(freq.keys()):
        ro_elements[qb] = linspace(freq[qb]-10e6, freq[qb]+10e6, 101)
    
    for idx, q in enumerate(freq):

        S21 = ds[f"y{2*idx}"] * cos(
                deg2rad(ds[f"y{2*idx+1}"])
            ) + 1j * ds[f"y{2*idx}"] * sin(deg2rad(ds[f"y{2*idx+1}"]))
        I, Q = real(S21), imag(S21)
        amp = sqrt(I**2+Q**2)
        fig, ax = plt.subplots()
        ax:plt.Axes
        ax.plot(ro_elements[q],amp)
        save_fig(fig,os.path.join("Modularize/under_test",f"test_{idx}.png"))