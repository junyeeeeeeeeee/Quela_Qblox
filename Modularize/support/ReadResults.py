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
    from numpy import moveaxis
    QD_agent = QDmanager('Modularize/QD_backup/2024_8_7/DR4#81_SumInfo.pkl')
    QD_agent.QD_loader()
    qs = ['q2']
    for q in qs:
        print(q,":")
        qubit = QD_agent.quantum_device.get_element(q)
        print(f"bare= {QD_agent.Notewriter.get_bareFreqFor(q)*1e-9} GHz")
        print(f"ROF = {qubit.clock_freqs.readout()*1e-9} GHz")
        print(f"XYF = {qubit.clock_freqs.f01()*1e-9} GHz")
        print(f"x = {(qubit.clock_freqs.readout()-QD_agent.Notewriter.get_bareFreqFor(q))*1e-6} MHz")
        print(f"g = {QD_agent.Notewriter.get_sweetGFor(q)*1e-6} MHz")
    # ds =  dh.to_gridded_dataset(open_dataset("Modularize/Meas_raw/2024_7_11/DR1MultiQ_PowerCavity_H15M20S44.nc"))
    # S21 = ds.y0 * cos(
    #             deg2rad(ds.y1)
    #         ) + 1j * ds.y0 * sin(deg2rad(ds.y1))
    # I, Q = real(S21), imag(S21)
    # amp = array(sqrt(I**2+Q**2))
    # print(amp.shape)

    # x = {"a":{"x":1}}
    # y = {"b":{"y":2}}
    # print(y|x)
    # from qcat.analysis.resonator.photon_dep.res_data import *
    # file = 'Modularize/Meas_raw/2024_7_16/Multiplex_CavityQuality_RTatte80dB_H17M8S0/DR4q2_CavitySpectro_H17M9S1.nc'
    # from scipy.io import loadmat
    # RT_atte = 0

    # a = loadmat(file)
    # I, Q = array(a['ZZI']), array(a['ZZQ'])
    # power = array(a['x']).reshape(-1)
    # freq = array(a['y']).reshape(-1)
    # resonator = PhotonDepResonator('q_test')
    # for power_idx, amp in enumerate(power):
    #     resonator.import_array(freq, I[power_idx]+1j*Q[power_idx], amp-RT_atte)
    # result = resonator.refined_analysis(os.path.split(file)[0])
    
    # csv = "/Users/ratiswu/Downloads/meas_res.csv"
    # import pandas as pd

    # ary = array(pd.read_csv(csv,sep=' '))
    # x = ["q0", "q1", "q2", "q3", "q4"]
    # dic = {}
    # for item in ary:
    #     val = item[0].split(",")
    #     dic[val[0]] = []
    #     for a_item in val[1:]:
    #         dic[val[0]].append(a_item)
    # fig, ax = plt.subplots(figsize=(8,5),dpi=120)
    # ax:plt.Axes
    # ax.errorbar(x, array(dic['T2']).astype(float), array(dic['T2_err']).astype(float),fmt='+-')
    # ax.set_ylim(0,20)
    # ax.xaxis.set_tick_params(labelsize=26)
    # ax.yaxis.set_tick_params(labelsize=26)
    # ax.set_ylabel("µs",fontsize=26)
    # ax.set_title("T2",fontsize=26)
    # plt.grid()
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig("/Users/ratiswu/Downloads/T2_distri.png")
    # plt.close()
    
    
    
    
    
    