# import xarray as xr
# import quantify_core.data.handling as dh
from Modularize.support.QDmanager import QDmanager
# from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
# from quantify_core.analysis.base_analysis import Basic2DAnalysis
# results_path = 'Modularize/Meas_raw/2024_3_14/DR2q1_2tone_H15M35S59.nc' 
# ds = xr.open_dataset(results_path)
# print(ds)
# ds = dh.to_gridded_dataset(ds)
# magnitude_ndarray = ds.y0.data
# phase_ndarray = ds.y1.data
# freq_ndarray = ds.x0.data
# power_ndarray = ds.x1.data


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
from numpy import array
QD_path = 'Modularize/QD_backup/2024_3_21/DR2#171_SumInfo.pkl'
QD_agent = QDmanager(QD_path)
QD_agent.QD_loader()
for i in ["q3"]:
    rof = QD_agent.quantum_device.get_element(i).clock_freqs.readout()
    xyf = QD_agent.quantum_device.get_element(i).clock_freqs.f01()
    
print(rof)
print(xyf)
# def aprx_fq(disper_MHz:float,bareF_MHz:float,g_MHz:float=45.0):
#     return bareF_MHz-(g_MHz**2)/disper_MHz


# print("aprx fq=",aprx_fq(x,bare))


