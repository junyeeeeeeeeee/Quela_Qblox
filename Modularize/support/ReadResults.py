import xarray as xr
import quantify_core.data.handling as dh
from Modularize.support import QDmanager
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
from quantify_core.analysis.base_analysis import Basic2DAnalysis
results_path = 'Modularize/Meas_raw/2024_3_14/DR2q1_2tone_H15M35S59.nc' 
ds = xr.open_dataset(results_path)
print(ds)
# ds = dh.to_gridded_dataset(ds)
# magnitude_ndarray = ds.y0.data
# phase_ndarray = ds.y1.data
# freq_ndarray = ds.x0.data
# power_ndarray = ds.x1.data


# from quantify_scheduler.helpers.collections import find_port_clock_path
# QD_agent = QDmanager('Modularize/QD_backup/2024_3_8/DR1#170_SumInfo.pkl')
# QD_agent.QD_loader()
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
# QD_path = 'Modularize/QD_backup/2024_3_6/DR2#171_SumInfo.pkl'
# QD_agent = QDmanager(QD_path)
# QD_agent.QD_loader()
# for i in ["q0"]:
#     q = QD_agent.quantum_device.get_element(i)
#     bare = QD_agent.Notewriter.get_bareFreqFor(i)*1e-6
#     x = q.clock_freqs.readout()*1e-6-bare
#     print("dispersive shift MHz = ",x) 

# def aprx_fq(disper_MHz:float,bareF_MHz:float,g_MHz:float=45.0):
#     return bareF_MHz-(g_MHz**2)/disper_MHz


# print("aprx fq=",aprx_fq(x,bare))
