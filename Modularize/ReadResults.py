# import xarray as xr
# from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
# from numpy import pi
# results_path = 'Modularize/Meas_raw/2024_2_26/q1_FluxCavSpec_H11M10S5.nc' 
# ds = xr.open_dataset(results_path)

# import quantify_core.data.handling as dh
# meas_datadir = '.data'
# dh.set_datadir(meas_datadir)
# print(ds.attrs["tuid"])
# ana = ResonatorFluxSpectroscopyAnalysis(tuid=ds.attrs["tuid"], dataset=ds).run()
# print(2*pi/ana.quantities_of_interest["frequency"].nominal_value)

from Modularize.support import QDmanager
QD_path = 'Modularize/QD_backup/2024_3_6/DR2#171_SumInfo.pkl'
QD_agent = QDmanager(QD_path)
QD_agent.QD_loader()
for i in ["q0"]:
    q = QD_agent.quantum_device.get_element(i)
    bare = QD_agent.Notewriter.get_bareFreqFor(i)*1e-6
    x = q.clock_freqs.readout()*1e-6-bare
    print("dispersive shift MHz = ",x) 

def aprx_fq(disper_MHz:float,bareF_MHz:float,g_MHz:float=45.0):
    return bareF_MHz-(g_MHz**2)/disper_MHz


print("aprx fq=",aprx_fq(x,bare))