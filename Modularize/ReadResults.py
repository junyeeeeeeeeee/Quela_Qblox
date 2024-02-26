import xarray as xr
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from numpy import pi
results_path = 'Modularize/Meas_raw/2024_2_26/q1_FluxCavSpec_H11M10S5.nc' 
ds = xr.open_dataset(results_path)

import quantify_core.data.handling as dh
meas_datadir = '.data'
dh.set_datadir(meas_datadir)
print(ds.attrs["tuid"])
ana = ResonatorFluxSpectroscopyAnalysis(tuid=ds.attrs["tuid"], dataset=ds).run()
print(2*pi/ana.quantities_of_interest["frequency"].nominal_value)