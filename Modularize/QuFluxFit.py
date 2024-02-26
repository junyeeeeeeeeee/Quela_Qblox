import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from numpy import flip, arange, argmin, argmax

"""
from scipy.ndimage import convolve1d as conv
from scipy.ndimage import gaussian_filter1d as g_filter
from numpy import  hanning

kern = hanning(15)   # a Hanning window with width 50
kern /= kern.sum()
mag = conv(mag,kern,0)
"""
results_path = 'Modularize/Meas_raw/2024_2_26/q1_FluxQ_H14M56S55.nc' 
ds = dh.to_gridded_dataset(xr.open_dataset(results_path))
ROF = flip(ds["x0"].to_numpy())
z = ds["x1"].to_numpy()
mag = flip(ds["y0"].to_numpy(),axis=0)
pha = ds["y1"].to_numpy()

# default the peak is at the minima
F_on_res = ROF[argmax(mag, axis=0)]


import plotly.express as pt
import plotly.graph_objects as go

data = [
    go.Heatmap(
        z=mag,
        x=z,
        y=ROF,
        colorscale="YlOrRd"
    ),
    go.Scatter(x=z,y=F_on_res,mode='markers',marker_color='rgb(0, 0, 255)',name="Dressed_fitted")
]
layout = go.Layout(
    title = "Cavity Flux Dependence with fitted Dressed cavity",
    xaxis=dict(title='Flux (V)'),
    yaxis=dict(title='Frequency (GHz)'),
    showlegend=True,
    legend={'x': .9, 'y': .0001}
)
fig = go.Figure(data=data, layout=layout)
fig.show()