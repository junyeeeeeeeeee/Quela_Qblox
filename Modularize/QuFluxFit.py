import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from numpy import flip, arange, argmin, argmax, diff, array, all, sqrt, std, mean, sort, cos, sin

from numpy import ndarray

def Z_sperate_del(datapoint:ndarray,flux_range:float):
    """
    Considering a given z_array after F_advan_del processed, if a neighboring z is seperated by the given flux_range (threshold), then remove the elements behind it.\n
    Return the remove starting index. If it's 0, there is unnecessary to remove.
    """
    z_ary = datapoint[:,0]
    break_idx = 0
    for i in range(z_ary.shape[0]):
        if i != z_ary.shape[0]-1:
            z_difference = z_ary[i+1] - z_ary[i]
            if z_difference >= flux_range:
                break_idx = i
                break

    if break_idx != 0:
        return array(datapoint[:break_idx])
    else:
        return datapoint


# Plotting
def plot_HeatScat(mag,x_heat_ary,y_heat_ary,x_scat_ary,y_scat_ary,fit_scat_ary):
    import plotly.graph_objects as go
    data = [
        go.Heatmap(
            z=mag,
            x=x_heat_ary,
            y=y_heat_ary,
            colorscale="YlOrRd"
        ),
        go.Scatter(x=x_scat_ary,y=y_scat_ary,mode='markers',marker_color='rgb(0, 0, 255)',name="XYF_maxima"),
        go.Scatter(x=x_scat_ary,y=fit_scat_ary,mode='lines',marker_color='rgb(0, 210, 0)',name="Fit line"),
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

def calc_g(fb,fq,coefA):
    return coefA*sqrt(fb*fq)/1000

def dispersive():
    pass

def FqEqn(x,Ec,coefA,a,b,d):
    return sqrt(8*coefA*Ec*sqrt(cos(a*(x-b))**2+d**2*sin(a*(x-b))**2))-Ec
     
def calc_Gcoef_inFbFqFd(bareF:float,Fq:float,dressF:float):
    return sqrt((1000**2)*(dressF-bareF)*(bareF-Fq)/(sqrt(bareF*Fq)**2))
worse = ''
better = 'Modularize/Meas_raw/2024_2_27/q2_FluxQ_H16M49S35.nc'

# Raw data extract
results_path = better
ds = dh.to_gridded_dataset(xr.open_dataset(results_path))
XYF = flip(ds["x0"].to_numpy())
z = ds["x1"].to_numpy()
mag = flip(ds["y0"].to_numpy(),axis=0)
pha = ds["y1"].to_numpy()


def filter_2D(raw_mag:ndarray,threshold:float=3.0):
    filtered_f_idx = []
    filtered_z_idx = []
    mu = mean(raw_mag.reshape(-1))
    sigma = std(raw_mag.reshape(-1))
    for i_z in range(raw_mag.shape[0]):
        for i_f in range(raw_mag.shape[1]):
            if raw_mag[i_f][i_z] > mu + threshold*sigma:
                filtered_f_idx.append(i_f)
                filtered_z_idx.append(i_z)
    return filtered_f_idx, filtered_z_idx


def sortAndDecora(raw_z:ndarray,raw_XYF:ndarray,raw_mag:ndarray,threshold:float=2.5):
    def takeLast(elem):
        return elem[-1]

    a, b = filter_2D(raw_mag,threshold)
    extracted = []
    for i in range(len(a)):
        extracted.append([raw_z[b[i]],raw_XYF[a[i]]])
    # filter out the freq when it's under the same bias
    filtered = []
    for point_i in range(len(extracted)):
        same_bias = []
        picked = extracted[point_i][0]
        same_bias.append(extracted[point_i])
        for point_j in range(point_i,len(extracted)):
            if extracted[point_j][0] == picked:
                same_bias.append(extracted[point_j])
        same_bias.sort(key=takeLast)
        if filtered != []:
            if filtered[-1][0] != same_bias[-1][0]:
                filtered.append(same_bias[-1])
        else:
            filtered.append(same_bias[-1])

    sort(array(filtered),axis=0)
    filtered = Z_sperate_del(array(filtered),0.04)

    return filtered




# Try using quadratic fit the symmetry axis
from numpy import polyfit
def quadra(x,a,b,c):
    return a*x**2+b*x+c
extracted = sortAndDecora(z,XYF,mag)
x_ary, y_ary = extracted[:,0], extracted[:,1]
coef = polyfit(x_ary,y_ary,deg=2)
fit_ary = quadra(x_ary,*coef)

sweetspot_z = x_ary[argmax(fit_ary)]
sweetspot_f = float(y_ary[argmax(fit_ary)])*1e-9
plot_HeatScat(mag=mag,x_heat_ary=z,x_scat_ary=x_ary,y_heat_ary=XYF,y_scat_ary=y_ary,fit_scat_ary=fit_ary)

# get bare and readout freq
from Modularize.support import QDmanager
QD_path = 'Modularize/QD_backup/2024_2_27/SumInfo.pkl'
QD_agent = QDmanager(QD_path)
QD_agent.QD_loader()
fb = float(QD_agent.Notewriter.get_bareFreqFor(target_q="q2"))*1e-9
fd = QD_agent.Fluxmanager.sin_for_cav(target_q="q2",bias_ary=sweetspot_z)*1e-9 # unit in GHz
print(f"fq = {sweetspot_f}")
print(f"bare f ={fb}")
print(f"dressed f ={fd}")
coefG = calc_Gcoef_inFbFqFd(fb,sweetspot_f,fd)
print(f"coefG = {coefG}")
g = calc_g(fb,sweetspot_f,coefG)
print(f"g (GHz)= {g}")



