import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from numpy import flip, arange, argmin, argmax, diff, array, all, sqrt, std, mean, sort, cos, sin
from Modularize.support import QDmanager
from Modularize.support.Pulse_schedule_library import IQ_data_dis
from numpy import ndarray, cos, sin, deg2rad, real, imag, transpose, delete, diff, where

def Z_sperate_del(datapoint:ndarray,flux_range:float):
    """
    Considering a given z_array after F_advan_del processed, if a neighboring z is seperated by the given flux_range (threshold), then remove the elements behind it.\n
    Return the remove starting index. If it's 0, there is unnecessary to remove.
    """
    z_ary = datapoint[:,0]
    suspicios_idx = []
    for i in range(z_ary.shape[0]):
        if i != z_ary.shape[0]-1:
            z_difference = abs(z_ary[i+1] - z_ary[i])
            if z_difference >= flux_range:
                suspicios_idx.append(i)
    
    if len(suspicios_idx) >= 2:
        max_gap_idx, = where(diff(suspicios_idx) == max(diff(suspicios_idx)))
        if max_gap_idx > datapoint.shape[0]/2:
            return array(datapoint[:max_gap_idx])
        else:
            return array(datapoint[max_gap_idx+1:])     
    elif len(suspicios_idx) == 1:
        if suspicios_idx[0] > datapoint.shape[0]/2:
            return array(datapoint[:suspicios_idx[0]])
        else:
            return array(datapoint[suspicios_idx[0]+1:])
    else:
        return datapoint

    
def F_seperate_del(datapoint:ndarray,freq_range:float=100e6):
    """
    Considering a given z_array after F_advan_del processed, if a neighboring frequency is seperated by the given freq_range (threshold), then remove the elements behind it.\n
    Return the remove starting index. If it's 0, there is unnecessary to remove.\n
    The default threshold is 100 MHz.
    """
    f_ary = datapoint[:,1]
    del_idx = []
    mark = False
    for i in range(f_ary.shape[0]):
        if i != f_ary.shape[0]-1:
            if not mark :
                f_difference = abs(f_ary[i+1] - f_ary[i])
                if f_difference >= freq_range:
                    del_idx.append(i+1)
                    mark = True
            else:
                f_difference = abs(f_ary[i+1] - f_ary[i-1])
                mark = False
                if f_difference >= freq_range:
                    del_idx.append(i+1)
                    mark = True

    if len(del_idx) != 0:
        return delete(datapoint,del_idx,axis=0)
    else:
        return datapoint


# Plotting
def plot_HeatScat(mag,x_heat_ary,y_heat_ary,x_scat_ary,y_scat_ary,fit_scat_ary:ndarray=[]):
    import plotly.graph_objects as go
    if fit_scat_ary.shape[0] != 0:
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
    else:
        data = [
            go.Heatmap(
                z=mag,
                x=x_heat_ary,
                y=y_heat_ary,
                colorscale="YlOrRd"
            ),
            go.Scatter(x=x_scat_ary,y=y_scat_ary,mode='markers',marker_color='rgb(0, 0, 255)',name="XYF_maxima")
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

def filter_2D(raw_mag:ndarray,threshold:float=3.0):
    filtered_f_idx = []
    filtered_z_idx = []
    mu = mean(raw_mag.reshape(-1))
    sigma = std(raw_mag.reshape(-1))
    for i_z in range(raw_mag.shape[0]):
        sofar_max = min(raw_mag.reshape(-1))# the maximum signal in the same bias
        z_idx_champion = 0
        f_idx_champion = 0
        for i_f in range(raw_mag.shape[1]):
            if raw_mag[i_f][i_z] > mu + threshold*sigma:
                if raw_mag[i_f][i_z] > sofar_max:
                    f_idx_champion = i_f
                    z_idx_champion = i_z 
                    sofar_max = raw_mag[i_f][i_z]
        if z_idx_champion != 0 or f_idx_champion != 0:
            filtered_z_idx.append(z_idx_champion)
            filtered_f_idx.append(f_idx_champion)
    
    return filtered_f_idx, filtered_z_idx


def sortAndDecora(raw_z:ndarray,raw_XYF:ndarray,raw_mag:ndarray,threshold:float=3):
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
    if len(filtered) == 0:
        raise ValueError("The given filter threshold is too heavy, it should be decreased!")

    filtered = F_seperate_del(array(filtered))
    filtered = Z_sperate_del(array(filtered),0.04)
    # filtered = F_seperate_del(array(filtered))
    
    return array(filtered)
def get_arrays_from_netCDF(netCDF_path:str,ref_IQ:list=[0,0]):
    dataset_processed = dh.to_gridded_dataset(xr.open_dataset(netCDF_path))
    XYF = flip(dataset_processed["x0"].to_numpy())
    z = dataset_processed["x1"].to_numpy()
    S21 = transpose(dataset_processed.y0.data * cos(
            deg2rad(dataset_processed.y1.data)
        ) + 1j * dataset_processed.y0.data * sin(
            deg2rad(dataset_processed.y1.data)
        )
    )
    I, Q = real(S21), imag(S21)
    displaced_magnitude = flip(transpose(IQ_data_dis(I_data=I,Q_data=Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])),axis=0)
    return XYF, z, displaced_magnitude


def quadra(x,a,b,c):
    return a*x**2+b*x+c

def calc_g(fb,fq,coefA):
    return coefA*sqrt(fb*fq)/1000

def dispersive():
    pass

def FqEqn(x,Ec,coefA,a,b,d):
    return sqrt(8*coefA*Ec*sqrt(cos(a*(x-b))**2+d**2*sin(a*(x-b))**2))-Ec
     
def calc_Gcoef_inFbFqFd(bareF:float,Fq:float,dressF:float):
    return sqrt((1000**2)*(dressF-bareF)*(bareF-Fq)/(sqrt(bareF*Fq)**2))

# predict fq
def calc_fq_g_excluded(coef_inG:float,fdress:float,fbare:float):
    """
    After we got the coef in coupling g, use dispersive relation to calc fq with the given args.
    """
    return fbare*((fdress-fbare)/((fbare*(coef_inG/1000)**2)+(fdress-fbare)))

def z_aranger(loaded_QDagent:QDmanager,target_q:str,artif_shift_inZ:float=0.0,period_devider:int=8):
    sweet = loaded_QDagent.Fluxmanager.get_sweetBiasFor(target_q)
    bottom = loaded_QDagent.Fluxmanager.get_tuneawayBiasFor(target_q)
    waist = (sweet+bottom)/2
    z_span = loaded_QDagent.Fluxmanager.get_PeriodFor(target_q)/period_devider
    
    return {"sweet":[sweet-z_span+artif_shift_inZ, sweet+z_span+artif_shift_inZ],
            "waist":[waist-z_span+artif_shift_inZ, waist+z_span+artif_shift_inZ],
            "bottom":[bottom-z_span+artif_shift_inZ, bottom+z_span+artif_shift_inZ]}

def rof_setter(loaded_QDagent:QDmanager,target_q:str='q0',bias_position:str='sweet'):
    """
    Return the readout freq according to the given z bias.\n
    args:\n
    bias_position: (1) 'sweet' for sweet spot. (2) 'waist' for the middle between sweet spot and bottom. (3) 'bottom' for the bottom position.
    """
    if bias_position == 'sweet':
        z = loaded_QDagent.Fluxmanager.get_sweetBiasFor(target_q)
    elif bias_position == 'bottom':
        z = loaded_QDagent.Fluxmanager.get_tuneawayBiasFor(target_q)
    else:
        z = (loaded_QDagent.Fluxmanager.get_sweetBiasFor(target_q)+loaded_QDagent.Fluxmanager.get_tuneawayBiasFor(target_q))/2
    
    rof = loaded_QDagent.Fluxmanager.sin_for_cav(target_q, bias_ary=array([z]))
    return rof


def plot_raw_data(path:str):
    pass

if __name__ == "__main__":
    worse = ''
    better = 'Modularize/Meas_raw/2024_3_18/DR2q2_Flux2tone_H14M22S59.nc'
    QD_path = 'Modularize/QD_backup/2024_3_18/DR2#171_SumInfo.pkl'
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()



    # Raw data extract
    XYF, z, mag = get_arrays_from_netCDF(netCDF_path=better,ref_IQ=QD_agent.refIQ['q2'])#
    

    # Try using quadratic fit the symmetry axis
    from numpy import polyfit

    extracted = sortAndDecora(z,XYF,mag,threshold=1)
    x_wanted, y_wanted = extracted[:,0], extracted[:,1]
    coef = polyfit(x_wanted, y_wanted,deg=2)
    fit_ary = quadra(x_wanted,*coef)

    plot_HeatScat(mag=mag,x_heat_ary=z,x_scat_ary=x_wanted,y_heat_ary=XYF,y_scat_ary=y_wanted,fit_scat_ary=fit_ary)

    # get bare and readout freq
    # QD_path = 'Modularize/QD_backup/2024_2_27/SumInfo.pkl'
    # QD_agent = QDmanager(QD_path)
    # QD_agent.QD_loader()
    # fb = float(QD_agent.Notewriter.get_bareFreqFor(target_q="q2"))*1e-9
    # fd = QD_agent.Fluxmanager.sin_for_cav(target_q="q2",bias_ary=sweetspot_z)*1e-9 # unit in GHz
    # coef_inG = calc_Gcoef_inFbFqFd(fb,sweetspot_f,fd)
    # fq = calc_fq_g_excluded(coef_inG,fd,fb)
