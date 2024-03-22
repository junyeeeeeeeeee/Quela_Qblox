import os, json
import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from typing import Callable
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from numpy import flip, pi, linspace, array, sqrt, std, mean, sort, diag
from Modularize.support import QDmanager, Data_manager
from Modularize.support.Pulse_schedule_library import IQ_data_dis
from numpy import ndarray, cos, sin, deg2rad, real, imag, transpose, delete, diff, where
from scipy.optimize import curve_fit

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
        max_gap_idx, = where(diff(suspicios_idx) == max(diff(suspicios_idx)))[0]
        if max_gap_idx > datapoint.shape[0]/2:
            return array(datapoint[:max_gap_idx])
        else:
            print(max_gap_idx+1)
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
def plot_HeatScat(mag,x_heat_ary,y_heat_ary,x_scat_ary,y_scat_ary,fit_scat_ary:ndarray=array([]),q:str=''):
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
        title = f"{q} Flux vs. Transition frequency",
        xaxis=dict(title='Flux (V)'),
        yaxis=dict(title='XY Frequency (GHz)'),
        showlegend=True,
        legend={'x': .9, 'y': .0001}
    )
    fig = go.Figure(data=data, layout=layout)
    fig.show()

def filter_2D(raw_mag:ndarray,threshold:float=3.0):
    filtered_f_idx = []
    filtered_z_idx = []
    filtered_mag   = []
    mu = mean(raw_mag.reshape(-1))
    sigma = std(raw_mag.reshape(-1))
    for i_z in range(raw_mag.shape[1]):
        sofar_max = min(raw_mag.reshape(-1))# the maximum signal in the same bias
        z_idx_champion = 0
        f_idx_champion = 0
        for i_f in range(raw_mag.shape[0]):
            if raw_mag[i_f][i_z] > mu + threshold*sigma:
                if raw_mag[i_f][i_z] > sofar_max:
                    f_idx_champion = i_f
                    z_idx_champion = i_z 
                    sofar_max = raw_mag[i_f][i_z]
        if z_idx_champion != 0 or f_idx_champion != 0:
            filtered_z_idx.append(z_idx_champion)
            filtered_f_idx.append(f_idx_champion)
            filtered_mag.append(sofar_max)
    
    return filtered_f_idx, filtered_z_idx, filtered_mag


def sortAndDecora(raw_z:ndarray,raw_XYF:ndarray,raw_mag:ndarray,threshold:float=3):
    def takeLast(elem):
        return elem[-1]

    while True:
        f_idx, z_idx, XL_mag = filter_2D(raw_mag,threshold)
        if len(f_idx)==0 and len(z_idx)==0:
            if threshold != 0:
                threshold-=0.5
            else:
                print(f"This interval can't find the trend : XYF={raw_XYF[0]}~{raw_XYF[-1]}")
        else:
            break

    extracted = []
    for i in range(len(f_idx)):
        extracted.append([raw_z[z_idx[i]],raw_XYF[f_idx[i]],XL_mag[i]])
    # filter out the freq when it's under the same bias
    filtered = []
    for point_i in range(len(extracted)):
        same_bias = []
        picked_z = extracted[point_i][0]
        same_bias.append(extracted[point_i])
        for point_j in range(point_i,len(extracted)):
            if extracted[point_j][0] == picked_z:
                same_bias.append(extracted[point_j])
        same_bias.sort(key=takeLast)
        if filtered != []:
            if filtered[-1][0] != same_bias[-1][0]:
                filtered.append(same_bias[-1])
        else:
            filtered.append(same_bias[-1])

    sort(array(filtered),axis=0)

    # filtered = F_seperate_del(array(filtered))
    # filtered = Z_sperate_del(array(filtered),0.04)
    # filtered = F_seperate_del(array(filtered))
    return array(filtered)

def convert_netCDF_2_arrays(CDF_path:str):
    """
    For Qblox system, give a netCDF file path to return some ndarrays.
    ## Return: XYF, z, I, Q
    """
    dataset_processed = dh.to_gridded_dataset(xr.open_dataset(CDF_path))
    XYF = flip(dataset_processed["x0"].to_numpy())
    z = dataset_processed["x1"].to_numpy()
    S21 = transpose(dataset_processed.y0.data * cos(
            deg2rad(dataset_processed.y1.data)
        ) + 1j * dataset_processed.y0.data * sin(
            deg2rad(dataset_processed.y1.data)
        )
    )
    I, Q = real(S21), imag(S21)
    return XYF, z, I, Q

def mag_repalce_origin(I_ary:ndarray,Q_ary:ndarray,ref_IQ,qblox_rotation:bool=False):
    """
    For Qblox system, data needs to be transpose and flip, turn on the `qblox_rotation`.
    """
    if qblox_rotation:
        displaced_magnitude = flip(transpose(IQ_data_dis(I_data=I_ary,Q_data=Q_ary,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])),axis=0)
    else:
        displaced_magnitude = IQ_data_dis(I_data=I_ary,Q_data=Q_ary,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])

    return displaced_magnitude


def quadra(x,a,b,c):
    return a*x**2+b*x+c

def calc_g(fb,fq,coefA):
    return coefA*sqrt(fb*fq)/1000

def FqEqn(x,a,b,Ec,coefA,d):
    """
    a ~ period, b ~ offset, 
    """
    return sqrt(8*coefA*Ec*sqrt(cos(a*(x-b))**2+d**2*sin(a*(x-b))**2))-Ec
     
def calc_Gcoef_inFbFqFd(bareF:float,Fq:float,dressF:float):
    return sqrt((1000**2)*(dressF-bareF)*(bareF-Fq)/(sqrt(bareF*Fq)**2))

# predict fq
def calc_fq_g_excluded(coef_inG:float,fdress:float,fbare:float):
    """
    After we got the coef in coupling g, use dispersive relation to calc fq with the given args.
    """
    return fbare*((fdress-fbare)/((fbare*(coef_inG/1000)**2)+(fdress-fbare)))

# def z_aranger(loaded_QDagent:QDmanager,target_q:str,artif_shift_inZ:float=0.0,period_devider:int=8):
#     sweet = loaded_QDagent.Fluxmanager.get_sweetBiasFor(target_q)
#     bottom = loaded_QDagent.Fluxmanager.get_tuneawayBiasFor(target_q)
#     waist = (sweet+bottom)/2
#     z_span = loaded_QDagent.Fluxmanager.get_PeriodFor(target_q)/period_devider
    
#     return {"sweet":[sweet-z_span+artif_shift_inZ, sweet+z_span+artif_shift_inZ],
#             "waist":[waist-z_span+artif_shift_inZ, waist+z_span+artif_shift_inZ],
#             "bottom":[bottom-z_span+artif_shift_inZ, bottom+z_span+artif_shift_inZ]}

# def rof_setter(loaded_QDagent:QDmanager,target_q:str='q0',bias_position:str='sweet'):
#     """
#     Return the readout freq according to the given z bias.\n
#     args:\n
#     bias_position: (1) 'sweet' for sweet spot. (2) 'waist' for the middle between sweet spot and bottom. (3) 'bottom' for the bottom position.
#     """
#     if bias_position == 'sweet':
#         z = loaded_QDagent.Fluxmanager.get_sweetBiasFor(target_q)
#     elif bias_position == 'bottom':
#         z = loaded_QDagent.Fluxmanager.get_tuneawayBiasFor(target_q)
#     else:
#         z = (loaded_QDagent.Fluxmanager.get_sweetBiasFor(target_q)+loaded_QDagent.Fluxmanager.get_tuneawayBiasFor(target_q))/2
    
#     rof = loaded_QDagent.Fluxmanager.sin_for_cav(target_q, bias_ary=array([z]))
#     return rof

def data2plot(XYF_array:ndarray,z_array:ndarray,I_array:ndarray,Q_array:ndarray,specified_refIQ:list,filter2D_threshold:float=1.0,qblox:bool=False,q:str=''):
    
    mag = mag_repalce_origin(I_array,Q_array,ref_IQ=specified_refIQ,qblox_rotation=qblox)
    extracted = sortAndDecora(z_array,XYF_array,mag,threshold=filter2D_threshold)
    x_wanted, y_wanted, mag_wanted = extracted[:,0], extracted[:,1], extracted[:,-1]
    plot_HeatScat(mag=mag,x_heat_ary=z_array,x_scat_ary=x_wanted,y_heat_ary=XYF_array,y_scat_ary=y_wanted,q=q)
    return x_wanted, y_wanted, mag_wanted


def read_fq_data(json_path:str):
    """
    Read a json file contains the fq vs. flux data.\n
    The key in this dict should be 'x' and 'y', which repesents `flux` and `fq`.\n
    Return flux array and fq array in GHz.
    """
    data2fit = {}
    with open(json_path) as json_file:
        data2fit = json.load(json_file)
    data = []
    for i in range(len(data2fit['x'])):
        data.append([data2fit['x'][i],data2fit['y'][i]])
    data.sort(key=lambda x:x[0])
    flux = array(data)[:,0]
    f01 = array(data)[:,1]*1e-9 # in GHz
    
    return flux, f01


def set_fitting_paras(period:float,offset:float,flux_array:ndarray,Ec_guess_GHz:float=0.21,Ej_sum_guess_GHz:float=20.0,squid_ratio_guess:float=0.5):
    """
    There are 5 paras in Fq eqn, give the initial guess and the fitting constrains for curve fit.\n
    Return guess, upper_bound, bottom_bound.
    """
    f = 2*pi/period
    b = offset/f
    guess = (f,b,Ec_guess_GHz,Ej_sum_guess_GHz,squid_ratio_guess) #[a, b, Ec, Ej_sum, d]
    wide_period = 5*period/4
    narrow_period = period/4
    upper_bound = [2*pi/narrow_period,max(flux_array),0.25,100,1] #[a, b, Ec, Ej_sum, d]
    bottom_bound = [2*pi/wide_period,min(flux_array),0.15,1,0]

    return guess, upper_bound, bottom_bound

def plot_fq_fit(flux_array:ndarray,f01_array:ndarray,target_q:str,popt:dict,plot:bool=False,fig_path:str=''):
    flux_extend = linspace(flux_array[0],flux_array[-1],1000)
    plt.scatter(flux_array,f01_array,label='exp data')
    plt.plot(flux_extend, FqEqn(flux_extend,*popt),'r-',label=' Ec=%5.3f, Ej_sum=%5.3f, d=%5.3f' % tuple(popt[2:])) #, 'r-',label='Ej_sum=%5.3f, d=%5.3f, Ec=%5.3f' % tuple(popt)
    plt.legend()
    plt.title(f"{target_q} Flux vs. Fq fitting")
    plt.xlabel("Flux (V)")
    plt.ylabel("XY Frequency (GHz)")
    if fig_path != '':
        plt.savefig(fig_path)
    if plot:
        plt.show()
    else:
        plt.close()

def FitErrorFilter(eqn:Callable,eqn_paras:dict,exp_x_ary:ndarray,exp_y_ary:ndarray,threshold:float=1.5):
    """
    Statics for the distances along the same x axis between exp_y_data and fitting curve. Throw away the larger deviation. 
    """
    fit_y_ary = eqn(exp_x_ary,*eqn_paras)
    distances  = []
    # compute distances along a same x axis between fit curve and exp data
    for idx in range(fit_y_ary.shape[0]):
        distances.append(abs(fit_y_ary[idx]-exp_y_ary[idx]))
    away_bound = mean(array(distances)) + threshold*std(array(distances))
    
    throwAwayCounts = 0
    pass_exp_x = []
    pass_exp_y = []
    for idx in range(exp_y_ary.shape[0]):
        if distances[idx] < away_bound:
            pass_exp_y.append(exp_y_ary[idx])
            pass_exp_x.append(exp_x_ary[idx])
        else:
            throwAwayCounts += 1

    return array(pass_exp_x), array(pass_exp_y), throwAwayCounts

def fq_fit(QD:QDmanager,data2fit_path:str,target_q:str,plot:bool=True,savefig_path:str='',saveParas:bool=False, FitFilter_threshold:float=1.0):
    
    period = QD.Fluxmanager.get_PeriodFor(target_q)
    offset = QD.Fluxmanager.get_sweetBiasFor(target_q)

    flux, f01 = read_fq_data(data2fit_path)
    guess, upper_bound, bottom_bound = set_fitting_paras(period,offset,flux)
    popt, pcov = curve_fit(FqEqn, flux, f01,p0=guess,bounds=(bottom_bound,upper_bound))

    # try filter and fit again
    previous_throwCounts = 1e16 # For initialize, larger is better
    while True:
        advan_flux, advan_f01, thrownCounts = FitErrorFilter(FqEqn, popt, flux, f01, FitFilter_threshold)
        advan_popt, advan_pcov = curve_fit(FqEqn, advan_flux, advan_f01,p0=guess,bounds=(bottom_bound,upper_bound))
        print(thrownCounts)
        if thrownCounts == 0:
            break    
        else:
            # if previous_throwCounts >= thrownCounts:
            #     """ previous time throw more data and this is normal. use this results to continue """
            #     previous_throwCounts = thrownCounts
            #     if FitFilter_threshold  >= 0.5:
            #         FitFilter_threshold -= 0.5
            #         flux = advan_flux
            #         f01 = advan_f01
            #         popt = advan_popt
            #     else:
            #         print("No data had been thrown away !")
            #         break 
            # else:
            #     """ this time throw morw data, we take previous fitting results """
            #     advan_flux = flux
            #     advan_f01 = f01
            #     advan_popt = popt
            #     break
            if mean(sqrt(diag(advan_pcov))) <= mean(sqrt(diag(pcov))):
                if FitFilter_threshold  >= 0.5:
                    FitFilter_threshold += 0.5
                    flux = advan_flux
                    f01 = advan_f01
                    popt = advan_popt
                    pcov = advan_pcov
                else:
                    print("FitFilter threshold zeroed !")
                    break 
            else:
                advan_flux = flux
                advan_f01 = f01
                advan_popt = popt
                advan_pcov = pcov


    if savefig_path != '':
        savefig_path = os.path.join(savefig_path,f"{QD.Identity.split('#')[0]}{target_q}_FluxFqFit.png")
    
    if saveParas:
        QD.Fluxmanager.save_qubFittingParas_for(target_q,*advan_popt)
    
    plot_fq_fit(advan_flux, advan_f01,target_q,advan_popt,plot,savefig_path)
    print("Fitting completed!")



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from numpy import array, pi, linspace
    qd_path = 'Modularize/QD_backup/2024_3_21/DR2#171_SumInfo.pkl'
    json_path = 'Modularize/Meas_raw/2024_3_21/DR2q0_FluxFqFIT_H20M43S51.json'
    
    q = json_path.split("/")[-1].split("_")[0][-2:]
    QD_agent = QDmanager(qd_path)
    QD_agent.QD_loader()
    pic_parentpath = os.path.join(Data_manager().get_today_picFolder())
    fq_fit(QD=QD_agent,data2fit_path=json_path,target_q=q,plot=False,savefig_path=pic_parentpath,saveParas=False)

    