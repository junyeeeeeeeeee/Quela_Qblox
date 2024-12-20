import os, json, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from typing import Callable
from numpy import flip, pi, linspace, array, sqrt, std, median, mean, sort, diag, sign, absolute, where, argmax
from qblox_drive_AS.support import QDmanager, Data_manager
from qblox_drive_AS.support.Pulse_schedule_library import IQ_data_dis
from numpy import ndarray, cos, sin, deg2rad, real, imag, transpose, abs, convolve, ones
from scipy.optimize import curve_fit
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

def smooth_filter(raw_data:ndarray)->ndarray:
    
    degree = 2
    poly_features = PolynomialFeatures(degree=degree)
    x_poly = poly_features.fit_transform(raw_data.reshape(-1, 1))

    # Fit the polynomial to the noisy data
    model = LinearRegression()
    model.fit(x_poly, raw_data)

    # Predict values and calculate residuals
    predictions = model.predict(x_poly)
    residuals = raw_data - predictions

    # Identify points that are likely part of the parabola signal
    threshold = mean(abs(residuals))+0.5*std(abs(residuals)) # Adjust based on noise level

    signal_indices = where(abs(residuals) < threshold)[0]
    extracted_signal = raw_data[signal_indices]
    print(extracted_signal)
    return signal_indices, extracted_signal

def remove_outliers_with_window(data, window_length, m=2):
    filtered_data = []
    original_indices = []

    # Split data into segments with the specified window length
    for i in range(0, len(data), window_length):
        window_data = data[i:i + window_length]
        
        # Compute mean and std within each window
        avg = median(window_data)
        std_dev = std(window_data)

        # Identify inliers in the window
        mask = (abs(window_data - avg) / std_dev) < m
        filtered_data.extend(window_data[mask])  # Add inliers to the result
        original_indices.extend(where(mask)[0] + i)  # Track original indices

    return array(filtered_data), array(original_indices)

def remove_outlier_after_fit(fit_func:callable,x_data:ndarray,y_data:ndarray,threshold:float=2)->tuple[ndarray,ndarray,tuple]:
    paras, _ = curve_fit(fit_func,x_data,y_data)
    fit_y = fit_func(x_data,*paras)
    fit_raw_distance = absolute(array(fit_y)-array(y_data))
    _, idxes = remove_outliers_with_window(fit_raw_distance,int(fit_raw_distance.shape[0]/2),threshold)
    filtered_x = []
    filtered_y = []
    for i in idxes:
        filtered_x.append(x_data[i])
        filtered_y.append(y_data[i])
    
    want_paras, _ = curve_fit(fit_func,array(filtered_x),array(filtered_y))  

    return array(filtered_x), array(filtered_y), want_paras 






def plot_QbFlux(Qmanager:QDmanager, nc_path:str, target_q:str):
    if target_q in Qmanager.refIQ:
        ref = Qmanager.refIQ[target_q]
    else:
        ref = [0,0]
    
    # plot flux-qubit 
    f,z,i,q = convert_netCDF_2_arrays(nc_path)
    amp = array(sqrt((i-array(ref)[0])**2+(q-array(ref)[1])**2)).transpose()
    print(amp.shape)
    fig, ax = plt.subplots(figsize=(12,8))
    ax:plt.Axes
    c = ax.pcolormesh(z, f*1e-9, amp, cmap='RdBu')
    ax.set_xlabel("Flux Pulse amp (V)", fontsize=20)
    ax.set_ylabel("Driving Frequency (GHz)", fontsize=20)
    fig.colorbar(c, ax=ax, label='Contrast (V)')
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    plt.tight_layout()
    plt.show()

def plot_QbFlux_multiVersn(Qmanager:QDmanager,qubit:str,ds:xr.Dataset,fit_func:callable=None, save_pic_path:str=None, filter_outlier:bool=False):
    fitting_packs = {}
    # plot flux-qubit 
    if qubit in ds.attrs["ref_z"].split("_"):
        ref_z = float(ds.attrs["ref_z"].split("_")[int(ds.attrs["ref_z"].split("_").index(qubit))+1])
    else:
        ref_z = 0
    f,z,i,q = convert_netCDF_2_arrays(ds)
    amp = array(sqrt((i-array(Qmanager.refIQ[qubit])[0])**2+(q-array(Qmanager.refIQ[qubit])[1])**2))
    
    fig, ax = plt.subplots(figsize=(13,9))
    ax:plt.Axes
    c = ax.pcolormesh(z, f*1e-9, amp.transpose())
    if fit_func is not None:
        # try:
        fit_z = []
        fit_f = []
        def parabola(x,a,b,c):
            return a*array(x)**2+b*array(x)+c
        for z_idx, a_z_data in enumerate(amp):
            if min(f) <= fit_func(a_z_data,f).attrs['f01_fit'] and fit_func(a_z_data,f).attrs['f01_fit'] <= max(f):
                fit_f.append(fit_func(a_z_data,f).attrs['f01_fit']*1e-9) # GHz
                fit_z.append(z[z_idx])
        

        if not filter_outlier:
            paras, _ = curve_fit(parabola,fit_z,fit_f)
            plt.scatter(fit_z,fit_f,marker='*',c='black')
            plt.plot(z,parabola(z,*paras),'red')
        else:
            
            filtered_z, filtered_f, paras = remove_outlier_after_fit(parabola,fit_z,fit_f)
            plt.plot(fit_z,parabola(fit_z,*paras),'red')
            plt.scatter(fit_z,fit_f,marker='*',c='black')
            plt.scatter(filtered_z,filtered_f,marker='*',c='red')

        if min(z)<=float(-paras[1]/(2*paras[0])) and float(-paras[1]/(2*paras[0])) <= max(z):
            plt.vlines(float(-paras[1]/(2*paras[0])),min(f)*1e-9,max(f)*1e-9,colors='red',linestyles='--')
        
        fitting_packs["sweet_bias"] = ref_z + float(-paras[1]/(2*paras[0])) # offset + z_pulse_amp
        fitting_packs["xyf"] = float(parabola(float(-paras[1]/(2*paras[0])),*paras))*1e9
        fitting_packs["parabola_paras"] = list(paras)
        # except:
        #     print("Passed a wrong fitting for flux-qubit spectroscopy")
    
    ax.set_xlabel("Flux Pulse amp (V)", fontsize=20)
    ax.set_ylabel("Driving Frequency (GHz)", fontsize=20)
    fig.colorbar(c, ax=ax, label='Contrast (V)')
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    # plt.tight_layout()
    if len(list(fitting_packs.keys())) != 0:
        plt.title(f"{qubit} XYF={round(fitting_packs['xyf']*1e-9,3)} GHz with z_pulse amp={round(float(-paras[1]/(2*paras[0])),3)} V")
    if save_pic_path is None:
        plt.show()
    else:
        plt.savefig(os.path.join(save_pic_path,f"{qubit}_fluxQubitSpectro_{ds.attrs['execution_time'] if 'execution_time' in list(ds.attrs) else Data_manager().get_time_now()}.png"))
        plt.close()

    return fitting_packs

def plot_QbFlux_iq(Qmanager:QDmanager, nc_path:str, target_q:str):
    if target_q in Qmanager.refIQ:
        ref = Qmanager.refIQ[target_q]
    else:
        ref = [0,0]
    # plot flux-qubit 
    f,z,i,q = convert_netCDF_2_arrays(nc_path)
    
    fig, ax = plt.subplots(1,2,figsize=(12,8))
    ax0:plt.Axes=ax[0]
    c = ax0.pcolormesh(z, f*1e-9, i.transpose(), cmap='RdBu')
    ax0.set_xlabel("Flux Pulse amp (V)", fontsize=20)
    ax0.set_ylabel("Frequency (GHz)", fontsize=20)
    fig.colorbar(c, ax=ax0, label='Contrast (V)')
    ax0.xaxis.set_tick_params(labelsize=18)
    ax0.yaxis.set_tick_params(labelsize=18)
    ax0.set_title("Inphase")
    ax1:plt.Axes=ax[1]
    d = ax1.pcolormesh(z, f*1e-9, q.transpose(), cmap='RdBu')
    ax1.set_xlabel("Flux Pulse amp (V)", fontsize=20)
    ax1.set_ylabel("Frequency (GHz)", fontsize=20)
    fig.colorbar(d, ax=ax1, label='Contrast (V)')
    ax1.xaxis.set_tick_params(labelsize=18)
    ax1.yaxis.set_tick_params(labelsize=18)
    ax1.set_title("Quadrature")
    plt.tight_layout()
    plt.show()


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

    # filter peaks idx by a threshold satisfies that peak > mean + threshold*std
    while True:
        f_idx, z_idx, XL_mag = filter_2D(raw_mag,threshold)
        if len(f_idx)==0 and len(z_idx)==0:
            if threshold != 0:
                threshold-=0.5
            else:
                print(f"This interval can't find the trend : XYF={raw_XYF[0]}~{raw_XYF[-1]}")
        else:
            break
    
    # extract the values by the filtered idx
    extracted = [] # [z, xyf, mag]
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

    return array(filtered)

def convert_netCDF_2_arrays(nc:(str|xr.Dataset)):
    """
    1.  If the arg `nc` was given with a `str` type, it will be treat as a nc file path.
    2.  The other type of the arg is `xr.Dataset`, it will be treat as a Dataset.
    this file should be magnitude and phase, not I and Q.
    """
    if isinstance(nc,str):
        print("Get the nc file path...")
        dataset_processed = dh.to_gridded_dataset(xr.open_dataset(nc))
    elif isinstance(nc,xr.Dataset):
        print("Get the nc Dataset...")
        dataset_processed = dh.to_gridded_dataset(nc)
    else:
        raise TypeError(f"The type of the given was not correct, it should be str or xr.Dataset")
    
    XYF = dataset_processed["x0"].to_numpy()
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
        displaced_magnitude = array(IQ_data_dis(I_data=I_ary,Q_data=Q_ary,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])).transpose()
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
    return sqrt((1000**2)*abs(dressF-bareF)*abs(bareF-Fq)/(sqrt(bareF*Fq)**2))

# predict fq
def calc_fq_g_excluded(coef_inG:float,fdress:float,fbare:float):
    """
    After we got the coef in coupling g, use dispersive relation to calc fq with the given args.
    """
    return fbare*((fdress-fbare)/((fbare*(coef_inG/1000)**2)+(fdress-fbare)))



def data2plot(XYF_array:ndarray,z_array:ndarray,I_array:ndarray,Q_array:ndarray,specified_refIQ:list,filter2D_threshold:float=1.0,qblox:bool=False,q:str='',plot_scatter:bool=False):
    
    mag = mag_repalce_origin(I_array,Q_array,ref_IQ=specified_refIQ,qblox_rotation=qblox)
    extracted = sortAndDecora(z_array,XYF_array,mag,threshold=filter2D_threshold)
    x_wanted, y_wanted, mag_wanted = extracted[:,0], extracted[:,1], extracted[:,-1]
    if plot_scatter:
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


def set_fitting_paras(period:float,offset:float,Ec_guess_GHz:float=0.21,Ej_sum_guess_GHz:float=20.0,squid_ratio_guess:float=0.5):
    """
    There are 5 paras in Fq eqn, give the initial guess and the fitting constrains for curve fit.\n
    Return guess, upper_bound, bottom_bound.
    """
    f = pi/period
    b = offset
    guess = (f,b,Ec_guess_GHz,Ej_sum_guess_GHz,squid_ratio_guess) #[a, b, Ec, Ej_sum, d]

    acceptable_loss = 30 / 100

    if f < 0:
        up_f = f * (1-acceptable_loss)
        lo_f = f * (1+acceptable_loss)
    else:
        up_f = f * (1+acceptable_loss)
        lo_f = f * (1-acceptable_loss)

    if b < 0:
        up_b = b * (1-acceptable_loss)
        lo_b = b * (1+acceptable_loss)
    else:
        up_b = b * (1+acceptable_loss)
        lo_b = b * (1-acceptable_loss)
    
    upper_bound =  [up_f, up_b, 0.25, 100, 1] #[a, b, Ec, Ej_sum, d]
    bottom_bound = [lo_f, lo_b, 0.19,  1,  0]

    return guess, upper_bound, bottom_bound

def plot_fq_fit(flux_array:ndarray,f01_array:ndarray,target_q:str,popt:dict,plot:bool=False,fig_path:str=''):
    plt.close()
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

    exam_dict = {"meas":{"z":list(flux_array), "xyf":list(f01_array)},"fit":{"z":list(flux_extend),"xyf":list(FqEqn(flux_extend,*popt))}}
    with open(os.path.join(os.path.split(fig_path)[0],"measANDfit_points.json"), "w") as record_file:
        json.dump(exam_dict,record_file)

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
    original_datapoints = flux.shape[0]
    guess, upper_bound, bottom_bound = set_fitting_paras(period,offset,0.22)
    popt, pcov = curve_fit(FqEqn, flux, f01,p0=guess,bounds=(bottom_bound,upper_bound))

    # try filter and fit again
    while True:
        advan_flux, advan_f01, thrownCounts = FitErrorFilter(FqEqn, popt, flux, f01, FitFilter_threshold)
        advan_popt, advan_pcov = curve_fit(FqEqn, advan_flux, advan_f01,p0=guess,bounds=(bottom_bound,upper_bound))
        if thrownCounts == 0:
            print("No data points can be thrown, break!")
            break    
        else:
            previous_err = mean(sqrt(diag(pcov))[:4])
            new_err = mean(sqrt(diag(advan_pcov))[:4])
            if new_err <= previous_err:
                if new_err > previous_err/4 :
                    flux = advan_flux
                    f01 = advan_f01
                    popt = advan_popt
                    pcov = advan_pcov
                    if flux.shape[0] < original_datapoints/4:
                        break
                else:
                    print("FitFilter_threshold touchecd bottom !")
                    break 
            else:
                advan_flux = flux
                advan_f01 = f01
                advan_popt = popt
                advan_pcov = pcov
                print("New fitting error had increased, break!")
                break


    if savefig_path != '':
        savefig_path = os.path.join(savefig_path,f"{QD.Identity.split('#')[0]}{target_q}_FluxFqFit.png")
    
    if saveParas:
        QD.Fluxmanager.save_qubFittingParas_for(target_q,*advan_popt)
    
    plot_fq_fit(advan_flux, advan_f01,target_q,advan_popt,plot,savefig_path)
    print("Fitting completed!")



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from numpy import array, pi, linspace
    # qd_path = 'Modularize/QD_backup/2024_5_15/DR3#13_SumInfo.pkl'
    # json_path = 'Modularize/Meas_raw/2024_5_15/DR3q0_FluxFqFIT_H17M21S59.json'
    
    # q = os.path.split(json_path)[-1].split("_")[0][-2:]
    # QD_agent = QDmanager(qd_path)
    # QD_agent.QD_loader()
    # print(QD_agent.Fluxmanager.get_bias_dict()["q1"])
    # pic_parentpath = os.path.join(Data_manager().get_today_picFolder())
    # fq_fit(QD=QD_agent,data2fit_path=json_path,target_q=q,plot=True,savefig_path='',saveParas=False,FitFilter_threshold=2.5)
    file = 'Modularize/Meas_raw/2024_8_12/DR4q0_Flux2tone_H12M32S5.nc'
    f, z, i, q = convert_netCDF_2_arrays(file)
    data2plot(f,z,i,q,[0,0],q='q0',qblox=True,plot_scatter = True)

