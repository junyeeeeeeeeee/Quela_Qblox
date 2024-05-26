from numpy import asarray, ndarray, array, sqrt, sort, mean, std, pi, cos, sin, diag, linspace
import sympy as sp
from scipy.optimize import curve_fit
import json, os
from typing import Callable
import matplotlib.pyplot as plt


def find_nearest(ary:ndarray, near_target:float):
    """ find the element  which is closest to the given near_target in the given array"""
    ary = asarray(ary)
    idx = (abs(ary - near_target)).argmin()
    return float(ary[idx])

def IQ_data_dis(I_data:ndarray,Q_data:ndarray,ref_I:float,ref_Q:float):
    Dis= sqrt((I_data-ref_I)**2+(Q_data-ref_Q)**2)
    return Dis   

def FqEqn(x,a,b,Ec,coefA,d):
    """
    a ~ period, b ~ offset, 
    """
    return sqrt(8*coefA*Ec*sqrt(cos(a*(x-b))**2+d**2*sin(a*(x-b))**2))-Ec

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

    upper_bound = [f*1.1,offset*1.1,0.32,100,1] #[a, b, Ec, Ej_sum, d]
    bottom_bound = [f*0.9,offset*0.9,0.28,1,0]

    return guess, upper_bound, bottom_bound


def mag_repalce_origin(I_ary:ndarray,Q_ary:ndarray,ref_IQ,qblox_rotation:bool=False):
    """
    For Qblox system, data needs to be transpose and flip, turn on the `qblox_rotation`.
    """
    if qblox_rotation:
        displaced_magnitude = array(IQ_data_dis(I_data=I_ary,Q_data=Q_ary,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])).transpose()
    else:
        displaced_magnitude = IQ_data_dis(I_data=I_ary,Q_data=Q_ary,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])

    return displaced_magnitude

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

    return array(filtered)


def data2plot(XYF_array:ndarray,z_array:ndarray,I_array:ndarray,Q_array:ndarray,specified_refIQ:list,filter2D_threshold:float=1.0,qblox:bool=False,plot_scatter:bool=False):
    
    mag = mag_repalce_origin(I_array,Q_array,ref_IQ=specified_refIQ,qblox_rotation=qblox)
    extracted = sortAndDecora(z_array,XYF_array,mag,threshold=filter2D_threshold)
    x_wanted, y_wanted, mag_wanted = extracted[:,0], extracted[:,1], extracted[:,-1]
    if plot_scatter:
        plot_HeatScat(mag=mag,x_heat_ary=z_array,x_scat_ary=x_wanted,y_heat_ary=XYF_array,y_scat_ary=y_wanted)
    return x_wanted, y_wanted, mag_wanted


### plot 
def plot_HeatScat(mag,x_heat_ary,y_heat_ary,x_scat_ary,y_scat_ary,fit_scat_ary:ndarray=array([])):
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
        title = "Flux vs. Transition frequency",
        xaxis=dict(title='Flux (V)'),
        yaxis=dict(title='XY Frequency (GHz)'),
        showlegend=True,
        legend={'x': .9, 'y': .0001}
    )
    fig = go.Figure(data=data, layout=layout)
    fig.show()


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


### filters

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

def mag_static_filter(x_array:ndarray,y_array:ndarray,mag_array:ndarray,threshold:float=3.0):
    x_choosen = []
    y_choosen = []
    while True:
        bottom_line = mean(mag_array) - threshold*std(mag_array)
        for idx in range(x_array.shape[0]):
            if mag_array[idx] > bottom_line:
                x_choosen.append(x_array[idx])
                y_choosen.append(y_array[idx])
        if len(x_choosen) > 0 :
            break
        else:
            if threshold != 0.5:
                threshold -= 0.5
            else:
                print(f"No peaks are found in this Z interval: {x_array[0]}~{x_array[-1]} V")
    
    return x_choosen, y_choosen

def FitErrorFilter(eqn:Callable,eqn_paras:dict,exp_x_ary:ndarray,exp_y_ary:ndarray,threshold:float=1.5):
    """
    Statics for the distances along the same x axis between exp_y_data and fitting curve. Throw away the larger deviation. 
    """
    fit_y_ary:ndarray = eqn(exp_x_ary,*eqn_paras)
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



### fitting flux vs collected xyf
def fq_fit(period:float,offset:float,data2fit_path:str,target_q:str,plot:bool=True,savefig_path:str='', FitFilter_threshold:float=1.0):

    flux, f01 = read_fq_data(data2fit_path)
    original_datapoints = flux.shape[0]
    guess, upper_bound, bottom_bound = set_fitting_paras(period,offset,flux,0.3)
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

    
    plot_fq_fit(advan_flux, advan_f01,target_q,advan_popt,plot,savefig_path)
    print("Fitting completed!")

    return advan_popt


### get wantted z by fitting results
def get_biasWithFq_from(fitting_popts:list,target_fq_Hz:float, flux_guard:float=0.4):
    """
    After we fit the tarnsition freq vs bias, we can get the bias according to the given `target_fq_Hz` for the `target_q`.\n
    ### The given `target_fq_Hz` should unit in Hz.\n
    ### The given `fitting_popts` should follow the order: [a, b, Ec, Ej_sum, d]
    Return the bias unit in V.
    """
    
    z = sp.Symbol('z',real=True)
    
    a,b,Ec,Ej_sum,d = fitting_popts[0], fitting_popts[1], fitting_popts[2], fitting_popts[3], fitting_popts[4]
    to_solve = sp.sqrt(8*Ej_sum*Ec*sp.sqrt(sp.cos(a*(z-b))**2+d**2*sp.sin(a*(z-b))**2))-Ec - target_fq_Hz*1e-9
    candidators = array(sp.solvers.solve(to_solve, z))
    
    if candidators.shape[0] == 0:
        answer = 'n'
        print(f"Can NOT find a bias makes the fq @ {target_fq_Hz*1e-9} GHz !")
    else:
        answer = find_nearest(candidators, 0) if find_nearest(candidators, 0) < flux_guard else flux_guard

    return answer



if __name__ == "__main__":
    # required inputs
    period = 0.0 # flux period 
    xyf = array([]) # sweep xyf array
    z = array([]) # sweep z array offset on sweet spot
    ref_z = 0 # sweet spot bias
    ii = array([]) # I chennel signal array
    qq = array([]) # Q chennel signal array
    ref_IQ = [0,0] # ground state center I,Q
    json_path = "" # path to save the filltered data points
    peak_threshold = 2.0

    i_want_fq_at = 4e9

    # program start here

    x2static, y2static, mag2static = data2plot(xyf,z+ref_z,ii,qq,specified_refIQ=ref_IQ,qblox=False,filter2D_threshold=peak_threshold,plot_scatter=False)
    

    x2fit, y2fit = mag_static_filter(array(x2static),array(y2static),array(mag2static))
    data2fit = {"x":x2fit,"y":y2fit}
    
    with open(json_path, "w") as json_file:
        json.dump(data2fit, json_file)

    fitting_paras = fq_fit(period, ref_z, json_path,target_q="q",savefig_path="",saveParas=True,plot=False) 

    # get z-bias according to the wantted fq
    wantted_z = get_biasWithFq_from(fitting_paras, i_want_fq_at)