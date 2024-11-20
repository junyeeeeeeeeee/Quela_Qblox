
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import numpy as np
import xarray as xr
from scipy import special
from scipy.integrate import quad
from scipy.signal import butter,sosfiltfilt
from lmfit import Model,Parameter 
from quantify_scheduler.enums import BinMode
from quantify_scheduler.backends.graph_compilation import SerialCompiler
from quantify_scheduler.schedules.schedule import Schedule
from quantify_scheduler.operations.acquisition_library import SSBIntegrationComplex,Trace
from quantify_scheduler.operations.pulse_library import (IdlePulse,SetClockFrequency,SquarePulse,DRAGPulse,GaussPulse,SoftSquarePulse,NumericalPulse)
from quantify_scheduler.device_under_test.quantum_device import QuantumDevice
from quantify_scheduler.operations.gate_library import Reset, Measure
from quantify_scheduler.resources import ClockResource, BasebandClockResource
from quantify_scheduler.helpers.collections import find_port_clock_path
from qblox_drive_AS.support.WaveformCtrl import XY_waveform, s_factor, half_pi_ratio

""" Global pulse settings """
electrical_delay:float = 280e-9


def FlatTopGaussianPulse(Du:float,amp:float,s_factor:int=8,sampling_rate:float=4e-9):
    t_samples = np.arange(0,Du,sampling_rate)
    g_u_t = np.arange(0,Du/s_factor,sampling_rate)
    g_d_t = np.arange((s_factor-1)*Du/s_factor,Du,sampling_rate)
    flat_t= np.arange(0,Du,sampling_rate)[int(t_samples.shape[0]/s_factor):int(t_samples.shape[0]*(s_factor-1)/s_factor)]
    gaussian_up =  gauss_func(g_u_t,np.mean(g_u_t),np.max(g_u_t)/s_factor,amp)[:int(g_u_t.shape[0]/2)]
    gaussian_dn =  gauss_func(g_d_t,np.mean(g_d_t),(np.max(g_d_t)-np.min(g_d_t))/s_factor,amp)[int(g_u_t.shape[0]/2):]
    flatten = np.array([amp for i in range(flat_t.shape[0])])

    env_sample = np.hstack([gaussian_up,flatten,gaussian_dn])
    time_sample = np.hstack([g_u_t[:int(g_u_t.shape[0]/2)],flat_t,g_d_t[int(g_u_t.shape[0]/2):]])

    return time_sample, env_sample
    
 

#%% IQ displacement
def IQ_data_dis(I_data:np.ndarray,Q_data:np.ndarray,ref_I:float,ref_Q:float):
    Dis= np.sqrt((I_data-ref_I)**2+(Q_data-ref_Q)**2)
    return Dis    
 
def XY_waveform_controller(pulse_sche:Schedule, amp:float, duration:float, phase:float, q:str, rel_time:float, ref_op:Schedule, waveform:str, ref_pt:str="start", sigma_factor:float=4,transition:str='.01',**kwargs):
    """
    waveform switch\n
    sigma_factor determines the relations between pulse sigma and pulse duration is `sigma = duration/sigma_factor`.
    ### * kwargs: *
    1. if use DRAG, there should ne a argument named 'drag_amp'. It specifies the derivative part amplitude in DRAG as default is the same as amp.
    """
    match waveform.lower():
        case 'drag':
    
            if 'drag_amp' in list(kwargs.keys()):
                diff_amp = kwargs['drag_amp']
            else:
                diff_amp = amp
            return pulse_sche.add(DRAGPulse(G_amp=amp, D_amp=diff_amp, duration= duration, phase=phase, port=q+":mw", clock=q+transition,sigma=duration/sigma_factor),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
        case 'gauss':
            
            return pulse_sche.add(GaussPulse(G_amp=amp, phase=phase,duration=duration, port=q+":mw", clock=q+transition,sigma=duration/sigma_factor),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
        case _:
            pass

#%% fit function
def fft_oscillation_guess(data: np.ndarray, t: np.ndarray):
    amp = np.fft.fft(data)[: len(data) // 2] #use positive frequency
    freq = np.fft.fftfreq(len(data), t[1] - t[0])[: len(amp)]
    amp[0] = 0  # Remove DC part 
    power = np.abs(amp)
    f_guess = abs(freq[power == max(power)][0])
    phase_guess = 2 * np.pi - (2 * np.pi * t[data == max(data)] * f_guess)[0]
    return f_guess, phase_guess
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
def rot(I,Q,angle):
    sin=np.sin(angle)
    cos=np.cos(angle)
    return I*cos+Q*sin, -I*sin+Q*cos

def Rabi_func(x,A,f,offset):
    return A*np.cos(2*np.pi*f*x+np.pi)+offset
def T1_func(D,A,T1,offset):
    return A*np.exp(-D/T1)+offset
def Ramsey_func(D,A,T2,f,phase,offset):
    return A*np.exp(-D/T2)*np.cos(2*np.pi*f*D+phase)+offset
def Cosine_func(D,A,f,phase,offset):
    return A*np.cos(2*np.pi*f*D+phase)+offset
def Loren_func(x,x0,gamma,A,base):
    return (A/np.pi)*((gamma/2)/((x-x0)**2+(gamma/2)**2))+base
def gauss_func(x,c,sigma,A):
    return A*np.exp(-(x-c)**2/2/sigma**2)    
def gauss2d_func(I,Q,c_I,c_Q,sigma,A):
    return A*np.exp(-(I-c_I)**2/2/sigma**2)*np.exp(-(Q-c_Q)**2/2/sigma**2)
def bigauss2d_func(I,Q,cg_I,cg_Q,sigma,Ag,ce_I,ce_Q,Ae):
    return gauss2d_func(I,Q,cg_I,cg_Q,sigma,Ag)+gauss2d_func(I,Q,ce_I,ce_Q,sigma,Ae)
def Relax_cal(inte_i,tf,T1):
    tau= tf - inte_i
    Relax= 1+(T1/tau)*(np.exp(-(inte_i+tau)/T1)-np.exp(-inte_i/T1))
    return Relax
Rabi_model = Model(Rabi_func)
cosine_model = Model(Cosine_func)
T1_func_model = Model(T1_func)
Ramsey_func_model = Model(Ramsey_func)
Loren_func_model = Model(Loren_func)
gauss2d_func_model = Model(gauss2d_func, independent_vars=['I', 'Q'])
bigauss2d_func_model = Model(bigauss2d_func, independent_vars=['I', 'Q'])

def cos_fit_analysis(data:np.ndarray,freeDu:np.ndarray):
    f_guess,phase_guess= fft_oscillation_guess(data,freeDu)
    f_guess_=Parameter(name='f', value=f_guess , min=0, max=f_guess*10)
    result = cosine_model.fit(data,D=freeDu,A=abs(min(data)+max(data))/2,f=f_guess_,phase=phase_guess, offset=np.mean(data))

    A_fit= result.best_values['A']
    f_fit= result.best_values['f']
    phase_fit= result.best_values['phase']
    offset_fit= result.best_values['offset']
    para_fit= np.linspace(freeDu.min(),freeDu.max(),50*len(data))

    fitting= Cosine_func(para_fit,A_fit,f_fit,phase_fit,offset_fit)
    paras = [A_fit,f_fit,phase_fit,offset_fit]
    return xr.Dataset(data_vars=dict(data=(['freeDu'],data),fitting=(['para_fit'],fitting)),coords=dict(freeDu=(['freeDu'],freeDu),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="xyfcali",f=f_fit,phase=phase_fit,coefs=paras))



def T1_fit_analysis(data:np.ndarray,freeDu:np.ndarray,T1_guess:float=10*1e-6,return_error:bool=False):
    offset_guess= data[-1]
    a_guess = np.max(data)-offset_guess if data[-1] > offset_guess else np.min(data)-offset_guess
    result = T1_func_model.fit(data,D=freeDu,A=np.max(data)-offset_guess,T1=T1_guess,offset=offset_guess)
    
    A_fit= result.best_values['A']
    T1_fit= result.best_values['T1']
    offset_fit= result.best_values['offset']
    para_fit= np.linspace(freeDu.min(),freeDu.max(),50*len(data))
    fitting= T1_func(para_fit,A_fit,T1_fit,offset_fit)
    if not return_error:
        return xr.Dataset(data_vars=dict(data=(['freeDu'],data),fitting=(['para_fit'],fitting)),coords=dict(freeDu=(['freeDu'],freeDu),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="T1",T1_fit=T1_fit))
    else:
        fit_error = float(result.covar[1][1])*1e6
        return xr.Dataset(data_vars=dict(data=(['freeDu'],data),fitting=(['para_fit'],fitting)),coords=dict(freeDu=(['freeDu'],freeDu),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="T1",T1_fit=T1_fit)), fit_error
def T2_fit_analysis(data:np.ndarray,freeDu:np.ndarray,T2_guess:float=10*1e-6,return_error:bool=False):
    
    f_guess,phase_guess= fft_oscillation_guess(data,freeDu)
    T2=Parameter(name='T2', value= T2_guess, min=0.1e-6, max=5*T2_guess) 
    up_lim_f= 30*1e6
    f_guess_=Parameter(name='f', value=f_guess , min=0, max=up_lim_f)
    result = Ramsey_func_model.fit(data,D=freeDu,A=abs(min(data)+max(data))/2,T2=T2,f=f_guess_,phase=phase_guess, offset=np.mean(data))
    if return_error:
        fit_error = float(result.covar[1][1])*1e6
    A_fit= result.best_values['A']
    f_fit= result.best_values['f']
    phase_fit= result.best_values['phase']
    
    T2_fit= result.best_values['T2']
    offset_fit= result.best_values['offset']
    para_fit= np.linspace(freeDu.min(),freeDu.max(),50*len(data))
    fitting= Ramsey_func(para_fit,A_fit,T2_fit,f_fit,phase_fit,offset_fit)
    if not return_error:
        return xr.Dataset(data_vars=dict(data=(['freeDu'],data),fitting=(['para_fit'],fitting)),coords=dict(freeDu=(['freeDu'],freeDu),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="T2",T2_fit=T2_fit,f=f_fit,phase=phase_fit))
    else:
        return xr.Dataset(data_vars=dict(data=(['freeDu'],data),fitting=(['para_fit'],fitting)),coords=dict(freeDu=(['freeDu'],freeDu),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="T2",T2_fit=T2_fit,f=f_fit,phase=phase_fit)), fit_error
def QS_fit_analysis(data:np.ndarray,f:np.ndarray):
    fmin = f.min()
    fmax = f.max()
    width_max = fmax-fmin
    delta_f = np.diff(f)  
    min_delta_f = delta_f[delta_f > 0].min()
    width_min = min_delta_f
    width_guess = np.sqrt(width_min*width_max) 
    A= np.pi * width_guess * (np.max(data)-np.mean(data))
    result = Loren_func_model.fit(data,x=f,x0= f[np.argmax(data)],gamma=width_guess,A=A,base=np.mean(data))
    f01_fit= result.best_values['x0']
    bw= result.best_values['gamma']
    A_fit= result.best_values['A']
    base_fit= result.best_values['base']
    para_fit= np.linspace(fmin,fmax,50*len(data))
    fitting= Loren_func(para_fit,f01_fit,bw,A_fit,base_fit)
    return xr.Dataset(data_vars=dict(data=(['f'],data),fitting=(['para_fit'],fitting)),coords=dict(f=(['f'],f),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="QS",f01_fit=f01_fit,bandwidth=bw))

def Rabi_fit_analysis(data:np.ndarray,samples:np.ndarray, Rabi_type:str):
    f_guess,phase_guess= fft_oscillation_guess(data,samples)
    result = Rabi_model.fit(data,x=samples,A=abs(min(data)+max(data))/2,f=f_guess, offset=np.mean(data))
    A_fit= result.best_values['A'] 
    f_fit= result.best_values['f']
    offset_fit= result.best_values['offset']
    pi_2= 1/(2*f_fit)
    para_fit= np.linspace(samples.min(),samples.max(),50*len(data))
    fitting= Rabi_func(para_fit,A_fit,f_fit,offset_fit)
    return xr.Dataset(data_vars=dict(data=(['samples'],data),fitting=(['para_fit'],fitting)),coords=dict(samples=(['samples'],samples),para_fit=(['para_fit'],para_fit)),attrs=dict(exper="Rabi",Rabi_type=Rabi_type,pi_2=pi_2))

def Single_shot_ref_fit_analysis(data:tuple):
    I,Q= np.array(data[0]),np.array(data[1]) 
    bins=401
    I_=np.linspace(I.min(),I.max(),bins)  
    Q_=np.linspace(Q.min(),Q.max(),bins)
    hist, xedges, yedges = np.histogram2d(I,Q, bins=(bins,bins), density=True)
    I_guess, Q_guess, sig_guess=np.mean(I),np.mean(Q),(np.std(I)+np.std(Q))/2
    c_I=Parameter(name='c_I', value= I_guess, min=0.9*I_guess, max=1.1*I_guess) 
    c_Q=Parameter(name='c_Q', value= Q_guess, min=0.9*Q_guess, max=1.1*Q_guess)
    sigma=Parameter(name='sigma', value= sig_guess, min=0.2*sig_guess, max=1.5*sig_guess)
    X,Y= np.meshgrid(I_,Q_)
    result= gauss2d_func_model.fit(hist.transpose(),I=X,Q=Y,c_I=c_I,c_Q=c_Q,sigma=sigma,A=Parameter(name='A',value=np.max(hist), min=0.1*np.max(hist)) )     
    c_I_fit=result.best_values['c_I']
    c_Q_fit=result.best_values['c_Q']
    sigma_fit=result.best_values['sigma']
    A_fit=result.best_values['A']
    fit_pack=[c_I_fit,c_Q_fit,sigma_fit]
    fitting= gauss2d_func(X,Y,c_I_fit,c_Q_fit,sigma_fit,A_fit)
    #print ('Total prob. =',np.sum(hist)*((max(xedges)-min(xedges))/bins*(max(yedges)-min(yedges))/bins))
    return dict(data=[I,Q],data_hist=hist,coords=[X,Y],fitting=fitting,fit_pack=fit_pack)

def Qubit_state_single_shot_fit_analysis(data:dict, T1:float,tau:float):
    Ig_data,Qg_data,Ie_data,Qe_data= np.array(data['g'][0]), np.array(data['g'][1]) ,np.array(data['e'][0]), np.array(data['e'][1]) 
    bins=401
    g_predict= Single_shot_ref_fit_analysis(data['g'])['fit_pack']
    Ig_guess, Qg_guess, sig_guess= g_predict[0],g_predict[1],g_predict[2]
    I_mixdata, Q_mixdata= np.hstack([Ig_data,Ie_data]), np.hstack([Qg_data,Qe_data])
    I_=np.linspace(I_mixdata.min(),I_mixdata.max(),bins)  
    Q_=np.linspace(Q_mixdata.min(),Q_mixdata.max(),bins)
    X,Y= np.meshgrid(I_,Q_)
    hist_Ie, edges_Ie= np.histogram(Ie_data, bins=bins, density=True)
    hist_Qe, edges_Qe= np.histogram(Qe_data, bins=bins, density=True)
    hist, xedges, yedges = np.histogram2d(I_mixdata, Q_mixdata, bins=(bins,bins), density=True)
    Ie_guess_idx, Qe_guess_idx=find_nearest(hist_Ie,np.max(hist_Ie)),find_nearest(hist_Qe,np.max(hist_Qe))
    Ie_guess, Qe_guess= edges_Ie[Ie_guess_idx],edges_Qe[Qe_guess_idx]
    #Parameter ini-guess
    cg_I=Parameter(name='cg_I', value= Ig_guess, min=0.9*Ig_guess, max=1.5*Ig_guess) 
    cg_Q=Parameter(name='cg_Q', value= Qg_guess, min=0.9*Qg_guess, max=1.5*Qg_guess)
    ce_I=Parameter(name='ce_I', value= Ie_guess, min=0.3*Ie_guess, max=3*Ie_guess) 
    ce_Q=Parameter(name='ce_Q', value= Qe_guess, min=0.3*Qe_guess, max=3*Qe_guess)
    sigma=Parameter(name='sigma', value= sig_guess, min=0.1*sig_guess, max=3*sig_guess)
    Ag=Parameter(name='Ag',value=np.max(hist), min=0.1*np.max(hist)) 
    Ae=Parameter(name='Ae',value=np.max(hist), min=0.1*np.max(hist)) 
    # mixed data fit
    result= bigauss2d_func_model.fit(hist.transpose(),I=X,Q=Y,cg_I=cg_I,cg_Q=cg_Q,sigma=sigma,Ag=Ag,ce_I=ce_I,ce_Q=ce_Q,Ae=Ae)
    cg_I_fit=result.best_values['cg_I']
    cg_Q_fit=result.best_values['cg_Q']
    ce_I_fit=result.best_values['ce_I']
    ce_Q_fit=result.best_values['ce_Q']
    sigma_fit=result.best_values['sigma']
    # displace + rotate
    angle= np.angle(ce_I_fit-cg_I_fit+(ce_Q_fit-cg_Q_fit)*1j)
    rot_e_center= rot(ce_I_fit-cg_I_fit,ce_Q_fit-cg_Q_fit,angle)
    rot_g_IQ= rot(Ig_data-cg_I_fit,Qg_data-cg_Q_fit,angle)
    rot_e_IQ= rot(Ie_data-cg_I_fit,Qe_data-cg_Q_fit,angle)
    Ig_data_new, Qg_data_new, Ie_data_new, Qe_data_new= rot_g_IQ[0],rot_g_IQ[1],rot_e_IQ[0],rot_e_IQ[1]
    # along single quadrature fit
    R= 10*sigma_fit #range_factor
    xmin, xmax = np.minimum(0,rot_e_center[0])-R,np.maximum(0,rot_e_center[0])+R
    ymin, ymax = np.minimum(0,rot_e_center[1])-R,np.maximum(0,rot_e_center[1])+R
    hist_r_g, xedges, yedges = np.histogram2d(Ig_data_new, Qg_data_new, bins=(bins,bins),range=[[xmin, xmax], [ymin, ymax]], density=True)
    hist_r_e, xedges, yedges = np.histogram2d(Ie_data_new, Qe_data_new, bins=(bins,bins),range=[[xmin, xmax], [ymin, ymax]], density=True)
    New_axe_g_hist= np.sum(hist_r_g.transpose(), axis=0)
    New_axe_e_hist= np.sum(hist_r_e.transpose(), axis=0)
    I_ro = np.linspace(xmin, xmax,bins)
    I_fit= np.linspace(xmin, xmax,bins*5)
    def Reduced_bimodal_func(x,Ag,Ae):
        return gauss_func(x,0,sigma_fit,Ag)+gauss_func(x,rot_e_center[0],sigma_fit,Ae)
    bimodal_func_model = Model(Reduced_bimodal_func)
    result_g= bimodal_func_model.fit(New_axe_g_hist,x=I_ro,Ag=Parameter(name='Ag',value=np.max(New_axe_g_hist), min=0.5*np.max(New_axe_g_hist)) ,Ae=Parameter(name='Ae',value=0.01*np.max(New_axe_g_hist), min=0))
    result_e= bimodal_func_model.fit(New_axe_e_hist,x=I_ro,Ag=Parameter(name='Ag',value=0.01*np.max(New_axe_e_hist), min=0),Ae=Parameter(name='Ae',value=np.max(New_axe_e_hist), min=0.5*np.max(New_axe_e_hist)))
    Agg_fit=result_g.best_values['Ag']
    Aeg_fit=result_g.best_values['Ae']
    Age_fit=result_e.best_values['Ag']
    Aee_fit=result_e.best_values['Ae']
    inter_point = rot_e_center[0]/2
    def G_gg(I):
        return Agg_fit*np.exp(-I**2/(2*sigma_fit**2))
    def G_eg(I):
        return Aeg_fit*np.exp(-(I-rot_e_center[0])**2/(2*sigma_fit**2))
    def G_ge(I):
        return Age_fit*np.exp(-I**2/(2*sigma_fit**2))
    def G_ee(I):
        return Aee_fit*np.exp(-(I-rot_e_center[0])**2/(2*sigma_fit**2))
    Ggg= quad(G_gg,-np.inf,np.inf)
    Geg= quad(G_eg,-np.inf,np.inf)
    Gge= quad(G_ge,-np.inf,np.inf)
    Gee= quad(G_ee,-np.inf,np.inf)
    overlap_gg= quad(G_gg,inter_point,np.inf)
    transi_eg= quad(G_eg,inter_point,np.inf)
    transi_ge= quad(G_ge,-np.inf,inter_point)
    overlap_ee= quad(G_ee,-np.inf,inter_point)
    overlap= ((overlap_gg[0]-overlap_gg[1])/(Ggg[0]-Ggg[1])+(overlap_ee[0]-overlap_ee[1])/(Gee[0]-Gee[1]))/2
    Peg= (transi_eg[0]-transi_eg[1])/(Ggg[0]-Ggg[1]+Geg[0]-Geg[1])+(overlap_gg[0]-overlap_gg[1])/(Ggg[0]-Ggg[1])
    Pge= (transi_ge[0]-transi_ge[1])/(Gee[0]-Gee[1]+Gge[0]-Gge[1])+(overlap_ee[0]-overlap_ee[1])/(Gee[0]-Gee[1])
    Thermal= (Geg[0]-Geg[1])/(Geg[0]-Geg[1]+Ggg[0]-Ggg[1])
    Relax= (Gge[0]-Gge[1])/(Gee[0]-Gee[1]+Gge[0]-Gge[1])-Thermal
    Relax_predict= Relax_cal(0,tau,T1)
    Pre_decay= Relax-Relax_predict
    D= rot_e_center[0]
    SNR= D/sigma_fit
    overlap_predict= (1/2)*(1-special.erf(np.sqrt(SNR**2/8)))
    F_s= 1-overlap
    F_g= 1-Peg
    F_e= 1-Pge
    F=1-(1/2)*(Peg+Pge)
    fit_pack= [rot_e_center,Agg_fit,Aeg_fit,Age_fit,Aee_fit,sigma_fit,New_axe_g_hist,New_axe_e_hist]
    error_pack= dict(D=D,sigma=sigma_fit,SNR=SNR,overlap=overlap,overlap_predict=overlap_predict,Peg=Peg,Pge=Pge,Thermal=Thermal,Relax=Relax,Relax_predict=Relax_predict,Pre_decay=Pre_decay, F_s=F_s,F_g=F_g,F_e=F_e,F=F)
    return dict(rot_IQdata=[rot_g_IQ,rot_e_IQ],I_ro=I_ro,I_fit=I_fit,fit_pack=fit_pack,error_pack=error_pack)

def Readout_F_opt_analysis(data:dict,f_samples:np.ndarray):
    mag_g,mag_e= np.array(data['g'][0]), np.array(data['e'][0]) 
    f_g,f_e= f_samples[find_nearest(mag_g,np.min(mag_g))],f_samples[find_nearest(mag_e,np.min(mag_e))]
    
    return dict(f_g=f_g, f_e=f_e,mag_g=mag_g,mag_e=mag_e,f_samples=f_samples)

def butter_lowpass(highcut, fs, order):
    sos= butter(order, highcut, fs=fs, btype='low',output='sos')
    return sos

def butter_lowpass_filter(data, highcut, fs, order):
    sos = butter_lowpass(highcut, fs, order=order)
    y = sosfiltfilt(sos, data)
    return y

def Trace_filtering(data:np.ndarray,fc:float):
    filter_data= butter_lowpass_filter(data,fc,1e9,40)
    return filter_data

def Digital_down_convert(IF_Idata:np.ndarray,IF_Qdata:np.ndarray,IF:float,time_array:np.ndarray):
    return IF_Idata*np.cos(2*np.pi*IF*time_array), IF_Qdata*np.sin(2*np.pi*IF*time_array)
    

#%% pulse library
def Spec_pulse(sche,amp,Du,q,ref_pulse_sche,freeDu, ref_point:str="start"):
    delay_c= -Du-freeDu
    return sche.add(SquarePulse(duration=Du,amp=amp, port=q+":mw", clock=q+".01"),rel_time=delay_c,ref_op=ref_pulse_sche,ref_pt=ref_point)

def X_theta(sche,amp,Du,q,ref_pulse_sche,freeDu):

    if Du!=0:
        delay_c= -Du-freeDu
        return XY_waveform_controller(sche,amp,Du,0,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt='start')
    else: pass

def Y_theta(sche,amp,Du,q,ref_pulse_sche,freeDu):
    if Du!=0:
        delay_c= -Du-freeDu
        return XY_waveform_controller(sche,amp,Du,90,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt='start')
    else: pass

def X_pi_2_p(sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu):
    amp= pi_amp[q]*half_pi_ratio
    delay_c= -pi_Du-freeDu
    return XY_waveform_controller(sche,amp,pi_Du,0,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt='start')

def Y_pi_2_p(sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu):
    amp= pi_amp[q]*half_pi_ratio
    delay_c= -pi_Du-freeDu
    return XY_waveform_controller(sche,amp,pi_Du,90,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt='start')

def X_pi_p(sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
    amp= pi_amp[q]
    delay_c= -pi_Du-freeDu
    return XY_waveform_controller(sche,amp,pi_Du,0,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt=ref_point)

def Y_pi_p(sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
    amp= pi_amp[q]
    delay_c= -pi_Du-freeDu
    return XY_waveform_controller(sche,amp,pi_Du,90,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt=ref_point)

def X_12pi_p(sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
    amp= pi_amp[q]
    delay_c= -pi_Du-freeDu
    return XY_waveform_controller(sche,amp,pi_Du,0,q,delay_c,ref_pulse_sche,XY_waveform,ref_pt=ref_point)

def Z(sche,Z_amp,Du,q,ref_pulse_sche,freeDu,ref_position='start'):
    if Du!=0:
        delay_z= -Du-freeDu
        return sche.add(SquarePulse(duration= Du,amp=Z_amp, port=q+":fl", clock="cl0.baseband"),rel_time=delay_z,ref_op=ref_pulse_sche,ref_pt=ref_position,)


def Zc(sche,Z_amp,Du,qc,ref_pulse_sche,freeDu):
    if Du!=0:
        delay_z= -Du-freeDu
        return sche.add(SquarePulse(duration= Du,amp=Z_amp, port=qc+":fl", clock="cl0.baseband"),rel_time=delay_z,ref_op=ref_pulse_sche,ref_pt="start",)
    else: pass

def Readout(sche,q,R_amp,R_duration,powerDep=False):
    if powerDep is True:
        amp= R_amp
        Du= R_duration[q]
    else:    
        amp= R_amp[q]
        Du= R_duration[q]

    # tim_sample, env_sample = FlatTopGaussianPulse(Du,amp)
    # return sche.add(NumericalPulse(samples=env_sample,t_samples=tim_sample,port="q:res",clock=q+".ro",t0=0e-9),)
    return sche.add(SquarePulse(duration=Du,amp=amp,port="q:res",clock=q+".ro",t0=4e-9))

def Multi_Readout(sche,q,ref_pulse_sche,R_amp,R_duration,powerDep=False,):
    if powerDep is True:
        amp= R_amp
        Du= R_duration[q]
    else:    
        amp= R_amp[q]
        Du= R_duration[q]

    # tim_sample, env_sample = FlatTopGaussianPulse(Du,amp)

    # return sche.add(NumericalPulse(samples=env_sample,t_samples=tim_sample,port="q:res",clock=q+".ro",t0=0),ref_pt="start",ref_op=ref_pulse_sche,)
    return sche.add(SquarePulse(duration=Du,amp=amp,port="q:res",clock=q+".ro",t0=4e-9),ref_pt="start",ref_op=ref_pulse_sche,)

    
def Integration(sche,q,R_inte_delay:float,R_inte_duration,ref_pulse_sche,acq_index,acq_channel:int=0,single_shot:bool=False,get_trace:bool=False,trace_recordlength:float=5*1e-6):
    if single_shot== False:     
        bin_mode=BinMode.AVERAGE
    else: bin_mode=BinMode.APPEND
    # Trace acquisition does not support APPEND bin mode !!!
    if get_trace==False:
        return sche.add(SSBIntegrationComplex(
            duration=R_inte_duration[q],
            port="q:res",
            clock=q+".ro",
            acq_index=acq_index,
            acq_channel=acq_channel,
            bin_mode=bin_mode,
            ),rel_time=R_inte_delay
            ,ref_op=ref_pulse_sche,ref_pt="start")
    else:  
        return sche.add(Trace(
                duration=trace_recordlength,
                port="q:res",
                clock=q+".ro",
                acq_index=acq_index,
                acq_channel=acq_channel,
                bin_mode=BinMode.AVERAGE,
                ),rel_time=R_inte_delay
                ,ref_op=ref_pulse_sche,ref_pt="start")
    
def pulse_preview(quantum_device:QuantumDevice,sche_func:Schedule, sche_kwargs:dict, **kwargs):
    import plotly.io as pio
    pio.renderers.default='browser'
    device_compiler = SerialCompiler("Device compiler", quantum_device)
    comp_sched = device_compiler.compile(
        sche_func(**sche_kwargs)
    )
    comp_sched.plot_pulse_diagram(plot_backend="plotly",**kwargs).show() 
    
#%% schedule function

def One_tone_sche(
    frequencies: np.ndarray,
    q:str,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    powerDep:bool,
    repetitions:int=1,    
) -> Schedule:
    
    sched = Schedule("One tone spectroscopy (NCO sweep)",repetitions=repetitions)
    sched.add_resource(ClockResource(name=q+ ".ro", freq=frequencies.flat[0]))
    
    for acq_idx, freq in enumerate(frequencies):
        
        sched.add(SetClockFrequency(clock= q+ ".ro", clock_freq_new=freq))
        sched.add(Reset(q))
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=powerDep)
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
     
    return sched


# ? Doing...

def One_tone_multi_sche(
    frequencies: dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    powerDep:bool,
    bias:any=0,
    bias_dura:float=0,
    repetitions:int=1,    
) -> Schedule:
    
    qubits2read = list(frequencies.keys())
    sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
    sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)


    for acq_idx in range(sameple_idx):    

        for qubit_idx, q in enumerate(qubits2read):
            freq = frequencies[q][acq_idx]
            if acq_idx == 0:
                sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
            
            sched.add(Reset(q))
            sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
            sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")

            
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=powerDep)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration,powerDep=powerDep)
            
            if bias != 0 and bias_dura != 0: 
                Z(sched,bias,bias_dura,q,spec_pulse,0,ref_position='end')

            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
          
    return sched

def RabiSplitting_multi_sche(
    frequencies: dict,
    bias_couplers:list,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    powerDep:bool,
    bias:any=0,
    bias_dura:float=0,
    repetitions:int=1,    
) -> Schedule:
    
    qubits2read = list(frequencies.keys())
    sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
    sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)


    for acq_idx in range(sameple_idx):    

        for qubit_idx, q in enumerate(qubits2read):
            freq = frequencies[q][acq_idx]
            if acq_idx == 0:
                sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
            
            sched.add(Reset(q))
            sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
            sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")

            
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=powerDep)
                if bias != 0 and bias_dura != 0 and len(list(bias_couplers)):
                    for qc in bias_couplers:
                        Z(sched,bias,bias_dura,qc,spec_pulse,0,ref_position='end')
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration,powerDep=powerDep)

            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
          
    return sched 

def Two_tone_sche(
    frequencies: np.ndarray,
    q:str,
    spec_amp:float,
    spec_Du:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,   
    ref_pt:str='end'
    
) -> Schedule:
    sched = Schedule("Two tone spectroscopy (NCO sweep)",repetitions=repetitions)
    sched.add_resource(ClockResource(name=q+".01", freq=frequencies.flat[0]))
    for acq_idx, freq in enumerate(frequencies):
        sched.add(SetClockFrequency(clock= q+".01", clock_freq_new=freq))
        sched.add(Reset(q))
        # sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)

        Spec_pulse(sched,spec_amp,spec_Du,q,spec_pulse,electrical_delay, ref_point=ref_pt)
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)

     
    return sched

def multi_Two_tone_sche(
    frequencies:dict,
    spec_amp:any,
    spec_Du:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,   
    ref_pt:str='end'
    
) -> Schedule:
    sched = Schedule("Two tone spectroscopy (NCO sweep)",repetitions=repetitions)
    qubits2read = list(frequencies.keys())
    sameple_idx = array(frequencies[qubits2read[0]]).shape[0]

    for acq_idx in range(sameple_idx):    
        for qubit_idx, q in enumerate(qubits2read):
            freq = frequencies[q][acq_idx]
            if acq_idx == 0:
                sched.add_resource(ClockResource(name=q+ ".01", freq=array(frequencies[q]).flat[0]))
        
            sched.add(SetClockFrequency(clock= q+".01", clock_freq_new=freq))
            sched.add(Reset(q))
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration,powerDep=False)

            Spec_pulse(sched,spec_amp,spec_Du,q,spec_pulse,electrical_delay, ref_point=ref_pt)

            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)

     
    return sched

# try to put a bias on RO part 03/20
def Z_gate_two_tone_sche(
    frequencies: np.ndarray,
    q:str,
    Z_amp:any,
    spec_amp:float,
    spec_Du:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,   
    Z_ro_amp:float =0 
) -> Schedule:
    sched = Schedule("Zgate_two_tone spectroscopy (NCO sweep)",repetitions=repetitions)
    sched.add_resource(ClockResource(name=q+".01", freq=frequencies.flat[0]))
    for acq_idx, freq in enumerate(frequencies):
        sched.add(SetClockFrequency(clock= q+ ".01", clock_freq_new=freq))
        sched.add(Reset(q))
        # sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        Spec_pulse(sched,spec_amp,spec_Du,q,spec_pulse,electrical_delay)
        Z(sched,Z_amp,spec_Du,q,spec_pulse,electrical_delay)
        if Z_ro_amp != 0:
            Z(sched,Z_ro_amp,R_integration[q],q,spec_pulse,0,'end')

        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
     
    return sched

#? Warning: Should avoid add the z-gate on the same channel at the same time
def multi_Z_gate_two_tone_sche(
    frequencies: dict,
    bias_qs:list,
    Z_amp:any,
    spec_amp:float,
    spec_Du:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,   
    
) -> Schedule:
    sched = Schedule("Zgate_two_tone spectroscopy (NCO sweep)",repetitions=repetitions)
    
    additional_bias_list = []
    qubits2read = list(frequencies.keys())
    sameple_idx = array(frequencies[qubits2read[0]]).shape[0]

    for q_need_bias in bias_qs:
        if q_need_bias not in qubits2read:
            additional_bias_list.append(q_need_bias)



    for acq_idx in range(sameple_idx):    
        for qubit_idx, q in enumerate(qubits2read):
            freq = frequencies[q][acq_idx]
            if acq_idx == 0:
                sched.add_resource(ClockResource(name=q+".01", freq=array(frequencies[q]).flat[0]))
   
            sched.add(SetClockFrequency(clock= q+ ".01", clock_freq_new=freq))
            sched.add(Reset(q))
            
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
            
            Spec_pulse(sched,spec_amp,spec_Du,q,spec_pulse,electrical_delay)

            if q in bias_qs:
                Z(sched,Z_amp,spec_Du,q,spec_pulse,electrical_delay)
            
            if len(additional_bias_list) != 0:
                for qs in additional_bias_list:
                    Z(sched,Z_amp,spec_Du,qs,spec_pulse,electrical_delay)
            

            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
     
    return sched

def Qubit_state_heterodyne_spec_sched_nco(
    frequencies: np.ndarray,
    q:str,
    ini_state:str,
    pi_amp: dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("One tone qubit-state spectroscopy (NCO sweep)",repetitions=repetitions)
    sched.add_resource(ClockResource(name=q+ ".ro", freq=frequencies.flat[0]))
    for acq_idx, freq in enumerate(frequencies):
        sched.add(SetClockFrequency(clock= q+ ".ro", clock_freq_new=freq))
        sched.add(Reset(q))
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        if ini_state=='e': 
            X_pi_p(sched,pi_amp,q,spec_pulse,freeDu=0)
            
        else: None
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        # Integration must after all the pulses
    return sched


def Rabi_sche(
    q:str,
    New_fxy:any,
    XY_amp: any,
    XY_duration:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    XY_theta:str,
    Rabi_type:str,
    repetitions:int=1,
    chevron:bool=False
) -> Schedule:

    sched = Schedule(Rabi_type,repetitions=repetitions)
    amps = np.asarray(XY_amp)
    amps = amps.reshape(amps.shape or (1,))
    durations = np.asarray(XY_duration)
    durations = durations.reshape(durations.shape or (1,))
    
    if Rabi_type=='TimeRabi':
       Para_XY_amp= amps*np.ones(np.shape(durations))
       Para_XY_Du= durations
       
    elif Rabi_type=='PowerRabi':
        Para_XY_amp =amps
        Para_XY_Du = XY_duration*np.ones(np.shape(amps))   
    else: raise KeyError ('Typing error: Rabi_type')
    
    
    
    for acq_idx, (amp, duration) in enumerate(zip(Para_XY_amp,Para_XY_Du)):
        
        
        sched.add(Reset(q))
        if chevron:
            sched.add(
                SetClockFrequency(clock=q+ ".01", clock_freq_new= New_fxy))
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        if XY_theta== 'X_theta':
            X_theta(sched,amp,duration,q,spec_pulse,freeDu=electrical_delay)
        elif XY_theta== 'Y_theta':
            Y_theta(sched,amp,duration,q,spec_pulse,freeDu=electrical_delay)
        else: raise KeyError ('Typing error: XY_theta')
       
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
    
    return sched

def multi_PowerRabi_sche(
    pi_amp:dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    XY_theta:str,
    repetitions:int=1,
    OS_or_not:bool=False
) -> Schedule:
    
    qubits2read = list(pi_amp.keys())
    sample_len = pi_amp[qubits2read[0]].shape[0]
    sched = Schedule("RabiOscillation",repetitions=repetitions)

    match XY_theta:
        case 'Y_theta':
            gate:callable = Y_theta
        case _:
            gate:callable = X_theta

    
    for acq_idx in range(sample_len):    
        for qubit_idx, q in enumerate(qubits2read):
            sched.add(Reset(q))
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
            
            
            gate(sched,pi_amp[q][acq_idx],pi_dura[q],q,spec_pulse,freeDu=electrical_delay)
           
        
            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=OS_or_not,get_trace=False,trace_recordlength=0)
    
    return sched

def multi_TimeRabi_sche(
    pi_amp:dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    XY_theta:str,
    repetitions:int=1,
    OS_or_not:bool=False
) -> Schedule:
    
    qubits2read = list(pi_dura.keys())
    sample_len = pi_dura[qubits2read[0]].shape[0]
    sched = Schedule("RabiOscillation",repetitions=repetitions)

    match XY_theta:
        case 'Y_theta':
            gate:callable = Y_theta
        case _:
            gate:callable = X_theta

    
    for acq_idx in range(sample_len):    
        for qubit_idx, q in enumerate(qubits2read):
            sched.add(Reset(q))
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
            
            
            gate(sched,pi_amp[q],pi_dura[q][acq_idx],q,spec_pulse,freeDu=electrical_delay)
           
        
            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=OS_or_not,get_trace=False,trace_recordlength=0)
    
    return sched

def multi_RabiSchevron_sche():
    pass

def Rabi_12_sche(
    q:str,
    XY_12_amp: any,
    XY_12_duration:any,
    pi_amp: dict,
    pi_dura:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    XY_theta:str,
    Rabi_type:str,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule(Rabi_type,repetitions=repetitions)
    amps = np.asarray(XY_12_amp)
    amps = amps.reshape(amps.shape or (1,))
    durations = np.asarray(XY_12_duration)
    durations = durations.reshape(durations.shape or (1,))
    
    if Rabi_type=='TimeRabi':
       Para_XY_amp= amps*np.ones(np.shape(durations))
       Para_XY_Du= durations
       
    elif Rabi_type=='PowerRabi':
        Para_XY_amp =amps
        Para_XY_Du = XY_12_duration*np.ones(np.shape(amps))   
    else: raise KeyError ('Typing error: Rabi_type')
    
    
    
    for acq_idx, (amp, duration) in enumerate(zip(Para_XY_amp,Para_XY_Du)):
        
        
        sched.add(Reset(q))
        
        # sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)

        if XY_theta== 'X_theta':
            X_theta(sched,pi_amp[q],pi_dura,q,spec_pulse,freeDu=electrical_delay)
        elif XY_theta== 'Y_theta':
            Y_theta(sched,pi_amp[q],pi_dura,q,spec_pulse,freeDu=electrical_delay)
        else: raise KeyError ('Typing error: XY_theta')
       
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
    
    return sched


def Zgate_Rabi_sche(
    q:str,
    XY_amp: any,
    XY_duration:any,
    Z_amp:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    XY_theta:str,
    Rabi_type:str,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule('Zgate_'+Rabi_type,repetitions=repetitions)
    amps = np.asarray(XY_amp)
    amps = amps.reshape(amps.shape or (1,))
    durations = np.asarray(XY_duration)
    durations = durations.reshape(durations.shape or (1,))
    
    if Rabi_type=='TimeRabi':
       Para_XY_amp= amps*np.ones(np.shape(durations))
       Para_XY_Du= durations
       
    elif Rabi_type=='PowerRabi':
        Para_XY_amp =amps
        Para_XY_Du = XY_duration*np.ones(np.shape(amps))   
    else: raise KeyError ('Typing error: Rabi_type')
    
    
    
    for acq_idx, (amp, duration) in enumerate(zip(Para_XY_amp,Para_XY_Du)):
        
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        if XY_theta== 'X_theta':
            X_theta(sched,amp,duration,q,spec_pulse,freeDu=0)
        elif XY_theta== 'Y_theta':
            Y_theta(sched,amp,duration,q,spec_pulse,freeDu=0)
        else: raise KeyError ('Typing error: XY_theta')
        
        Z(sched,Z_amp,duration,q,spec_pulse,0)
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
    
    return sched

def T1_sche(
    q:str,
    pi_amp: dict,
    pi_dura:float,
    freeduration:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("T1", repetitions=repetitions)
    
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
    
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        
        
        X_pi_p(sched,pi_amp,q,pi_dura,spec_pulse,freeDu+electrical_delay)
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def multi_T1_sche(
    freeduration:dict,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,
) -> Schedule:
    qubits2read = list(freeduration.keys())
    sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
    sched = Schedule("T1", repetitions=repetitions)
    
    for acq_idx in range(sameple_idx):
        for qubit_idx, q in enumerate(qubits2read):
            freeDu = freeduration[q][acq_idx]
            sched.add(Reset(q))
        
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
        
        
            X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu+electrical_delay)
            
            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def mix_T1_sche(
    q:str,
    pi_amp: dict,
    freeduration:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:


    sched = Schedule("T1", repetitions=repetitions)
    
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
    
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        
        Spec_pulse(sched,pi_amp[q],10e-6,q,spec_pulse,freeDu) # drive to mix state
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def XY_Z_timing_sche(
    q:str,
    pi_amp: dict,
    freeduration:float,
    Z_delay:any,
    Z_amp:float,
    Z_Du:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("XY_Z_timing", repetitions=repetitions)
    
    for acq_idx, delay in enumerate(Z_delay):
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        
        X_pi_p(sched,pi_amp,q,spec_pulse,freeDu=freeduration)
        
        Z(sched,Z_amp,Z_Du,q,spec_pulse,freeDu=freeduration)
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def Zgate_T1_sche(
    q:str,
    pi_amp: dict,
    pi_dura:dict,
    freeduration:any,
    Z_amp:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Zgate_T1", repetitions=repetitions)
    
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")

    
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        
        X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu+electrical_delay)
        
        Z(sched,Z_amp,freeDu,q,spec_pulse,freeDu=electrical_delay)
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def multi_Zgate_T1_sche(
    freeduration:dict,
    Z_amp:any,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,
) -> Schedule:
    qubits2read = list(freeduration.keys())
    sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
    sched = Schedule("Zgate-T1", repetitions=repetitions)

    for acq_idx in range(sameple_idx):
        for qubit_idx, q in enumerate(qubits2read):
            freeDu = freeduration[q][acq_idx]
            sched.add(Reset(q))
        
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
        
        
            X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu+electrical_delay)
            Z(sched,Z_amp,freeDu,q,spec_pulse,freeDu=electrical_delay)
            
            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched


def Ramsey_sche(
    q:str,
    pi_amp: dict,
    New_fxy:any,
    freeduration:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    pi_dura:float=20e-9,
    repetitions:int=1,
    echo_pi_num:int = 0,
    second_pulse_phase:str='x',
) -> Schedule:

    sched = Schedule("Ramsey", repetitions=repetitions)
    
    pi_Du= pi_dura
    
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(
            SetClockFrequency(clock=q+ ".01", clock_freq_new= New_fxy))
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        # we start construction from readout
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        if second_pulse_phase.lower() == 'x':
            first_half_pi = X_pi_2_p(sched,pi_amp,q,pi_Du,spec_pulse,freeDu=electrical_delay)
        else:
            first_half_pi = Y_pi_2_p(sched,pi_amp,q,pi_Du,spec_pulse,freeDu=electrical_delay)
        
        if echo_pi_num != 0:
            a_separate_free_Du = freeDu / echo_pi_num
            for pi_idx in range(echo_pi_num):
                if pi_idx == 0 :
                    pi = X_pi_p(sched,pi_amp,q,pi_Du,first_half_pi,0.5*a_separate_free_Du)
                else:
                    pi = X_pi_p(sched,pi_amp,q,pi_Du,pi,1*a_separate_free_Du)
            X_pi_2_p(sched,pi_amp,q,pi_Du,pi,0.5*a_separate_free_Du)
        else:
            X_pi_2_p(sched,pi_amp,q,pi_Du,first_half_pi,freeDu)

        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def multi_ramsey_sche(
    freeduration:dict,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    echo_pi_num:dict,
    repetitions:int=1,
    second_pulse_phase:str='x',
) -> Schedule:
    qubits2read = list(freeduration.keys())
    sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
    sched = Schedule("Ramsey", repetitions=repetitions)
    
    for acq_idx in range(sameple_idx):   
        for qubit_idx, q in enumerate(qubits2read):
            freeDu = freeduration[q][acq_idx]
            sched.add(Reset(q))

            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)

            if second_pulse_phase.lower() == 'x':
                first_half_pi = X_pi_2_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
            else:
                first_half_pi = Y_pi_2_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
        
            if echo_pi_num[q] != 0:
                a_separate_free_Du = freeDu / echo_pi_num[q]
                for pi_idx in range(echo_pi_num[q]):
                    if pi_idx == 0 :
                        pi = X_pi_p(sched,pi_amp,q,pi_dura[q],first_half_pi,0.5*a_separate_free_Du)
                    else:
                        pi = X_pi_p(sched,pi_amp,q,pi_dura[q],pi,1*a_separate_free_Du)
                X_pi_2_p(sched,pi_amp,q,pi_dura[q],pi,0.5*a_separate_free_Du)
            else:
                X_pi_2_p(sched,pi_amp,q,pi_dura[q],first_half_pi,freeDu)

            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def Ramsey_readOther_sche(
    q:str,
    read_q:str,
    pi_amp: dict,
    New_fxy:float,
    freeduration:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    pi_dura:float=20e-9,
    repetitions:int=1,
    echo_pi_num:int = 0,
    second_pulse_phase:str='x',
) -> Schedule:

    sched = Schedule("Ramsey", repetitions=repetitions)
    
    pi_Du= pi_dura
    print(f"second pulse phase: {second_pulse_phase}")
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(
            SetClockFrequency(clock=q+ ".01", clock_freq_new= New_fxy))
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        # we start construction from readout
        spec_pulse = Readout(sched,read_q,R_amp,R_duration,powerDep=False)
        if second_pulse_phase.lower() == 'x':
            first_half_pi = X_pi_2_p(sched,pi_amp,q,pi_Du,spec_pulse,freeDu=electrical_delay)
        else:
            first_half_pi = Y_pi_2_p(sched,pi_amp,q,pi_Du,spec_pulse,freeDu=electrical_delay)
        
        if echo_pi_num != 0:
            a_separate_free_Du = freeDu / echo_pi_num
            for pi_idx in range(echo_pi_num):
                if pi_idx == 0 :
                    pi = X_pi_p(sched,pi_amp,q,pi_Du,first_half_pi,0.5*a_separate_free_Du)
                else:
                    pi = X_pi_p(sched,pi_amp,q,pi_Du,pi,1*a_separate_free_Du)
            X_pi_2_p(sched,pi_amp,q,pi_Du,pi,0.5*a_separate_free_Du)
        else:
            X_pi_2_p(sched,pi_amp,q,pi_Du,first_half_pi,freeDu)

        Integration(sched,read_q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def Zgate_Ramsey_sche(
    q:str,
    pi_amp: dict,
    New_fxy:float,
    freeduration:any,
    Z_amp:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Zgate_Ramsey", repetitions=repetitions)
    
    pi_Du= 20*1e-9


   
    for acq_idx, freeDu_ in enumerate(freeduration):
        
        sched.add(
            SetClockFrequency(clock=q+ ".01", clock_freq_new= New_fxy))
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        
        X_pi_2_p(sched,pi_amp,q,spec_pulse,freeDu=freeDu_+pi_Du+electrical_delay)
        
        X_pi_2_p(sched,pi_amp,q,spec_pulse,freeDu=0+electrical_delay)
        
        Z(sched,Z_amp,freeDu_,q,spec_pulse,freeDu=pi_Du+electrical_delay)   
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)

    return sched

def Zz_Interaction(
    q:str,
    excite_qubit:str,
    pi_amp: dict,
    excite_pi_amp: dict,
    New_fxy:float,
    freeduration:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    pi_dura:float=20e-9,
    excite_pi_dura:float=20e-9,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Ramsey", repetitions=repetitions)
    
    pi_Du= pi_dura

    excite_pi_Du= excite_pi_dura
   
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(
            SetClockFrequency(clock=q+ ".01", clock_freq_new= New_fxy))
        
        sched.add(Reset(q))

        sched.add(Reset(excite_qubit))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)

        X_pi_p(sched,excite_pi_amp,excite_qubit,excite_pi_Du,spec_pulse,freeDu=freeDu+pi_Du*2)

        X_pi_2_p(sched,pi_amp,q,pi_Du,spec_pulse,freeDu=freeDu+pi_Du)
        
        X_pi_2_p(sched,pi_amp,q,pi_Du,spec_pulse,freeDu=0)

        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
    return sched

def Qubit_SS_sche(
    q:str,
    ini_state:str,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Single shot", repetitions=repetitions)
    
    sched.add(Reset(q))
    
    sched.add(IdlePulse(duration=5000*1e-9))
    
    spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
    
    if ini_state=='e': 
        X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
        
    else: None
    
    Integration(sched,q,R_inte_delay,R_integration,spec_pulse,0,single_shot=True,get_trace=False,trace_recordlength=0)

    return sched

def multi_Qubit_SS_sche(
    ini_state:str,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Single shot", repetitions=repetitions)

    for qubit_idx, q in enumerate(R_integration):

        sched.add(Reset(q))
        if qubit_idx == 0:
            spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        else:
            Multi_Readout(sched,q,spec_pulse,R_amp,R_duration,powerDep=False)
    
        if ini_state=='e': 
            X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
        
    
        Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=0,acq_channel=qubit_idx,single_shot=True,get_trace=False,trace_recordlength=0)

    return sched

#? Calibrations :
def ROF_Cali_sche(
    q:str,
    ro_freq:np.ndarray,
    ini_state:str,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Single shot", repetitions=repetitions)
    sched.add_resource(ClockResource(name=q+ ".ro", freq=ro_freq.flat[0]))
    for acq_idx, rof in enumerate(ro_freq):
        sched.add(SetClockFrequency(clock= q+ ".ro", clock_freq_new=rof))
        sched.add(Reset(q))
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)

        if ini_state=='e': 
            X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
            
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)

    return sched

def multi_ROF_Cali_sche(
    ro_freq:dict,
    ini_state:str,
    pi_amp: dict,
    pi_dura:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,
) -> Schedule:
    qubits2read = list(ro_freq.keys())
    sameple_idx = array(ro_freq[qubits2read[0]]).shape[0]
    sched = Schedule("ROF calibration", repetitions=repetitions)
    for acq_idx in range(sameple_idx):    

        for qubit_idx, q in enumerate(qubits2read):
            freq = ro_freq[q][acq_idx]
            if acq_idx == 0:
                sched.add_resource(ClockResource(name=q+ ".ro", freq=array(ro_freq[q]).flat[0]))

            sched.add(SetClockFrequency(clock= q+ ".ro", clock_freq_new=freq))
            sched.add(Reset(q))
            sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")

            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)

            if ini_state=='e': 
                X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
                
            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
          
    return sched

def multi_ROL_Cali_sche(
    R_amp_coefs: dict,
    ini_state:str,
    pi_amp: dict,
    pi_dura:dict,
    R_duration: dict,
    R_amp: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,
) -> Schedule:
    qubits2read = list(R_amp_coefs.keys())
    sameple_idx = array(R_amp_coefs[qubits2read[0]]).shape[0]
    sched = Schedule("ROL calibration", repetitions=repetitions)
    
    for acq_idx in range(sameple_idx):    

        for qubit_idx, q in enumerate(qubits2read):
            rol_coef = R_amp_coefs[q][acq_idx]
            
            sched.add(Reset(q))
    
            if qubit_idx == 0:
                spec_pulse = Readout(sched,q,{q:R_amp[q]*rol_coef},R_duration)
            else:
                Multi_Readout(sched,q,spec_pulse,{q:R_amp[q]*rol_coef},R_duration)

            if ini_state=='e': 
                X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
                
            Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
          
    return sched

def PI_amp_cali_sche(
    q:str,
    XY_amp: dict,
    pi_amp_coefs: np.ndarray,
    pi_pair_num:int,
    XY_duration:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
    )-> Schedule:


    sched = Schedule("Pi amp modification", repetitions=repetitions)
    for acq_idx, amp_coef in enumerate(np.asarray(pi_amp_coefs)):
        
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
    
        read_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        
        for pi_num in range(pi_pair_num):
            for pi_idx in range(2):
                spec_pulse = X_pi_p(sched,{str(q):float(XY_amp[q])*amp_coef},q,XY_duration,read_pulse if (pi_num == 0 and pi_idx == 0) else spec_pulse, freeDu=electrical_delay if (pi_num == 0 and pi_idx == 0) else 0)
                
        Integration(sched,q,R_inte_delay,R_integration,read_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)

    return sched

def multi_PI_amp_cali_sche(
    pi_amp_coefs: dict,
    pi_pair_num:int,
    XY_amp: dict,
    XY_duration:dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:dict,
    repetitions:int=1,
    )-> Schedule:
    qubits2read = list(pi_amp_coefs.keys())
    sameple_idx = array(pi_amp_coefs[qubits2read[0]]).shape[0]
    sched = Schedule("Pi amp modification", repetitions=repetitions)

    if np.median(pi_amp_coefs[qubits2read[0]]) > 0.9 : # should be 1 when pi amp cali
        pairs:int = 2
    else:                          # should be 0.5 when half-pi amp cali
        pairs:int = 4
    
    for acq_idx in range(sameple_idx):
        for qubit_idx, q in enumerate(qubits2read):
            amp_coef = pi_amp_coefs[q][acq_idx]
            
            
            sched.add(Reset(q))
            
            if qubit_idx == 0:
                read_pulse = Readout(sched,q,R_amp,R_duration)
            else:
                Multi_Readout(sched,q,read_pulse,R_amp,R_duration)
            
            for pi_num in range(pi_pair_num):
                for pi_idx in range(pairs):
                    spec_pulse = X_pi_p(sched,{str(q):float(XY_amp[q])*amp_coef},q,XY_duration[q],read_pulse if (pi_num == 0 and pi_idx == 0) else spec_pulse, freeDu=electrical_delay if (pi_num == 0 and pi_idx == 0) else 0)
                    
            Integration(sched,q,R_inte_delay[q],R_integration,read_pulse,acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)

    return sched

def pi_half_cali_sche(
    q:str,
    pi_amp: dict,
    pi_half_coefs: np.ndarray,
    half_pi_quadruple_num:int,
    XY_duration:float,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
    )-> Schedule:
    
    sched = Schedule("Pi amp modification", repetitions=repetitions)
    for acq_idx, amp_coef in enumerate(np.asarray(pi_half_coefs)):
        sched.add(Reset(q))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
    
        read_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
        for pi_num in range(half_pi_quadruple_num):
            for pi_idx in range(4):
                spec_pulse = X_pi_p(sched,{str(q):float(pi_amp[q])*amp_coef},q,XY_duration,read_pulse if (pi_num == 0 and pi_idx == 0) else spec_pulse,freeDu=electrical_delay if (pi_num == 0 and pi_idx == 0) else 0)
                
        Integration(sched,q,R_inte_delay,R_integration,read_pulse,acq_idx,single_shot=False,get_trace=False,trace_recordlength=0)

    return sched

def drag_coef_cali():
    """ (X,Y/2) and (Y,X/2)"""
    pass

def StarkShift_cali():
    """ N * (X,-X)"""
    pass

def Qubit_amp_SS_sche(
    q:str,
    ini_state:str,
    pi_amp: dict,
    R_amp: any,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Single shot", repetitions=repetitions)

    sched.add(Reset(q))
    
    sched.add(IdlePulse(duration=5000*1e-9))
    
    spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=True)
    
    if ini_state=='e': 
        X_pi_p(sched,pi_amp,q,spec_pulse,freeDu=0)
        
    else: None
    Integration(sched,q,R_inte_delay,R_integration,spec_pulse,0,single_shot=True,get_trace=False,trace_recordlength=0)

    return sched

def Trace_sche(
    q:str,
    ini_state:str,
    pi_amp: dict,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    trace_recordlength:float,
    repetitions:int=1,
) -> Schedule:
   
    sched = Schedule("Trace", repetitions=repetitions)
    
    sched.add(Reset(q))
    
    sched.add(IdlePulse(duration=5000*1e-9))
    
    spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
    
    if ini_state=='e': 
        X_pi_p(sched,pi_amp,q,spec_pulse,freeDu=0)
        
    else: None
    
    Integration(sched,q,0,R_integration,spec_pulse,acq_index=0,single_shot=False,get_trace=True,trace_recordlength=trace_recordlength)


    return sched

def iSwap_sche(
    q:list,
    pi_amp: dict,
    pi_dura:dict,
    target_pi:str,
    freeduration:any,
    target_Z1:str,
    target_Z2:str,
    target_Zc:str,
    Z1_amp:any,
    Z2_amp:any,
    Zc_amp:any,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    Z_on_off:dict,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("iSwap", repetitions=repetitions)
    
    for acq_idx, freeDu in enumerate(freeduration):
        
        sched.add(Reset(q[0]))
        
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")

    
        spec_pulse = Readout(sched,q[0],R_amp,R_duration,powerDep=False)
        
        X_pi_p(sched,pi_amp,target_pi,pi_dura,spec_pulse,freeDu)
        
        for i in range(1,len(q)):
            Multi_Readout(sched,q[i],spec_pulse,R_amp,R_duration,powerDep=False)    
            Integration(sched,q[i],R_inte_delay,R_integration,spec_pulse,acq_index=acq_idx,acq_channel=i,single_shot=False,get_trace=False,trace_recordlength=0)
            
        if Z_on_off['Z1_on_off']== True:
            Z(sched,Z1_amp,freeDu,target_Z1,spec_pulse,freeDu=0)
        else:
            None
        if Z_on_off['Z2_on_off']== True:
            Z(sched,Z2_amp,freeDu,target_Z2,spec_pulse,freeDu=0)
        else:
            None    
        if Z_on_off['Zc_on_off']== True:
            Zc(sched,Zc_amp,freeDu,target_Zc,spec_pulse,freeDu=0)
        else:
            None

        
        Integration(sched,q[0],R_inte_delay,R_integration,spec_pulse,acq_index=acq_idx,acq_channel=0,single_shot=False,get_trace=False,trace_recordlength=0)
        
        
    return sched

#%% plot

def Readout_F_opt_Plot(quantum_device:QuantumDevice, results:dict):

    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250)
    f= results['f_samples']
    f_g= results['f_g']
    f_e= results['f_e']
    mag_g= results['mag_g']
    mag_e= results['mag_e']
    f_eff_bare= (f_g+f_e)/2
    ax.plot(f, mag_g*1000,'b',alpha=0.8,label=r"$0$",lw=2)
    ax.plot(f, mag_e*1000,'r',alpha=0.8,label=r"$\pi$",lw=2)
    ax.axvline(f_g, color = "blue", ls = "--")
    ax.axvline(f_e, color = "red", ls = "--")
    ax.axvline(f_eff_bare, color = "k", ls = "--")
    ax.legend(fontsize=10)
    x_label= 'Frequency'+' [GHz]'
    y_label= r"$\vert S_{21}\vert$"+' [mV]'
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    fig.suptitle(f"Resonator spectroscopy, {quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()


def Amp_phase_plot(quantum_device:QuantumDevice, results:dict,title:str):
    fig, ax = plt.subplots(1,2,figsize=plt.figaspect(1/2), sharey = False)
    dh.to_gridded_dataset(results.dataset).y0.plot(ax = ax[0],lw=5)
    dh.to_gridded_dataset(results.dataset).y1.plot(ax = ax[1],lw=5)
    ax[0].tick_params('both',labelsize='12', length=4, width=2)
    ax[1].tick_params('both',labelsize='12', length=4, width=2)
    ax[0].xaxis.get_label().set_fontsize(20)
    ax[0].yaxis.get_label().set_fontsize(20)
    ax[1].xaxis.get_label().set_fontsize(20)
    ax[1].yaxis.get_label().set_fontsize(20)
    fig.suptitle(title+f"{quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()
    
def Amp_plot(quantum_device:QuantumDevice, results:dict,title:str):
    fig, ax = plt.subplots(1,1,figsize=plt.figaspect(1), sharey = False)
    dh.to_gridded_dataset(results.dataset).y0.plot(ax = ax)
    
    fig.suptitle(title+f"{quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()
    
def hist_plot(q:str,data:dict,title:str, save_path:str='', show:bool=True):
    fig, ax = plt.subplots(nrows =1,figsize =(4,3),dpi =250) 
    m, bins, patches = ax.hist(np.array(data[q]), bins='auto', density=False)
    ax.axvline(np.median(np.array(data[q])), color = "k", ls = "--",lw=1)
    ax.set_ylabel('Counts')
    if title.lower() in ['t1','t2','t2*']:
        ax.set_xlabel(f"{title} (µs)")
        plt.title(f"{title} = {round(np.median(np.array(data[q])),1)}  $\pm$ {round(np.std(np.array(data[q])),1)} µs")
    else:
        if title.lower() in ["thermalpop"]:
            ax.set_xlabel(f"{title} (%)")
            plt.title(f"{title} = {round(np.median(np.array(data[q])),1)}  $\pm$ {round(np.std(np.array(data[q])),1)} %")
        else:
            ax.set_xlabel(f"{title} (mK)")
            plt.title(f"{title} = {round(np.median(np.array(data[q])),1)}  $\pm$ {round(np.std(np.array(data[q])),1)} mK")
        
    plt.tight_layout()
    if save_path != '':
        fig.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close()
    
def Z_bias_error_bar_plot(q:str,data:dict,title:str):
    times, Z_bias= data['plot_parameters'][0],data['plot_parameters'][1]
    T1_array= np.array(data[q]).reshape(times, len(Z_bias))
    mean, sigma= T1_array.mean(axis=0),T1_array.std(axis=0)
    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250) 
    ax:plt.Axes
    ax.errorbar(Z_bias, mean, yerr=2*sigma, fmt='o', color='blue',
              ecolor='blue', elinewidth=2, capsize=5)
    ax.set_ylabel(title)
    ax.set_xlabel(r"$Z\ bias$")
    ax.set_title(r"$times= %.0f $" %(times))
    fig.tight_layout()
    plt.show()
    return T1_array

def Ramsey_F_Z_bias_error_bar_plot(q:str,data:dict):
    times, Z_bias= data['plot_parameters'][0],data['plot_parameters'][1]
    Ramsey_F_array= np.array(data[q]).reshape(times, len(Z_bias))
    Ramsey_F_mean= Ramsey_F_array.mean(axis=0)
    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250) 
    ax.plot(Z_bias,Ramsey_F_mean*1e-6,'o', color="blue", alpha=0.5, ms=10)
    ax.set_ylabel(r"$detuning\ $[MHz]")
    ax.set_xlabel(r"$Z\ bias$")
    ax.set_title(r"$times= %.0f $" %(times))
    fig.tight_layout()
    
    return Ramsey_F_array
    
def dataset_to_array(dataset:xr.core.dataset.Dataset,dims:float):
    if dims==1:
        I= dataset.y0.data
        Q= dataset.y1.data
    elif dims==2:
        gridded_dataset = dh.to_gridded_dataset(dataset)
        I=gridded_dataset.y0.data
        Q=gridded_dataset.y1.data
    else: raise KeyError ('dims is not 1 or 2')  
    
    return I,Q

def Multi_dataset_to_array(dataset:xr.core.dataset.Dataset,dims:float,Q:list):
    I_data={}
    Q_data={}
    if dims==1:
        for i in range(len(Q)):
            I_data[Q[i]]= dataset['y'+str(i)].data
            Q_data[Q[i]]= dataset['y'+str(1+i)].data
    elif dims==2:
        gridded_dataset = dh.to_gridded_dataset(dataset)
        for i in range(len(Q)):
            I_data[Q[i]]= gridded_dataset['y'+str(i)].data
            Q_data[Q[i]]= gridded_dataset['y'+str(1+i)].data         
        
    else: raise KeyError ('dims is not 1 or 2')  
    
    return I_data,Q_data

def plot_textbox(ax,text, **kw):
    box_props = dict(boxstyle="round", pad=0.4, facecolor="white", alpha=0.5)
    new_kw_with_defaults = dict(
        x=1.05,
        y=0.95,
        transform=ax.transAxes,
        bbox=box_props,
        verticalalignment="top",
        s=text,
    )
    new_kw_with_defaults.update(kw)
    t_obj = ax.text(**new_kw_with_defaults)
    return t_obj


def Qubit_state_single_shot_plot(results:dict,Plot_type:str,y_scale:str):
    ce_I,ce_Q,sig=1000*results['fit_pack'][0][0],1000*results['fit_pack'][0][1],1000*results['fit_pack'][5]
    Inte_g_data,Inte_e_data= results['fit_pack'][6],results['fit_pack'][7]
    Ig,Qg= results['rot_IQdata'][0][0],results['rot_IQdata'][0][1]
    Ie,Qe= results['rot_IQdata'][1][0],results['rot_IQdata'][1][1]
    I_ro,I_fit= 1000*results['I_ro'],1000*results['I_fit']
    Mgg= gauss_func(I_fit,0,sig,results['fit_pack'][1])
    Meg= gauss_func(I_fit,ce_I,sig,results['fit_pack'][2])
    Mge= gauss_func(I_fit,0,sig,results['fit_pack'][3])
    Mee= gauss_func(I_fit,ce_I,sig,results['fit_pack'][4])
    
    fig, ax = plt.subplots(nrows =1,figsize =(3,3),dpi =200)
    fig1,ax1= plt.subplots(nrows =1,figsize =(6,4),dpi =200)
    if Plot_type== 'g_only':
        ax.scatter(1000*Ig, 1000*Qg, color="blue", alpha=0.5, s=5)  
        ax1.plot(I_ro, Inte_g_data,'bo',alpha=0.5,ms=8)
        ax1.plot(I_fit, Mgg,'b',alpha=0.8,lw=2)
        ax1.plot(I_fit, Meg,'--b',alpha=0.8,lw=2)
        ax1.plot(I_fit, Mgg+Meg,'--k',alpha=1,lw=1)
    elif Plot_type== 'e_only':
        ax.scatter(1000*Ie, 1000*Qe, color="red", alpha=0.5, s=5)
        ax1.plot(I_ro, Inte_e_data,'ro',alpha=0.5,ms=8)
        ax1.plot(I_fit, Mee,'r',alpha=0.8,lw=2)
        ax1.plot(I_fit, Mge,'--r',alpha=0.8,lw=2)
        ax1.plot(I_fit, Mge+Mee,'--k',alpha=1,lw=1)
    elif Plot_type== 'both': 
        ax.scatter(1000*Ig, 1000*Qg, color="blue", alpha=0.5, s=0.3)
        ax.scatter(1000*Ie, 1000*Qe, color="red", alpha=0.5, s=0.3)
        ax1.plot(I_ro, Inte_g_data,'bo',alpha=0.5,ms=8)
        ax1.plot(I_fit, Mgg,'b',alpha=0.8,lw=2)
        ax1.plot(I_fit, Meg,'--b',alpha=0.8,lw=2)
        ax1.plot(I_fit, Mgg+Meg,'--k',alpha=1,lw=1)
        ax1.plot(I_ro, Inte_e_data,'ro',alpha=0.5,ms=8)
        ax1.plot(I_fit, Mee,'r',alpha=0.8,lw=2)
        ax1.plot(I_fit, Mge,'--r',alpha=0.8,lw=2)
        ax1.plot(I_fit, Mge+Mee,'--k',alpha=1,lw=1)
    else: raise KeyError ('Incorrect statement of Plot_type')
    if y_scale=='log':
        ax1.set_yscale('log')
        ax1.set_ylim(np.max(Inte_g_data)*1e-3,np.max(Inte_g_data)*10)
    elif y_scale=='linear': 
        pass
    else: raise KeyError ('Incorrect statement of y_scale')
    ax.scatter(0,0,c='k',s=15)
    ax.scatter(ce_I,ce_Q,c='k',s=15)
    ax.add_patch(Ellipse(xy=[0,0],width=sig*4,height=sig*4,fill=False, alpha=0.8, facecolor= None, edgecolor="k", linewidth=0.8, linestyle='--',angle=0))
    ax.add_patch(Ellipse(xy=[ce_I,ce_Q],width=sig*4,height=sig*4,fill=False, alpha=0.8, facecolor= None, edgecolor="k", linewidth=0.8, linestyle='--',angle=0))
    ax.set_xlim(1000*min(Ig),1000*max(Ie))
    ax.set_ylim(1000*np.minimum(min(Qg),min(Qe)),1000*np.maximum(max(Qg),max(Qe)))
    ax.set_xlabel(r"$I\ $(mV)",size ='15')
    ax.set_ylabel(r"$Q\ $(mV)",size ='15')
    ax.set_xlim(-ce_I*4,ce_I*4)
    ax.set_ylim(-ce_I*4,ce_I*4)
    ax.set_title('Single shot rotated data')
    fig.tight_layout()
    ax1.set_xlabel(r"$I^{'}\ $(mV)",size ='15')
    ax1.set_ylabel(r'$PDF$',size ='15')
    ax1.set_title('Single shot rotated data')
    ax1.set_xlim(-4*sig,ce_I+4*sig)
    fig1.tight_layout()
    plt.show()


def Single_shot_fit_plot(results:dict,title_qubit:str=None):
    c_I,c_Q,sig=1000*results['fit_pack'][0],1000*results['fit_pack'][1],1000*results['fit_pack'][2]
    I,Q= results['data'][0],results['data'][1]
    
    fig, ax = plt.subplots(nrows =1,figsize =(3,3),dpi =200)
    ax.scatter(1000*I, 1000*Q, color="blue", alpha=0.5, s=5)       
    ax.scatter(c_I,c_Q,c='k',s=15)
    ax.add_patch(Ellipse(xy=[c_I,c_Q],width=sig*4,height=sig*4,fill=False, alpha=0.8, facecolor= None, edgecolor="k", linewidth=0.8, linestyle='--',angle=0))
    ax.set_xlabel(r"$I\ $(mV)",size ='15')
    ax.set_ylabel(r"$Q\ $(mV)",size ='15')
    ax.set_title(f'{title_qubit if title_qubit is not None else ""} Single shot raw data')
    ax.set_xlim(1000*np.minimum(min(I),min(Q)),1000*np.maximum(max(I),max(Q)))
    ax.set_ylim(1000*np.minimum(min(I),min(Q)),1000*np.maximum(max(I),max(Q)))
    fig.tight_layout()
    
    fig, ax = plt.subplots(nrows =1,figsize =(5,4),dpi =200)
    cmap = plt.get_cmap('jet')
    vmax= np.max(results['data_hist'])
    bounds = np.linspace(0,vmax,256)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    pcm = ax.pcolormesh(1000*results['coords'][0],1000*results['coords'][1], results['data_hist'].transpose(), norm=norm, cmap=cmap,shading='auto')
    cbar =fig.colorbar(pcm, ax=ax, extend='both', orientation='vertical')
    cbar.ax.tick_params(labelsize=10)
    ax.set_xlabel(r"$I\ $(mV)",size ='15')
    ax.set_ylabel(r"$Q\ $(mV)",size ='15')
    cbar.set_label(r'$PDF$',size ='10')
    ax.set_title('Single shot data histogram')
    ax.scatter(c_I,c_Q,c='k',s=15)
    ax.add_patch(Ellipse(xy=[c_I,c_Q],width=sig*4,height=sig*4,fill=False, alpha=0.8, facecolor= None, edgecolor="w", linewidth=0.8, linestyle='--',angle=0))
    fig.tight_layout()
    
    fig, ax = plt.subplots(nrows =1,figsize =(5,4),dpi =200)
    cmap = plt.get_cmap('jet')
    vmax= np.max(results['fitting'])
    bounds = np.linspace(0,vmax,256)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    pcm = ax.pcolormesh(1000*results['coords'][0],1000*results['coords'][1], results['fitting'], norm=norm, cmap=cmap,shading='auto')
    cbar =fig.colorbar(pcm, ax=ax, extend='both', orientation='vertical')
    cbar.ax.tick_params(labelsize=10)
    ax.set_xlabel(r"$I\ $(mV)",size ='15')
    ax.set_ylabel(r"$Q\ $(mV)",size ='15')
    cbar.set_label(r'$PDF$',size ='10')
    ax.set_title('Single shot fitting')
    ax.scatter(c_I,c_Q,c='k',s=15)
    ax.add_patch(Ellipse(xy=[c_I,c_Q],width=sig*4,height=sig*4,fill=False, alpha=0.8, facecolor= None, edgecolor="w", linewidth=0.8, linestyle='--',angle=0))
    fig.tight_layout()
    plt.show()
    
def Qubit_state_Avgtimetrace_plot(results:dict,fc:float,Digital_downconvert:bool,IF:float):
    fc=fc
    trace_recordlength= results['trace_recordlength']
    offset_Ig,offset_Qg= np.mean(results['g'][0][-100:-1]), np.mean(results['g'][1][-100:-1])
    offset_Ie,offset_Qe= np.mean(results['e'][0][-100:-1]), np.mean(results['e'][1][-100:-1])
    raw_Ig,raw_Qg= results['g'][0]-offset_Ig, results['g'][1]-offset_Qg
    raw_Ie,raw_Qe= results['e'][0]-offset_Ie, results['e'][1]  -offset_Qe 
    time_array= np.linspace(0,trace_recordlength,int(trace_recordlength*1e9))
    
    
    if Digital_downconvert:
        dw_Ig,dw_Qg= Digital_down_convert(raw_Ig,raw_Qg,IF=IF,time_array=time_array)
        dw_Ie,dw_Qe= Digital_down_convert(raw_Ie,raw_Qe,IF=IF,time_array=time_array)
        Ig,Qg= Trace_filtering(dw_Ig,fc), Trace_filtering(dw_Qg,fc)
        Ie,Qe= Trace_filtering(dw_Ie,fc), Trace_filtering(dw_Qe,fc)
    else:
         Ig,Qg= raw_Ig, raw_Qg
         Ie,Qe= raw_Ie, raw_Qe

    trace= np.linspace(0,trace_recordlength*1e9,int(trace_recordlength*1e9))
    
    fig,ax= plt.subplots(nrows =2,figsize =(6,4),dpi =200)
    ax[0].plot(trace/1000, Ig*1000,'b',alpha=0.8,lw=2)
    ax[0].plot(trace/1000, Ie*1000,'r',alpha=0.8,lw=2)
    ax[0].set_ylabel(r"$I\ $(mV)",size ='15')
    ax[0].set_title('Avg_IQ_timetrace')

    ax[1].plot(trace/1000, Qg*1000,'b',alpha=0.8,lw=2)
    ax[1].plot(trace/1000, Qe*1000,'r',alpha=0.8,lw=2)
    ax[1].set_xlabel(r"$t\ (\mu$s)",size ='15')
    ax[1].set_ylabel(r"$Q\ $(mV)",size ='15')
    fig.tight_layout()
    plt.show()
    
    return dict(Ig=Ig,Qg=Qg,Ie=Ie,Qe=Qe)
    
def Fit_analysis_plot(results:xr.core.dataset.Dataset, P_rescale:bool, Dis:any, save_path:str='', spin_echo=False, q:str='', fq_MHz:int=None):
    if P_rescale is not True:
        Nor_f=1/1000
        y_label= 'Contrast'+' [mV]'
    elif P_rescale is True:
        Nor_f= Dis
        y_label= r"$P_{1}\ $"
    else: raise KeyError ('P_rescale is not bool') 
    plot_fit:bool = True
    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250)
    text_msg = "Fit results\n"
    if results.attrs['exper'] == 'QS':
        title= 'Two tone spectroscopy'
        x_label= 'Frequency'+' [GHz]'
        x= results.coords['f']*1e-9
        x_fit= results.coords['para_fit']*1e-9
        text_msg += r"$f_{01}= %.4f $"%(results.attrs['f01_fit']*1e-9) +' GHz\n'
        text_msg += r"$BW= %.2f $"%(results.attrs['bandwidth']*1e-6) +' MHz\n'
        
    elif results.attrs['exper'] == 'T1': 
        if spin_echo:
            title = "Spin echo" if q =="" else f"{q} Spin echo"
            text_msg += r"$T_{2}= %.3f $"%(results.attrs['T1_fit']*1e6) +r"$\ [\mu$s]"+'\n'
        else:
            title= 'T1 relaxation' if q =="" else f"{q} T1 relaxation"
            text_msg += r"$T_{1}= %.3f $"%(results.attrs['T1_fit']*1e6) +r"$\ [\mu$s]"+'\n'
        x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
        x= results.coords['freeDu']*1e6
        x_fit= results.coords['para_fit']*1e6
        
        if np.array(results.data_vars['fitting'])[-1] < np.array(results.data_vars['fitting'])[0]:
            ans = np.exp(-1)*(np.array(results.data_vars['fitting']).max()-np.array(results.data_vars['fitting']).min())/Nor_f+np.array(results.data_vars['fitting']).min()/Nor_f
        else:
            ans = np.array(results.data_vars['fitting']).max()/Nor_f-np.exp(-1)*(np.array(results.data_vars['fitting']).max()-np.array(results.data_vars['fitting']).min())/Nor_f
        ax.axhline(y=ans,linestyle='--',xmin=np.array(x).min(),xmax=np.array(x).max(),color="#DCDCDC")
    elif results.attrs['exper'] == 'T2':  
        title= 'Ramsey' if q =="" else f"{q} Ramsey"
        if fq_MHz is not None: title += f" @ fq = {int(fq_MHz)} MHz"
        x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
        x= results.coords['freeDu']*1e6
        x_fit= results.coords['para_fit']*1e6
        text_msg += r"$T_{2}^{*}= %.3f $"%(results.attrs['T2_fit']*1e6) +r"$\ [\mu$s]"+'\n'        
        text_msg += r"$detuning= %.3f $"%(results.attrs['f']*1e-6) +' MHz\n'
        ans = np.exp(-1)*(np.array(results.data_vars['fitting'])[0]-np.array(results.data_vars['fitting'])[-1])/Nor_f+np.array(results.data_vars['fitting'])[-1]/Nor_f
        ax.axhline(y=ans,linestyle='--',xmin=np.array(x).min(),xmax=np.array(x).max(),color="#DCDCDC")
    
    elif results.attrs['exper'] == 'xyfcali':  
        title= 'XYF-calibration' if q =="" else f"{q} XYF-calibration"
        if fq_MHz is not None: title += f" @ fq = {int(fq_MHz)} MHz"
        x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
        x= results.coords['freeDu']*1e6
        x_fit= results.coords['para_fit']*1e6   
        text_msg += r"$detuning= %.3f $"%(results.attrs['f']*1e-6) +' MHz\n'
        
    elif results.attrs['exper'] == 'Rabi': 
        title= results.attrs['Rabi_type']
        
        if title=='PowerRabi':
            if q != "": title = f"{q} {title}"
            pi_2= results.attrs['pi_2']
            x_unit= r"$\ [V]$"
            x_label= r"$XY\ amp$"+x_unit
            x= results.coords['samples']
            x_fit= results.coords['para_fit']
            
        elif title=='TimeRabi':
            if q != "": title = f"{q} {title}"
            pi_2= results.attrs['pi_2']*1e9
            x_unit= r"$\ [ns]$"
            x_label= r"$XY\ duration$"+x_unit
            x= results.coords['samples']*1e9
            x_fit= results.coords['para_fit']*1e9

        if abs(float(pi_2)) > max(x):
            plot_fit = False

        text_msg += r"$\pi= %.3f $"%(pi_2) +x_unit       
        if plot_fit: 
            ax.axvline(x=pi_2, ymin=0, ymax= np.array(results.data_vars['data']/Nor_f).max(),color='r',linestyle='dashed', alpha=0.8,lw=1)
        
    ax.plot(x,results.data_vars['data']/Nor_f,'-', color="blue",label=r"$data$", alpha=0.5, ms=4)
    if plot_fit:
        ax.plot(x_fit,results.data_vars['fitting']/Nor_f,'-', color="red",label=r"$fit$", alpha=0.8, lw=1)       
    ax.set_xlabel(x_label)
    ax.set_title(title)
    ax.set_ylabel(y_label)
    plot_textbox(ax,text_msg)
    fig.tight_layout()
    if save_path != "":
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()

def Fit_T2_cali_analysis_plot(all_ramsey_results:list, P_rescale:bool, Dis:any):
    if P_rescale is not True:
        Nor_f=1/1000
        y_label= 'Contrast'+' [mV]'
    elif P_rescale is True:
        Nor_f= Dis
        y_label= r"$P_{1}\ $"
    else: raise KeyError ('P_rescale is not bool') 
    
    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250)
    text_msg = "Fit results\n"

    if all_ramsey_results[0].attrs['exper'] == 'T2':  
        title= 'Ramsey'
        x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
        for i in all_ramsey_results:
            x= i.coords['freeDu']*1e6
            x_fit= i.coords['para_fit']*1e6
            text_msg += r"$T_{2}= %.3f $"%(i.attrs['T2_fit']*1e6) +r"$\ [\mu$s]"+'\n'        
            text_msg += r"$detuning= %.3f $"%(i.attrs['f']*1e-6) +' MHz\n'
            ax.plot(x,i.data_vars['data']/Nor_f,'-',label=r"$data$", alpha=0.5, ms=4)
            ax.plot(x_fit,i.data_vars['fitting']/Nor_f,'-',label=r"$fit$", alpha=0.8, lw=1) 

    ax.set_xlabel(x_label)
    ax.set_title(title)
    ax.set_ylabel(y_label)
    plot_textbox(ax,text_msg)
    fig.tight_layout()
    plt.show()

def Fit_T2_cali_analysis_plot(all_ramsey_results:list, P_rescale:bool, Dis:any):
    if P_rescale is not True:
        Nor_f=1/1000
        y_label= 'Contrast'+' [mV]'
    elif P_rescale is True:
        Nor_f= Dis
        y_label= r"$P_{1}\ $"
    else: raise KeyError ('P_rescale is not bool') 
    
    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250)
    text_msg = "Fit results\n"

    if all_ramsey_results[0].attrs['exper'] == 'T2':  
        title= 'Ramsey'
        x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
        for i in all_ramsey_results:
            x= i.coords['freeDu']*1e6
            x_fit= i.coords['para_fit']*1e6
            text_msg += r"$T_{2}= %.3f $"%(i.attrs['T2_fit']*1e6) +r"$\ [\mu$s]"+'\n'        
            text_msg += r"$detuning= %.3f $"%(i.attrs['f']*1e-6) +' MHz\n'
            ax.plot(x,i.data_vars['data']/Nor_f,'-',label=r"$data$", alpha=0.5, ms=4)
            ax.plot(x_fit,i.data_vars['fitting']/Nor_f,'-',label=r"$fit$", alpha=0.8, lw=1) 

    ax.set_xlabel(x_label)
    ax.set_title(title)
    ax.set_ylabel(y_label)
    plot_textbox(ax,text_msg)
    fig.tight_layout()
    plt.show()

from numpy import array, max, mean, min
def twotone_comp_plot(results:xr.core.dataset.Dataset,substrate_backgroung:any=[], turn_on:bool=False, save_path='', title_q:str=''):
    fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250)
    ax.axvline(x=results.attrs['f01_fit']*1e-9, ymin=0, ymax= max(results.data_vars['data'].to_numpy())*1000,color='green',linestyle='dashed', alpha=0.8,lw=1)
    text_msg = "Fit results\n"
    title= 'Two tone spectroscopy'
    x_label= 'Frequency'+' [GHz]'
    y_label= 'Contrast'+' [mV]'

    x= results.coords['f']*1e-9
    x_fit= results.coords['para_fit']*1e-9

    text_msg += r"$f_{01}= %.4f $"%(results.attrs['f01_fit']*1e-9) +' GHz\n'
    text_msg += r"$BW= %.2f $"%(results.attrs['bandwidth']*1e-6) +' MHz\n'
    ax.plot(x,results.data_vars['data']*1000,'-', color="red",label=r"$data$", alpha=0.75, ms=4)
    cond_1 = max(array(results.data_vars['fitting'])) < max(results.data_vars['data'])+0.5*(max(results.data_vars['data']) - mean(array(results.data_vars['data'])))
    cond_2 = min(array(results.data_vars['fitting'])) > min(results.data_vars['data'])-0.5*(mean(results.data_vars['data']) - min(array(results.data_vars['data'])))
    if cond_1 and cond_2:
        ax.plot(x_fit,results.data_vars['fitting']*1000,'-', color="brown",label=r"$fit$", alpha=1, lw=1)  
         
    if array(substrate_backgroung).shape[0] != 0 :
        ax.plot(x,substrate_backgroung*1000,'-', color="blue",label=r"$background$", alpha=0.5, ms=4)
    
    
    ax.set_xlabel(x_label)
    ax.set_title(title)
    ax.set_ylabel(y_label)
    plot_textbox(ax,text_msg)
    if title_q != "":
        plt.title(f" {title_q} 2tone results")
    fig.tight_layout()
    plt.legend()
    if save_path != '':
        plt.savefig(save_path)
    if turn_on:
        plt.show()
    else:
        plt.close()

#%%    
def set_LO_frequency(quantum_device:QuantumDevice,q:str,module_type:str,LO_frequency:float):
    
    qubit= quantum_device.get_element(q)
    if module_type== 'drive':
        clock=qubit.name + ".01"
        port=qubit.ports.microwave()
        
    elif module_type== 'readout':
        clock=qubit.name + ".ro"
        port= "q:res"#qubit.ports.readout()
        
    else: raise KeyError ('module_type is not drive or readout')  
    
    hw_config = quantum_device.hardware_config()
    
    output_path = find_port_clock_path(
        hw_config, port=port, clock= clock)
    
    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    
    module = hw_config[cluster_key][module_key]
    module[output_key]["lo_freq"] = LO_frequency
    
    if "interm_freq" in module[output_key]["portclock_configs"][0]:
        del module[output_key]["portclock_configs"][0]["interm_freq"]
    else: pass

        
    quantum_device.hardware_config(hw_config)








