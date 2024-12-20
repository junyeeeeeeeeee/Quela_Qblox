from numpy import array, mean, median, argmax, linspace, arange, moveaxis, empty_like, std, average, transpose, where, arctan2, sort, polyfit, delete, degrees
from numpy import sqrt, pi
from numpy import ndarray
from xarray import Dataset, DataArray, open_dataset
from qcat.analysis.base import QCATAna
import matplotlib.pyplot as plt
import os, sys, traceback
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import quantify_core.data.handling as dh
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support.QDmanager import QDmanager, Data_manager
from qblox_drive_AS.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis, T2_fit_analysis, Fit_analysis_plot, T1_fit_analysis, cos_fit_analysis, IQ_data_dis, twotone_comp_plot, gate_phase_fit_analysis
from qblox_drive_AS.support.QuFluxFit import remove_outlier_after_fit
from scipy.optimize import curve_fit
from qblox_drive_AS.support.QuFluxFit import remove_outliers_with_window
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.analysis.state_discrimination import p01_to_Teff
from qcat.analysis.state_discrimination.discriminator import get_proj_distance
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from qcat.analysis.resonator.photon_dep.res_data import ResonatorData
from qblox_drive_AS.support import rotate_onto_Inphase, rotate_data
from qcat.visualization.qubit_relaxation import plot_qubit_relaxation
from qcat.analysis.qubit.relaxation import qubit_relaxation_fitting
from datetime import datetime
from matplotlib.figure import Figure

def parabola(x,a,b,c):
    return a*array(x)**2+b*array(x)+c


def find_minima(f, phase, start, stop):
    """ from the given astart and stop, which is counted for the minimua number  """
    minima = []
    for n in range(start, stop):
        x_min = ((2 * n + 1) * pi - phase) / (2 * pi * f)
        minima.append(x_min)
    return minima

def plot_FreqBiasIQ(f:ndarray,z:ndarray,I:ndarray,Q:ndarray,refIQ:list=[], ax:plt.Axes=None)->plt.Axes:
    """
        I and Q shape in (z, freq)
    """
    if len(refIQ) != 2:
        refIQ = [0,0]
        
    data = sqrt((I-refIQ[0])**2+(Q-refIQ[1])**2)
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,8))
        ax:plt.Axes
    c = ax.pcolormesh(z, f*1e-9, data.transpose(), cmap='RdBu')
    ax.set_xlabel("Flux Pulse amp (V)", fontsize=20)
    ax.set_ylabel("Frequency (GHz)", fontsize=20)
    fig.colorbar(c, ax=ax, label='Contrast (V)')
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    
    return ax

def zgate_T1_fitting(time_array:DataArray,IQ_array:DataArray, ref_IQ:list, fit:bool=True)->tuple[list,list]:
    I_array = list(IQ_array[0]) # shape in (bias, time)
    Q_array = list(IQ_array[1]) # shape in (bias, time)

    T1s = []
    signals = []
    if len(ref_IQ) == 1:
        data = moveaxis(IQ_array,0,1)
        for bias_idx, bias_data in enumerate(data): 
            signals.append(rotate_data(bias_data,ref_IQ[0])[0])
            if fit:
                T1s.append(qubit_relaxation_fitting(array(time_array),signals[-1]).params["tau"].value)
    else:
        for bias_idx, bias_data in enumerate(I_array):
            signals.append(sqrt((array(bias_data)-ref_IQ[0])**2+(array(Q_array[bias_idx])-ref_IQ[1])**2))
            if fit:
                T1s.append(qubit_relaxation_fitting(array(time_array),signals[-1]).params["tau"].value)

    return signals, T1s

def build_result_pic_path(dir_path:str,folder_name:str="")->str:
    parent = os.path.split(dir_path)[0]
    new_path = os.path.join(parent,"ZgateT1_pic" if folder_name == "" else folder_name)
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    return new_path

def anal_rof_cali(I_e:ndarray,Q_e:ndarray,I_g:ndarray,Q_g:ndarray,dis_diff:ndarray,ro_f_samples:ndarray):
    max_dif_idx = where(dis_diff==max(dis_diff))[0][0]
    optimized_rof = ro_f_samples[max_dif_idx]
    fig, ax = plt.subplots(3,1)
    amp_e = sqrt(I_e**2+Q_e**2)
    amp_g = sqrt(I_g**2+Q_g**2)
    ax[0].plot(ro_f_samples,amp_e,label='|1>')
    ax[0].plot(ro_f_samples,amp_g,label='|0>')
    ax[0].set_xlabel('ROF (Hz)')
    ax[0].set_ylabel("Amplitude (V)")
    ax[0].legend()

    pha_e = arctan2(Q_e,I_e)
    pha_g = arctan2(Q_g,I_g)
    ax[1].plot(ro_f_samples,pha_e,label='|1>')
    ax[1].plot(ro_f_samples,pha_g,label='|0>')
    ax[1].set_xlabel('ROF (Hz)')
    ax[1].set_ylabel("phase")
    ax[1].legend()

    ax[2].plot(ro_f_samples,dis_diff,label='diff')
    ax[2].vlines(x=array([optimized_rof]),ymin=min(dis_diff),ymax=max(dis_diff),colors='black',linestyles='--',label='optimal')
    ax[2].vlines(x=array([mean(ro_f_samples)]),ymin=min(dis_diff),ymax=max(dis_diff),colors='#DCDCDC',linestyles='--',label='original')
    ax[2].set_xlabel('ROF (Hz)')
    ax[2].set_ylabel("diff (V)")
    ax[2].legend()
    plt.show()

def sort_timeLabel(files):
    sorted_files = sorted(files, key=lambda x: datetime.strptime(x.split("_")[-1].split(".")[0], "%Y%m%d%H%M%S"))
    return sorted_files



class Artist():
    def __init__(self,pic_title:str,save_pic_path:str=None):
        self.title:str = pic_title
        self.pic_save_path:str = save_pic_path
        self.title_fontsize:int = 22
        self.axis_subtitle_fontsize:int = 19
        self.axis_label_fontsize:int = 19
        self.tick_number_fontsize:int = 16
        self.color_cndidates:list = ["#1E90FF","#D2691E","#FFA500","#BDB76B","#FF0000","#0000FF","#9400D3"]
        self.axs:list = []

    def export_results(self):
        import warnings
        warnings.filterwarnings("ignore", message="No artists with labels found to put in legend")
        plt.legend()
        plt.tight_layout()
        if self.pic_save_path is not None:
            slightly_print(f"pic saved located:\n{self.pic_save_path}")
            plt.savefig(self.pic_save_path)
            plt.close()
        else:
            plt.show()

    def includes_axes(self,axes:list):
        for ax in axes:
            self.axs.append(ax)
    
    def set_LabelandSubtitles(self,info:list):
        """ 
        The arg `info` is a list contains multiples of dict. each dict must have three key name: 'subtitle', 'xlabel' and 'ylabel'.\n
        The order is compatible with the arg `axes` of method `self.includes_axes`. So, the first dict in `info` will be fill into the first ax in axes.
        """
        for idx, element in enumerate(info):
            ax:plt.Axes = self.axs[idx]
            if element['subtitle'] != "": ax.set_title(element['subtitle'],fontsize=self.axis_subtitle_fontsize)
            if element['xlabel'] != "": ax.set_xlabel(element['xlabel'],fontsize=self.axis_label_fontsize)
            if element['ylabel'] != "": ax.set_ylabel(element['ylabel'],fontsize=self.axis_label_fontsize)
            ax.legend()
 
    def build_up_plot_frame(self,subplots_alignment:tuple=(3,1),fig_size:tuple=None)->tuple[Figure,list]:
        fig, axs = plt.subplots(subplots_alignment[0],subplots_alignment[1],figsize=(subplots_alignment[0]*9,subplots_alignment[1]*6) if fig_size is None else fig_size)
        plt.title(self.title,fontsize=self.title_fontsize)
        if subplots_alignment[0] == 1 and subplots_alignment[1] == 1:
            axs = [axs]
        for ax in axs:
            ax:plt.Axes
            ax.grid()
            ax.xaxis.set_tick_params(labelsize=self.tick_number_fontsize)
            ax.yaxis.set_tick_params(labelsize=self.tick_number_fontsize)
            ax.xaxis.label.set_size(self.axis_label_fontsize) 
            ax.yaxis.label.set_size(self.axis_label_fontsize)
        return fig, axs
    
    def add_colormesh_on_ax(self,x_data:ndarray,y_data:ndarray,z_data:ndarray,fig:Figure,ax:plt.Axes,colorbar_name:str=None)->plt.Axes:
        im = ax.pcolormesh(x_data,y_data,transpose(z_data),shading='nearest')
        fig.colorbar(im, ax=ax, label="" if colorbar_name is None else colorbar_name)
        return ax

    def add_scatter_on_ax(self,x_data:ndarray,y_data:ndarray,ax:plt.Axes,**kwargs)->plt.Axes:
        ax.scatter(x_data,y_data,**kwargs)
        return ax
    
    def add_plot_on_ax(self,x_data:ndarray,y_data:ndarray,ax:plt.Axes,**kwargs)->plt.Axes:
        ax.plot(x_data,y_data,**kwargs)
        return ax
    
    def add_verline_on_ax(self,x:float,y_data:ndarray,ax:plt.Axes,**kwargs)->plt.Axes:
        ax.vlines(x,min(y_data),max(y_data),**kwargs)
        return ax
    
    def add_horline_on_ax(self,x_data:float,y:float,ax:plt.Axes,**kwargs)->plt.Axes:
        ax.hlines(y,min(x_data),max(x_data),**kwargs)
        return ax


#? Analyzer use tools to analyze data
################################
####     Analysis tools    ####
################################
class analysis_tools():
    def __init__(self):
        self.ds:Dataset = None
        self.fit_packs = {}
    
    def import_dataset(self,ds:Dataset):
        self.ds = ds  # inheritance 


    
    def fluxCoupler_ana(self, var_name:str, refIQ:list=[]):
        self.qubit = var_name
        self.freqs = {}
        self.freqs[var_name] = array(self.ds.data_vars[f"{var_name}_freq"])[0][0]
        self.bias = array(self.ds.coords["bias"])
        self.refIQ = refIQ
        
    def fluxCoupler_plot(self, save_pic_path:str=None):
        I_data = array(self.ds.data_vars[self.qubit])[0]
        Q_data = array(self.ds.data_vars[self.qubit])[1]
        ax = plot_FreqBiasIQ(self.freqs[self.qubit],self.bias,I_data,Q_data,self.refIQ)
        ax.set_title(f"{self.qubit} Readout",fontsize=20)
        ax.set_xlabel(f"couplers {self.ds.attrs['cntrl_couplers'].replace('_', ' & ')} bias amplitude (V)")
        plt.grid()
        plt.tight_layout()
        if save_pic_path is None:
            plt.show()
        else:
            pic_path = os.path.join(save_pic_path,f"{self.qubit}_FluxCoupler_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png")
            slightly_print(f"pic saved located:\n{pic_path}")
            plt.savefig(pic_path)
            plt.close()

    def fluxCavity_ana(self, var_name:str):

        self.qubit = var_name
        self.freqs = array(self.ds.data_vars[f"{var_name}_freq"])[0][0]
        self.bias = array(self.ds.coords["bias"])
        S21 = array(self.ds.data_vars[self.qubit])[0] + 1j * array(self.ds.data_vars[self.qubit])[1]
        self.mags:ndarray = abs(S21)
        self.collected_freq = []
        self.collected_flux = []
        try:
            freq_fit = []
            fit_err = []
            for idx, _ in enumerate(S21):
                res_er = ResonatorData(freq=self.freqs,zdata=array(S21[idx]))
                result, _, _ = res_er.fit()
                freq_fit.append(result['fr'])
                fit_err.append(result['chi_square'])

            _, indexs = remove_outliers_with_window(array(fit_err),int(len(fit_err)/3),m=1)
            
            for i in indexs:
                self.collected_freq.append(freq_fit[i])
                self.collected_flux.append(list(self.bias)[i])
            self.fit_results = cos_fit_analysis(array(self.collected_freq),array(self.collected_flux))
            paras = array(self.fit_results.attrs['coefs'])

            from scipy.optimize import minimize
            from qblox_drive_AS.support.Pulse_schedule_library import Cosine_func
            def cosine(x):
                return -Cosine_func(x,*paras)

            sweet_spot = minimize(cosine,x0=0,bounds=[(min(self.bias),max(self.bias))])


            sweet_flux = float(sweet_spot.x[0])
            sweet_freq = -sweet_spot.fun

            self.fit_packs = {"A":paras[0],"f":paras[1],"phi":paras[2],"offset":paras[3],"sweet_freq":sweet_freq,"sweet_flux":sweet_flux}
        except Exception as err:
            print(f"While fitting got error = {err}")
            traceback.print_exc()
            self.fit_packs = {}

    def fluxCavity_plot(self,save_pic_folder):
        Plotter = Artist(pic_title=f"FluxCavity-{self.qubit}")
        fig, axs = Plotter.build_up_plot_frame([1,1])
        ax = Plotter.add_colormesh_on_ax(self.bias,self.freqs,self.mags,fig,axs[0])
        ax.set_xlim(min(self.bias), max(self.bias))
        ax.set_ylim(min(self.freqs),max(self.freqs))

        if len(list(self.fit_packs.keys())) != 0:
            fit_x = self.fit_results.coords['para_fit']
            fit_y = self.fit_results.data_vars['fitting']
            ax = Plotter.add_scatter_on_ax(self.fit_packs["sweet_flux"],self.fit_packs["sweet_freq"],ax,c='red',marker="*",s=300)
            ax = Plotter.add_plot_on_ax(fit_x,fit_y,ax,c='red')
        
        if len(self.collected_flux) != 0 and len(self.collected_freq) != 0:
            ax = Plotter.add_scatter_on_ax(self.collected_flux,self.collected_freq,ax,c='white',marker="o",s=30)

        Plotter.includes_axes([ax])
        Plotter.set_LabelandSubtitles([{"subtitle":"","xlabel":f"{self.qubit} flux (V)","ylabel":"RO freqs (Hz)"}])
        Plotter.pic_save_path = os.path.join(save_pic_folder,f"FluxCavity_{self.qubit}.png")
        Plotter.export_results()

    def conti2tone_ana(self, var, fit_func:callable=None, ref:list=[]):
        self.target_q = var
        self.xyf = array(self.ds[f"{self.target_q}_freq"])[0][0]
        self.xyl = array(self.ds.coords["xy_amp"])
        ii = array(self.ds[self.target_q])[0]
        qq = array(self.ds[self.target_q])[1]
        if len(ref) == 2:
            self.contrast = sqrt((ii-ref[0])**2+(qq-ref[1])**2).reshape(self.xyl.shape[0],self.xyf.shape[0])
        else:
            self.contrast = rotate_data(array(self.ds[self.target_q]),ref[0])[0]
        self.fit_f01s = []
        self.fif_amps = []
        if fit_func is not None and self.xyl.shape[0] != 1:
            for xyl_idx, an_amp_data in enumerate(self.contrast):
                if self.xyl[xyl_idx] != 0:
                    if min(self.xyf) <= fit_func(an_amp_data,self.xyf).attrs['f01_fit'] and fit_func(an_amp_data,self.xyf).attrs['f01_fit'] <= max(self.xyf):
                        self.fit_f01s.append(fit_func(an_amp_data,self.xyf).attrs['f01_fit']) # Hz
                        self.fif_amps.append(self.xyl[xyl_idx])
        
        if self.xyl.shape[0] == 1:
            self.fit_packs[self.target_q] = {"xyf_data":self.xyf,"contrast":self.contrast[0]}
    
    def conti2tone_plot(self, save_pic_path:str=None):
        Plotter = Artist(pic_title=f"2Tone-{self.target_q}")
        fig, axs = Plotter.build_up_plot_frame([1,1])

        if self.xyl.shape[0] != 1:
            ax = Plotter.add_colormesh_on_ax(self.xyf,self.xyl,transpose(self.contrast),fig,axs[0])
            if len(self.fit_f01s) != 0 :
                ax = Plotter.add_scatter_on_ax(self.fit_f01s,self.fif_amps,ax,marker="*",c='red')
            ax.set_xlim(min(self.xyf), max(self.xyf))
            ax.set_ylim(min(self.xyl),max(self.xyl))
            Plotter.includes_axes([ax])
            Plotter.set_LabelandSubtitles([{"subtitle":"","xlabel":f"{self.target_q} XYF (Hz)","ylabel":"XY Power (V)"}])
        else:
            res = QS_fit_analysis(self.fit_packs[self.target_q]["contrast"],f=self.fit_packs[self.target_q]["xyf_data"])
            ax = Plotter.add_verline_on_ax(x=res.attrs['f01_fit']*1e-9, y_data=res.data_vars['data'].to_numpy()*1000, ax=axs[0], color='green',linestyle='dashed', alpha=0.8,lw=1)
            ax = Plotter.add_plot_on_ax(res.coords['f']*1e-9,res.data_vars['data']*1000,ax,linestyle='-', color="red", alpha=0.75, ms=4)
            Plotter.includes_axes([ax])
            Plotter.set_LabelandSubtitles([{"subtitle":"","xlabel":f"XYF (GHz)","ylabel":"Contrast (mV)"}])
        
    
        Plotter.pic_save_path = os.path.join(save_pic_path,f"{self.target_q}_Conti2tone_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png")
        Plotter.export_results()

    def fluxQb_ana(self, var, fit_func:callable=None, refIQ:list=[], filter_outlier:bool=True):
        self.qubit = var
        self.filtered_z, self.filtered_f = [], []
        self.ref_z = float(self.ds.attrs[f"{self.qubit}_z_ref"])
        self.f = array(self.ds[f"{self.qubit}_freq"])[0][0]
        self.z = array(self.ds.coords["bias"])
        IQarray = array(self.ds[f"{self.qubit}"])


        if len(refIQ) == 2:
            self.contrast = array(sqrt((IQarray[0]-refIQ[0])**2+(IQarray[1]-refIQ[1])**2))
        else:
            self.contrast = rotate_data(IQarray,refIQ[0])[0]
        if fit_func is not None:
            # try:
            self.fit_z = []
            self.fit_f = []
            def parabola(x,a,b,c):
                return a*array(x)**2+b*array(x)+c
            for z_idx, a_z_data in enumerate(self.contrast):
                if min(self.f) <= fit_func(a_z_data,self.f).attrs['f01_fit'] and fit_func(a_z_data,self.f).attrs['f01_fit'] <= max(self.f):
                    self.fit_f.append(fit_func(a_z_data,self.f).attrs['f01_fit']*1e-9) # GHz
                    self.fit_z.append(self.z[z_idx])
            
            if not filter_outlier:
                self.paras, _ = curve_fit(parabola,self.fit_z,self.fit_f)
            else:
                self.filtered_z, self.filtered_f, self.paras = remove_outlier_after_fit(parabola,self.fit_z,self.fit_f)
            
            self.fit_packs["sweet_bias"] = self.ref_z + float(-self.paras[1]/(2*self.paras[0])) # offset + z_pulse_amp
            self.fit_packs["xyf"] = float(parabola(float(-self.paras[1]/(2*self.paras[0])),*self.paras))*1e9
            self.fit_packs["parabola_paras"] = list(self.paras)
    
    def fluxQb_plot(self, save_pic_path:str=None):
        fig, ax = plt.subplots(figsize=(13,9))
        ax:plt.Axes
        c = ax.pcolormesh(self.z, self.f*1e-9, self.contrast.transpose())

        plt.scatter(self.fit_z,self.fit_f,marker='*',c='black')
        if array(self.filtered_z).shape[0] != 0 and array(self.filtered_f).shape[0] != 0:
            plt.scatter(self.filtered_z,self.filtered_f,marker='*',c='red')
        plt.plot(self.z,parabola(self.z,*self.paras),'red')
        if min(self.z)<=float(-self.paras[1]/(2*self.paras[0])) and float(-self.paras[1]/(2*self.paras[0])) <= max(self.z):
            plt.vlines(float(-self.paras[1]/(2*self.paras[0])),min(self.f)*1e-9,max(self.f)*1e-9,colors='red',linestyles='--')

        ax.set_xlabel("Flux Pulse amp (V)", fontsize=20)
        ax.set_ylabel("Driving Frequency (GHz)", fontsize=20)
        fig.colorbar(c, ax=ax, label='Contrast (V)')
        ax.xaxis.set_tick_params(labelsize=18)
        ax.yaxis.set_tick_params(labelsize=18)
        ax.set_xlim(min(self.z), max(self.z))
        ax.set_ylim(min(self.f)*1e-9,max(self.f)*1e-9)

        if len(list(self.fit_packs.keys())) != 0:
            plt.title(f"{self.qubit} XYF={round(self.fit_packs['xyf']*1e-9,3)} GHz with z_pulse amp={round(float(-self.paras[1]/(2*self.paras[0])),3)} V")
        if save_pic_path is None:
            plt.show()
        else:
            pic_path = os.path.join(save_pic_path,f"{self.qubit}_fluxQubitSpectro_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png")
            slightly_print(f"pic saved located:\n{pic_path}")
            plt.savefig(pic_path)
            plt.close()

    def rabi_ana(self,var:str):
        self.qubit = var
        self.rabi_type = self.ds.attrs["rabi_type"]
        i_data = array(self.ds[var])[0]
        q_data = array(self.ds[var])[1]
        if not self.ds.attrs["OS_mode"]:
            if len(self.refIQ) == 2:
                data_to_fit = sqrt((i_data-self.refIQ[0])**2+(q_data-self.refIQ[1])**2)
            else:
                data_to_fit = rotate_data(array(self.ds[var]),float(self.refIQ[0]))[0]
        else:
            # wait GMM
            pass

        if self.rabi_type.lower() == 'powerrabi':
            x_data = array(self.ds[f"{var}_piamp"])[0]
            fixed = self.ds.attrs[f'{var}_pidura']
        else:
            x_data = array(self.ds[f"{var}_pidura"])[0]
            fixed = self.ds.attrs[f'{var}_piamp']
    
        self.fit_packs = Rabi_fit_analysis(data_to_fit,x_data,self.rabi_type)
        self.fit_packs.attrs["fix_variable"] = fixed

    def rabi_plot(self, save_pic_path:str=None):
        save_pic_path = os.path.join(save_pic_path,f"{self.qubit}_{self.rabi_type}_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png") if save_pic_path is not None else ""
        if save_pic_path != "": slightly_print(f"pic saved located:\n{save_pic_path}")
        Fit_analysis_plot(self.fit_packs,save_path=save_pic_path,P_rescale=None,Dis=None,q=self.qubit)

    def oneshot_ana(self,data:DataArray,tansition_freq_Hz:float=None):
        self.fq = tansition_freq_Hz
        self.gmm2d_fidelity = GMMROFidelity()
        self.gmm2d_fidelity._import_data(data)
        self.gmm2d_fidelity._start_analysis()
        g1d_fidelity = self.gmm2d_fidelity.export_G1DROFidelity()
        
        p00 = g1d_fidelity.g1d_dist[0][0][0]
        self.thermal_populations = g1d_fidelity.g1d_dist[0][0][1]
        p11 = g1d_fidelity.g1d_dist[1][0][1]
        if self.fq is not None:
            self.effT_mK = p01_to_Teff(self.thermal_populations, self.fq)*1000
        else:
            self.effT_mK = 0
        self.RO_fidelity_percentage = (p00+p11)*100/2


        _, self.RO_rotate_angle = rotate_onto_Inphase(self.gmm2d_fidelity.mapped_centers[0],self.gmm2d_fidelity.mapped_centers[1])
        z = moveaxis(array(data),0,1) # (IQ, state, shots) -> (state, IQ, shots)
        self.rotated_data = empty_like(array(data))
        for state_idx, state_data in enumerate(z):
            self.rotated_data[state_idx] = rotate_data(state_data,self.RO_rotate_angle)
        
        self.fit_packs = {"effT_mK":self.effT_mK,"thermal_population":self.thermal_populations,"RO_fidelity":self.RO_fidelity_percentage,"RO_rotation_angle":self.RO_rotate_angle}

    def oneshot_plot(self,save_pic_path:str=None):
        da = DataArray(moveaxis(self.rotated_data,0,1), coords= [("mixer",["I","Q"]), ("prepared_state",[0,1]), ("index",arange(array(self.rotated_data).shape[2]))] )
        self.gmm2d_fidelity._import_data(da)
        self.gmm2d_fidelity._start_analysis()
        g1d_fidelity = self.gmm2d_fidelity.export_G1DROFidelity()

        plot_readout_fidelity(da, self.gmm2d_fidelity, g1d_fidelity,self.fq,save_pic_path,plot=True if save_pic_path is None else False)
        plt.close()

    def T2_ana(self,var:str,ref:list):
        raw_data = self.ds[var]
        time_samples = array(self.ds[f"{var}_x"])[0][0]
        reshaped = moveaxis(array(raw_data),0,1)  # (repeat, IQ, idx)
        self.qubit = var


        self.T2_fit = []
        for idx, data in enumerate(reshaped):
            self.echo:bool=False if raw_data.attrs["spin_num"] == 0 else True
            if len(ref) == 1:
                self.data_n = rotate_data(data,ref[0])[0]
            else:
                self.data_n = sqrt((data[0]-ref[0])**2+(data[1]-ref[1])**2)
            self.plot_item = {"data":self.data_n*1000,"time":array(time_samples)*1e6}
            if not self.echo:
                self.ans = T2_fit_analysis(self.plot_item["data"]/1000,self.plot_item["time"]/1e6)
                self.T2_fit.append(self.ans.attrs["T2_fit"]*1e6)
                self.fit_packs["freq"] = self.ans.attrs['f']
            else:
                self.ans = qubit_relaxation_fitting(self.plot_item["time"],self.plot_item["data"])
                self.T2_fit.append(self.ans.params["tau"].value)
                self.fit_packs["freq"] = 0

        
        self.fit_packs["median_T2"] = median(array(self.T2_fit))
        self.fit_packs["mean_T2"] = mean(array(self.T2_fit))
        self.fit_packs["std_T2"] = std(array(self.T2_fit))
    
    def XYF_cali_ana(self,var:str,ref:list):
        raw_data = self.ds[f'{var}']
        time_samples =  array(self.ds[f'{var}_x'])[0][0]
        reshaped = moveaxis(array(raw_data),0,1)  # (repeat, IQ, idx)
        self.qubit = var
        for idx, data in enumerate(reshaped):
            self.echo:bool=False if raw_data.attrs["spin_num"] == 0 else True
            if len(ref) == 1:
                data = rotate_data(data,ref[0])[0] # I
            else:
                data = sqrt((data[0]-ref[0])**2+(data[1]-ref[1])**2)
            
            self.ans = cos_fit_analysis(data,array(time_samples))
            p_deg = degrees(self.ans.attrs['coefs'][2])  # Converts radians to degrees
            # Normalize to range [0, 360)
            self.fit_packs['phase'] = (p_deg % 360 + 360) % 360
            self.fit_packs["freq"] = self.ans.attrs['f']
            eyeson_print(f"phase fit = {round(self.fit_packs['phase'],2)} deg")
        
    def T2_plot(self,save_pic_path:str=None):
        save_pic_path = os.path.join(save_pic_path,f"{self.qubit}_{'Echo' if self.echo else 'Ramsey'}_{self.ds.attrs['execution_time'].replace(' ', '_')}.png") if save_pic_path is not None else ""
        if save_pic_path != "" : slightly_print(f"pic saved located:\n{save_pic_path}")
        if not self.echo:
            Fit_analysis_plot(self.ans,P_rescale=False,Dis=None,spin_echo=self.echo,save_path=save_pic_path,q=self.qubit)
        else:
            fig, ax = plt.subplots()
            ax = plot_qubit_relaxation(self.plot_item["time"],self.plot_item["data"], ax, self.ans)
            ax.set_title(f"{self.qubit} T2 = {round(self.ans.params['tau'].value,1)} µs" )
            if save_pic_path != "" : 
                plt.savefig(save_pic_path)
                plt.close()
            else:
                plt.show()

        if len(self.T2_fit) > 1:
            Data_manager().save_histo_pic(None,{str(self.qubit):self.T2_fit},self.qubit,mode=f"{'t2' if self.echo else 't2*'}",pic_folder=os.path.split(save_pic_path)[0])

    def XYF_cali_plot(self,save_pic_path:str=None,fq_MHz:int=None):
        save_pic_path = os.path.join(save_pic_path,f"{self.qubit}_XYFcalibration_{Data_manager().get_time_now()}.png") if save_pic_path is not None else ""
        if save_pic_path != "" : slightly_print(f"pic saved located:\n{save_pic_path}")
        Fit_analysis_plot(self.ans,P_rescale=False,Dis=None,save_path=save_pic_path,q=self.qubit,fq_MHz=fq_MHz)

    def T1_ana(self,var:str,ref:list):
        raw_data = self.ds[var]
        time_samples = array(self.ds[f"{var}_x"])[0][0]
        reshaped = moveaxis(array(raw_data),0,1)  # (repeat, IQ, idx)
        self.qubit = var
        self.plot_item = {"time":array(time_samples)*1e6}
        self.T1_fit = []
        for idx, data in enumerate(reshaped):
            if len(ref) == 1:
                self.plot_item["data"] = rotate_data(data,ref[0])[0]*1000 # I
            else:
                self.plot_item["data"] = sqrt((data[0]-ref[0])**2+(data[1]-ref[1])**2)*1000

            self.ans = qubit_relaxation_fitting(self.plot_item["time"],self.plot_item["data"])
            self.T1_fit.append(self.ans.params["tau"].value)
        self.fit_packs["median_T1"] = median(array(self.T1_fit))
        self.fit_packs["mean_T1"] = mean(array(self.T1_fit))
        self.fit_packs["std_T1"] = std(array(self.T1_fit))
    
    def T1_plot(self,save_pic_path:str=None):
        save_pic_path = os.path.join(save_pic_path,f"{self.qubit}_T1_{self.ds.attrs['execution_time'].replace(' ', '_')}.png") if save_pic_path is not None else ""
        fig, ax = plt.subplots()
        ax = plot_qubit_relaxation(self.plot_item["time"],self.plot_item["data"], ax, self.ans)
        ax.set_title(f"{self.qubit} T1 = {round(self.ans.params['tau'].value,1)} µs" )
        if save_pic_path != "" : 
            slightly_print(f"pic saved located:\n{save_pic_path}")
            plt.savefig(save_pic_path)
            plt.close()
        else:
            plt.show()
        if len(self.T1_fit) > 1:
            Data_manager().save_histo_pic(None,{str(self.qubit):self.T1_fit},self.qubit,mode="t1",pic_folder=os.path.split(save_pic_path)[0])

    def ZgateT1_ana(self,time_sort:bool=False):
        self.time_trace_mode = time_sort
        self.prepare_excited = self.ds.attrs["prepare_excited"]
        self.qubit = [var for var in self.ds.data_vars if var.split("_")[-1] != "time"][0]
        end_times = array(self.ds.coords["end_time"])
        self.evotime_array = array(self.ds[f"{self.qubit}_time"])[0][0][0]*1e6
        z_pulse_amplitudes = array(self.ds.coords["z_voltage"])
        if not time_sort:
            self.T1_per_time = []
            self.I_chennel_per_time = []  
            for time_dep_data in array(self.ds[self.qubit]): # time_dep_data shape in  (mixer, bias, evo-time)
                signals, T1s = zgate_T1_fitting(self.evotime_array,time_dep_data,self.refIQ,self.prepare_excited)
                if self.prepare_excited:
                    self.T1_per_time.append(T1s)
            self.I_chennel_per_time.append(signals)


            self.avg_I_data = average(array(self.I_chennel_per_time),axis=0)
            
            if self.prepare_excited:
                self.avg_T1 = average(array(self.T1_per_time),axis=0)
                self.std_T1_percent = array([round(i,1) for i in std(array(self.T1_per_time),axis=0)*100/self.avg_T1])

        else:
            self.T1_rec = {}
            for end_time_idx, time_dep_data in enumerate(array(self.ds[self.qubit])):
                signals, T1s = zgate_T1_fitting(self.evotime_array,time_dep_data,self.refIQ,self.prepare_excited)
                self.T1_rec[end_times[end_time_idx]] = {}
                self.T1_rec[end_times[end_time_idx]]["T1s"] = T1s
            self.T1_rec = dict(sorted(self.T1_rec.items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S")))
        
        self.z = z_pulse_amplitudes+self.ds.attrs["z_offset"]
    
    def ZgateT1_plot(self,save_pic_folder:str=None):
        Plotter = Artist(pic_title=f"zgate-T1-{self.qubit}")
        if not self.time_trace_mode:
            fig, axs = Plotter.build_up_plot_frame([1,1])
            ax = Plotter.add_colormesh_on_ax(self.z,self.evotime_array,self.avg_I_data,fig,axs[0])
            if self.prepare_excited:
                for idx, T1_per_dataset in enumerate(self.T1_per_time):
                    if idx == 0:
                        ax = Plotter.add_scatter_on_ax(x_data=self.z,y_data=T1_per_dataset,ax=ax,c='pink',s=5,label="raw data")
                    else:
                        ax = Plotter.add_scatter_on_ax(x_data=self.z,y_data=T1_per_dataset,ax=ax,c='pink',s=5)
                ax = Plotter.add_scatter_on_ax(x_data=self.z,y_data=self.avg_T1,ax=ax,c='red',marker='*',s=30,label="avg data")
            ax.set_ylim(0,max(self.evotime_array))
            Plotter.includes_axes([ax])
            Plotter.set_LabelandSubtitles([{"subtitle":"","xlabel":f"{self.qubit} flux (V)","ylabel":"Free evolution time (µs)"}])

        else:
            fig, axs = Plotter.build_up_plot_frame([1,1])
            z = self.z
            ans_dict = self.T1_rec
            times = [datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S") for time_str in ans_dict.keys()]
            time_diffs = [(t - times[0]).total_seconds() / 60 for t in times]  # Convert to minutes

            # Get values
            T1_data = [entry["T1s"] for entry in ans_dict.values()]
            T1_data = array(T1_data).reshape(z.shape[0],len(times))
            ax = Plotter.add_colormesh_on_ax(self.z,time_diffs,T1_data,fig,axs[0],"T1 (µs)")
            Plotter.includes_axes([ax])
            Plotter.set_LabelandSubtitles([{"subtitle":"","xlabel":f"{self.qubit} flux (V)","ylabel":"Time past (min)"}])

        Plotter.pic_save_path = os.path.join(build_result_pic_path(save_pic_folder),f"ZgateT1_{'zAVG' if  not self.time_trace_mode else 'TimeTrace'}_poster.png")
        Plotter.export_results()
    
    def ROF_cali_ana(self,var:str):
        self.qubit = var
        IQ_data = moveaxis(array(self.ds[var]),1,0) # shape (mixer, state, rof) -> (state, mixer, rof)
        self.rof = moveaxis(array(self.ds[f"{var}_rof"]),0,1)[0][0]
        self.I_g, self.I_e = IQ_data[0][0], IQ_data[1][0]
        self.Q_g, self.Q_e = IQ_data[0][1], IQ_data[1][1]
        I_diff = self.I_e-self.I_g
        Q_diff = self.Q_e-self.Q_g
        self.dis_diff = sqrt((I_diff)**2+(Q_diff)**2)
        self.fit_packs[var] = {"optimal_rof":self.rof[where(self.dis_diff==max(self.dis_diff))[0][0]]}

    def ROF_cali_plot(self,save_pic_folder:str=None):
        
        if save_pic_folder is not None: save_pic_folder = os.path.join(save_pic_folder,f"{self.qubit}_ROF_cali_{self.ds.attrs['execution_time']}.png")
        Plotter = Artist(pic_title=f"{self.qubit}_ROF_calibration",save_pic_path=save_pic_folder)
        fig, axs = Plotter.build_up_plot_frame((3,1),fig_size=(12,9))
        ax0:plt.Axes = axs[0]
        Plotter.add_plot_on_ax(self.rof,sqrt(self.I_g**2+self.Q_g**2),ax0,label="|0>")
        Plotter.add_plot_on_ax(self.rof,sqrt(self.I_e**2+self.Q_e**2),ax0,label="|1>")
        ax1:plt.Axes = axs[1]
        Plotter.add_plot_on_ax(self.rof,arctan2(self.Q_g,self.I_g)/pi,ax1,label="|0>")
        Plotter.add_plot_on_ax(self.rof,arctan2(self.Q_e,self.I_e)/pi,ax1,label="|1>")
        ax2:plt.Axes = axs[2]
        Plotter.add_plot_on_ax(self.rof,self.dis_diff,ax2,label='diff')
        Plotter.add_verline_on_ax(self.fit_packs[self.qubit]["optimal_rof"],self.dis_diff,ax2,label="optimal",colors='black',linestyles='--')
        Plotter.add_verline_on_ax(float(self.ds.attrs[f"{self.qubit}_ori_rof"]),self.dis_diff,ax2,label="original",colors='#DCDCDC',linestyles='--')
        
        Plotter.includes_axes([ax0,ax1,ax2])
        Plotter.set_LabelandSubtitles(
            [{'subtitle':"", 'xlabel':"ROF (Hz)", 'ylabel':"Magnitude (V)"},
                {'subtitle':"", 'xlabel':"ROF (Hz)", 'ylabel':"Phase (π)"},
                {'subtitle':"", 'xlabel':"ROF (Hz)", 'ylabel':"Diff (V)"}]
        )
        Plotter.export_results()

    def ROL_cali_ana(self,var:str):
        self.qubit = var
        IQ_data = moveaxis(array(self.ds[var]),1,0) # shape (mixer, state, rof) -> (state, mixer, rof)
        self.rol = moveaxis(array(self.ds[f"{var}_rol"]),0,1)[0][0]
        self.I_g, self.I_e = IQ_data[0][0], IQ_data[1][0]
        self.Q_g, self.Q_e = IQ_data[0][1], IQ_data[1][1]
        I_diff = self.I_e-self.I_g
        Q_diff = self.Q_e-self.Q_g
        self.dis_diff = sqrt((I_diff)**2+(Q_diff)**2)

    def ROL_cali_plot(self,save_pic_folder:str=None):
        if save_pic_folder is not None: save_pic_folder = os.path.join(save_pic_folder,f"{self.qubit}_ROL_cali_{self.ds.attrs['execution_time']}.png")
        Plotter = Artist(pic_title=f"{self.qubit}_ROL_calibration",save_pic_path=save_pic_folder)
        fig, axs = Plotter.build_up_plot_frame((2,1),fig_size=(12,9))
        ax0:plt.Axes = axs[0]
        Plotter.add_plot_on_ax(self.rol,self.dis_diff,ax0,label="|1>-|0>")
        ax1:plt.Axes = axs[1]
        Plotter.add_plot_on_ax(self.rol,arctan2(self.Q_e,self.I_e)/pi-arctan2(self.Q_g,self.I_g)/pi,ax1,label="|1>-|0>")
        # Plotter.add_verline_on_ax(self.fit_packs[q]["optimal_rol"],self.dis_diff,ax2,label="optimal",colors='black',linestyles='--')
        # Plotter.add_verline_on_ax(mean(self.rol),self.dis_diff,ax2,label="original",colors='#DCDCDC',linestyles='--')
        
        Plotter.includes_axes([ax0,ax1])
        Plotter.set_LabelandSubtitles(
            [{'subtitle':"", 'xlabel':"ROL-coef", 'ylabel':"Magnitude (V)"},
                {'subtitle':"", 'xlabel':"ROL-coef", 'ylabel':"Phase (π)"},]
        )
        Plotter.export_results()

    def piamp_cali_ana(self,var:str):
        self.qubit = var
        self.pi_amp_coef =  moveaxis(array(self.ds[f"{var}_PIcoef"]),1,0)[0][0]
        self.pi_pair_num = array(self.ds.coords["PiPairNum"])
        data = moveaxis(array(self.ds[var]),1,0)
        refined_data_folder = []
        for PiPairNum_dep_data in data:
            if len(self.refIQ) == 2:
                refined_data = IQ_data_dis(PiPairNum_dep_data[0],PiPairNum_dep_data[1],self.refIQ[0],self.refIQ[1])
            else:
                refined_data = rotate_data(PiPairNum_dep_data,self.refIQ[0])[0]
            refined_data_folder.append(cos_fit_analysis(refined_data,self.pi_amp_coef))
        
        # attrs["coefs"] = [A_fit,f_fit,phase_fit,offset_fit]
        
        candidates = []
        self.ans = []
        for idx, PiPairNum_dep_fitting_ds in enumerate(refined_data_folder[:2]):
            # Find minima in the range n = -3 to n = 3
            mini = list(find_minima(PiPairNum_dep_fitting_ds.attrs["coefs"][1], PiPairNum_dep_fitting_ds.attrs["coefs"][2], -3, 4))
            minimas = []
            
            for i in mini:
                if i <= max(self.pi_amp_coef) and i >= min(self.pi_amp_coef):
                    if idx == 0:
                        
                        minimas.append(i)
                    else:
                        if len(candidates) > 0:
                            for candit in candidates:
                                if abs(candit-i) <= (max(self.pi_amp_coef)-min(self.pi_amp_coef))/10:
                                    self.ans.append(i)
                                    candidates.remove(candit)
                        else:
                            break
            
            if idx == 0:
                candidates = minimas
        if len(self.ans) == 0:
            self.ans.append(1)
        
        self.fit_packs = {var:refined_data_folder,"ans":self.ans[0]}
    
    def piamp_cali_plot(self,save_pic_folder:str=None):
        q = list(self.fit_packs.keys())[0]
        if save_pic_folder is not None: save_pic_folder = os.path.join(save_pic_folder,f"{q}_PIamp_cali_{self.ds.attrs['execution_time']}.png")
        Plotter = Artist(pic_title=f"{q} PI-pulse amp coef calibration",save_pic_path=save_pic_folder)
        fig, axs = Plotter.build_up_plot_frame((1,1),fig_size=(12,9))
        ax:plt.Axes = axs[0]
        for idx, refined_data in enumerate(self.fit_packs[q]):
            x = refined_data.coords['freeDu']
            x_fit = refined_data.coords['para_fit']  
            Plotter.add_plot_on_ax(x,refined_data.data_vars['data'],ax,linestyle='--',label=f"{self.pi_pair_num[idx]} PI pairs", alpha=0.8, ms=4)
            Plotter.add_plot_on_ax(x_fit,refined_data.data_vars['fitting'],ax,linestyle='-', alpha=1, lw=2)    
        
        if len(self.ans) != 0:
            for an_ans in self.ans:
                ax = Plotter.add_verline_on_ax(an_ans, self.fit_packs[q][0].data_vars['data'], ax, colors='red',linestyles='--',label='Optimal')
            
        Plotter.includes_axes([ax])
        Plotter.set_LabelandSubtitles(
            [{'subtitle':"", 'xlabel':"PI-amp coef.", 'ylabel':"I signals (V)"}]
        )
        Plotter.export_results()
    
    def halfpiamp_cali_ana(self,var:str):
        self.qubit = var
        self.pi_amp_coef =  moveaxis(array(self.ds[f"{var}_HalfPIcoef"]),1,0)[0][0]
        self.pi_pair_num = array(self.ds.coords["PiPairNum"])
        data = moveaxis(array(self.ds[var]),1,0)
        refined_data_folder = []
        for PiPairNum_dep_data in data:
            if len(self.refIQ) == 2:
                refined_data = IQ_data_dis(PiPairNum_dep_data[0],PiPairNum_dep_data[1],self.refIQ[0],self.refIQ[1])
            else:
                refined_data = rotate_data(PiPairNum_dep_data,self.refIQ[0])[0]
            refined_data_folder.append(cos_fit_analysis(refined_data,self.pi_amp_coef))

        candidates = []
        self.ans = []
        for idx, PiPairNum_dep_fitting_ds in enumerate(refined_data_folder[:2]):
            # Find minima in the range n = -3 to n = 3
            mini = list(find_minima(PiPairNum_dep_fitting_ds.attrs["coefs"][1], PiPairNum_dep_fitting_ds.attrs["coefs"][2], -3, 4))
            minimas = []
            
            for i in mini:
                if i <= max(self.pi_amp_coef) and i >= min(self.pi_amp_coef):
                    if idx == 0:
                        
                        minimas.append(i)
                    else:
                        if len(candidates) > 0:
                            for candit in candidates:
                                if abs(candit-i) <= (max(self.pi_amp_coef)-min(self.pi_amp_coef))/10:
                                    self.ans.append(i)
                                    candidates.remove(candit)
                        else:
                            break
            
            if idx == 0:
                candidates = minimas
        
        if len(self.ans) == 0:
            self.ans.append(1)
        self.fit_packs = {var:refined_data_folder,"ans":self.ans[0]}
    
    def halfpiamp_cali_plot(self,save_pic_folder:str=None):
        q = list(self.fit_packs.keys())[0]
        if save_pic_folder is not None: save_pic_folder = os.path.join(save_pic_folder,f"{q}_halfPIamp_cali_{self.ds.attrs['execution_time']}.png")
        Plotter = Artist(pic_title=f"{q} half PI-pulse amp coef calibration",save_pic_path=save_pic_folder)
        fig, axs = Plotter.build_up_plot_frame((1,1),fig_size=(12,9))
        ax:plt.Axes = axs[0]
        for idx, refined_data in enumerate(self.fit_packs[q]):
            x = refined_data.coords['freeDu']
            x_fit = refined_data.coords['para_fit']  
            Plotter.add_plot_on_ax(x,refined_data.data_vars['data'],ax,linestyle='--',label=f"{self.pi_pair_num[idx]} PI quadruples", alpha=0.8, ms=4)
            Plotter.add_plot_on_ax(x_fit,refined_data.data_vars['fitting'],ax,linestyle='-', alpha=1, lw=2)    
        if len(self.ans) != 0:
            for an_ans in self.ans:
                ax = Plotter.add_verline_on_ax(an_ans, self.fit_packs[q][0].data_vars['data'], ax, colors='red',linestyles='--')
        Plotter.includes_axes([ax])
        Plotter.set_LabelandSubtitles(
            [{'subtitle':"", 'xlabel':"half PI-amp coef.", 'ylabel':"I signals (V)"}]
        )
        Plotter.export_results()

    def dragCali_ana(self,var:str):
        self.qubit = var
        self.drag_coef =  moveaxis(array(self.ds[f"{var}_dragcoef"]),1,0)[0][0]
        self.operations = array(self.ds.coords["operations"])
        data = moveaxis(array(self.ds[var]),1,0)
        self.plot_item = {}

        for idx, operation_dep_data in enumerate(data):
            self.plot_item[self.operations[idx]] = {}
            if len(self.refIQ) == 2:
                refined_data = IQ_data_dis(operation_dep_data[0],operation_dep_data[1],self.refIQ[0],self.refIQ[1])
            else:
                refined_data = rotate_data(operation_dep_data,self.refIQ[0])[0]
            
            self.plot_item[self.operations[idx]]['data'] = refined_data
            coefficients = polyfit(self.drag_coef, refined_data, deg=1)  
            self.plot_item[self.operations[idx]]["fit_para"] = coefficients  # [slope, intercept]

        if self.plot_item[self.operations[0]]["fit_para"][0] != self.plot_item[self.operations[1]]["fit_para"][1]:  # Ensure the lines are not parallel
            x_intersect = (self.plot_item[self.operations[1]]["fit_para"][1] - self.plot_item[self.operations[0]]["fit_para"][1]) / (self.plot_item[self.operations[0]]["fit_para"][0] - self.plot_item[self.operations[1]]["fit_para"][0])
            self.y_intersect = self.plot_item[self.operations[0]]["fit_para"][0] * x_intersect + self.plot_item[self.operations[0]]["fit_para"][1]
            self.fit_packs = {"optimal_drag_coef": x_intersect}
        else:
            print("The lines are parallel and do not intersect.")

    def dragCali_plot(self, save_pic_folder:str=None):
        def linear(x,a,b):
            return a*x+b

        if save_pic_folder is not None: save_pic_folder = os.path.join(save_pic_folder,f"{self.qubit}_DragCali_{self.ds.attrs['execution_time']}.png")
        Plotter = Artist(pic_title=f"{self.qubit} DRAG coef calibration",save_pic_path=save_pic_folder)
        fig, axs = Plotter.build_up_plot_frame((1,1),fig_size=(12,9))
        ax:plt.Axes = axs[0]
        for ope in self.plot_item:
            ax = Plotter.add_plot_on_ax(self.drag_coef,self.plot_item[ope]['data'],ax,linestyle='--',label=ope, alpha=0.8, ms=4)
            ax = Plotter.add_plot_on_ax(self.drag_coef,linear(self.drag_coef,*self.plot_item[ope]['fit_para']),ax,linestyle='-', alpha=1, lw=2)

        if len(list(self.fit_packs.keys())) != 0:
            ax = Plotter.add_scatter_on_ax(self.fit_packs["optimal_drag_coef"],self.y_intersect,ax,marker='*',c='red',s=200,label=f"Optimal={round(self.fit_packs['optimal_drag_coef'],2)}")
        Plotter.includes_axes([ax])
        Plotter.set_LabelandSubtitles(
            [{'subtitle':"", 'xlabel':"DRAG coef.", 'ylabel':"I signals (V)"}]
        )
        Plotter.export_results()

    def gateError_ana(self,var:str,tansition_freq_Hz:float=None):
        self.qubit = var
        datas = moveaxis(array(self.ds[var]),0,1)*1000 # shape (pulse_num, IQ, shots) 
        gate_num = array(self.ds.coords["pulse_num"])
        p0_data = datas[0]
        p1_data = datas[1] 
        
        gates, p_rec = [], []
        self.fq = tansition_freq_Hz  
        self.md = GMMROFidelity()
        self.train_set = DataArray(moveaxis(array([p0_data,p1_data]),0,1), coords= [("mixer",["I","Q"]), ("prepared_state",[0,1]), ("index",arange(p0_data.shape[-1]))] )
        self.md._import_data(self.train_set)
        self.md._start_analysis()
        self.g1d_fidelity = self.md.export_G1DROFidelity()
    
        lis = []
        for idx, data in enumerate(datas):
            if idx > 0: # idx = 0 is prepare ground
                lis.append(list(moveaxis(array([data]),0,1))) # gate num, mixer, prepared_state, index
                gates.append(gate_num[idx])
             
        # GMM1D (OKAY)
        # centers_2d, centers1d, sigmas = self.md.discriminator._export_1D_paras()
        # da = DataArray(moveaxis(moveaxis(array(lis),0,1),-1,1), coords= [("mixer",array(["I","Q"])), ("index",arange(p1_data.shape[-1])), ("gate_num",array(gates)) , ("prepared_state",array([0]))] )
        # train_data_proj = get_proj_distance(centers_2d.transpose(), da.transpose(*tuple(da.coords)).values)
        # dataset_proj = DataArray(train_data_proj,coords=[(i,array(da.coords[i])) for i in list(da.coords)[1:]])  
        # self.g1d_fidelity.discriminator._import_data(dataset_proj)
        # self.g1d_fidelity.discriminator._start_analysis()
        # ans = self.g1d_fidelity.discriminator.result.transpose() # (prepared_state, gate_num, index)
        
        # GMM
        da = DataArray(moveaxis(array(lis),0,1), coords= [("mixer",array(["I","Q"])), ("gate_num",array(gates)) , ("prepared_state",array([0])), ("index",arange(p1_data.shape[-1]))] )
        self.md.discriminator._import_data(da)
        self.md.discriminator._start_analysis()
        ans = self.md.discriminator._export_result()

        for dim_1_data in ans:
            for dim_2_data in dim_1_data:
                p = list(dim_2_data).count(1)/len(list(dim_2_data))
                p_rec.append(p)
        
        print(array(p_rec).shape)
        print(array(gates).shape)
        # Fit the data
        self.params = gate_phase_fit_analysis(array(p_rec),array(gates))
        self.fit_packs = {"gate_num":self.params.coords['freeDu'].values, "f":self.params.attrs['f']*1000, "tau":self.params.attrs['T2_fit']}

    def gateError_plot(self,save_pic_path:str=None):
        
        self.md._import_data(self.train_set)
        self.md._start_analysis()
        g1d_fidelity = self.md.export_G1DROFidelity()
        plot_readout_fidelity(self.train_set, self.md, g1d_fidelity, self.fq, save_pic_path if save_pic_path is not None else None, plot=True if save_pic_path is None else False)
        plt.close()

        data = self.params['data'].values
        gate_num = self.params.coords['freeDu'].values
        fitting = self.params['fitting'].values
        fit_x = self.params.coords['para_fit'].values

        
        Plotter = Artist(pic_title=f"{self.qubit} Gate Error Test",save_pic_path=save_pic_path+"_Error.png")
        fig, axs = Plotter.build_up_plot_frame((1,1))
        ax:plt.Axes = axs[0]
        ax = Plotter.add_scatter_on_ax(gate_num, data, ax)
        # Extract fitted parameters
        
        ax = Plotter.add_plot_on_ax(fit_x,fitting,ax,c='red',label=f"f={round(self.params.attrs['f']*1000,3)} mrad, tau={round(self.params.attrs['T2_fit'],2)}")
        Plotter.includes_axes([ax])
        Plotter.set_LabelandSubtitles(
            [{'subtitle':"", 'xlabel':"gate num", 'ylabel':"|1> populations"}]
        )
        Plotter.export_results()


################################
####   Analysis Interface   ####
################################


class Multiplex_analyzer(QCATAna,analysis_tools):
    def __init__(self, analyze_for_what_exp:str):
        QCATAna.__init__(self)
        analysis_tools.__init__(self)
        self.exp_name = analyze_for_what_exp.lower()
        self.fit_packs = {}


    def _import_data( self, data:Dataset|DataArray, var_dimension:int, refIQ:list=[], fit_func:callable=None, fq_Hz:float=None):
        self.ds = data
        self.dim = var_dimension
        self.refIQ = refIQ if len(refIQ) != 0 else [0,0]
        self.fit_func:callable = fit_func
        self.transition_freq = fq_Hz

    def _start_analysis(self,**kwargs):
        match self.exp_name:
            case 'm5':
                self.fluxCoupler_ana(kwargs["var_name"],self.refIQ)
            case 'm6':
                self.fluxCavity_ana(kwargs["var_name"])
            case 'm8': 
                self.conti2tone_ana(kwargs["var_name"],self.fit_func,self.refIQ)
            case 'm9': 
                self.fluxQb_ana(kwargs["var_name"],self.fit_func,self.refIQ)
            case 'm11': 
                self.rabi_ana(kwargs["var_name"])
            case 'm14': 
                self.oneshot_ana(self.ds,self.transition_freq)
            case 'm12':
                self.T2_ana(kwargs["var_name"],self.refIQ)
            case 'm13':
                self.T1_ana(kwargs["var_name"],self.refIQ)
            case 'c1':
                self.ROF_cali_ana(kwargs["var_name"])
            case 'c2':
                self.XYF_cali_ana(kwargs["var_name"],self.refIQ)
            case 'c3':
                self.piamp_cali_ana(kwargs["var_name"])
            case 'c4':
                self.halfpiamp_cali_ana(kwargs["var_name"])
            case 'c5':
                self.ROL_cali_ana(kwargs["var_name"])
            case 'c6':
                self.dragCali_ana(kwargs["var_name"])
            case 'auxa':
                self.ZgateT1_ana(**kwargs)
            case 't1':
                self.gateError_ana(kwargs["var_name"],self.transition_freq)
            case _:
                raise KeyError(f"Unknown measurement = {self.exp_name} was given !")

    def _export_result( self, pic_save_folder=None):
        match self.exp_name:
            case 'm5':
                self.fluxCoupler_plot(pic_save_folder)
            case 'm6':
                self.fluxCavity_plot(pic_save_folder)
            case 'm8':
                self.conti2tone_plot(pic_save_folder)
            case 'm9':
                self.fluxQb_plot(pic_save_folder)
            case 'm11': 
                self.rabi_plot(pic_save_folder)
            case 'm14': 
                self.oneshot_plot(pic_save_folder)
            case 'm12':
                self.T2_plot(pic_save_folder)
            case 'm13':
                self.T1_plot(pic_save_folder)
            case 'c1':
                self.ROF_cali_plot(pic_save_folder)
            case 'c2':
                self.XYF_cali_plot(pic_save_folder,round(self.transition_freq*1e-6) if self.transition_freq is not None else None)
            case 'c3':
                self.piamp_cali_plot(pic_save_folder)
            case 'c4':
                self.halfpiamp_cali_plot(pic_save_folder)
            case 'c5':
                self.ROL_cali_plot(pic_save_folder)
            case 'c6':
                self.dragCali_plot(pic_save_folder)
            case 'auxa':
                self.ZgateT1_plot(pic_save_folder)
            case 't1':
                self.gateError_plot(pic_save_folder)




