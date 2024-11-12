import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import xarray as xr
import matplotlib.pyplot as plt 
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
from qblox_drive_AS.support.QDmanager import QDmanager
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support.Pulse_schedule_library import IQ_data_dis, T1_fit_analysis
from numpy import array, mean, median, std, average, round, max, min, transpose, abs, sqrt, cos, sin, pi, linspace, arange,ndarray, log10, ndarray, asarray
from datetime import datetime
from scipy.optimize import curve_fit 


def build_result_pic_path(dir_path:str,folder_name:str="")->str:
    parent = os.path.split(dir_path)[0]
    new_path = os.path.join(parent,"ZgateT1_pic" if folder_name == "" else folder_name)
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    return new_path

def zgate_T1_fitting(dataset:xr.Dataset, ref_IQ:list, T1_guess:float, fit:bool=True):
    
    time = dataset.coords["time"].values
    flux = dataset.coords["z_voltage"].values

    T1s = []
    signals = []
    
    for ro_name, data in dataset.data_vars.items():
        I = data.values[0]
        Q = data.values[1]
        
        for z_idx in range(array(data.values[0]).shape[0]):
            data = IQ_data_dis(I[z_idx],Q[z_idx],ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
            signals.append(data)
            if fit:
                try:
                    data_fit = T1_fit_analysis(data=data,freeDu=array(time),T1_guess=T1_guess)
                    # Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)
                    if data_fit.attrs['T1_fit'] < 10*T1_guess and data_fit.attrs['T1_fit'] > 0:
                        T1s.append(data_fit.attrs['T1_fit']*1e6)
                    else:
                        T1s.append(0)
                except:
                    T1s.append(0e-6)

    return time*1e6, flux, T1s, signals


class ZT1_Analyzer():
    def __init__(self,refIQ:list,bias_ref:float,T1_guess:float,excited_prepared:bool=True):
        self.refIQ = refIQ
        self.ref_bias = bias_ref
        self.prepare_excited = excited_prepared
        self.T1_guess = T1_guess

    def import_data(self,folder_path:str):
        self.folder_path = folder_path
        self.datasets = []
        # Iterate directory
        for path in os.listdir(dir_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(dir_path, path)):
                file_path = os.path.join(dir_path,path)
                if file_path.split(".")[-1] == 'nc':
                    self.datasets.append(xr.open_dataset(file_path))
    def start_analysis(self,time_sort:bool=False):
        """ 
        If time_sort: 
            self.T1_rec = {"%Y-%m-%d %H:%M:%S":{"T1s":[],"mean_T1":0, "median_T1":0, "std_T1":0}.
        else:
            
        """
        if not time_sort:
            self.T1_per_dataset = []
            self.I_chennel_per_dataset = []   
            for dataset in self.datasets :
                # for dataset in sub_set:
                self.qubit = list(dataset.data_vars)[0]
                self.time, bias, T1s, Isignal = zgate_T1_fitting(dataset,self.refIQ,self.T1_guess,fit=self.prepare_excited)
                if self.prepare_excited:
                    self.T1_per_dataset.append(T1s)
                self.I_chennel_per_dataset.append(Isignal)
            self.avg_I_data = average(array(self.I_chennel_per_dataset),axis=0)
            self.z = bias+self.ref_bias
            if self.prepare_excited:
                self.avg_T1 = average(array(self.T1_per_dataset),axis=0)
                self.std_T1_percent = round(std(array(self.T1_per_dataset),axis=0)*100/self.avg_T1,1)

        else:
            self.T1_rec = {}
            for dataset in self.datasets :
                # for dataset in sub_set:
                self.qubit = list(dataset.data_vars)[0]
                self.time, bias, T1s, Isignal = zgate_T1_fitting(dataset,self.refIQ,self.T1_guess)
                self.T1_rec[dataset.attrs["end_time"]] = {}
                self.T1_rec[dataset.attrs["end_time"]]["T1s"] = T1s
            self.T1_rec = dict(sorted(self.T1_rec.items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S")))
            self.z = bias+self.ref_bias

class Drawer():
    def __init__(self,pic_title:str,save_pic_path:str=None):
        self.title = pic_title
        self.pic_save_path = save_pic_path
        self.title_fontsize:int = 20

    def export_results(self):
        plt.title(self.title,fontsize=self.title_fontsize)
        plt.legend()
        plt.grid()
        plt.tight_layout()
        if self.pic_save_path is not None:
            plt.savefig(self.pic_save_path)
            plt.close()
        else:
            plt.show()
 
    
    def build_up_plot_frame(self,subplots_alignment:list=[3,1])->tuple[Figure,list]:
        fig, axs = plt.subplots(subplots_alignment[0],subplots_alignment[1],figsize=(subplots_alignment[0]*9,subplots_alignment[1]*6))
        if subplots_alignment[0] == 1 and subplots_alignment[1] == 1:
            axs = [axs]
        return fig, axs
    
    def add_colormesh_on_ax(self,x_data:ndarray,y_data:ndarray,z_data:ndarray,fig:Figure,ax:plt.Axes)->plt.Axes:
        im = ax.pcolormesh(x_data,y_data,transpose(z_data),shading="nearest")
        fig.colorbar(im, ax=ax)
        return ax

    def add_scatter_on_ax(self,x_data:ndarray,y_data:ndarray,ax:plt.Axes,**kwargs)->plt.Axes:
        ax.scatter(x_data,y_data,**kwargs)
        return ax
    
    def add_plot_on_ax(self,x_data:ndarray,y_data:ndarray,ax:plt.Axes,**kwargs)->plt.Axes:
        ax.plot(x_data,y_data,**kwargs)
        return ax
    
    def set_xaxis_number_size(self,axs:list,fontsize:int):
        for ax in axs:
            ax:plt.Axes
            ax.xaxis.set_tick_params(labelsize=fontsize)
    
    def set_yaxis_number_size(self,axs:list,fontsize:int):
        for ax in axs:
            ax:plt.Axes
            ax.yaxis.set_tick_params(labelsize=fontsize)



def zT1_poster(dir_path:str,refIQ:list,bias_ref:float,T1:float,prepare_excite:bool=True,time_trace_mode:bool=False,fq:ndarray=None):
    ZT1_ana = ZT1_Analyzer(refIQ=refIQ,bias_ref=bias_ref,T1_guess=T1,excited_prepared=prepare_excite)
    ZT1_ana.import_data(dir_path)
    ZT1_ana.start_analysis(time_trace_mode)
    Plotter = Drawer(pic_title=f"zgate-T1-{ZT1_ana.qubit}")
    x_label = f"{ZT1_ana.qubit} flux (V)"
    if fq is not None:
        ZT1_ana.z = fq
        x_label = f"{ZT1_ana.qubit} tran. Frequency (GHz)"
    if not time_trace_mode:
        fig, axs = Plotter.build_up_plot_frame([1,1])
        ax = Plotter.add_colormesh_on_ax(ZT1_ana.z,ZT1_ana.time,ZT1_ana.avg_I_data,fig,axs[0])
        if prepare_excite:
            for idx, T1_per_dataset in enumerate(ZT1_ana.T1_per_dataset):
                if idx == 0:
                    ax = Plotter.add_scatter_on_ax(x_data=ZT1_ana.z,y_data=T1_per_dataset,ax=ax,c='red',marker='*',s=50,label="raw data")
                else:
                    ax = Plotter.add_scatter_on_ax(x_data=ZT1_ana.z,y_data=T1_per_dataset,ax=ax,c='red',marker='*',s=50)
            Plotter.ax = Plotter.add_scatter_on_ax(x_data=ZT1_ana.z,y_data=ZT1_ana.avg_T1,ax=ax,c='pink',s=10,label="avg data")
        Plotter.ax.set_xlabel(x_label,fontsize=Plotter.title_fontsize)
        Plotter.ax.set_ylabel("Free evolution time (µs)",fontsize=Plotter.title_fontsize)
        Plotter.ax.set_ylim(0,max(ZT1_ana.time))
        

    else:
        
        fig, axs = Plotter.build_up_plot_frame([1,1])
        z = ZT1_ana.z
        ans_dict = ZT1_ana.T1_rec
        times = [datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S") for time_str in ans_dict.keys()]
        time_diffs = [(t - times[0]).total_seconds() / 60 for t in times]  # Convert to minutes

        # Get values
        T1_data = [entry["T1s"] for entry in ans_dict.values()]
        T1_data = array(T1_data).reshape(z.shape[0],len(times))
        Plotter.ax = Plotter.add_colormesh_on_ax(ZT1_ana.z,time_diffs,T1_data,fig,axs[0])
        Plotter.ax.set_xlabel(x_label,fontsize=Plotter.title_fontsize)
        Plotter.ax.set_ylabel("Time past (min)",fontsize=Plotter.title_fontsize)

    Plotter.pic_save_path = os.path.join(build_result_pic_path(ZT1_ana.folder_path),f"ZgateT1_{'zAVG' if  not time_trace_mode else 'TimeTrace'}_poster.png")
    slightly_print(f"pic saved located:\n{Plotter.pic_save_path}")
    Plotter.export_results()



if __name__ == "__main__":
    # fq, p, z, rate, std_T1_percent = plot_z_gateT1_poster(dir_path,z_fq_map["sweet"]["Z_v"],ref_IQ)
    # plot_background("Modularize/Meas_raw/z_gate_T1_test/z_gate_T1_pi_False/ToQM",ref_IQ,0)
    # plot_purcell_compa(fq, p, z, z_fq_map["sweet"]["Z_v"], rate, std_T1_percent, kappa)
    target_q = 'q2'
    t1_guess = 1e-6
    background_dir_path = "" # This folder contains all the ZgateT1 BACKGROUND nc files
    dir_path = "/Users/ratiswu/Downloads/ZgateT1_q2_H22M41S14/ToQM" # This folder contains all the ZgateT1 nc files
    QD_file = "Modularize/QD_backup/20241104/DRKE#242_SumInfo.pkl"
    QD_agent = QDmanager(QD_file)
    QD_agent.QD_loader()

    sweet_fq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()*1e-9
    sweet_bias=QD_agent.Fluxmanager.get_proper_zbiasFor(target_q)

    ZT1_ana = ZT1_Analyzer(refIQ=QD_agent.refIQ[target_q],bias_ref=sweet_bias,T1_guess=1e-6,excited_prepared=True)
    ZT1_ana.import_data(dir_path)
    ZT1_ana.start_analysis(False)

    
    def set_fit_paras(): # **** manually set
        Ec = 0.2 #GHz
        Ej_sum = 25
        init = (Ec,Ej_sum)
        up_b = (0.22,50)
        lo_b = (0.2,10)

        return init, lo_b, up_b

    def FqEqn(x,Ec,coefA):
        """
        a ~ period, b ~ offset, 
        """
        a = pi/QD_agent.Fluxmanager.get_PeriodFor(target_q)
        b = sweet_bias
        d = 0
        return sqrt(8*coefA*Ec*sqrt(cos(a*(x-b))**2+d**2*sin(a*(x-b))**2))-Ec

 
    init_, lo_b, up_b = set_fit_paras()
    p, e = curve_fit(FqEqn,array([sweet_bias]),array([sweet_fq]),p0=init_,bounds=(lo_b,up_b))
    fq = FqEqn(ZT1_ana.z,*p)

    zT1_poster(dir_path,QD_agent.refIQ[target_q],QD_agent.Fluxmanager.get_sweetBiasFor(target_q),t1_guess,time_trace_mode=True, fq=fq)
    


    
