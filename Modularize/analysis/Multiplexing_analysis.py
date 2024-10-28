from numpy import array, mean, median, argmax, linspace, arange
from numpy import sqrt
from numpy import ndarray
from xarray import Dataset
from qcat.analysis.base import QCATAna
import matplotlib.pyplot as plt
import os, sys 
import quantify_core.data.handling as dh
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from Modularize.support.QDmanager import QDmanager, Data_manager
from Modularize.analysis.raw_data_demolisher import Conti2tone_dataReducer, fluxQub_dataReductor, Rabi_dataReducer
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis, Fit_analysis_plot
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays, remove_outlier_after_fit
from scipy.optimize import curve_fit
def parabola(x,a,b,c):
    return a*array(x)**2+b*array(x)+c


#? Analyzer use tools to analyze data
################################
####     Analysis tools    ####
################################


class analysis_tools():
    def __init__(self):
        self.ds:Dataset = None
    
    def import_dataset(self,ds:Dataset):
        self.ds = ds  # inheritance 

    def conti2tone_ana(self, fit_func:callable=None, refIQ:list=[]):
        self.fit_packs = {}
        dataset_processed = dh.to_gridded_dataset(self.ds)
        self.target_q = self.ds.attrs['target_q']
        self.xyf = array(dataset_processed['x0'])
        self.xyl = array(dataset_processed['x1'])
        ii = array(self.ds["y0"])
        qq = array(self.ds["y1"])

        self.contrast = sqrt((ii-refIQ[0])**2+(qq-refIQ[1])**2).reshape(self.xyl.shape[0],self.xyf.shape[0])
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
    
    def conti2tone_plot(self):
        if self.xyl.shape[0] != 1:
            plt.pcolormesh(self.xyf,self.xyl,self.contrast)
            if len(self.fit_f01s) != 0 :
                plt.scatter(self.fit_f01s,self.fif_amps,marker="*",c='red')
            plt.title(f"{self.target_q} power-dep 2tone")
            plt.ylabel("XY power (V)")
            
        else:
            plt.plot(self.xyf,self.contrast[0])
            plt.ylabel("Contrast")
            plt.title(f"{self.target_q} 2tone with XY power {self.xyl[0]} V")
        plt.xlabel("XY frequency (Hz)")
        plt.grid()
        plt.show()

    def fluxQb_ana(self, fit_func:callable=None, refIQ:list=[], filter_outlier:bool=True):
        qubit = self.ds.attrs['target_q']
        self.fit_packs = {}
        self.filtered_z, self.filtered_f = [], []
        if qubit in self.ds.attrs["ref_z"].split("_"):
            self.ref_z = float(self.ds.attrs["ref_z"].split("_")[int(self.ds.attrs["ref_z"].split("_").index(qubit))+1])
        else:
            self.ref_z = 0

        self.f,self.z,i,q = convert_netCDF_2_arrays(self.ds)
        self.contrast = array(sqrt((i-refIQ[0])**2+(q-refIQ[1])**2))
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

        if len(list(self.fit_packs.keys())) != 0:
            plt.title(f"{self.ds.attrs['target_q']} XYF={round(self.fit_packs['xyf']*1e-9,3)} GHz with z_pulse amp={round(float(-self.paras[1]/(2*self.paras[0])),3)} V")
        if save_pic_path is None:
            plt.show()
        else:
            plt.savefig(os.path.join(save_pic_path,f"{self.ds.attrs['target_q']}_fluxQubitSpectro_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png"))
            plt.close()

    def rabi_ana(self, refIQ:list=[]):
        x_data = array(self.ds['x0'])
        title = self.ds['x0'].attrs["long_name"]
        x_axis_name = self.ds['x0'].attrs["name"]
        x_axis_unit = self.ds['x0'].attrs["unit"]

        i_data = array(self.ds['y0'])
        q_data = array(self.ds['y1'])

        contrast = sqrt((i_data-refIQ[0])**2+(q_data-refIQ[1])**2)
        self.fit_packs = Rabi_fit_analysis(contrast,x_data,title)

    def rabi_plot(self):
        Fit_analysis_plot(self.fit_packs,P_rescale=None,Dis=None,q=self.ds.attrs['target_q'])




################################
####   Analysis Interface   ####
################################


class Multiplex_analyzer(QCATAna,analysis_tools):
    def __init__(self, analyze_for_what_exp:str):
        QCATAna.__init__(self)
        analysis_tools.__init__(self)
        self.exp_name = analyze_for_what_exp.lower()


    def _import_data( self, data:Dataset, var_dimension:int, refIQ:list=[], fit_func:callable=None):
        self.ds = data
        self.dim = var_dimension
        self.refIQ = refIQ if len(refIQ) != 0 else [0,0]
        self.fit_func:callable = fit_func

    def _start_analysis(self):
        match self.exp_name:
            case 'm88': 
                self.conti2tone_ana(self.fit_func,self.refIQ)
            case 'm99': 
                self.fluxQb_ana(self.fit_func,self.refIQ)
            case 'm1111': 
                self.rabi_ana(self.refIQ)
            case _:
                raise KeyError(f"Unknown measurement = {self.exp_name} was given !")

    def _export_result( self, plot:bool=False):
        if plot:
            match self.exp_name:
                case 'm88':
                    self.conti2tone_plot()
                case 'm99':
                    self.fluxQb_plot()
                case 'm1111': 
                    self.rabi_plot()
        return self.fit_packs



if __name__ == "__main__":

    """ Continuous wave 2 tone """
    # m88_file = "Modularize/Meas_raw/20241026/DR2multiQ_2tone_H13M01S43.nc"
    # QD_file = "Modularize/QD_backup/20241026/DR2#10_SumInfo.pkl"
    # QD_agent = QDmanager(QD_file)
    # QD_agent.QD_loader()

    # dss = Conti2tone_dataReducer(m88_file)  
    # ANA = Multiplex_analyzer("m88")      
    # for q in dss:
    #     ANA._import_data(dss[q],2,QD_agent.refIQ[q],QS_fit_analysis)
    #     ANA._start_analysis()
    #     ans = ANA._export_result(plot=True)


    """ Flux Qubit Spectroscopy """
    # m99_file = "Modularize/Meas_raw/20241027/DR2multiQ_Flux2tone_H18M50S29.nc"
    # QD_file = "Modularize/QD_backup/20241027/DR2#10_SumInfo.pkl"
    # QD_agent = QDmanager(QD_file)
    # QD_agent.QD_loader()

    # dss = Conti2tone_dataReducer(m99_file)  
    # ANA = Multiplex_analyzer("m99")      
    # for q in dss:
    #     ANA._import_data(dss[q],2,QD_agent.refIQ[q],QS_fit_analysis)
    #     ANA._start_analysis()
    #     ans = ANA._export_result(plot=True)

    """ Rabi Oscillation """
    m1111_file = "Modularize/Meas_raw/20241028/DR2multiQ_Rabi_H19M51S07.nc"
    QD_file = "Modularize/QD_backup/20241027/DR2#10_SumInfo.pkl"
    QD_agent = QDmanager(QD_file)
    QD_agent.QD_loader()

    dss = Rabi_dataReducer(m1111_file)  
    ANA = Multiplex_analyzer("m1111")      
    for q in dss:
        ANA._import_data(dss[q],2,QD_agent.refIQ[q])
        ANA._start_analysis()
        ans = ANA._export_result(plot=True)

    
