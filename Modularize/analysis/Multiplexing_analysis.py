from numpy import array, mean, median, argmax, linspace, arange, column_stack, moveaxis, empty_like
from numpy import sqrt
from numpy import ndarray
from xarray import Dataset, DataArray
from qcat.analysis.base import QCATAna
import matplotlib.pyplot as plt
import os, sys 
import quantify_core.data.handling as dh
from Modularize.support.UserFriend import *
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from Modularize.support.QDmanager import QDmanager, Data_manager
from Modularize.analysis.raw_data_demolisher import Conti2tone_dataReducer, fluxQub_dataReductor, Rabi_dataReducer, OneShot_dataReducer
from Modularize.support.Pulse_schedule_library import QS_fit_analysis, Rabi_fit_analysis, Fit_analysis_plot
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays, remove_outlier_after_fit
from scipy.optimize import curve_fit
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.analysis.state_discrimination import p01_to_Teff
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from Modularize.support import rotate_onto_Inphase, rotate_data


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
    
    def conti2tone_plot(self, save_pic_path:str=None):
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

        if save_pic_path is None:
            plt.show()
        else:
            pic_path = os.path.join(save_pic_path,f"{self.target_q}_Conti2tone_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png")
            slightly_print(f"pic saved located:\n{pic_path}")
            plt.savefig(pic_path)
            plt.close()


        plt.show()

    def fluxQb_ana(self, fit_func:callable=None, refIQ:list=[], filter_outlier:bool=True):
        self.qubit = self.ds.attrs['target_q']
        self.fit_packs = {}
        self.filtered_z, self.filtered_f = [], []
        if self.qubit in self.ds.attrs["ref_z"].split("_"):
            self.ref_z = float(self.ds.attrs["ref_z"].split("_")[int(self.ds.attrs["ref_z"].split("_").index(self.qubit))+1])
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
            plt.title(f"{self.qubit} XYF={round(self.fit_packs['xyf']*1e-9,3)} GHz with z_pulse amp={round(float(-self.paras[1]/(2*self.paras[0])),3)} V")
        if save_pic_path is None:
            plt.show()
        else:
            pic_path = os.path.join(save_pic_path,f"{self.qubit}_fluxQubitSpectro_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png")
            slightly_print(f"pic saved located:\n{pic_path}")
            plt.savefig(pic_path)
            plt.close()

    def rabi_ana(self, ref:list=[]):
        self.qubit = self.ds.attrs["target_q"]
        x_data = array(self.ds['x0'])
        title = self.ds['x0'].attrs["long_name"]

        i_data = array(self.ds['y0'])
        q_data = array(self.ds['y1'])

        if not self.ds.attrs["OS_mode"]:
            if len(ref) == 2:
                data_to_fit = sqrt((i_data-ref[0])**2+(q_data-ref[1])**2)
            else:
                iq_data = column_stack((i_data,q_data)).T
                data_to_fit = rotate_data(iq_data,float(ref[0]))[0]
                plt.plot(data_to_fit)
                plt.show()
        else:
            # wait GMM
            pass

        self.fit_packs = Rabi_fit_analysis(data_to_fit,x_data,title)
        

    def rabi_plot(self, save_pic_path:str=None):
        save_pic_path = os.path.join(save_pic_path,f"{self.qubit}_Rabi_{self.ds.attrs['execution_time'] if 'execution_time' in list(self.ds.attrs) else Data_manager().get_time_now()}.png") if save_pic_path is not None else ""
        if save_pic_path is not None: slightly_print(f"pic saved located:\n{save_pic_path}")
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


        _, self.RO_rotate_angle = rotate_onto_Inphase(self.gmm2d_fidelity.centers[0],self.gmm2d_fidelity.centers[1])
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


################################
####   Analysis Interface   ####
################################


class Multiplex_analyzer(QCATAna,analysis_tools):
    def __init__(self, analyze_for_what_exp:str):
        QCATAna.__init__(self)
        analysis_tools.__init__(self)
        self.exp_name = analyze_for_what_exp.lower()


    def _import_data( self, data:Dataset|DataArray, var_dimension:int, refIQ:list=[], fit_func:callable=None, fq_Hz:float=None):
        self.ds = data
        self.dim = var_dimension
        self.refIQ = refIQ if len(refIQ) != 0 else [0,0]
        self.fit_func:callable = fit_func
        self.transition_freq = fq_Hz

    def _start_analysis(self):
        match self.exp_name:
            case 'm88': 
                self.conti2tone_ana(self.fit_func,self.refIQ)
            case 'm99': 
                self.fluxQb_ana(self.fit_func,self.refIQ)
            case 'm1111': 
                self.rabi_ana(self.refIQ)
            case 'm1414': 
                self.oneshot_ana(self.ds,self.transition_freq)
            case _:
                raise KeyError(f"Unknown measurement = {self.exp_name} was given !")

    def _export_result( self, pic_save_folder=None):
        match self.exp_name:
            case 'm88':
                self.conti2tone_plot(pic_save_folder)
            case 'm99':
                self.fluxQb_plot(pic_save_folder)
            case 'm1111': 
                self.rabi_plot(pic_save_folder)
            case 'm1414': 
                self.oneshot_plot(pic_save_folder)



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
    #     ans = ANA._export_result()


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
    #     ans = ANA._export_result()

    """ Rabi Oscillation """
    # m1111_file = "Modularize/Meas_raw/20241029/DR2multiQ_Rabi_H15M46S09.nc" # power rabi
    # # m1111_file = "Modularize/Meas_raw/20241029/DR2multiQ_Rabi_H11M03S46.nc" # time  rabi
    # QD_file = "Modularize/QD_backup/20241029/DR2#10_SumInfo.pkl"
    # QD_agent = QDmanager(QD_file)
    # QD_agent.QD_loader()

    # dss = Rabi_dataReducer(m1111_file)  
    # ANA = Multiplex_analyzer("m1111")      
    # for q in dss:
    #     ANA._import_data(dss[q],1,QD_agent.refIQ[q])
    #     ANA._start_analysis()
    #     ans = ANA._export_result()

    """ Single shot """
    m1414_file = "Modularize/Meas_raw/20241029/DR2multiQ_SingleShot(0)_H20M03S35.nc"
    QD_file = "Modularize/QD_backup/20241029/DR2#10_SumInfo.pkl"
    QD_agent = QDmanager(QD_file)
    QD_agent.QD_loader()

    ds = OneShot_dataReducer(m1414_file)
    for var in ds.data_vars:
        ANA = Multiplex_analyzer("m1414")
        ANA._import_data(ds[var],var_dimension=0,fq_Hz=QD_agent.quantum_device.get_element(var).clock_freqs.f01())
        ANA._start_analysis()
        pic_path = None #os.path.join(Data_manager().get_today_picFolder(),f"{var}_SingleShot_{ds.attrs['execution_time']}")
        ANA._export_result(pic_path)
    
