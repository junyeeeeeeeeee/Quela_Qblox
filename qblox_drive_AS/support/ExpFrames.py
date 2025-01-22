from numpy import ndarray
from abc import ABC
import os
from datetime import datetime
from xarray import Dataset
from qblox_drive_AS.support.QDmanager import QDmanager, Data_manager
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer, sort_timeLabel
from qblox_drive_AS.support.UserFriend import *
from xarray import open_dataset
from numpy import array, linspace, arange, logspace, mean, median, std, sort
from abc import abstractmethod
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, advise_where_fq
from qblox_drive_AS.support.Pulse_schedule_library import set_LO_frequency, QS_fit_analysis
from quantify_scheduler.helpers.collections import find_port_clock_path
from qblox_drive_AS.analysis.raw_data_demolisher import ZgateT1_dataReducer
from qblox_drive_AS.analysis.TimeTraceAna import time_monitor_data_ana


class ExpGovernment(ABC):
    def __init__(self):
        self.QD_path:str = ""
    
    @abstractmethod
    def SetParameters(self,*args,**kwargs):
        pass

    @abstractmethod
    def PrepareHardware(self,*args,**kwargs):
        pass

    @abstractmethod
    def RunMeasurement(self,*args,**kwargs):
        pass

    @abstractmethod
    def RunAnalysis(self,*args,**kwargs):
        pass

    @abstractmethod
    def CloseMeasurement(self,*args,**kwargs):
        pass

    @abstractmethod
    def WorkFlow(self):
        pass



class BroadBand_CavitySearching(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, target_qs:list,freq_start:float, freq_end:float, freq_pts:int):
        self.counter:int = len(target_qs)
        self.target_qs = target_qs
        self.freq_start = freq_start
        self.freq_end = freq_end
        self.freq_pts = freq_pts
    
    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[len(self.target_qs)-self.counter], 'ro'))
        # Readout select
        qrmRF_slot_idx = int(find_port_clock_path(self.QD_agent.quantum_device.hardware_config(),"q:res",f"{self.target_qs[int(len(self.target_qs)-self.counter)]}.ro")[1].split("_")[-1][6:])
        self.readout_module = self.cluster.modules[qrmRF_slot_idx-1]
    
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.wideCS import wideCS
        dataset = wideCS(self.readout_module,self.freq_start,self.freq_end,self.freq_pts)
        if self.save_dir is not None:
            self.save_path = os.path.join(self.save_dir,f"BroadBandCS_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
            self.__raw_data_location = self.save_path + ".nc"
            dataset.to_netcdf(self.__raw_data_location)
            self.save_fig_path = self.save_path+".png"
        else:
            self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)
        self.counter -= 1


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.wideCS import plot_S21
        if new_QD_path is None:
            QD_file = self.QD_path
        else:
            QD_file = new_QD_path

        if new_file_path is None:
            file_path = self.__raw_data_location
            fig_path = self.save_fig_path
        else:
            file_path = new_file_path
            fig_path = os.path.join(os.path.split(new_file_path)[0],"S21.png")

        QD_savior = QDmanager(QD_file)
        QD_savior.QD_loader()

        ds = open_dataset(file_path)

        plot_S21(ds,fig_path)
        ds.close()
        QD_savior.QD_keeper()


    def WorkFlow(self):
        while self.counter > 0 :
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()

class Zoom_CavitySearching(ExpGovernment):
    """ Helps you get the **BARE** cavities. """
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_range:dict, freq_pts:int=100, avg_n:int=100, execution:bool=True):
        """ freq_range: {"q0":[freq_start, freq_end], ...}, sampling function use linspace """
        self.freq_range = {}
        for q in freq_range:
            self.freq_range[q] = linspace(freq_range[q][0], freq_range[q][1], freq_pts)

        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(self.freq_range.keys())


    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro'))
        
        
    
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.CavitySpec import QD_RO_init, CavitySearch
        QD_RO_init(self.QD_agent,self.freq_range)
        meas = CavitySearch()
        meas.ro_elements = self.freq_range
        meas.execution = self.execution
        meas.n_avg = self.avg_n
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.run()
        dataset = meas.dataset
        
        
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"zoomCS_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.CavitySpec import CS_ana
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)

            CS_ana(QD_savior,ds,fig_path)
            ds.close()
            QD_savior.QD_keeper()


    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()
        
        self.CloseMeasurement()
        
class PowerCavity(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_span_range:dict, roamp_range:list, roamp_sampling_func:str, freq_pts:int=100, avg_n:int=100, execution:bool=True):
        """ ### Args:
            * freq_span_range: {"q0":[freq_span_start, freq_span_end], ...}, sampling function use linspace\n
            * roamp_range: [amp_start, amp_end, pts]\n
            * roamp_sampling_func (str): 'linspace', 'arange', 'logspace'
        """
        self.freq_range = {}
        self.tempor_freq:list = [freq_span_range,freq_pts] # After QD loaded, use it to set self.freq_range

        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(freq_span_range.keys())
        if roamp_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(roamp_sampling_func)
        else:
            sampling_func:callable = linspace
        self.roamp_samples = sampling_func(*roamp_range)

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.PowCavSpec import PowerDep_spec
        from qblox_drive_AS.SOP.CavitySpec import QD_RO_init
        
        # set self.freq_range
        for q in self.tempor_freq[0]:
            rof = self.QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            self.freq_range[q] = linspace(rof+self.tempor_freq[0][q][0],rof+self.tempor_freq[0][q][1],self.tempor_freq[1])
        QD_RO_init(self.QD_agent,self.freq_range)
        dataset = PowerDep_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.roamp_samples,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"PowerCavity_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.PowCavSpec import plot_powerCavity_S21
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)

            plot_powerCavity_S21(ds,QD_savior,fig_path)
            ds.close()
        # QD_savior.QD_keeper()


    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()      

class Dressed_CavitySearching(ExpGovernment):
    """ Helps you get the **Dressed** cavities. """
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_range:dict, ro_amp:dict, freq_pts:int=100, avg_n:int=100, execution:bool=True):
        """ 
        ### Args:\n
        * freq_range: {"q0":[freq_start, freq_end], ...}, sampling function use linspace\n
        * ro_amp: {"q0":0.1, "q2":.... }
        """
        self.freq_range = {}
        for q in freq_range:
            self.freq_range[q] = linspace(freq_range[q][0], freq_range[q][1], freq_pts)
        self.ro_amp = ro_amp
        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(self.freq_range.keys())


    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro')) 
        
    
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.CavitySpec import QD_RO_init, CavitySearch
        QD_RO_init(self.QD_agent,self.freq_range)
        for q in self.ro_amp:
            self.QD_agent.quantum_device.get_element(q).measure.pulse_amp(self.ro_amp[q])
        
        meas = CavitySearch()
        meas.ro_elements = self.freq_range
        meas.execution = self.execution
        meas.n_avg = self.avg_n
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.run()
        dataset = meas.dataset
        
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"dressedCS_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.CavitySpec import CS_ana
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            for q in self.ro_amp:
                QD_savior.quantum_device.get_element(q).measure.pulse_amp(self.ro_amp[q])
            CS_ana(QD_savior,ds,fig_path,keep_bare=False)
            ds.close()
            QD_savior.QD_keeper()


    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 
        
class FluxCoupler(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_span_range:dict, bias_elements:list, flux_range:list, flux_sampling_func:str, freq_pts:int=100, avg_n:int=100, execution:bool=True):
        """ ### Args:
            * freq_span_range: {"q0":[freq_span_start, freq_span_end], ...}, sampling function use linspace\n
            * bias_elements (list): ["c0", "c1",... ]\n
            * flux_range: [amp_start, amp_end, pts]\n
            * flux_sampling_func (str): 'linspace', 'arange', 'logspace'
        """
        self.freq_range = {}
        self.tempor_freq:list = [freq_span_range,freq_pts] # After QD loaded, use it to set self.freq_range
        self.bias_targets = bias_elements
        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(freq_span_range.keys())
        if flux_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(flux_sampling_func)
        else:
            sampling_func:callable = linspace
        self.flux_samples = sampling_func(*flux_range)

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro'))

        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.CouplerFluxSpec import fluxCoupler_spec
        from qblox_drive_AS.SOP.CavitySpec import QD_RO_init
        # set self.freq_range
        for q in self.tempor_freq[0]:
            rof = self.QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            self.freq_range[q] = linspace(rof+self.tempor_freq[0][q][0],rof+self.tempor_freq[0][q][1],self.tempor_freq[1])
        QD_RO_init(self.QD_agent,self.freq_range)
        dataset = fluxCoupler_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.bias_targets,self.flux_samples,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"FluxCoupler_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            for var in ds.data_vars:
                ANA = Multiplex_analyzer("m5")
                if var.split("_")[-1] != 'freq':
                    ANA._import_data(ds,2)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)
            ds.close()



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()          

class FluxCavity(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_span_range:dict, flux_range:list, flux_sampling_func:str, freq_pts:int=100, avg_n:int=100, execution:bool=True):
        """ ### Args:
            * freq_span_range: {"q0":[freq_span_start, freq_span_end], ...}, sampling function use linspace\n
            * flux_range: [amp_start, amp_end, pts]\n
            * flux_sampling_func (str): 'linspace', 'arange', 'logspace'
        """
        self.freq_range = {}
        self.tempor_freq:list = [freq_span_range,freq_pts] # After QD loaded, use it to set self.freq_range
        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(freq_span_range.keys())
        if flux_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(flux_sampling_func)
        else:
            sampling_func:callable = linspace
        self.flux_samples = sampling_func(*flux_range)

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro'))
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.FluxCavSpec import FluxCav_spec
        from qblox_drive_AS.SOP.CavitySpec import QD_RO_init
        # set self.freq_range
        for q in self.tempor_freq[0]:
            rof = self.QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            self.freq_range[q] = linspace(rof+self.tempor_freq[0][q][0],rof+self.tempor_freq[0][q][1],self.tempor_freq[1])
            
        QD_RO_init(self.QD_agent,self.freq_range)
        dataset = FluxCav_spec(self.QD_agent,self.meas_ctrl,self.Fctrl,self.freq_range,self.flux_samples,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"FluxCavity_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.FluxCavSpec import update_flux_info_in_results_for
        if self.execution:
            
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'freq':
                    ANA = Multiplex_analyzer("m6")
                    ANA._import_data(ds,2)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)
                    answer[var] = ANA.fit_packs
            ds.close()
            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                update_flux_info_in_results_for(QD_savior,permi,answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in list(answer.keys()):
                    update_flux_info_in_results_for(QD_savior,q,answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   

class IQ_references(ExpGovernment):
    """ Helps you get the **Dressed** cavities. """
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, ro_amp_factor:dict, shots:int=100, execution:bool=True):
        """ 
        ### Args:\n
        * ro_amp_factor: {"q0":1.2, "q2":.... }, new ro amp = ro_amp*ro_amp_factor
        """
        self.ask_save:bool = False
        self.ro_amp = ro_amp_factor
        self.avg_n = shots
        self.execution = execution
        self.target_qs = list(self.ro_amp.keys())


    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro')) 
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        for q in self.ro_amp:
            self.Fctrl[q](float(self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q)))
    
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.RefIQ import Single_shot_ref_spec
       

        dataset = Single_shot_ref_spec(self.QD_agent,self.ro_amp,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"IQref_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RefIQ import IQ_ref_ana
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            
            answer = {}
            for q in ds.data_vars:
                answer[q] = IQ_ref_ana(ds,q,fig_path)
            ds.close()
            
            QD_savior.memo_refIQ(answer)
            QD_savior.QD_keeper()


    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

class PowerConti2tone(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_range:dict, xyl_range:list, xyl_sampling_func:str, freq_pts:int=100, avg_n:int=100, ro_xy_overlap:bool=False, execution:bool=True):
        """ ### Args:
            * freq_range: {"q0":[freq_start, freq_end], ...}, sampling function use linspace\n
                * if someone is 0 like {"q0":[0]}, system will calculate an advised value.
            * xyl_range: [amp_start, amp_end, pts], if only one value inside, we only use that value. \n
            * xyl_sampling_func (str): 'linspace', 'arange', 'logspace'
        """
        self.freq_range = {}
        self.overlap:bool = ro_xy_overlap
        self.f_pts = freq_pts
        for q in freq_range:
            if len(freq_range[q]) == 1 and freq_range[q][0] == 0:
                self.freq_range[q] = freq_range[q][0]
            else:
                self.freq_range[q] = linspace(freq_range[q][0],freq_range[q][1],freq_pts)
        
        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(freq_range.keys())
        if xyl_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(xyl_sampling_func)
        else:
            sampling_func:callable = linspace
        if len(xyl_range) != 1:
            self.xyl_samples = list(sampling_func(*xyl_range))
        else:
            self.xyl_samples = list(xyl_range)

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro'), xy_out_att=0)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # set driving LO and offset bias
        for q in self.freq_range:
            self.Fctrl[q](float(self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q)))
            if isinstance(self.freq_range[q],ndarray):
                print(f"{q} LO @ {max(self.freq_range[q])}")
                set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=max(self.freq_range[q]))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.Cnti2Tone import Two_tone_spec
        # set self.freq_range
        for q in self.freq_range:
            if not isinstance(self.freq_range[q],ndarray):
                advised_fq = advise_where_fq(self.QD_agent,q,self.QD_agent.Notewriter.get_sweetGFor(q)) 
                eyeson_print(f"fq advice for {q} @ {round(advised_fq*1e-9,4)} GHz")
                IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
                if advised_fq-IF_minus < 2e9:
                    raise ValueError(f"Attempting to set {q} driving LO @ {round((advised_fq-IF_minus)*1e-9,1)} GHz")
                set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=advised_fq-IF_minus)
                self.freq_range[q] = linspace(advised_fq-IF_minus-500e6,advised_fq-IF_minus,self.f_pts)

        dataset = Two_tone_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.xyl_samples,self.avg_n,self.execution,self.overlap)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"PowerCnti2tone_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)

    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.Cnti2Tone import update_2toneResults_for
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'freq':
                    ANA = Multiplex_analyzer("m8")     
                    ANA._import_data(ds,2,QD_savior.refIQ[var],QS_fit_analysis)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)
                    if ANA.fit_packs != {}:
                        analysis_result = QS_fit_analysis(ANA.fit_packs[var]["contrast"],f=ANA.fit_packs[var]["xyf_data"])
                        update_2toneResults_for(QD_savior,var,{str(var):analysis_result},ANA.xyl[0])
            ds.close()
            QD_savior.QD_keeper()

    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

class FluxQubit(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_span_range:dict, bias_targets:list,z_amp_range:list, z_amp_sampling_func:str, freq_pts:int=100, avg_n:int=100, execution:bool=True):
        """ ### Args:
            * freq_span_range: {"q0":[freq_span_start, freq_span_end], ...}, sampling function use linspace\n
            * bias_targets: list, what qubit need to be bias, like ['q0', 'q1', ...]\n
            * z_amp_range: [amp_start, amp_end, pts]\n
            * z_amp_sampling_func (str): 'linspace', 'arange', 'logspace'
        """
        self.freq_range = {}
        self.bias_elements = bias_targets
        self.tempor_freq:list = [freq_span_range,freq_pts] # After QD loaded, use it to set self.freq_range
        self.avg_n = avg_n
        self.execution = execution
        self.target_qs = list(freq_span_range.keys())
        if z_amp_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(z_amp_sampling_func)
        else:
            sampling_func:callable = linspace
        self.z_amp_samples = sampling_func(*z_amp_range)

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # Set the system attenuations
        init_system_atte(self.QD_agent.quantum_device,self.target_qs,ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(self.target_qs[0], 'ro'),xy_out_att=0)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias
        self.z_ref = {}
        for q in self.target_qs:
            self.z_ref[q] = float(self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            self.Fctrl[q](self.z_ref[q])
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.FluxQubit import Zgate_two_tone_spec
        
        # set self.freq_range
        for q in self.target_qs:
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            if abs(self.tempor_freq[0][q][0]-self.tempor_freq[0][q][1]) >500e6:
                raise ValueError(f"Attempting to span over 500 MHz for driving on {q}")
            self.freq_range[q] = linspace(xyf+self.tempor_freq[0][q][0],xyf+self.tempor_freq[0][q][1],self.tempor_freq[1])
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=max(self.freq_range[q]))
            
        dataset = Zgate_two_tone_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.bias_elements,self.z_amp_samples,self.avg_n,self.execution)
        if self.execution:
            for q in self.z_ref:
                dataset.attrs[f"{q}_z_ref"] = self.z_ref[q]
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"FluxQubit_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.FluxQubit import update_by_fluxQubit
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'freq':
                    ANA = Multiplex_analyzer("m9") 
                    ANA._import_data(ds,2,QD_savior.refIQ[var],QS_fit_analysis)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)
                    if len(list(ANA.fit_packs.keys())) != 0: answer[var] = ANA.fit_packs        
            ds.close()
            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no :").lower()
            if permi in list(answer.keys()):
                update_by_fluxQubit(QD_savior,answer[q],q)
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in list(answer.keys()):
                    update_by_fluxQubit(QD_savior,answer[q],q)
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   

class PowerRabiOsci(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, pi_amp:dict, pi_amp_sampling_func:str, pi_amp_pts_or_step:float=100, avg_n:int=100, execution:bool=True, OSmode:bool=False):
        """ ### Args:
            * pi_amp: {"q0":[pi_amp_start, pi_amp_end], ...}\n
            * pi_amp_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * pi_amp_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        if pi_amp_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(pi_amp_sampling_func)
        else:
            sampling_func:callable = linspace
        
        self.pi_amp_samples = {}
        for q in pi_amp:
            self.pi_amp_samples[q] = sampling_func(*pi_amp[q],pi_amp_pts_or_step)
            
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(pi_amp.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and driving atte
        self.pi_dura = {}
        
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
            self.pi_dura[q] = self.QD_agent.quantum_device.get_element(q).rxy.duration()
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.RabiOsci import PowerRabi
    
        dataset = PowerRabi(self.QD_agent,self.meas_ctrl,self.pi_amp_samples,self.pi_dura,self.avg_n,self.execution,self.OSmode)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"PowerRabi_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RabiOsci import conditional_update_qubitInfo
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'piamp':
                    ANA = Multiplex_analyzer("m11")      
                    ANA._import_data(ds,1,QD_savior.refIQ[var] if QD_savior.rotate_angle[var][0] == 0 else QD_savior.rotate_angle[var])
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)
                    conditional_update_qubitInfo(QD_savior,ANA.fit_packs,var)  

            ds.close()
            QD_savior.QD_keeper()
            



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   

class TimeRabiOsci(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, pi_dura:dict, pi_dura_sampling_func:str, pi_dura_pts_or_step:float=100, avg_n:int=100, execution:bool=True, OSmode:bool=False):
        """ ### Args:
            * pi_amp: {"q0": pi-amp in V, ...}\n
            * pi_dura: {"q0":[pi_duration_start, pi_duration_end], ...}\n
            * pi_dura_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * pi_dura_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if pi_dura_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(pi_dura_sampling_func)
        else:
            sampling_func:callable = linspace
        
        self.pi_dura_samples = {}
        for q in pi_dura:
            self.pi_dura_samples[q] = sort_elements_2_multiples_of(sampling_func(*pi_dura[q],pi_dura_pts_or_step)*1e9,4)*1e-9
    
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(pi_dura.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        self.pi_amp = {}
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
            self.pi_amp[q] = self.QD_agent.quantum_device.get_element(q).rxy.amp180()
            
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.RabiOsci import TimeRabi
    
        dataset = TimeRabi(self.QD_agent,self.meas_ctrl,self.pi_amp,self.pi_dura_samples,self.avg_n,self.execution,self.OSmode)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"TimeRabi_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RabiOsci import conditional_update_qubitInfo
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'pidura':
                    ANA = Multiplex_analyzer("m11")      
                    ANA._import_data(ds,1,QD_savior.refIQ[var] if QD_savior.rotate_angle[var][0] == 0 else QD_savior.rotate_angle[var])
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)
                    conditional_update_qubitInfo(QD_savior,ANA.fit_packs,var)  
                    
            ds.close()
            # QD_savior.QD_keeper(new_QD_dir)
            



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   

class SingleShot(ExpGovernment):
    """ Helps you get the **Dressed** cavities. """
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location
    
    @RawDataPath.setter
    def RawDataPath(self,path:str):
        self.__raw_data_location = path

    def SetParameters(self, target_qs:list, histo_counts:int=1, shots:int=10000, execution:bool=True):
        """ 
        ### Args:\n
        * target_qs: list, like ["q0", "q1", ...]
        """
        self.use_time_label:bool = False
        self.avg_n = shots
        self.execution = execution
        self.target_qs = target_qs
        self.histos = histo_counts
        self.counter:int = 0
        if self.histos > 0:
            self.use_time_label = True



    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and driving atte
        for q in self.target_qs:
            # self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.SingleShot import Qubit_state_single_shot
       

        dataset = Qubit_state_single_shot(self.QD_agent,self.target_qs,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"SingleShot_{datetime.now().strftime('%Y%m%d%H%M%S') if (self.JOBID is None or self.use_time_label) else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None,histo_ana:bool=False):
        """ if histo_ana, it will check all the data in the same folder with the given new_file_path """
    
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            if not histo_ana:
                ds = open_dataset(file_path)
                for var in ds.data_vars:
                    ANA = Multiplex_analyzer("m14")
                    # QD_savior.StateDiscriminator = ANA.gmm2d_fidelity # will be in the future
                    ANA._import_data(ds[var]*1000,var_dimension=0,fq_Hz=QD_savior.quantum_device.get_element(var).clock_freqs.f01())
                    ANA._start_analysis()
                    pic_path = os.path.join(fig_path,f"{var}_SingleShot_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                    ANA._export_result(pic_path)
                    highlight_print(f"{var} rotate angle = {round(ANA.fit_packs['RO_rotation_angle'],2)} in degree.")
                    QD_savior.rotate_angle[var] = [ANA.fit_packs["RO_rotation_angle"]]
                ds.close()
                
                QD_savior.QD_keeper()
            else:

                eff_T, thermal_pop = {}, {}
                files = sort_timeLabel([os.path.join(fig_path,name) for name in os.listdir(fig_path) if (os.path.isfile(os.path.join(fig_path,name)) and name.split(".")[-1]=='nc')])
                for nc_idx, nc_file in enumerate(files):
                    ds = open_dataset(nc_file)
                    for var in ds.data_vars:
                        if nc_idx == 0: eff_T[var], thermal_pop[var] = [], []
                        ANA = Multiplex_analyzer("m14")
                        ANA._import_data(ds[var]*1000,var_dimension=0,fq_Hz=QD_savior.quantum_device.get_element(var).clock_freqs.f01())
                        ANA._start_analysis()
                        eff_T[var].append(ANA.fit_packs["effT_mK"])
                        thermal_pop[var].append(ANA.fit_packs["thermal_population"]*100)
                
                for qubit in eff_T:
                    highlight_print(f"{qubit}: {round(median(array(eff_T[qubit])),1)} +/- {round(std(array(eff_T[qubit])),1)} mK")
                    Data_manager().save_histo_pic(QD_savior,eff_T,qubit,mode="ss",pic_folder=fig_path)
                    Data_manager().save_histo_pic(QD_savior,thermal_pop,qubit,mode="pop",pic_folder=fig_path)

    def WorkFlow(self):
        for i in range(self.histos):
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement() 

            self.counter += 1

class Ramsey(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID
        self.histos:int = 0

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, time_range:dict, time_sampling_func:str, time_pts_or_step:int|float=100,histo_counts:int=1, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * time_range: {"q0":[time_start, time_end], ...}\n
            * histo_counts: int, larger than 100 use while loop.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(*time_range[q],time_pts_or_step)*1e9,4)*1e-9
        else:
            for q in time_range:
                self.time_samples[q] = sampling_func(*time_range[q],time_pts_or_step)

        self.avg_n = avg_n

        if histo_counts <= 100:
            self.want_while = False
            self.histos = histo_counts
        else:
            self.want_while = True
            self.histos = 1
        
        self.execution = execution
        self.OSmode = OSmode
        self.spin_num = {}
        self.target_qs = list(time_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            slightly_print(f"{q} arti-detune = {round(self.QD_agent.Notewriter.get_artiT2DetuneFor(q)*1e-6,2)} MHz")
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()+self.QD_agent.Notewriter.get_artiT2DetuneFor(q)
            self.QD_agent.quantum_device.get_element(q).clock_freqs.f01(xyf)
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.T2 import Ramsey
    
        dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"Ramsey_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
        
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    ANA = Multiplex_analyzer("m12")
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA._import_data(ds,var_dimension=2,refIQ= ref)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_T2_for(ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            QD_savior.QD_keeper()
            

    def WorkFlow(self,freq_detune_Hz:float=None):
        while True:
            self.PrepareHardware()

            if freq_detune_Hz is not None:
                for q in self.target_qs:
                    self.QD_agent.quantum_device.get_element(q).clock_freqs.f01(self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()+freq_detune_Hz)

            self.RunMeasurement()

            self.CloseMeasurement()   
            if not self.want_while:
                break

class SpinEcho(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, time_range:dict, time_sampling_func:str, time_pts_or_step:int|float=100,histo_counts:int=1, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * time_range: {"q0":[time_start, time_end], ...}\n
            * histo_counts: int, larger than 100 use while loop.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        self.spin_num = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(*time_range[q],time_pts_or_step)*1e9,8)*1e-9
        else:
            for q in time_range:
                self.time_samples[q] = sampling_func(*time_range[q],time_pts_or_step)

        self.avg_n = avg_n

        if histo_counts <= 100:
            self.want_while = False
            self.histos = histo_counts
        else:
            self.want_while = True
            self.histos = 1
        
        self.execution = execution
        self.OSmode = OSmode
        
        self.target_qs = list(time_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.spin_num[q] = 1
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.T2 import Ramsey
    
        dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"SpinEcho_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
        
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    ANA = Multiplex_analyzer("m12")
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA._import_data(ds,var_dimension=2,refIQ=ref)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_echoT2_for(ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            QD_savior.QD_keeper()
            

    def WorkFlow(self):
        while True:
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()   
            if not self.want_while:
                break

class CPMG(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, time_range:dict, pi_num:dict, time_sampling_func:str, time_pts_or_step:int|float=100,histo_counts:int=1, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * time_range: {"q0":[time_start, time_end], ...}\n
            * pi_num: {"q0":1, "q1":2, ...}
            * histo_counts: int, larger than 100 use while loop.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        self.spin_num = pi_num
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(*time_range[q],time_pts_or_step)*1e9,(1+int(pi_num[q]))*4)*1e-9
        else:
            for q in time_range:
                self.time_samples[q] = sampling_func(*time_range[q],time_pts_or_step)
        
        self.avg_n = avg_n

        if histo_counts <= 100:
            self.want_while = False
            self.histos = histo_counts
        else:
            self.want_while = True
            self.histos = 1
        
        self.execution = execution
        self.OSmode = OSmode
        
        self.target_qs = list(time_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.T2 import Ramsey
    
        dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"CPMG_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
        
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    ANA = Multiplex_analyzer("m12")
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA._import_data(ds,var_dimension=2,refIQ=ref)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_echoT2_for(ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            QD_savior.QD_keeper()
            

    def WorkFlow(self, freq_detune_Hz:float=None):
        while True:
            self.PrepareHardware()
            
            if freq_detune_Hz is not None:
                for q in self.target_qs:
                    self.QD_agent.quantum_device.get_element(q).clock_freqs.f01(self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()+freq_detune_Hz)

            self.RunMeasurement()

            self.CloseMeasurement()   
            if not self.want_while:
                break

class EnergyRelaxation(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, time_range:dict, time_sampling_func:str, time_pts_or_step:int|float=100,histo_counts:int=1, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * time_range: {"q0":[time_start, time_end], ...}\n
            * histo_counts: int, larger than 100 use while loop.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(*time_range[q],time_pts_or_step)*1e9,4)*1e-9
        else:
            for q in time_range:
                self.time_samples[q] = sampling_func(*time_range[q],time_pts_or_step)

        self.avg_n = avg_n

        if histo_counts <= 100:
            self.want_while = False
            self.histos = histo_counts
        else:
            self.want_while = True
            self.histos = 1
        
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(time_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.T1 import T1, EnergyRelaxPS
        meas = EnergyRelaxPS()
        meas.set_time_samples = self.time_samples
        meas.set_os_mode = self.OSmode
        meas.set_n_avg = self.avg_n
        meas.set_repeat = self.histos
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        
        meas.run()
        dataset = meas.dataset


        # dataset = T1(self.QD_agent,self.meas_ctrl,self.time_samples,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"T1_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
        
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    ANA = Multiplex_analyzer("m13")
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA._import_data(ds,var_dimension=2,refIQ=ref)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_T1_for(ANA.fit_packs["median_T1"],var)

            ds.close()
            QD_savior.QD_keeper()
            

    def WorkFlow(self):
        while True:
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()   
            if not self.want_while:
                break

class XYFcali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, target_qs:list, evo_time:float=0.5e-6, detu:float=0, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * target_qs: ["q0", "q1", ...]\n
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
    
        self.time_samples = {}
        for q in target_qs:
            self.time_samples[q] = sort_elements_2_multiples_of(linspace(40e-9,evo_time,100)*1e9,4)*1e-9
        self.detu = detu
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.spin_num = {}
        self.target_qs = target_qs
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()+self.detu
            self.QD_agent.quantum_device.get_element(q).clock_freqs.f01(xyf)
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.T2 import Ramsey
    
        dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,1,self.avg_n,self.execution,second_phase='y')
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"XYFcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    ANA = Multiplex_analyzer("c2")
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    
                    ANA._import_data(ds,var_dimension=2,refIQ=ref)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)

                    if ANA.fit_packs['phase'] > 180:
                        sign = -1
                    else:
                        sign = 1
                    
                    answer[var] = sign*ANA.fit_packs['freq']
                    highlight_print(f"{var}: actual detune = {round(answer[var]*1e-6,4)} MHz")
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).clock_freqs.f01(QD_savior.quantum_device.get_element(permi).clock_freqs.f01()-answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).clock_freqs.f01(QD_savior.quantum_device.get_element(q).clock_freqs.f01()-answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   

class ROFcali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_span_range:dict, freq_pts:int=100, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * freq_span_range: {"q0":[freq_span_start, freq_span_end], ....}\n
        """
        self.freq_samples = {}
        self.tempor_freq = [freq_span_range, freq_pts]

        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(freq_span_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
            rof = self.QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            self.freq_samples[q] = linspace(rof+self.tempor_freq[0][q][0],rof+self.tempor_freq[0][q][1],self.tempor_freq[1])
        
    def RunMeasurement(self):
        from qblox_drive_AS.Calibration_exp.RofCali import rofCali
    
        dataset = rofCali(self.QD_agent,self.meas_ctrl,self.freq_samples,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"ROFcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'rof':
                    ANA = Multiplex_analyzer("c1")
                    ANA._import_data(ds,var_dimension=1)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(fig_path)
                    answer[var] = ANA.fit_packs[var]["optimal_rof"]
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).clock_freqs.readout(answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).clock_freqs.readout(answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")


    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

class PiAcali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, piamp_coef_range:dict, amp_sampling_funct:str, coef_ptsORstep:int=100, pi_pair_num:list=[2,3], avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * piamp_coef_range: {"q0":[0.9, 1.1], "q1":[0.95, 1.15], ...]\n
            * amp_sampling_funct: str, `linspace` or `arange`.\n
            * pi_pair_num: list, like [2, 3] will try 2 exp, the first uses 2\*2 pi-pulse, and the second exp uses 3*2 pi-pulse
        """
        if amp_sampling_funct in ['linspace','logspace','arange']:
            sampling_func:callable = eval(amp_sampling_funct)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {amp_sampling_funct}")
        
        self.amp_coef_samples = {}
        for q in piamp_coef_range:
           self.amp_coef_samples[q] = sampling_func(*piamp_coef_range[q],coef_ptsORstep)
        self.pi_pair_num = pi_pair_num
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(piamp_coef_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))

        
    def RunMeasurement(self):
        from qblox_drive_AS.Calibration_exp.PI_ampCali import pi_amp_cali
    
        dataset = pi_amp_cali(self.QD_agent,self.meas_ctrl,self.amp_coef_samples,self.pi_pair_num,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"PIampcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'PIcoef':
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA = Multiplex_analyzer("c3")
                    ANA._import_data(ds,var_dimension=1,refIQ=ref)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(fig_path)
                    answer[var] = ANA.fit_packs["ans"]
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).rxy.amp180(QD_savior.quantum_device.get_element(permi).rxy.amp180()*answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).rxy.amp180(QD_savior.quantum_device.get_element(q).rxy.amp180()*answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

class hPiAcali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, piamp_coef_range:dict, amp_sampling_funct:str, coef_ptsORstep:int=100, halfPi_pair_num:list=[3,5], avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * piamp_coef_range: {"q0":[0.9, 1.1], "q1":[...], ...], this is the coef tiles to QDmanager.Waveformer.__xylog[q]["halfPI_ratio"].\n
            * amp_sampling_funct: str, `linspace` or `arange`.\n
            * pi_pair_num: list, like [3, 5] will try 2 exp, the first uses 3\*4 half pi-pulse, and the second exp uses 5*4 half pi-pulse
        """
        if amp_sampling_funct in ['linspace','logspace','arange']:
            sampling_func:callable = eval(amp_sampling_funct)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {amp_sampling_funct}")
        
        self.amp_coef_samples = {}
        for q in piamp_coef_range:
           self.amp_coef_samples[q] = sampling_func(*piamp_coef_range[q],coef_ptsORstep)
        self.halfPi_pair_num = halfPi_pair_num
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(piamp_coef_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))

        
    def RunMeasurement(self):
        from qblox_drive_AS.Calibration_exp.halfPI_ampCali import half_pi_amp_cali
    
        dataset = half_pi_amp_cali(self.QD_agent,self.meas_ctrl,self.amp_coef_samples,self.halfPi_pair_num,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"halfPIampcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'HalfPIcoef':
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA = Multiplex_analyzer("c4")
                    ANA._import_data(ds,var_dimension=1,refIQ=ref)
                    ANA._start_analysis(var_name = var) 
                    ANA._export_result(fig_path)
                    answer[var] = ANA.fit_packs["ans"]

            ds.close()
            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.Waveformer.set_halfPIratio_for(permi, QD_savior.Waveformer.get_halfPIratio_for(permi)*answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.Waveformer.set_halfPIratio_for(q, QD_savior.Waveformer.get_halfPIratio_for(q)*answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")
    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

class ROLcali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, roamp_coef_range:dict, coef_sampling_func:str, ro_coef_ptsORstep:int=100, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * roamp_coef_range: {"q0":[0.85, 1.3], ... }, rule:  q_name:[roamp_coef_start, roamp_coef_end]. exp with ro-amp *=  roamp_coef\n
            * coef_sampling_func: str, `'linspace'` or `'arange'`.
        """
        if coef_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(coef_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {coef_sampling_func}")
        
        self.amp_coef_samples = {}
        for q in roamp_coef_range:
           self.amp_coef_samples[q] = sampling_func(*roamp_coef_range[q],ro_coef_ptsORstep)

        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(roamp_coef_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))

        
    def RunMeasurement(self):
        from qblox_drive_AS.Calibration_exp.RO_ampCali import rolCali
    
        dataset = rolCali(self.QD_agent,self.meas_ctrl,self.amp_coef_samples,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"ROLcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            
            for var in ds.data_vars:
                if var.split("_")[-1] != 'rol':
                    ANA = Multiplex_analyzer("c5")
                    ANA._import_data(ds,var_dimension=1)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(fig_path)
                    
            ds.close()


    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

class ZgateEnergyRelaxation(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, time_range:dict, time_sampling_func:str, bias_range:list, prepare_excited:bool=True, bias_sample_func:str='linspace', time_pts_or_step:int|float=100,Whileloop:bool=False, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * time_range: {"q0":[time_start, time_end], ...}\n
            * whileloop: bool, use while loop or not.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step.\n
            * bias_range:list, [z_amp_start, z_amp_end, pts/step]\n
            * bias_sampling_func:str, `linspace`(default) or `arange`. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(*time_range[q],time_pts_or_step)*1e9,4)*1e-9
        else:
            for q in time_range:
                self.time_samples[q] = sampling_func(*time_range[q],time_pts_or_step)

        self.avg_n = avg_n
        if bias_sample_func in ['linspace', 'arange']:
            self.bias_samples = eval(bias_sample_func)(*bias_range)
        else:
            raise ValueError(f"bias sampling function must be 'linspace' or 'arange' !")
        self.want_while = Whileloop
        self.prepare_1 = prepare_excited
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(time_range.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.aux_measurement.ZgateT1 import Zgate_T1
    
        dataset = Zgate_T1(self.QD_agent,self.meas_ctrl,self.time_samples,self.bias_samples,self.avg_n,self.execution,no_pi_pulse= not self.prepare_1)
        if self.execution:
            if self.save_dir is not None:
                if self.want_while:
                    self.JOBID = None
                self.save_path = os.path.join(self.save_dir,f"zgateT1_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self, new_QD_path:str=None,new_file_path:str=None, time_dep_plot:bool=False):
        """ If new file path was given, check all the data in that folder. """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()


            nc_paths = ZgateT1_dataReducer(fig_path)
            for q in nc_paths:
                if QD_savior.rotate_angle[q][0] != 0:
                    ref = QD_savior.rotate_angle[q]
                else:
                    eyeson_print(f"{q} rotation angle is 0, use contrast to analyze.")
                    ref = QD_savior.refIQ[q]

                ds = open_dataset(nc_paths[q])
                ANA = Multiplex_analyzer("auxA")
                ANA._import_data(ds,var_dimension=2,refIQ=ref)
                ANA._start_analysis(time_sort=time_dep_plot)
                ANA._export_result(nc_paths[q])
                ds.close()
            

    def WorkFlow(self, histo_counts:int=None):
        idx = 1
        start_time = datetime.now()
        while True:
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()  

            slightly_print(f"It's the {idx}-th measurement, about {round((datetime.now() - start_time).total_seconds()/60,1)} mins recorded.")
            
            if histo_counts is not None:
                # ensure the histo_counts you set is truly a number
                try: 
                    a = int(histo_counts)/100
                    self.want_while = True
                except:
                    raise TypeError(f"The arg `histo_counts` you set is not a number! We see it's {type(histo_counts)}...")
                if histo_counts == idx:
                    break
            idx += 1
            if not self.want_while:
                break
            
                
class QubitMonitor():
    def __init__(self, QD_path:str, save_dir:str, execution:bool=True):
        self.QD_path = QD_path
        self.save_dir = save_dir
        self.Execution = execution
        self.T1_time_range:dict = {}
        self.T2_time_range:dict = {}
        self.OS_target_qs:list = []
        self.echo_pi_num:int = 0
        self.a_little_detune_Hz:float = 0.1e6
        self.time_sampling_func = 'linspace'
        self.time_ptsORstep:int|float = 100
        self.OS_shots:int = 10000
        self.AVG:int = 300
        self.idx = 0

    def StartMonitoring(self):
        start_time = datetime.now()
        pi_num_dict = {}
        if self.T2_time_range is not None:
            for q in self.T2_time_range:
                pi_num_dict[q] = self.echo_pi_num
        while True:
            if self.T1_time_range is not None:
                if len(list(self.T1_time_range.keys())) != 0:
                    EXP = EnergyRelaxation(QD_path=self.QD_path,data_folder=self.save_dir)
                    EXP.SetParameters(self.T1_time_range,self.time_sampling_func,self.time_ptsORstep,1,self.AVG,self.Execution)
                    EXP.WorkFlow()

            if self.T2_time_range is not None:
                if len(list(self.T2_time_range.keys())) != 0:
                    EXP = CPMG(QD_path=self.QD_path,data_folder=self.save_dir)
                    EXP.SetParameters(self.T2_time_range,pi_num_dict,self.time_sampling_func,self.time_ptsORstep,1,self.AVG,self.Execution)
                    EXP.WorkFlow(freq_detune_Hz=self.a_little_detune_Hz)

            if self.OS_target_qs is not None:
                if  self.OS_shots != 0:
                    if len(self.OS_target_qs) == 0:
                        self.OS_target_qs = list(set(list(self.T1_time_range.keys())+list(self.T2_time_range.keys())))
                    EXP = SingleShot(QD_path=self.QD_path,data_folder=self.save_dir)
                    EXP.SetParameters(self.OS_target_qs,1,self.OS_shots,self.Execution)
                    EXP.WorkFlow()
            slightly_print(f"It's the {self.idx}-th measurement, about {round((datetime.now() - start_time).total_seconds()/3600,2)} hrs recorded.")
            self.idx += 1

    def TimeMonitor_analysis(self,New_QD_path:str=None,New_data_file:str=None,save_all_fit_fig:bool=False):
        if New_QD_path is not None:
            self.QD_path = New_QD_path
        if New_data_file is not None:
            self.save_dir = os.path.split(New_data_file)[0]

        QD_agent = QDmanager(self.QD_path)
        QD_agent.QD_loader()
        time_monitor_data_ana(QD_agent,self.save_dir,save_all_fit_fig)


class DragCali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, drag_coef_range:dict, coef_sampling_funct:str, coef_ptsORstep:int=100, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * piamp_coef_range: {"q0":[0.9, 1.1], "q1":[0.95, 1.15], ...]\n
            * amp_sampling_funct: str, `linspace` or `arange`.\n
            * pi_pair_num: list, like [2, 3] will try 2 exp, the first uses 2\*2 pi-pulse, and the second exp uses 3*2 pi-pulse
        """
        if coef_sampling_funct in ['linspace','logspace','arange']:
            sampling_func:callable = eval(coef_sampling_funct)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {coef_sampling_funct}")
        
        self.drag_coef_samples = {}
        for q in drag_coef_range:
           self.drag_coef_samples[q] = sampling_func(*drag_coef_range[q],coef_ptsORstep)
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(self.drag_coef_samples.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))

        
    def RunMeasurement(self):
        from qblox_drive_AS.Calibration_exp.DRAGcali import drag_cali 
    
        dataset = drag_cali(self.QD_agent,self.meas_ctrl,self.drag_coef_samples,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"DragCali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'dragcoef':
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    ANA = Multiplex_analyzer("c6")
                    ANA._import_data(ds,var_dimension=1,refIQ=ref)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(fig_path)
                    answer[var] = ANA.fit_packs
            
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.Waveformer.set_dragRatio_for(permi,answer[permi]["optimal_drag_coef"])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.Waveformer.set_dragRatio_for(q, answer[q]["optimal_drag_coef"])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 


class XGateErrorTest(ExpGovernment):
    """ Helps you get the **Dressed** cavities. """
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location
    
    @RawDataPath.setter
    def RawDataPath(self,path:str):
        self.__raw_data_location = path

    def SetParameters(self, target_qs:list, shots:int=10000, MaxGate_num:int=300, execution:bool=True, use_untrained_wf:bool=False):
        """ 
        ### Args:\n
        * target_qs: list, like ["q0", "q1", ...]
        """
        self.use_time_label:bool = False
        self.avg_n = shots
        self.Max_Gate_num = MaxGate_num
        self.execution = execution
        self.use_de4t_wf = use_untrained_wf
        self.target_qs = target_qs
        



    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and driving atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    
    def RunMeasurement(self):
        from qblox_drive_AS.aux_measurement.GateErrorTest import XGateError_single_shot
       

        dataset = XGateError_single_shot(self.QD_agent,self.target_qs,self.Max_Gate_num,self.avg_n,self.use_de4t_wf,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"XGateErrorTest_{datetime.now().strftime('%Y%m%d%H%M%S') if (self.JOBID is None or self.use_time_label) else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
        """ if histo_ana, it will check all the data in the same folder with the given new_file_path """
    
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            
            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                ANA = Multiplex_analyzer("t1")
                ANA._import_data(ds,var_dimension=0,fq_Hz=QD_savior.quantum_device.get_element(var).clock_freqs.f01())
                ANA._start_analysis(var_name=var)
                pic_path = os.path.join(fig_path,f"{var}_XGateErrorTest_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                ANA._export_result(pic_path)

                answer[var] = ANA.fit_packs 
                highlight_print(f"{var} X-gate phase error ~ {round(answer[var]['f'], 3)} mrad")
                
            ds.close()

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 




if __name__ == "__main__":
    EXP = Ramsey(QD_path="")
    EXP.execution = 1
    EXP.RunAnalysis(new_QD_path="/Users/ratiswu/Desktop/FTP-ASqcMeas/ramsey_updates/qblox_ExpConfigs_20250117192846/DR1#11_SumInfo.pkl",new_file_path="/Users/ratiswu/Desktop/FTP-ASqcMeas/ramsey_updates/Ramsey_20250117192846.nc")
    