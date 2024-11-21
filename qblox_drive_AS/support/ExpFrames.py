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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.wideCS import plot_S21
        ds = open_dataset(self.__raw_data_location)

        QD_savior = QDmanager(self.QD_path)
        QD_savior.QD_loader()
        if new_QD_dir is None:
            new_QD_dir = self.QD_path
        else:
            new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

        plot_S21(ds,self.save_fig_path)
        ds.close()
        QD_savior.QD_keeper(new_QD_dir)


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
        from qblox_drive_AS.SOP.CavitySpec import QD_RO_init, Cavity_spec
        QD_RO_init(self.QD_agent,self.freq_range)
        dataset = Cavity_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"zoomCS_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.CavitySpec import CS_ana
        if self.execution:
            ds = open_dataset(self.__raw_data_location)

            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            CS_ana(QD_savior,ds,self.save_dir)
            ds.close()
            QD_savior.QD_keeper(new_QD_dir)


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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.PowCavSpec import plot_powerCavity_S21
        if self.execution:
            ds = open_dataset(self.__raw_data_location)

            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            plot_powerCavity_S21(ds,QD_savior,self.save_dir)
            ds.close()
        # QD_savior.QD_keeper(new_QD_dir)


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
        from qblox_drive_AS.SOP.CavitySpec import Cavity_spec, QD_RO_init
        QD_RO_init(self.QD_agent,self.freq_range)
        for q in self.ro_amp:
            self.QD_agent.quantum_device.get_element(q).measure.pulse_amp(self.ro_amp[q])
        dataset = Cavity_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"dressedCS_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.CavitySpec import CS_ana
        if self.execution:
            ds = open_dataset(self.__raw_data_location)

            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])
            for q in self.ro_amp:
                QD_savior.quantum_device.get_element(q).measure.pulse_amp(self.ro_amp[q])
            CS_ana(QD_savior,ds,self.save_dir,keep_bare=False)
            ds.close()
            QD_savior.QD_keeper(new_QD_dir)


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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        if self.execution:
            ds = open_dataset(self.__raw_data_location)
            for var in ds.data_vars:
                ANA = Multiplex_analyzer("m5")
                if var.split("_")[-1] != 'freq':
                    ANA._import_data(ds,2)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(self.save_dir)
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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.FluxCavSpec import update_flux_info_in_results_for
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            answer = {}
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'freq':
                    ANA = Multiplex_analyzer("m6")
                    ANA._import_data(ds,2)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(self.save_dir)
                    answer[var] = ANA.fit_packs
            ds.close()
            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                update_flux_info_in_results_for(QD_savior,permi,answer[permi])
                QD_savior.QD_keeper(new_QD_dir)
            elif permi in ["all",'y','yes']:
                for q in list(answer.keys()):
                    update_flux_info_in_results_for(QD_savior,q,answer[q])
                QD_savior.QD_keeper(new_QD_dir)
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
        for i in self.ro_amp:
            if self.ro_amp[i] != 1:
                self.ask_save = True


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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RefIQ import IQ_ref_ana
        if self.execution:
            ds = open_dataset(self.__raw_data_location)

            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])
            
            answer = {}
            for q in ds.data_vars:
                answer[q] = IQ_ref_ana(ds,q,self.save_dir)
            ds.close()
            if self.ask_save:
                permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
                if permi in list(answer.keys()):
                    QD_savior.memo_refIQ({permi:answer[permi]})
                    QD_savior.QD_keeper(new_QD_dir)
                elif permi in ["all",'y','yes']:
                    QD_savior.memo_refIQ(answer)
                    QD_savior.QD_keeper(new_QD_dir)
                else:
                    print("Updating got denied ~")
            else:
                QD_savior.QD_keeper(new_QD_dir)


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

    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.Cnti2Tone import update_2toneResults_for
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'freq':
                    ANA = Multiplex_analyzer("m8")     
                    ANA._import_data(ds,2,self.QD_agent.refIQ[var] if self.QD_agent.rotate_angle[var] == 0 else [self.QD_agent.rotate_angle[var]],QS_fit_analysis)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(self.save_dir)
                    if ANA.fit_packs != {}:
                        analysis_result = QS_fit_analysis(ANA.fit_packs[var]["contrast"],f=ANA.fit_packs[var]["xyf_data"])
                        update_2toneResults_for(QD_savior,var,{str(var):analysis_result},ANA.xyl[0])
            ds.close()
            QD_savior.QD_keeper(new_QD_dir)

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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.FluxQubit import update_by_fluxQubit
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            answer = {}
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'freq':
                    ANA = Multiplex_analyzer("m9") 
                    ANA._import_data(ds,2,QD_savior.refIQ[var],QS_fit_analysis)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(self.save_dir)
                    if len(list(ANA.fit_packs.keys())) != 0: answer[var] = ANA.fit_packs        
            ds.close()
            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no :").lower()
            if permi in list(answer.keys()):
                update_by_fluxQubit(QD_savior,answer[q],q)
                QD_savior.QD_keeper(new_QD_dir)
            elif permi in ["all",'y','yes']:
                for q in list(answer.keys()):
                    update_by_fluxQubit(QD_savior,answer[q],q)
                QD_savior.QD_keeper(new_QD_dir)
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

    def SetParameters(self, pi_amp:dict, pi_dura:dict, pi_amp_sampling_func:str, pi_amp_pts_or_step:float=100, avg_n:int=100, execution:bool=True, OSmode:bool=False):
        """ ### Args:
            * pi_amp: {"q0":[pi_amp_start, pi_amp_end], ...}\n
            * pi_dura: {"q0":pi_duration_in_seconds, ...}\n
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
            
        self.pi_dura = pi_dura
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(pi_amp.keys())
        

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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RabiOsci import conditional_update_qubitInfo
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'piamp':
                    ANA = Multiplex_analyzer("m11")      
                    ANA._import_data(ds,1,self.QD_agent.refIQ[var] if self.QD_agent.rotate_angle[var][0] == 0 else self.QD_agent.rotate_angle[var])
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(self.save_dir)
                    conditional_update_qubitInfo(QD_savior,ANA.fit_packs,var)  

            ds.close()
            QD_savior.QD_keeper(new_QD_dir)
            



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

    def SetParameters(self, pi_dura:dict, pi_amp:dict, pi_dura_sampling_func:str, pi_dura_pts_or_step:float=100, avg_n:int=100, execution:bool=True, OSmode:bool=False):
        """ ### Args:
            * pi_amp: {"q0": pi-amp in V, ...}\n
            * pi_dura: {"q0":[pi_duration_start, pi_duration_end], ...}\n
            * pi_dura_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * pi_dura_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
        if pi_dura_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(pi_dura_sampling_func)
        else:
            sampling_func:callable = linspace
        
        self.pi_dura_samples = {}
        for q in pi_dura:
            time_sample_sec = list(set([round_to_nearest_multiple_of_multipleNS(x) for x in sampling_func(*pi_dura[q],pi_dura_pts_or_step)*1e9]))
            time_sample_sec.sort()
            self.pi_dura_samples[q] = array(time_sample_sec)*1e-9
    

        self.pi_amp = pi_amp
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(pi_dura.keys())
        

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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RabiOsci import conditional_update_qubitInfo
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'pidura':
                    ANA = Multiplex_analyzer("m11")      
                    ANA._import_data(ds,1,self.QD_agent.refIQ[var] if self.QD_agent.rotate_angle[var][0] == 0 else self.QD_agent.rotate_angle[var])
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(self.save_dir)
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
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
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


    def RunAnalysis(self,new_QD_dir:str=None,histo_ana:bool=False):
        """ User callable analysis function pack """
    
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])
            if not histo_ana:
                ds = open_dataset(self.__raw_data_location)
                for var in ds.data_vars:
                    ANA = Multiplex_analyzer("m14")
                    ANA._import_data(ds[var]*1000,var_dimension=0,fq_Hz=QD_savior.quantum_device.get_element(var).clock_freqs.f01())
                    ANA._start_analysis()
                    pic_path = os.path.join(self.save_dir,f"{var}_SingleShot_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                    ANA._export_result(pic_path)
                    highlight_print(f"{var} rotate angle = {round(ANA.fit_packs['RO_rotation_angle'],2)} in degree.")
                    QD_savior.rotate_angle[var] = [ANA.fit_packs["RO_rotation_angle"]]
                ds.close()
                
                QD_savior.QD_keeper(new_QD_dir)
            else:
                eff_T, thermal_pop = {}, {}
                files = sort_timeLabel([os.path.join(self.save_dir,name) for name in os.listdir(self.save_dir) if (os.path.isfile(os.path.join(self.save_dir,name)) and name.split(".")[-1]=='nc')])
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
                    Data_manager().save_histo_pic(QD_savior,eff_T,qubit,mode="ss",pic_folder=self.save_dir)
                    Data_manager().save_histo_pic(QD_savior,thermal_pop,qubit,mode="pop",pic_folder=self.save_dir)

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
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort(array(list(set([round_to_nearest_multiple_of_multipleNS(x,4) for x in sampling_func(*time_range[q],time_pts_or_step)*1e9])))*1e-9)
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
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
        
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
                    ANA._export_result(self.save_dir)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_T2_for(ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            QD_savior.QD_keeper(new_QD_dir)
            

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
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        self.spin_num = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort(array(list(set([round_to_nearest_multiple_of_multipleNS(x,8) for x in sampling_func(*time_range[q],time_pts_or_step)*1e9])))*1e-9)
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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
        
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
                    ANA._export_result(self.save_dir)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_echoT2_for(ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            QD_savior.QD_keeper(new_QD_dir)
            

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
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        self.spin_num = pi_num
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort(array(list(set([round_to_nearest_multiple_of_multipleNS(x,(1+int(pi_num[q]))*4) for x in sampling_func(*time_range[q],time_pts_or_step)*1e9])))*1e-9)
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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
        
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
                    ANA._export_result(self.save_dir)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_echoT2_for(ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            QD_savior.QD_keeper(new_QD_dir)
            

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
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort(array(list(set([round_to_nearest_multiple_of_multipleNS(x,4) for x in sampling_func(*time_range[q],time_pts_or_step)*1e9])))*1e-9)
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
        from qblox_drive_AS.SOP.T1 import T1
    
        dataset = T1(self.QD_agent,self.meas_ctrl,self.time_samples,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"T1_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
        
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
                    ANA._export_result(self.save_dir)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_T1_for(ANA.fit_packs["median_T1"],var)

            ds.close()
            QD_savior.QD_keeper(new_QD_dir)
            

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

    def SetParameters(self, target_qs:list, evo_time:float=0.5e-6, detune:float=10e6, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * target_qs: ["q0", "q1", ...]\n
        """
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
    
        self.time_samples = {}
        self.detune = detune
        for q in target_qs:
            self.time_samples[q] = sort(array(list(set([round_to_nearest_multiple_of_multipleNS(x,4) for x in  linspace(0,evo_time,100)*1e9])))*1e-9)

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
            self.QD_agent.quantum_device.get_element(q).clock_freqs.f01(self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()+self.detune)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.T2 import Ramsey
    
        dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,1,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"XYFcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
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
                    ANA._export_result(self.save_dir)
                     
                    answer[var] = self.detune-ANA.fit_packs['freq']
                    highlight_print(f"{var}: actual detune = {round(answer[var]*1e-6,4)} MHz")
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).clock_freqs.f01(QD_savior.quantum_device.get_element(permi).clock_freqs.f01()+answer[permi])
                QD_savior.QD_keeper(new_QD_dir)
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).clock_freqs.f01(QD_savior.quantum_device.get_element(q).clock_freqs.f01()+answer[q])
                QD_savior.QD_keeper(new_QD_dir)
            else:
                print("Updating got denied ~")

            QD_savior.QD_keeper(new_QD_dir)
            

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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'rof':
                    ANA = Multiplex_analyzer("c1")
                    ANA._import_data(ds,var_dimension=1)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(self.save_dir)
                    answer[var] = ANA.fit_packs[var]["optimal_rof"]
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).clock_freqs.readout(answer[permi])
                QD_savior.QD_keeper(new_QD_dir)
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).clock_freqs.readout(answer[q])
                QD_savior.QD_keeper(new_QD_dir)
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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            
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
                    fit_pic_folder = Data_manager().get_today_picFolder()
                    ANA._export_result(fit_pic_folder)
                    fit_packs = ANA.fit_packs
            ds.close()

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
            * piamp_coef_range: {"q0":[0.4, 0.6], "q1":[0.45, 0.55], ...], must be close to 0.5 because it's a half pi-pulse.\n
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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            
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
                    ANA._export_result(self.save_dir)
                    fit_packs = ANA.fit_packs

            ds.close()

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


    def RunAnalysis(self,new_QD_dir:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            ds = open_dataset(self.__raw_data_location)
            
            for var in ds.data_vars:
                if var.split("_")[-1] != 'rol':
                    ANA = Multiplex_analyzer("c5")
                    ANA._import_data(ds,var_dimension=1)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(self.save_dir)
                    
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
        from qblox_drive_AS.SOP.RabiOsci import round_to_nearest_multiple_of_multipleNS
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in time_range:
                self.time_samples[q] = sort(array(list(set([round_to_nearest_multiple_of_multipleNS(x,4) for x in sampling_func(*time_range[q],time_pts_or_step)*1e9])))*1e-9)
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
                self.save_path = os.path.join(self.save_dir,f"zgateT1_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self, new_QD_dir:str=None, time_dep_plot:bool=False):
        """ User callable analysis function pack """
        
        if self.execution:
            QD_savior = QDmanager(self.QD_path)
            QD_savior.QD_loader()
            if new_QD_dir is None:
                new_QD_dir = self.QD_path
            else:
                new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

            nc_paths = ZgateT1_dataReducer(self.save_dir)
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
            

    def WorkFlow(self):
        idx = 1
        start_time = datetime.now()
        while True:
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()  

            slightly_print(f"It's the {idx}-th measurement, about {round((datetime.now() - start_time).total_seconds()/60,1)} mins recorded.")
            
            idx += 1
            if not self.want_while:
                break
                
