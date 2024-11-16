from numpy import ndarray
from abc import ABC
import os
from datetime import datetime
from xarray import Dataset
from qblox_drive_AS.support.QDmanager import QDmanager
from xarray import open_dataset
from numpy import array, linspace, arange, logspace
from abc import abstractmethod
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down
from quantify_scheduler.helpers.collections import find_port_clock_path


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
        qrmRF_slot_idx = int(find_port_clock_path(self.QD_agent.quantum_device.hardware_config(),"q:res",f"{self.target_qs[int(len(self.target_qs)-self.counter)]}.ro")[1][-1])
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
        ds = open_dataset(self.__raw_data_location)

        QD_savior = QDmanager(self.QD_path)
        QD_savior.QD_loader()
        if new_QD_dir is None:
            new_QD_dir = self.QD_path
        else:
            new_QD_dir = os.path.join(new_QD_dir,os.path.split(self.QD_path)[-1])

        CS_ana(ds,self.save_dir)
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
        # set self.freq_range
        for q in self.tempor_freq[0]:
            rof = self.QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            self.freq_range[q] = linspace(rof+self.tempor_freq[0][q][0],rof+self.tempor_freq[0][q][1],self.tempor_freq[1])

        
    def RunMeasurement(self):
        from qblox_drive_AS.SOP.PowCavSpec import PowerDep_spec
        dataset = PowerDep_spec(self.QD_agent,self.meas_ctrl,self.freq_range,self.roamp_samples,self.avg_n,self.execution)
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

        
        
    
