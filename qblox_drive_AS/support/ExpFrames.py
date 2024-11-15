from numpy import ndarray
from abc import ABC
import os
from datetime import datetime
from xarray import Dataset
from qblox_drive_AS.support.QDmanager import QDmanager
from xarray import open_dataset
from numpy import array, linspace, arange
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
        print("### initialized")
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


        
        

        
        
    
