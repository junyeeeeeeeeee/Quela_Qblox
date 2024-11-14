from numpy import ndarray
from numpy import array, linspace, arange
from abc import ABC, abstractmethod
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down
from quantify_scheduler.helpers.collections import find_port_clock_path
import os
from datetime import datetime

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
        self.JOBID = JOBID

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
        self.dataset = wideCS(self.readout_module,self.freq_start,self.freq_end,self.freq_pts)
        if self.save_dir is not None:
            path = os.path.join(self.save_dir,f"BroadBandCS_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
            self.dataset.to_netcdf(path+".nc")
            self.save_fig_path = path+".png"
        else:
            self.save_fig_path = None
        

    def RunAnalysis(self):
        from qblox_drive_AS.SOP.wideCS import plot_S21
        plot_S21(self.dataset,self.save_fig_path)

    
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)
        self.counter -= 1

    def WorkFlow(self):
        while self.counter > 0 :
            self.PrepareHardware()

            self.RunMeasurement()

            self.RunAnalysis()

            self.CloseMeasurement()


        
        

        
        
    
