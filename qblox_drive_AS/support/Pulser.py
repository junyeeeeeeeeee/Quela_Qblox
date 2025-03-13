""" Connection between ExpFrames and pulse-schedule """
from abc import ABC, abstractmethod
from qblox_drive_AS.support import QDmanager
from quantify_core.measurement.control import MeasurementControl
from xarray import Dataset
from quantify_scheduler.schedules.schedule import Schedule



class ScheduleConductor(ABC):
    """ Some adjustable parameters please name started with "_" like: "self._pi_dura" """
    def __init__(self):
        self.QD_agent = QDmanager()
        self.meas_ctrl =  MeasurementControl("dummy")
        self._execution:bool = True
        self.__dataset:Dataset = {}
        self.__sche = Schedule("dummy")
    
    @property
    def schedule(self):
        return self.__sche
    @schedule.setter
    def schedule(self,ps:Schedule):
        self.__sche = ps
    @property
    def dataset(self):
        return self.__dataset
    @dataset.setter
    def dataset(self,ds:Dataset):
        self.__dataset = ds
    @property
    def execution( self ):
        return self._execution
    @execution.setter
    def execution( self, exec:bool):
        self._execution = exec
    
    @abstractmethod
    def __PulseSchedule__(self,*args,**kwargs):
        """ Write your pulse schedule here """
        pass

    @abstractmethod
    def __SetParameters__(self,*args,**kwargs):
        """ Something your pulse schedule will use """
        pass

    @abstractmethod
    def __Compose__(self,*args,**kwargs):
        """ Combine parameters and schedule together """
        pass

    @abstractmethod
    def __RunAndGet__(self,*args,**kwargs):
        """ Run schedule and get the dataset """
        pass
    
    def get_adjsutable_paras(self,display:bool=False)->list:
        """ Return what parameters can be set manually """
    
        can_be_set_paras = [
            (attr, getattr(self, attr)) for attr in dir(self) if  attr.startswith("_") 
            and not attr.startswith("__") 
            and not callable(getattr(self, attr)) 
            and not isinstance(getattr(type(self), attr, None), property)
            and not attr.startswith("_abc_")
            and not attr.startswith("_ScheduleConductor_")
        ]
        if display:
            for par in can_be_set_paras:
                print(f"{par[0]}: Type {type(par[1])}")
        
        return can_be_set_paras
    
    def run(self):
        self.__SetParameters__()
        self.__Compose__()
        self.__RunAndGet__()
        