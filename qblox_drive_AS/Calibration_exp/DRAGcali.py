import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from numpy import array, arange, ndarray
from qblox_drive_AS.support import Data_manager
from quantify_scheduler.gettables import ScheduleGettable

from qblox_drive_AS.support import check_acq_channels
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import BinMode, Schedule, pulse_preview, X, Y, Y90, X90, Reset, electrical_delay
from quantify_scheduler.operations.gate_library import Measure
from quantify_scheduler.operations.pulse_library import IdlePulse

class DRAGcalibrationPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._ro_elements:dict = {}
        self._avg_n:int = 100
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles

    @property
    def n_avg(self):
        return self._avg_n
    @n_avg.setter
    def n_avg(self, avg:int):
        self._avg_n = avg

    @property
    def execution(self):
        return self._execution
    @execution.setter
    def execution(self, execu:bool):
        self._execution = execu

    def __PulseSchedule__(self, 
        drag_samples:ndarray,
        qubits2read:list,
        operation:str,
        repetitions:int=1,    
    ) -> Schedule:
        

        sched = Schedule("Motzoi Calibration",repetitions=repetitions)

        for acq_idx, motzoi in enumerate(drag_samples):    
            align_pulse = sched.add(IdlePulse(4e-9))
            for q in qubits2read:
                reset = sched.add(Reset(q), ref_op=align_pulse)
                match operation:
                    case "(X,Y/2)":
                        sched.add(X(q, motzoi=motzoi), ref_op=reset)
                        final_pulse = sched.add(Y90(q, motzoi=motzoi))
                    case "(Y,X/2)":
                        sched.add(Y(q, motzoi=motzoi), ref_op=reset)
                        final_pulse = sched.add(X90(q, motzoi=motzoi))
                    
            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE), rel_time=electrical_delay, ref_op=final_pulse)
                
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(0,len(list(list(self._ro_elements.values())[0])))
        self.__motzoi_samples = array(list(self._ro_elements.values())[0])
        self.__operations = array(["(X,Y/2)", "(Y,X/2)"])
        
        for q in self._ro_elements:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
           
        self.__motzoi = ManualParameter(name="motzoi", unit="", label="motzoi")
        self.__motzoi.batched = True
        self.__op =  ManualParameter(name="operation", unit="", label="operation")
        self.__op.batched = False


        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))
        self.__spec_sched_kwargs = dict(   
        drag_samples=self.__motzoi_samples,
        qubits2read=list(self._ro_elements.keys()),
        operation=self.__op
        )

    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
            self.QD_agent.quantum_device,
            schedule_function=self.__PulseSchedule__, 
            schedule_kwargs=self.__spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(self._ro_elements.keys())),
            )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables([self.__motzoi, self.__op])
            self.meas_ctrl.setpoints_grid([self.__motzoi_samples, self.__operations])
        
        else:
            preview_para = array([self.__motzoi_samples[1],self.__motzoi_samples[-2]])
            self.__spec_sched_kwargs['drag_samples']= preview_para
            self.__spec_sched_kwargs['operation']= array(["(X,Y/2)"])
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rs_ds = self.meas_ctrl.run("DRAG calibration")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                coefs = 2*2*list(self._ro_elements[q])
                i_data = array(rs_ds[f'y{2*q_idx}']).reshape(self.__operations.shape[0],array(self._ro_elements[q]).shape[0])
                q_data = array(rs_ds[f'y{2*q_idx+1}']).reshape(self.__operations.shape[0],array(self._ro_elements[q]).shape[0])
                dict_[q] = (["mixer", "operations", "dragCoef"],array([i_data,q_data]))
                dict_[f'{q}_dragcoef'] = (["mixer", "operations", "dragCoef"],array(coefs).reshape(2,self.__operations.shape[0],self.__motzoi_samples.shape[0]))

            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"operations":self.__operations,"dragCoef":self.__motzoi_samples})
            
            ds.attrs["execution_time"] = Data_manager().get_time_now()
            ds.attrs["method"] = "Average"
            ds.attrs["system"] = "qblox"
            self.dataset = ds
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)
