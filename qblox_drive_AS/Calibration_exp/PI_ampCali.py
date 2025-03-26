"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from numpy import array, arange
from qblox_drive_AS.support import Data_manager, BasicTransmonElement
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import check_acq_channels
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import BinMode, Schedule, pulse_preview, X, Reset, IdlePulse, Measure, electrical_delay

class PiAcalibrationPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._ro_elements:dict = {}
        self._pi_pairs:list = []
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles
    @property
    def pi_pairs(self):
        return self._pi_pairs
    @pi_pairs.setter
    def pi_pairs(self, pairs:list):
        if not isinstance(pairs, list):
            raise TypeError("pairs must be a list !")
        self._pi_pairs = pairs
    
    

    def __PulseSchedule__(self, 
        new_pi_amp:dict,
        pi_pair_num:any,
        repetitions:int=1,
    )-> Schedule:
        qubits2read = list(new_pi_amp.keys())
        sched = Schedule("Pi amp modification", repetitions=repetitions)

        for acq_idx in range(array(new_pi_amp[qubits2read[0]]).shape[0]):
            align_pulse = sched.add(IdlePulse(4e-9))
            for q in qubits2read:
                new_amp180 = new_pi_amp[q][acq_idx]
                
                sched.add(Reset(q), ref_op=align_pulse)
                
                for pi_num in range(pi_pair_num):
                    for pi_idx in range(2):
                        pi = sched.add(X(q, amp180=new_amp180))

            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol="SSBIntegrationComplex", bin_mode=BinMode.AVERAGE), rel_time=electrical_delay)
                        
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(len(list(list(self._ro_elements.values())[0])))
        self.__amp180_samples = {}
        
        for q in self._ro_elements:
            qubit_info:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
            self.__amp180_samples[q] = qubit_info.rxy.amp180() * self._ro_elements[q] 
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
           
        self.__amp180 = ManualParameter(name="amp180", unit="V", label="amplitude")
        self.__amp180.batched = True
        self.__pi_pair =  ManualParameter(name="pi_pairs", unit="", label="number")
        self.__pi_pair.batched = False


        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))
        self.__spec_sched_kwargs = dict(   
        new_pi_amp=self.__amp180_samples,
        pi_pair_num=self.__pi_pair,
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
            self.meas_ctrl.settables([self.__amp180, self.__pi_pair])
            self.meas_ctrl.setpoints_grid([self.__datapoint_idx, self._pi_pairs])
        
        else:
            preview_para = {}
            for q in self.__amp180_samples:
                preview_para[q] = array([self.__amp180_samples[q][1], self.__amp180_samples[q][-2]])
            self.__spec_sched_kwargs['new_pi_amp']= preview_para
            self.__spec_sched_kwargs['pi_pair_num']= self._pi_pairs[0]
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rs_ds = self.meas_ctrl.run("pi-amp calibration")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                coefs = 2*len(self._pi_pairs)*list(self._ro_elements[q])
                i_data = array(rs_ds[f'y{2*q_idx}']).reshape(len(self._pi_pairs),self.__datapoint_idx.shape[0])
                q_data = array(rs_ds[f'y{2*q_idx+1}']).reshape(len(self._pi_pairs),self.__datapoint_idx.shape[0])
                dict_[q] = (["mixer", "PiPairNum", "PiAmpCoef"],array([i_data,q_data]))
                dict_[f'{q}_PIcoef'] = (["mixer", "PiPairNum", "PiAmpCoef"],array(coefs).reshape(2,len(self._pi_pairs),self.__datapoint_idx.shape[0]))

            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"PiPairNum":array(self._pi_pairs),"PiAmpCoef":self.__datapoint_idx})
            
            ds.attrs["execution_time"] = Data_manager().get_time_now()
            ds.attrs["method"] = "Average"
            ds.attrs["system"] = "qblox"
            self.dataset = ds
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)


    

