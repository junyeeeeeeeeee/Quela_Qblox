
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from quantify_scheduler.gettables import ScheduleGettable
from numpy import array, arange
from qblox_drive_AS.support import Data_manager, check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import Measure, Schedule, pulse_preview, BinMode, Reset, X, IdlePulse
from qblox_drive_AS.support.Pulser import ScheduleConductor
 

class ROLcalibrationPS(ScheduleConductor):
    def __init__(self):
        super().__init__() 
        self._avg_n:int = 100
        self._power_samples:dict = []
    
    
    @property
    def power_samples(self):
        return self._power_samples
    @power_samples.setter
    def power_samples(self,powers:dict):
        self._power_samples = powers

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
        Ramps: dict,
        state:int,
        repetitions:int=1,    
    ) -> Schedule:
        
        qubits2read = list(Ramps.keys())
        sameple_idx = array(Ramps[qubits2read[0]]).shape[0]
        sched = Schedule("ROL calibration",repetitions=repetitions)

        for acq_idx in range(sameple_idx):    
            reset = sched.add(Reset(*qubits2read))
            for qubit_idx, q in enumerate(qubits2read):
                R_amp = Ramps[q][acq_idx]
                
                if state:
                    pi_pulse = sched.add(X(q), ref_op=reset)
            
                sched.add(Measure(q, acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE, pulse_amp=R_amp), ref_op=pi_pulse if state else reset)
                
        self.schedule =  sched  
        return sched
    

    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(0,len(list(list(self._power_samples.values())[0])))
        self.__r_amps = {}
        quantum_device = self.QD_agent.quantum_device
        for q in self._power_samples:
            self.__r_amps[q] = quantum_device.get_element(q).measure.pulse_amp()*self._power_samples[q]
            
        self.QD_agent = check_acq_channels(self.QD_agent, list(self._power_samples.keys()))
        
        self.__ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
        self.__ro_pulse_amp.batched = True
        self.__state =  ManualParameter(name="prepared_state", unit="", label="state")
        self.__state.batched = False
        self.__prepared_states = array([0, 1])

        self.__spec_sched_kwargs = dict(   
        Ramps=self.__r_amps,
        state=self.__state
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
            self.QD_agent.quantum_device,
            schedule_function=self.__PulseSchedule__, 
            schedule_kwargs=self.__spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(self._power_samples.keys())),
            )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables([self.__ro_pulse_amp, self.__state])
            self.meas_ctrl.setpoints_grid([self.__datapoint_idx, self.__prepared_states])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._power_samples:
                preview_para[q] = self.__r_amps[q][:n_s]
             
            self.__spec_sched_kwargs['state'] = array([1])
            self.__spec_sched_kwargs['Ramps']= preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rp_ds = self.meas_ctrl.run("PowerCavity")
            dict_ = {}
            for q_idx, q in enumerate(list(self.__r_amps.keys())):
                rol = 2*self.__prepared_states.shape[0]*list(self._power_samples[q])
                i_data = array(rp_ds[f'y{2*q_idx}']).reshape(self.__prepared_states.shape[0],self.__r_amps[q].shape[0])
                q_data = array(rp_ds[f'y{2*q_idx+1}']).reshape(self.__prepared_states.shape[0],self.__r_amps[q].shape[0])
                dict_[q] = (["mixer","state","rol"],array([i_data,q_data]))
                dict_[f'{q}_rol'] = (["mixer","state","rol"],array(rol).reshape(2,self.__prepared_states.shape[0],self.__r_amps[q].shape[0]))
            
            dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"state":array(["g","e"]),"rol":self.__datapoint_idx})
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["method"] = "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)

