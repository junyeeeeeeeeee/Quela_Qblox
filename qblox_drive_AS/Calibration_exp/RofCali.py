
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from quantify_scheduler.gettables import ScheduleGettable
from numpy import array, arange
from numpy import nan as NaN
from qblox_drive_AS.support import Data_manager, check_acq_channels

from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import BinMode, Schedule, pulse_preview, X, Reset, Measure, IdlePulse, SetClockFrequency, ClockResource

class ROFcalibrationPS(ScheduleConductor):
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
        frequencies: dict,
        state:int,
        repetitions:int=1,    
    ) -> Schedule:
        
        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
        sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)

        for acq_idx in range(sameple_idx):    
            
            reset = sched.add(Reset(*qubits2read))
            for qubit_idx, q in enumerate(qubits2read):
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
                
                sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
                sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")
                if state:
                    pi_pulse = sched.add(X(q), ref_op=reset)

            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE), ref_op=pi_pulse if state else reset)
                
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(0,len(list(list(self._ro_elements.values())[0])))
        self.__prepared_states = array([0, 1])
        
        self.__ro_f_origin = {}
        for q in self._ro_elements:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            self.__ro_f_origin[q] = qubit_info.clock_freqs.readout()
            qubit_info.clock_freqs.readout(NaN)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
           
        self.__freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
        self.__freq.batched = True
        self.__state =  ManualParameter(name="prepared_state", unit="", label="state")
        self.__state.batched = False


        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))
        self.__spec_sched_kwargs = dict(   
        frequencies=self._ro_elements,
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
            num_channels=len(list(self._ro_elements.keys())),
            )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables([self.__freq, self.__state])
            self.meas_ctrl.setpoints_grid([self.__datapoint_idx, self.__prepared_states])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
        
            self.__spec_sched_kwargs['frequencies']= preview_para
            self.__spec_sched_kwargs['state']= array([1])
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rs_ds = self.meas_ctrl.run("ROF calibration")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                freq_values = 2*self.__prepared_states.shape[0]*list(self._ro_elements[q])
                i_data = array(rs_ds[f'y{2*q_idx}']).reshape(self.__prepared_states.shape[0],array(self._ro_elements[q]).shape[0])
                q_data = array(rs_ds[f'y{2*q_idx+1}']).reshape(self.__prepared_states.shape[0],array(self._ro_elements[q]).shape[0])
                dict_[q] = (["mixer","state","rof"],array([i_data,q_data]))
                dict_[f'{q}_rof'] = (["mixer","state","rof"],array(freq_values).reshape(2,self.__prepared_states.shape[0],array(self._ro_elements[q]).shape[0]))

            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"state":self.__prepared_states,"rof":self.__datapoint_idx})
            
            ds.attrs["execution_time"] = Data_manager().get_time_now()
            ds.attrs["method"] = "Average"
            ds.attrs["system"] = "qblox"
            for q in self.__ro_f_origin:
                ds.attrs[f"{q}_ori_rof"] = self.__ro_f_origin[q]
            self.dataset = ds
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)
