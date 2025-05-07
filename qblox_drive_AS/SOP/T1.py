import time
from xarray import Dataset
from numpy import array, arange
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, BinMode, Measure, Reset, X, ConditionalReset, IdlePulse
from quantify_scheduler.operations.gate_library import Reset


class EnergyRelaxPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._time_samples:dict = {}
        self._repeat:int = 1
       

    @property
    def time_samples( self ):
        return self._time_samples
    @time_samples.setter
    def set_time_samples(self, time_samples:dict):
        if not isinstance(time_samples,dict):
            raise TypeError("Arg 'time_samples' must be a dict !")
        self._time_samples = time_samples
    @property
    def repeat( self ):
        return self._repeat
    @repeat.setter
    def set_repeat(self, repeat_num:int):
        if not isinstance(repeat_num,(int,float)):
            raise TypeError("Ard 'repeat_num' must be a int or float !")
        self._repeat = int(repeat_num)



    def __PulseSchedule__(self, 
        freeduration:dict,
        repetitions:int=1,
        singleshot:bool=False,
        activeReset:bool=False
        )->Schedule:

        qubits2read = list(freeduration.keys())
        free_Dus = array(freeduration[qubits2read[0]])
        sameple_idx = free_Dus.shape[0]
        sched = Schedule("T1", repetitions=repetitions)

        for acq_idx in range(sameple_idx):
            align_pulse = sched.add(IdlePulse(4e-9))
            for qubit_idx, q in enumerate(qubits2read):
                if activeReset:
                    ring_down_wait = sched.add(IdlePulse(1e-6), ref_op=align_pulse)
                    reset = sched.add(
                        ConditionalReset(q, acq_index=acq_idx,acq_channel=qubit_idx+len(qubits2read)),
                        label=f"Reset {q} {acq_idx}",
                        ref_op=ring_down_wait
                    )
                else:
                    reset = sched.add(Reset(q), ref_op=align_pulse)

                sched.add(X(q), ref_op=reset)
            sched.add(Measure(*qubits2read,acq_index=acq_idx, acq_protocol='SSBIntegrationComplex' if not activeReset else 'ThresholdedAcquisition', bin_mode=BinMode.APPEND if singleshot else BinMode.AVERAGE), rel_time=free_Dus[acq_idx])
            
        self.schedule = sched
        return sched

    def __SetParameters__(self, *args, **kwargs):
        
        self.__time_data_idx = arange(len(list(self._time_samples.values())[0]))
        for q in self._time_samples:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            if max(self._time_samples[q]) >= qubit_info.reset.duration():
                qubit_info.reset.duration(qubit_info.reset.duration()+50e-6)
            if self.QD_agent.activeReset:
                eyeson_print("Active reset ON !")
            else:
                eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__repeat_data_idx = arange(self._repeat)

        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")

        self.QD_agent = check_acq_channels(self.QD_agent, list(self._time_samples.keys()))

        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        singleshot=self._os_mode,
        activeReset=self.QD_agent.activeReset
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=2*len(list(self._time_samples.keys())) if self.QD_agent.activeReset else len(list(self._time_samples.keys())),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            if not self._os_mode:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,self.__repeat_data_idx))
            else:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__one_shot_para,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,arange(self._avg_n),self.__repeat_data_idx))
        else:
            preview_para = {}
            for q in self._time_samples:
                preview_para[q] = array([self._time_samples[q][0],self._time_samples[q][-1]])
            self.__sched_kwargs['freeduration']= preview_para
    
    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run('T1')
            dict_ = {}
            if not self.QD_agent.activeReset:
                if not self._os_mode:
                    for q_idx, q in enumerate(self._time_samples):
                        i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._time_samples[q].shape[0])
                        q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._time_samples[q].shape[0])
                        dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
                        time_values = list(self._time_samples[q])*2*self._repeat
                        dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,self._repeat,self._time_samples[q].shape[0]))
                    
                    dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":self.__repeat_data_idx,"idx":self.__time_data_idx})
                else:
                    dict_ = {}
                    for q_idx, q in enumerate(self._time_samples):
                        i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._avg_n,self._time_samples[q].shape[0])
                        q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._avg_n,self._time_samples[q].shape[0])
                        dict_[q] = (["mixer","prepared_state","repeat","index","time_idx"],array([[i_data],[q_data]]))
                        time_values = list(self._time_samples[q])*2*self._repeat*self._avg_n
                        dict_[f"{q}_x"] = (["mixer","prepared_state","repeat","index","time_idx"],array(time_values).reshape(2,1,self._repeat,self._avg_n,self._time_samples[q].shape[0]))

                    dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":self.__repeat_data_idx,"prepared_state":array([1]),"index":arange(self._avg_n),"time_idx":self.__time_data_idx})
                    
                    
                dataset.attrs["execution_time"] = Data_manager().get_time_now()
                dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
                dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
                dataset.attrs["system"] = "qblox"
                self.dataset = dataset
            else:
                self.dataset = ds
                for attr in ds.attrs:
                    if isinstance(ds.attrs[attr], bool):
                        self.dataset.attrs[attr] = int(ds.attrs[attr])

                for var in ds.data_vars:
                    for attr in ds.data_vars[var].attrs:
                        if isinstance(ds.data_vars[var].attrs[attr], bool):
                            self.dataset.data_vars[var].attrs[attr] = int(ds.data_vars[var].attrs[attr])

                for var in ds.coords:
                    for attr in ds.coords[var].attrs:
                        if isinstance(ds.coords[var].attrs[attr], bool):
                            self.dataset.coords[var].attrs[attr] = int(ds.coords[var].attrs[attr])

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)



if __name__ == "__main__":
    x = False
    print(int(x))