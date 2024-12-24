import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from numpy import array, arange
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_T1_sche, pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Readout, Multi_Readout, Integration, electrical_delay
from quantify_scheduler.operations.gate_library import Reset


def T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,time_samples:dict,repeat:int=1,n_avg:int=300,run:bool=True):
    sche_func= multi_T1_sche

    for q in time_samples:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        time_data_idx = arange(time_samples[q].shape[0])
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
    Para_repeat.batched = False
    repeat_data_idx = arange(repeat)

    sched_kwargs = dict(
        freeduration=time_samples,
        pi_amp=compose_para_for_multiplexing(QD_agent,time_samples,'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,time_samples,'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,time_samples,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,time_samples,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,time_samples,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,time_samples,'r2'),
        )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(time_samples.keys())),
            )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([Para_free_Du,Para_repeat])
        meas_ctrl.setpoints_grid((time_data_idx,repeat_data_idx))

        ds = meas_ctrl.run('T1')
        dict_ = {}
        for q_idx, q in enumerate(time_samples):
            i_data = array(ds[f'y{2*q_idx}']).reshape(repeat,time_samples[q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(repeat,time_samples[q].shape[0])
            dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
            time_values = list(time_samples[q])*2*repeat
            dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,repeat,time_samples[q].shape[0]))
        
        T1_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":repeat_data_idx,"idx":time_data_idx})
        T1_ds.attrs["execution_time"] = Data_manager().get_time_now()
        T1_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

    else:
        preview_para = {}
        for q in time_samples:
            preview_para[q] = array([time_samples[q][0],time_samples[q][-1]])
        sched_kwargs['freeduration']= preview_para
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        T1_ds = ''
    
    return T1_ds


class EnergyRelaxPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._time_samples:dict = {}
        self._repeat:int = 1
        self._avg_n:int = 300


    def __PulseSchedule__(self, 
        freeduration:dict,
        pi_amp: dict,
        pi_dura:dict,
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        repetitions:int=1,
        )->Schedule:

        qubits2read = list(freeduration.keys())
        sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
        sched = Schedule("T1", repetitions=repetitions)
        for acq_idx in range(sameple_idx):
            for qubit_idx, q in enumerate(qubits2read):
                freeDu = freeduration[q][acq_idx]
                sched.add(Reset(q))
            
                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
            
            
                self.QD_agent.Waveformer.X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu+electrical_delay)
                
                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
        self.schedule = sched
        return sched

    def __SetParameters__(self, *args, **kwargs):
        self.__time_data_idx = arange(array(self._time_samples.values())[0].shape[0])
        for q in self._time_samples:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__repeat_data_idx = arange(self._repeat)

        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        pi_amp=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'d1'),
        pi_dura=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'d3'),
        R_amp=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r1'),
        R_duration=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r3'),
        R_integration=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r4'),
        R_inte_delay=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r2'),
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(self._time_samples.keys())),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables([self.__Para_free_Du,self.__Para_repeat])
            self.meas_ctrl.setpoints_grid((self.__time_data_idx,self.__repeat_data_idx))
        else:
            preview_para = {}
            for q in self._time_samples:
                preview_para[q] = array([self._time_samples[q][0],self._time_samples[q][-1]])
            self.__sched_kwargs['freeduration']= preview_para
    
    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run('T1')
            dict_ = {}
            for q_idx, q in enumerate(self._time_samples):
                i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._time_samples[q].shape[0])
                q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._time_samples[q].shape[0])
                dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
                time_values = list(self._time_samples[q])*2*self._repeat
                dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,self._repeat,self._time_samples[q].shape[0]))
            
            dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":self.__repeat_data_idx,"idx":self.__time_data_idx})
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)


if __name__ == "__main__":
    t1 = EnergyRelaxPS()
    ps = t1.get_adjsutable_paras(display=True)