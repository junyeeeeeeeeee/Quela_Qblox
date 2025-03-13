import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from numpy import array, moveaxis, arange, ndarray
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import Data_manager, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Readout, Multi_Readout, Integration, electrical_delay
from quantify_scheduler.operations.gate_library import Reset


class ReadoutFidelityPS( ScheduleConductor ):
    def __init__( self ):
        super().__init__()
        self._target_qs:list = []
        self._shots:int = 5001
        self._repeat:int = 1

    @property
    def target_qs( self ):
        return self._target_qs
    @target_qs.setter
    def set_target_qs( self, target_qs:list):
        if not isinstance(target_qs, list):
            raise TypeError("Target_qs must be given as a list !")
        else:
            self._target_qs = target_qs
    @property
    def shots( self ):
        return self._shots
    @shots.setter
    def set_shots( self, shot_num:int):
        self._shots = shot_num
    @property
    def repeat( self ):
        return self._repeat
    @repeat.setter
    def set_repeat(self, repeat_num:int):
        if not isinstance(repeat_num,(int,float)):
            raise TypeError("Arg 'repeat_num' must be a int or float !")
        self._repeat = int(repeat_num)

    def __PulseSchedule__(self,
        pi_amp: dict,
        pi_dura:dict,
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        states:ndarray = array([0, 1]),
        repetitions:int=1,
        ) -> Schedule:
        
        sameple_idx = states.shape[0]
        sched = Schedule("Single shot", repetitions=repetitions)
        for acq_idx in range(sameple_idx):
        
            for qubit_idx, q in enumerate(R_integration):

                sched.add(Reset(q))
                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=False)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration,powerDep=False)
            
                xyl = {q:pi_amp[q]*states[acq_idx]}
                self.QD_agent.Waveformer.X_pi_p(sched,xyl,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=True,get_trace=False,trace_recordlength=0)
        
        self.schedule = sched
        return sched
    
    def __SetParameters__(self, *args, **kwargs):
        for q in self._target_qs:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
            eyeson_print(f"{q} Integration time: {round(qubit_info.measure.integration_time()*1e6,1)} µs")
        
        self.__Para_state = ManualParameter(name="state")
        self.__Para_state.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__one_shot_para =  ManualParameter(name="Shot")

        self.__repeat_data_idx = arange(self._repeat)
        self.__state_idx = array([0, 1])
        
        
        self.__sched_kwargs = dict(   
            pi_amp=compose_para_for_multiplexing(self.QD_agent,self._target_qs,'d1'),
            pi_dura=compose_para_for_multiplexing(self.QD_agent,self._target_qs,'d3'),
            R_amp=compose_para_for_multiplexing(self.QD_agent,self._target_qs,'r1'),
            R_duration=compose_para_for_multiplexing(self.QD_agent,self._target_qs,'r3'),
            R_integration=compose_para_for_multiplexing(self.QD_agent,self._target_qs,'r4'),
            R_inte_delay=compose_para_for_multiplexing(self.QD_agent,self._target_qs,'r2'),
        )

    def __Compose__(self, *args, **kwargs):
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(self._target_qs),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._shots)
            self.meas_ctrl.gettables(self.__gettable)

            self.meas_ctrl.settables([self.__Para_state,self.__one_shot_para,self.__Para_repeat])
            self.meas_ctrl.setpoints_grid((self.__state_idx,arange(self._shots),self.__repeat_data_idx))
    
    def __RunAndGet__(self, *args, **kwargs):
        if self._execution:
            ds = self.meas_ctrl.run('Readout_Fidelity')
            dict_ = {}
            for q_idx, q in enumerate(self._target_qs):
                i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._shots,2)
                q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._shots,2)
                raw_data = moveaxis(moveaxis(array([i_data,q_data]),0,1),2,-1)  # (mixer, repeat, index, prepared_state) -> (repeat, mixer, prepared_state, index)
                print(raw_data.shape)
                dict_[q] = (["repeat","mixer","prepared_state","index"],raw_data)
        
            dataset = Dataset(dict_,coords={"repeat":self.__repeat_data_idx,"mixer":array(["I","Q"]),"prepared_state":array([0, 1]),"index":arange(self._shots)})
                
                
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["method"] = "Shot"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)
        

        
    
