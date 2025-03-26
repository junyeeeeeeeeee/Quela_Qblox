import time
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from numpy import array, moveaxis, arange, ndarray
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import Data_manager, check_acq_channels
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Rxy, IdlePulse, Measure, BinMode, Reset, pulse_preview, electrical_delay



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
        meas_qs: list,
        states:ndarray = array([0,1]),
        repetitions:int=1,
        activeReset:bool=False
        ) -> Schedule:
        
        sameple_idx = states.shape[0]
        sched = Schedule("Single shot", repetitions=repetitions)
        for acq_idx in range(sameple_idx):
            align_pulse = sched.add(IdlePulse(4e-9))
            for qubit_idx, q in enumerate(meas_qs):

                reset = sched.add(Reset(q), ref_op=align_pulse)
                
                if acq_idx == 1:
                    pi_pulse = sched.add(Rxy(qubit=q, theta=180, phi=0), ref_op=reset)
                
            sched.add(Measure(*meas_qs,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.APPEND), rel_time=electrical_delay, ref_op=pi_pulse if acq_idx==1 else reset)
            
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
        
        self.QD_agent = check_acq_channels(self.QD_agent, self._target_qs)
        
        self.__sched_kwargs = dict(   
            meas_qs=self._target_qs,
            states=self.__state_idx
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
                dict_[q] = (["repeat","mixer","prepared_state","index"],raw_data)
        
            dataset = Dataset(dict_,coords={"repeat":self.__repeat_data_idx,"mixer":array(["I","Q"]),"prepared_state":array([0, 1]),"index":arange(self._shots)})
                
                
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["method"] = "Shot"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)
        

        
    
