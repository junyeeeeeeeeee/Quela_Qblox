import time
from qcodes.parameters import ManualParameter
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from quantify_scheduler.gettables import ScheduleGettable
from numpy import arange, array, arange, rad2deg
from qblox_drive_AS.support import Data_manager, check_acq_channels
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, X90, Y90, Measure, IdlePulse, Reset, X, ConditionalReset, BinMode, Rxy, pulse_preview
from numpy import pi as PI

class RamseyT2PS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._time_samples:dict = {}
        self._use_arti_detune:bool = False
        self._repeat:int = 1
        self._spin_num:dict = {}
        self._second_phase:str = 'x'

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
    @property
    def arti_detune_status(self):
        return self._use_arti_detune
    @arti_detune_status.setter
    def enable_arti_detune( self, switch:bool):
        if not isinstance(switch,bool):
            if switch in [0,1]:
                pass
            else:
                raise TypeError("Arg 'os_mode' must be a bool, 0 or 1 !")
        
        self._use_arti_detune = switch

    @property
    def spin_num( self ):
        return self._spin_num
    @spin_num.setter
    def set_spin_num(self, pi_num:dict):
        self._spin_num = pi_num
    @property
    def second_phase( self ):
        return self._second_phase
    @second_phase.setter
    def set_second_phase(self, phase:str):
        if phase.lower() in ["x", "y"]:
            self._second_phase = phase.lower()
        else:
            print(f"Second phase should be 'x' or 'y', but '{phase}' was given. We choose to use 'x' instead.")
            self._second_phase = 'x'


    def __PulseSchedule__(self,
        freeduration:dict,
        echo_pi_num:dict,
        repetitions:int=1,
        second_pulse_phase:str='x',
        arti_detune:dict={},
        singleshot:bool=False,
        activeReset:bool=False
    ) -> Schedule:
        qubits2read = list(echo_pi_num.keys())
        sched = Schedule("Ramsey", repetitions=repetitions)
        
        for acq_idx in range(array(freeduration[qubits2read[0]]).shape[0]): 
            if not activeReset:
                reset = sched.add(Reset(*qubits2read))
            else:
                for idx, q in enumerate(qubits2read):
                    if idx == 0:
                        reset = sched.add(
                            ConditionalReset(q, acq_index=acq_idx),
                            label=f"Reset {acq_idx}",
                        )
                    else:
                        sched.add(
                            ConditionalReset(q, acq_index=acq_idx),
                            label=f"Reset {acq_idx}",
                            ref_op=reset,
                            ref_pt="start"
                        )
            
            for qubit_idx, q in enumerate(qubits2read):
                
                freeDu = freeduration[q][acq_idx]
                
                # first pi/2
                first_pulse = sched.add(X90(q), ref_op=reset)
                
                
                if echo_pi_num[q] != 0:
                    a_separate_free_Du = freeDu / echo_pi_num[q]
                    # middle pulses
                    for pi_idx in range(echo_pi_num[q]):
                        if pi_idx == 0 :
                            pi = sched.add(X(q), ref_op=first_pulse, rel_time=0.5*a_separate_free_Du)
                        else:
                            pi = sched.add(X(q), ref_op=pi, rel_time=1*a_separate_free_Du)
                    
                    # The last pulse
                    sched.add(X90(q) if second_pulse_phase.lower()=='x' else Y90(q), ref_op=pi, rel_time=0.5*a_separate_free_Du)  
                else:
                    if q in arti_detune:
                        recovery_phase = rad2deg(2 * PI * arti_detune[q] * freeDu)
                        sched.add(
                            Rxy(theta=90, phi=recovery_phase, qubit=q), ref_op=first_pulse, rel_time=freeDu
                        )
                    else:
                        sched.add(X90(q) if second_pulse_phase.lower()=='x' else Y90(q), ref_op=first_pulse, rel_time=freeDu)  
                
            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex' if not activeReset else 'ThresholdedAcquisition', bin_mode=BinMode.APPEND if singleshot else BinMode.AVERAGE))
        
        return sched
        

    def __SetParameters__(self, *args, **kwargs):
        
        self.__time_data_idx = arange(len(list(self._time_samples.values())[0]))
        for q in self._time_samples:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            if q not in list(self._spin_num.keys()):
                self._spin_num[q] = 0
            qubit_info.reset.duration(qubit_info.reset.duration()+max(self._time_samples[q]))
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__repeat_data_idx = arange(self._repeat)

        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")

        self.QD_agent = check_acq_channels(self.QD_agent, list(self._time_samples.keys()))
        arti_detunes = {}
        if self._use_arti_detune:
            for q in list(self._time_samples.keys()):
                arti_detunes[q] = self.QD_agent.Notewriter.get_artiT2DetuneFor(q)

        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        singleshot=self._os_mode,
        echo_pi_num=self._spin_num,
        second_pulse_phase=self._second_phase,
        activeReset=self.QD_agent.activeReset,
        arti_detune=arti_detunes
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
            ds = self.meas_ctrl.run('RamseyT2')
            dict_ = {}
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
                
                
            for var in [i for i in dataset.data_vars if "_" not in i]:
                dataset[var].attrs["spin_num"] = self._spin_num[var]
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)










        
    
    
        
            
        


    