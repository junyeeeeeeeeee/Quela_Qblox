"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from numpy import array, arange, ndarray, round, full, concatenate
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import compose_para_for_multiplexing
from xarray import Dataset
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Readout, Multi_Readout, Integration, electrical_delay
from quantify_scheduler.operations.gate_library import Reset

#? The way to merge two dict a and b : c = {**a,**b}

def sort_elements_2_multiples_of(x:ndarray, specific_multiple:int=None):
    """
    Adjusts all numbers in a 1D NumPy array to the nearest multiples of 4
    and ensures the array has exactly `target_size` elements.
    
    Parameters:
        array (numpy.ndarray): The input 1D array of numbers.
        target_size (int): The desired number of elements in the adjusted array.
    
    Returns:
        numpy.ndarray: The adjusted array with all numbers as multiples of 4 and
                       exactly `target_size` elements.
    """
    if specific_multiple is None: specific_multiple = 4
    if not isinstance(x, ndarray):
        raise ValueError("Input must be a NumPy array.")
    
    target_size = x.shape[0]

    # Step 1: Adjust numbers to the nearest multiples of 4
    adjusted_array = round(x / specific_multiple) * specific_multiple
    
    # Step 2: Adjust the array size to match the target size
    current_size = adjusted_array.size
    if current_size < target_size:
        # Extend the array by repeating the largest multiple of 4
        max_value = max(adjusted_array)
        additional_values = full(target_size - current_size, max_value)
        adjusted_array = concatenate([adjusted_array, additional_values])
    elif current_size > target_size:
        # Trim the array to the target size
        adjusted_array = adjusted_array[:target_size]
    
    return adjusted_array.astype(int)


def conditional_update_qubitInfo(QD_agent:QDmanager,fit_results:Dataset,target_q:str):
    if fit_results.attrs['pi_2'] >= min(fit_results.coords['samples']) and fit_results.attrs['pi_2'] <= max(fit_results.coords['samples']) :
        qubit = QD_agent.quantum_device.get_element(target_q)
        match str(fit_results.attrs['Rabi_type']).lower():
            case 'powerrabi':
                qubit.rxy.amp180(fit_results.attrs['pi_2'])
    
            case 'timerabi':
                qubit.rxy.duration(fit_results.attrs['pi_2'])
    else:
        warning_print(f"Results for {target_q} didn't satisfy the update condition !")


class RabiPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._RabiType:str = "power"
        self._variables:dict = {}
        self._os_mode:bool = False 
        self._avg_n:int = 300

    @property
    def RabiType( self ):
        return self._RabiType
    @RabiType.setter
    def set_RabiType( self, type:str):
        if type.lower() in ["power", 'time']:
            self._RabiType = type.lower()
        else:
            raise ValueError("Arg 'type' must be given as 'time' or 'power' !")

    @property
    def samples( self ):
        return self._variables
    @samples.setter
    def set_samples( self, samples:dict):
        if not isinstance(samples, dict):
            raise TypeError("Arg 'samples' must be a dict !")
        self._variables = samples
    @property
    def os_mode( self ):
        return self._os_mode
    @os_mode.setter
    def set_os_mode( self, mode:bool):
        if not isinstance(mode,bool):
            if mode in [0,1]:
                pass
            else:
                raise TypeError("Arg 'mode' must be a bool, 0 or 1 !")
        
        self._os_mode = mode
    @property
    def n_avg( self ):
        return self._avg_n
    @n_avg.setter
    def set_n_avg(self, avg_num:int):
        if not isinstance(avg_num,(int,float)):
            raise TypeError("Ard 'avg_num' must be a int or float !")
        self._avg_n = int(avg_num)

    def __PulseSchedule__(self,
        pi_amp:dict,
        pi_dura:dict,
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        XY_theta:str,
        repetitions:int=1,
        OS_or_not:bool=False
    ) -> Schedule:

        qubits2read = list(pi_amp.keys()) if self._RabiType.lower() == 'power' else list(pi_dura.keys())
        sample_len = pi_amp[qubits2read[0]].shape[0] if self._RabiType.lower() == 'power' else pi_dura[qubits2read[0]].shape[0]
        sched = Schedule("RabiOscillation",repetitions=repetitions)

        match XY_theta:
            case 'Y_theta':
                gate:callable = self.QD_agent.Waveformer.Y_pi_p
            case _:
                gate:callable = self.QD_agent.Waveformer.X_pi_p

        
        for acq_idx in range(sample_len):    
            for qubit_idx, q in enumerate(qubits2read):
                sched.add(Reset(q))
                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
                
                if self._RabiType.lower() == 'power':
                    gate(sched,{q:pi_amp[q][acq_idx]},q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
                else:
                    gate(sched,{q:pi_amp[q]},q,pi_dura[q][acq_idx],spec_pulse,freeDu=electrical_delay)
            
            
                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=OS_or_not,get_trace=False,trace_recordlength=0)
        self.schedule = sched
        return sched

    def __SetParameters__(self, *args, **kwargs):
        self.amps, self.duras = {}, {}
        for q in self._variables:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
            eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
            self.__variable_idx = arange(0,self._variables[q].shape[0])

        
        if self._RabiType == 'power':
            self.__Sweep_para = ManualParameter(name="PowerRabi", unit="V", label="amplitude")
            self.amps = self._variables
            for q in self._variables:
                self.duras[q] = self.QD_agent.quantum_device.get_element(q).rxy.duration()
        else:
            self.__Sweep_para = ManualParameter(name="TimeRabi", unit="sec", label="time")
            self.duras = self._variables
            for q in self._variables:
                self.amps[q] = self.QD_agent.quantum_device.get_element(q).rxy.amp180()
       
        self.__Sweep_para.batched = True
        
        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")
        
        self.__sched_kwargs = dict(
            pi_amp = self.amps,
            pi_dura = self.duras,
            R_amp=compose_para_for_multiplexing(self.QD_agent,self._variables,'r1'),
            R_duration=compose_para_for_multiplexing(self.QD_agent,self._variables,'r3'),
            R_integration=compose_para_for_multiplexing(self.QD_agent,self._variables,'r4'),
            R_inte_delay=compose_para_for_multiplexing(self.QD_agent,self._variables,'r2'),
            XY_theta='X_theta',
            OS_or_not=self._os_mode
            )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(self._variables.keys())),
            )
        
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)

            if not self._os_mode:
                self.meas_ctrl.settables(self.__Sweep_para)
                self.meas_ctrl.setpoints(self.__variable_idx)
            else:
                self.meas_ctrl.settables([self.__Sweep_para, self.__one_shot_para])
                self.meas_ctrl.setpoints_grid((self.__variable_idx, arange(self._avg_n)))
        
        else:
            preview_para = {}
            for q in self._variables:
                preview_para[q] = array([self._variables[q][0],self._variables[q][-1]])
            
            if self._RabiType == "power":
                self.__sched_kwargs['pi_amp'] = preview_para
            else:
                self.__sched_kwargs['pi_dura'] = preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run("RabiOscillation")
            dict_ = {}

            if not self._os_mode:
                for idx, q in enumerate(list(self._variables.keys())):
                    I = array(ds[f'y{2*idx}'])
                    Q = array(ds[f'y{2*idx+1}'])
                    dict_[q] = (['mixer','var_idx'],array([I,Q]))
                    dict_[f"{q}_variable"] = (['mixer','var_idx'],array(2*list(self._variables[q])).reshape(2,self._variables[q].shape[0]))
                
                dataset = Dataset(dict_, coords={"mixer":array(["I","Q"]),"var_idx":self.__variable_idx})
                
            else:
                for idx, q in enumerate(list(self._variables.keys())):
                    I = array(ds[f'y{2*idx}']).reshape(self._avg_n,self._variables[q].shape[0])
                    Q = array(ds[f'y{2*idx+1}']).reshape(self._avg_n,self._variables[q].shape[0])
                    dict_[q] = (["mixer","prepared_state","index","var_idx"],array([[I],[Q]]))
                    variables =  list(self._variables[q])*2*self._avg_n
                    dict_[f"{q}_variable"] = (["mixer","prepared_state","index","var_idx"],array(variables).reshape(2,1,self._avg_n,self._variables[q].shape[0]))
                
                dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"prepared_state":array([1]),"index":arange(self._avg_n),"var_idx":self.__variable_idx})
                    
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["rabi_type"] = self._RabiType
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"            
            if self._RabiType == 'time':
                for q in self.amps:
                    dataset.attrs[f"{q}_piamp"] = self.amps[q]
            else:
                for q in self.duras:
                    dataset.attrs[f"{q}_pidura"] = self.duras[q]
            
            self.dataset = dataset