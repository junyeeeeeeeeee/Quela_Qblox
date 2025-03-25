
from numpy import nan as NaN
from numpy import ndarray
from xarray import Dataset
from numpy import array, arange, sqrt
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, check_acq_channels
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support.Pulse_schedule_library import Measure, Schedule, BinMode, pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from quantify_scheduler.operations.gate_library import Reset 
from quantify_scheduler.operations.pulse_library import IdlePulse,SetClockFrequency, SquarePulse
from quantify_scheduler.resources import ClockResource


z_pulse_amp_OVER_const_z = sqrt(2)/2.5

def update_by_fluxQubit(QD_agent:QDmanager,correct_results:dict,target_q:str):
    """
    correct_results dict in the form: {"xyf":float,"sweet_bias":float}
    """
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.clock_freqs.f01(correct_results["xyf"])
    QD_agent.Fluxmanager.check_offset_and_correctFor(target_q=target_q,new_offset=correct_results["sweet_bias"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=target_q,bias=correct_results["sweet_bias"])

class FluxDepQubitPS(ScheduleConductor):
    def __init__(self):
        super().__init__() 
        self._ro_elements:dict = {}
        self._avg_n:int = 100
        self._flux_samples:ndarray = []
        self._bias_elements:list = []
        self._ROoverlapDriving:bool=False
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles
    @property
    def bias_elements(self):
        return self._bias_elements
    @bias_elements.setter
    def bias_elements(self, bias_eles:dict):
        self._bias_elements = bias_eles
    
    @property
    def flux_samples(self):
        return self._flux_samples
    @flux_samples.setter
    def flux_samples(self,fluxs:ndarray):
        self._flux_samples = fluxs

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
        bias_qs:list,
        Z_amp:any,
        spec_amp:float,
        spec_Du:float,
        repetitions:int=1,     
    ) -> Schedule:
        
        sched = Schedule("Zgate_two_tone spectroscopy (NCO sweep)",repetitions=repetitions)

        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]

        for acq_idx in range(sameple_idx):    
            align_pulse = sched.add(IdlePulse(4e-9))
            for qubit_idx, q in enumerate(qubits2read):
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+".01", freq=array(frequencies[q]).flat[0]))
    
                sched.add(SetClockFrequency(clock= q+ ".01", clock_freq_new=freq))
                sched.add(IdlePulse(4e-9))
                reset = sched.add(Reset(q), ref_op=align_pulse)
                
                if qubit_idx == 0:
                    for qb in bias_qs:
                        sched.add(SquarePulse(amp=Z_amp, duration=spec_Du, port=qb+":fl", clock="cl0.baseband"), ref_op=reset)
                
                sched.add(SquarePulse(amp=spec_amp,duration=spec_Du,port=f"{q}:mw",clock=f"{q}.01"), ref_op=reset)
                

            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE),rel_time=280e-9)
                
     
        return sched
    

    def __SetParameters__(self, *args, **kwargs):
        self.__drive_pulse_length = 10e-6
        original_xyfs = {}
        for q in self._ro_elements:
            self.__freq_datapoint_idx = arange(0,self._ro_elements[q].shape[0])
            qubit_info = self.QD_agent.quantum_device.get_element(q)  
            original_xyfs[q] = qubit_info.clock_freqs.f01()
            qubit_info.clock_freqs.f01(NaN)
            self.__spec_pulse_amp = self.QD_agent.Notewriter.get_2tone_piampFor(q)
    

        self.__freq = ManualParameter(name="XYfreq", unit="Hz", label="Frequency")
        self.__freq.batched = True
        self.__flux = ManualParameter(name="flux", unit="V", label="flux")
        self.__flux.batched = False

        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))

        self.__spec_sched_kwargs = dict(   
            frequencies=self.ro_elements,
            bias_qs=self._bias_elements,
            Z_amp=self.__flux,
            spec_amp=self.__spec_pulse_amp,
            spec_Du=self.__drive_pulse_length,
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
            self.meas_ctrl.settables([self.__freq, self.__flux])
            self.meas_ctrl.setpoints_grid([self.__freq_datapoint_idx, self._flux_samples*z_pulse_amp_OVER_const_z])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
             
            self.__spec_sched_kwargs['Z_amp']= self._flux_samples[-1]
            self.__spec_sched_kwargs['frequencies']= preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            qs_ds = self.meas_ctrl.run("Zgate-two-tone")
        
            dict_ = {}
            for idx, q in enumerate(self._ro_elements):
                freq_values = 2*self._flux_samples.shape[0]*list(self._ro_elements[q])
                i_data = array(qs_ds[f'y{2*idx}']).reshape(self._flux_samples.shape[0],array(self._ro_elements[q]).shape[0])
                q_data = array(qs_ds[f'y{2*idx+1}']).reshape(self._flux_samples.shape[0],array(self._ro_elements[q]).shape[0])
                dict_[q] = (["mixer","bias","freq"],array([i_data,q_data]))
        
                dict_[f"{q}_freq"] = (["mixer","bias","freq"],array(freq_values).reshape(2,self._flux_samples.shape[0],array(self._ro_elements[q]).shape[0]))

            rfs_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":self._flux_samples,"freq":self.__freq_datapoint_idx})
            rfs_ds.attrs["execution_time"] = Data_manager().get_time_now()
            rfs_ds.attrs["method"] = "Average"
            rfs_ds.attrs["system"] = "qblox"
            self.dataset = rfs_ds
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)


