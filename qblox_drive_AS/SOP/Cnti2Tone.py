
from numpy import nan as NaN
from numpy import array, arange, ndarray
from xarray import Dataset
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support.QuFluxFit import calc_Gcoef_inFbFqFd, calc_g
from qblox_drive_AS.support import check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import Measure, Schedule, pulse_preview, BinMode

from qblox_drive_AS.support.Pulser import ScheduleConductor
from quantify_scheduler.operations.gate_library import Reset
from quantify_scheduler.operations.pulse_library import SetClockFrequency,SquarePulse, IdlePulse
from quantify_scheduler.resources import ClockResource

def update_2toneResults_for(QD_agent:QDmanager,qb:str,QS_results:dict,XYL:float):
    qubit = QD_agent.quantum_device.get_element(qb)
    Revised_f01 = QS_results[qb].attrs['f01_fit']
    fb = float(QD_agent.Notewriter.get_bareFreqFor(target_q=qb))*1e-6
    fd = QD_agent.quantum_device.get_element(qb).clock_freqs.readout()*1e-6
    A = calc_Gcoef_inFbFqFd(fb,Revised_f01*1e-6,fd)
    sweet_g = calc_g(fb,Revised_f01*1e-6,A)
    qubit.clock_freqs.f01(Revised_f01)
    QD_agent.Notewriter.save_2tone_piamp_for(qb,XYL)
    QD_agent.Notewriter.save_CoefInG_for(target_q=qb,A=A)
    QD_agent.Notewriter.save_sweetG_for(target_q=qb,g_Hz=sweet_g*1e6)

def tune_away_setup(QD_agent:QDmanager,Fctrl:dict,bias_setup:dict,ro_qs:list,zero_mode:bool=False):
    for q in ro_qs:
        if zero_mode:
            Fctrl[q](0.0)
        else:
            offset = QD_agent.Fluxmanager.get_proper_zbiasFor(q)
            if q not in list(bias_setup.keys()):
                bias_setup[q] = 0
            want_bias = offset+bias_setup[q]
            if bias_setup[q] != 0:
                rof = QD_agent.Fluxmanager.sin_for_cav(q,array([want_bias]))[0]
                QD_agent.quantum_device.get_element(q).clock_freqs.readout(rof)
                QD_agent.Fluxmanager.save_tuneawayBias_for('manual',q,want_bias)
                warning_print(f"meas under flux = {round(offset,3)}+{round(bias_setup[q],3)} V")
            Fctrl[q](want_bias) 



class PowerDepQubitPS(ScheduleConductor):
    def __init__(self):
        super().__init__() 
        self._ro_elements:dict = {}
        self._power_samples:ndarray = []
        self._ROoverlapDriving:bool=False
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles

    @property
    def overlap(self):
        return self._ROoverlapDriving
    @overlap.setter
    def overlap(self, overlapping:bool|int):
        if overlapping in [0, 1]:
            self._ROoverlapDriving = overlapping
        elif isinstance(overlapping, bool):
            self._ROoverlapDriving = overlapping
        else:
            raise TypeError("Overlap setting must be a bool, 0 or 1 !")
    
    @property
    def power_samples(self):
        return self._power_samples
    @power_samples.setter
    def power_samples(self,powers:ndarray):
        self._power_samples = powers


    def __PulseSchedule__(self, 
        frequencies: dict,
        mixdrive_amp: any,
        mixdrive_dura:float,
        overlapping:bool=False,
        repetitions:int=1,    
    ) -> Schedule:
        
        sched = Schedule("Two tone spectroscopy (NCO sweep)",repetitions=repetitions)
        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]

        for acq_idx in range(sameple_idx):   
            align_pulse = sched.add(IdlePulse(4e-9)) 
            for q in qubits2read:
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+ ".01", freq=array(frequencies[q]).flat[0]))
            
                sched.add(SetClockFrequency(clock= q+".01", clock_freq_new=freq))
                sched.add(Reset(q))
                sched.add(SquarePulse(amp=mixdrive_amp, duration=mixdrive_dura, port=f'{q}:mw', clock=f'{q}.01'), ref_op=align_pulse)

            
            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE), ref_pt="start" if overlapping else "end", rel_time=280e-9 if not overlapping else 0)
                   
        return sched
    

    def __SetParameters__(self, *args, **kwargs):
        self.__drive_pulse_length = 100e-6
        original_qubit_info = {}
        for q in self._ro_elements:
            original_qubit_info[q] = {}
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            original_qubit_info[q]["ROW"] = qubit_info.measure.pulse_duration()
            original_qubit_info[q]["ITW"] = qubit_info.measure.integration_time()
            original_qubit_info[q]["XYF"] = qubit_info.clock_freqs.f01()
            self.__freq_datapoint_idx = arange(0,self._ro_elements[q].shape[0])
                
            if self._ROoverlapDriving:
                qubit_info.measure.pulse_duration(self.__drive_pulse_length)
                qubit_info.measure.integration_time(self.__drive_pulse_length)

            
            qubit_info.reset.duration(4e-9) 
            qubit_info.clock_freqs.f01(NaN)
            eyeson_print(f"Inte_time= {round(qubit_info.measure.integration_time()*1e6,1)} µs")
            eyeson_print(f"Reset_time= {round(qubit_info.reset.duration()*1e6,1)} µs")

        self.__freq = ManualParameter(name="XYfreq", unit="Hz", label="Frequency")
        self.__freq.batched = True
        self.__xy_amp = ManualParameter(name="XYamp", unit="V", label="Amplitude")
        self.__xy_amp.batched = False

        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))

        self.__spec_sched_kwargs = dict(   
            frequencies=self.ro_elements,
            mixdrive_amp=self.__xy_amp,
            mixdrive_dura=self.__drive_pulse_length,
            overlapping=self._ROoverlapDriving
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
            self.meas_ctrl.settables([self.__freq, self.__xy_amp])
            self.meas_ctrl.setpoints_grid([self.__freq_datapoint_idx, self._power_samples])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
             
            self.__spec_sched_kwargs['mixdrive_amp']= self._power_samples[-1]
            self.__spec_sched_kwargs['frequencies']= preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run("Two-tone")
        
            dict_ = {}
            for idx, q in enumerate(self._ro_elements):
                freq_values = 2*self._power_samples.shape[0]*list(self._ro_elements[q])
                i_data = array(ds[f'y{2*idx}']).reshape(self._power_samples.shape[0],array(self._ro_elements[q]).shape[0])
                q_data = array(ds[f'y{2*idx+1}']).reshape(self._power_samples.shape[0],array(self._ro_elements[q]).shape[0])
                dict_[q] = (["mixer","xy_amp","freq"],array([i_data,q_data]))
                dict_[f"{q}_freq"] = (["mixer","xy_amp","freq"],array(freq_values).reshape(2,self._power_samples.shape[0],array(self._ro_elements[q]).shape[0]))

            
            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"xy_amp":self._power_samples,"freq":self.__freq_datapoint_idx})
            ds.attrs["execution_time"] = Data_manager().get_time_now()
            ds.attrs["method"] = "Average"
            ds.attrs["system"] = "qblox" 
            self.dataset = ds
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)


