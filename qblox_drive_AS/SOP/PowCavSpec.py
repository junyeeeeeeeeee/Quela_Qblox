""" 
Base on a BARE cavity observe a dispersive shift in RO-freq with the variable RO-amp.  
"""
import os
from numpy import array, sqrt, mean, arange, ndarray
from xarray import Dataset
from numpy import nan as NaN
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager
from quantify_scheduler.gettables import ScheduleGettable
import matplotlib.pyplot as plt
from qblox_drive_AS.support import check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import Measure, Schedule, pulse_preview, BinMode
from qblox_drive_AS.support.Pulser import ScheduleConductor
from quantify_scheduler.operations.gate_library import Reset
from quantify_scheduler.operations.pulse_library import IdlePulse,SetClockFrequency
from quantify_scheduler.resources import ClockResource



def plot_powerCavity_S21(ds:Dataset, QD_agent:QDmanager=None, save_fig_folder:str=None):
    """
    Plot |S21| from a given power cavity nc file and save it in the pic folder within the same day.
    """
    power = array(ds.coords['ro_amp'])
    # ds.x0 = freq. ; ds.x1 = power
    for idx, q in enumerate(ds.data_vars):
        if str(q).split("_")[-1] != 'freq':
            amp:ndarray = sqrt(array(ds[q])[0]**2+array(ds[q])[1]**2)
            s21 = []
            for i in range(amp.shape[0]):
                if power[i] != 0:
                    s21.append(list(array(amp[i])/power[i]))
                else:
                    s21.append(list(amp[i]))
            s21 = array(s21)
            freq = array(ds[f"{q}_freq"])[0][0]

            fig, ax = plt.subplots()
            ax:plt.Axes
            d = ax.pcolormesh(freq*1e-9, power, s21, shading='gouraud',cmap='RdBu')
           
            ax.vlines(mean(freq)*1e-9 if QD_agent is None else QD_agent.quantum_device.get_element(q).clock_freqs.readout()*1e-9,ymin=min(power),ymax=max(power),linestyles='--',colors='#FF00FF',label='window_shift_baseline')
            fig.colorbar(d, ax=ax)
            plt.xlabel("frequency (GHz)")
            plt.ylabel("Power (V)")
            plt.minorticks_on()
            if QD_agent is not None:
                ro_atte:str = f"_{int(QD_agent.Notewriter.get_DigiAtteFor(q,'ro'))}dB"
            else:
                ro_atte:str = "" 

            plt.title(f"PowerCavity_{q}{ro_atte}")
            plt.grid()
            plt.tight_layout()
            plt.legend()

            if save_fig_folder is not None:
                plt.savefig(os.path.join(save_fig_folder,f"PowerCavity_{q}{ro_atte}.png"))
                plt.close()
            else:
                plt.show()


class PowerDepCavityPS(ScheduleConductor):
    def __init__(self):
        super().__init__() 
        self._ro_elements:dict = {}
        self._power_samples:ndarray = []
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles
    
    @property
    def power_samples(self):
        return self._power_samples
    @power_samples.setter
    def power_samples(self,powers:ndarray):
        self._power_samples = powers


    def __PulseSchedule__(self, 
        frequencies: dict,
        R_amp: any,
        repetitions:int=1,    
    ) -> Schedule:
        
        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
        sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)


        for acq_idx in range(sameple_idx):    

            for qubit_idx, q in enumerate(qubits2read):
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
                
                sched.add(Reset(q))
                sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
                sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")
            
            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE, pulse_amp=R_amp))
                
        self.schedule =  sched  
        return sched
    

    def __SetParameters__(self, *args, **kwargs):
         
        self.__freq_datapoint_idx = arange(0,len(list(list(self._ro_elements.values())[0])))
        original_rof = {}
        quantum_device = self.QD_agent.quantum_device
        for q in self._ro_elements:
            original_rof[q] = quantum_device.get_element(q).clock_freqs.readout()
            quantum_device.get_element(q).clock_freqs.readout(NaN) # avoid cluster clock warning
            
        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))
        self.__freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
        self.__freq.batched = True
        self.__ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
        self.__ro_pulse_amp.batched = False

        self.__spec_sched_kwargs = dict(   
        frequencies=self._ro_elements,
        R_amp=self.__ro_pulse_amp,
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
            self.meas_ctrl.settables([self.__freq, self.__ro_pulse_amp])
            self.meas_ctrl.setpoints_grid([self.__freq_datapoint_idx, self._power_samples])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
             
            self.__spec_sched_kwargs['R_amp']= self._power_samples[-1]
            self.__spec_sched_kwargs['frequencies']= preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rp_ds = self.meas_ctrl.run("PowerCavity")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                i_data = array(rp_ds[f'y{2*q_idx}']).reshape(self._power_samples.shape[0],self._ro_elements[q].shape[0])
                q_data = array(rp_ds[f'y{2*q_idx+1}']).reshape(self._power_samples.shape[0],self._ro_elements[q].shape[0])
                dict_[q] = (["mixer","ro_amp","freq"],array([i_data,q_data]))
                dict_[f'{q}_freq'] = (["mixer","ro_amp","freq"],array(list(self._ro_elements[q])*2*self._power_samples.shape[0]).reshape(2,self._power_samples.shape[0],self._ro_elements[q].shape[0]))
            
            dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"ro_amp":self._power_samples,"freq":self.__freq_datapoint_idx})
            dataset.attrs["method"] = "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)


