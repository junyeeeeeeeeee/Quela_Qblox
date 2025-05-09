""" 
Base on a BARE cavity observe a dispersive shift in RO-freq with the variable RO-amp.  
"""
import os
from numpy import array, sqrt, mean, arange, ndarray, moveaxis, argmin, average
from xarray import Dataset
from numpy import nan as NaN
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager
from quantify_scheduler.gettables import ScheduleGettable
import matplotlib.pyplot as plt
from qblox_drive_AS.support.UserFriend import slightly_print
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
    rof = {}
    # ds.x0 = freq. ; ds.x1 = power
    for idx, q in enumerate(ds.data_vars):
        if str(q).split("_")[-1] != 'freq':
            freq = array(ds[f"{q}_freq"])[0][0]
            amp:ndarray = sqrt(array(ds[q])[0]**2+array(ds[q])[1]**2)
            s21 = []
            deep = []
            for i in range(amp.shape[0]):
                norm_s21 = amp[i] / max(abs(amp[i]))
                s21.append(norm_s21)
                deep.append(freq[argmin(norm_s21)])
            s21 = array(s21)
            

            fig, ax = plt.subplots()
            ax:plt.Axes
            d = ax.pcolormesh(freq*1e-9, power, s21, shading='gouraud',cmap='RdBu')
            threshold = mean(array(deep))
            dress = []
            bare = []
            if QD_agent.Notewriter.get_bareFreqFor(q) > threshold: # bare > dress
                for jdx, fr in enumerate(deep):
                    if fr > threshold:
                        bare.append([fr, power[jdx]])
                    else:
                        dress.append([fr, power[jdx]])
            else:                                                  # bare <= dress
                for jdx, fr in enumerate(deep):
                    if fr <= threshold:
                        bare.append([fr, power[jdx]])
                    else:
                        dress.append([fr, power[jdx]])
            
            
            ax.vlines(average(array(bare)[:,0]*1e-9), min(array(bare)[:,1]),  min(array(dress)[:,1]), colors='green', linestyles="--", label='Bare')
            ax.vlines(average(array(dress)[:,0]*1e-9), min(array(dress)[:,1]), max(array(dress)[:,1]), colors='yellow', linestyles="--", label='Dress')
            rof[q] = average(array(dress)[:,0])
            fig.colorbar(d, ax=ax)
            plt.xlabel("frequency (GHz)")
            plt.ylabel("Atte. (dB)")
            plt.minorticks_on()

            plt.title(f"PowerCavity_{q}")
            plt.grid()
            plt.tight_layout()
            plt.legend()

            if save_fig_folder is not None:
                pt = os.path.join(save_fig_folder,f"PowerCavity_{q}.png")
                plt.savefig(pt)
                slightly_print(f"pic saved at {pt}")
                plt.close()

            else:
                plt.show()

    return rof    


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
                
                # sched.add(Reset(q))
                sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
                sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")

            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE) )
                
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

        self.__spec_sched_kwargs = dict(   
        frequencies=self._ro_elements,
        )
    
    def __Compose__(self, ro_atte:int):
        
        if self._execution:
            print(f"ro_atte: {ro_atte} dB")
            if ro_atte%2 != 0:
                raise ValueError("ro_atte must be even")
            
            for q in self.QD_agent.quantum_device.elements():
                self.QD_agent.quantum_device.hardware_config()["hardware_options"]["output_att"][f'{q}:res-{q}.ro'] = ro_atte
            
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
            self.meas_ctrl.settables(self.__freq)
            self.meas_ctrl.setpoints(self.__freq_datapoint_idx)
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
             
            self.__spec_sched_kwargs['frequencies']= preview_para

    def __RunAndGet__(self, *args, **kwargs):
        if self._execution:
            s21_folder = {}
            print(self._power_samples)
            for atte in self._power_samples:
                self.__Compose__(atte)
                rs_ds = self.meas_ctrl.run("1D-cavity-sweep")
            
                for q_idx, q in enumerate(list(self._ro_elements.keys())):
                    if q not in s21_folder:
                        s21_folder[q] = []
                    i_data = array(rs_ds[f'y{2*q_idx}'])
                    q_data = array(rs_ds[f'y{2*q_idx+1}'])
                    s21_folder[q].append([i_data,q_data])
            
            dict_ = {}
            for  q in s21_folder: # shape of s21_folder[q] is (self._power_samples.shape[0],2,self._ro_elements[q].shape[0]) right is (atte, mixer, freq)
                
                dict_[q] = (["mixer","ro_amp","freq"], moveaxis(array(s21_folder[q]),1,0)) # reshape to (mixer, atte, freq)
                dict_[f'{q}_freq'] = (["mixer","ro_amp","freq"],array(list(self._ro_elements[q])*2*self._power_samples.shape[0]).reshape(2,self._power_samples.shape[0],self._ro_elements[q].shape[0]))
            
            dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"ro_amp":self._power_samples,"freq":self.__freq_datapoint_idx})
            dataset.attrs["method"] = "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)




 
    
        