
from numpy import array, arange, ndarray
from numpy import nan as NaN
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, check_acq_channels
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support.Pulse_schedule_library import pulse_preview, BinMode
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Measure, Schedule, pulse_preview, SquarePulse
from quantify_scheduler.operations.pulse_library import IdlePulse,SetClockFrequency
from quantify_scheduler.resources import ClockResource
from qblox_drive_AS.SOP.FluxQubit import z_pulse_amp_OVER_const_z
from xarray import Dataset


def update_flux_info_in_results_for(QD_agent:QDmanager,qb:str,FD_results:dict):
    qubit = QD_agent.quantum_device.get_element(qb)
    qubit.clock_freqs.readout(FD_results["sweet_freq"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=FD_results["sweet_flux"])
    QD_agent.Fluxmanager.save_period_for(target_q=qb, period=1/FD_results['f'])
    QD_agent.Fluxmanager.save_tuneawayBias_for(target_q=qb,mode='auto')
    QD_agent.Fluxmanager.save_cavFittingParas_for(target_q=qb,
        f=FD_results['f'],
        amp=FD_results['A'],
        phi=FD_results['phi'],
        offset=FD_results['offset']
    )



class FluxDepCavityPS(ScheduleConductor):
    def __init__(self):
        super().__init__() 
        self._ro_elements:dict = {}
        self._flux_samples:ndarray = []
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles
    
    @property
    def flux_samples(self):
        return self._flux_samples
    @flux_samples.setter
    def flux_samples(self,fluxs:ndarray):
        self._flux_samples = fluxs

    def __PulseSchedule__(self, 
        frequencies: dict,
        bias_dura:float,
        bias:any,
        repetitions:int=1,    
    ) -> Schedule:
        
        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
        sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)
        

        for acq_idx in range(sameple_idx):    
            align_pulse = sched.add(IdlePulse(duration=4e-9), label=f"buffer {acq_idx}")
            for qubit_idx, q in enumerate(qubits2read):
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
                      
                sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
                sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")
                sched.add(SquarePulse(bias, bias_dura, port=q+":fl", clock="cl0.baseband"),ref_op=align_pulse)
                
            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE), ref_op=align_pulse, ref_pt='end')
            

        self.schedule =  sched  
        return sched
    

    def __SetParameters__(self, *args, **kwargs):
         
        self.__freq_datapoint_idx = arange(0,len(list(list(self._ro_elements.values())[0])))
        original_rof = {}
        for q in self._ro_elements:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            flux_dura = qubit_info.measure.integration_time()
            original_rof[q] = qubit_info.clock_freqs.readout()
            qubit_info.clock_freqs.readout(NaN)
            
        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))

        self.__freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
        self.__freq.batched = True
        self.__bias = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
        self.__bias.batched = False

        self.__spec_sched_kwargs = dict(   
            frequencies=self._ro_elements,
            bias_dura=flux_dura,
            bias=self.__bias
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
            self.meas_ctrl.settables([self.__freq, self.__bias])
            self.meas_ctrl.setpoints_grid([self.__freq_datapoint_idx, self._flux_samples*z_pulse_amp_OVER_const_z])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
             
            self.__spec_sched_kwargs['bias']= self._flux_samples[-1]*z_pulse_amp_OVER_const_z
            self.__spec_sched_kwargs['frequencies']= preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run("One-tone-Flux")
            dict_ = {}
            for idx, q in enumerate(self._ro_elements):
                freq_values = 2*self._flux_samples.shape[0]*list(self._ro_elements[q])
                i_data = array(ds[f'y{2*idx}']).reshape(self._flux_samples.shape[0],array(self._ro_elements[q]).shape[0])
                q_data = array(ds[f'y{2*idx+1}']).reshape(self._flux_samples.shape[0],array(self._ro_elements[q]).shape[0])
                dict_[q] = (["mixer","bias","freq"],array([i_data,q_data]))
                dict_[f"{q}_freq"] = (["mixer","bias","freq"],array(freq_values).reshape(2,self._flux_samples.shape[0],array(self._ro_elements[q]).shape[0]))

            
            rfs_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":self._flux_samples,"freq":self.__freq_datapoint_idx})
            rfs_ds.attrs["method"] = "Average"
            rfs_ds.attrs["system"] = "qblox"
            self.dataset = rfs_ds
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)