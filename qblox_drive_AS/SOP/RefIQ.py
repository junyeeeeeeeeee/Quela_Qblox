import os
from numpy import array, arange, minimum, maximum, cos, sin, pi, linspace
from numpy import nan as NaN
from xarray import Dataset
import matplotlib.pyplot as plt
from qblox_drive_AS.support import check_acq_channels
from qblox_drive_AS.support.UserFriend import *
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support.Pulse_schedule_library import Single_shot_ref_fit_analysis, pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import BinMode, Schedule
from quantify_scheduler.operations.gate_library import Measure, Reset



def Single_shot_fit_plot(results:dict,title_qubit:str=None,save_pic_folder:str=None):
    c_I,c_Q,sig=1000*results['fit_pack'][0],1000*results['fit_pack'][1],1000*results['fit_pack'][2]
    I,Q= results['data'][0],results['data'][1]
    radius = sig*2
    theta = linspace(0, 2 * pi, 720)
    x = c_I + radius * cos(theta)
    y = c_Q + radius * sin(theta)

    fig, ax = plt.subplots(nrows =1,figsize =(8,8),dpi =200)
    ax.scatter(1000*I, 1000*Q, color="blue", alpha=0.5, s=5)       
    ax.scatter(c_I,c_Q,c='k',s=15)
    ax.plot(x, y, linewidth=0.8, linestyle='--', c='red')
    ax.set_xlabel(r"$I\ $(mV)",size ='15')
    ax.set_ylabel(r"$Q\ $(mV)",size ='15')
    ax.set_title(f'{title_qubit if title_qubit is not None else ""} Single shot raw data')
    ax.set_xlim(1000*minimum(min(I),min(Q)),1000*maximum(max(I),max(Q)))
    ax.set_ylim(1000*minimum(min(I),min(Q)),1000*maximum(max(I),max(Q)))
    fig.tight_layout()
    plt.grid()
    if save_pic_folder is not None:
        plt.savefig(os.path.join(save_pic_folder,f"RefIQ_{title_qubit}.png"))
        plt.close()
    else:
        plt.show()

def IQ_ref_ana(ds:Dataset, q:str, save_pic_folder:str=None):
    analysis_result = Single_shot_ref_fit_analysis(ds[q])
    ref_IQ = array([analysis_result['fit_pack'][0],analysis_result['fit_pack'][1]])
    Single_shot_fit_plot(analysis_result,q,save_pic_folder)
    return ref_IQ



class RefIQPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._ro_elements:dict = {}
        self._avg_n:int = 100
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles

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
        meas_qs:list,
        repetitions:int=1,    
    ) -> Schedule:
        
        
        sched = Schedule("Ref IQ positioning Sche",repetitions=repetitions)

        for q in meas_qs:
            sched.add(Reset(q))
            
        sched.add(Measure(*meas_qs,  acq_index=0, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.APPEND) )
        
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        for q in self._ro_elements:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
        
            eyeson_print(f"Inte_time= {round(qubit_info.measure.integration_time()*1e6,1)} µs")
            eyeson_print(f"Reset_time= {round(qubit_info.reset.duration()*1e6,1)} µs")
        
            qubit_info.measure.pulse_amp(self._ro_elements[q]*float(qubit_info.measure.pulse_amp()))
            if qubit_info.rxy.amp180() is NaN:
                qubit_info.rxy.amp180(0)
            if qubit_info.rxy.duration() is NaN:
                qubit_info.rxy.duration(0)

        
        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))

        self.__spec_sched_kwargs = dict(   
        meas_qs=list(self._ro_elements.keys())
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
            
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            iq_tuples = self.__gettable.get()
            
            dict_ = {}
            for q_idx, q in enumerate(self._ro_elements):
                IQ_array = array([iq_tuples[2*q_idx],iq_tuples[2*q_idx+1]])
                dict_[q] = (["mixer","shots"],IQ_array)
            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"shots":arange(self._avg_n)})
            ds.attrs["system"] = "qblox"
            ds.attrs["method"] = "Average"
            
            self.dataset = ds
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)

  
