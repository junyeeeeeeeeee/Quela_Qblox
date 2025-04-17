"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from numpy import array, arange
from qblox_drive_AS.support import Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import check_acq_channels, BasicTransmonElement
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import BinMode, Schedule, pulse_preview, Rxy, Reset, IdlePulse, Measure

class hPiAcalibrationPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._ro_elements:dict = {}
        self._hpi_quads:list = []
    
    @property
    def ro_elements(self):
        return self._ro_elements
    @ro_elements.setter
    def ro_elements(self, ro_eles:dict):
        self._ro_elements = ro_eles
    @property
    def hpi_quads(self):
        return self._hpi_quads
    @hpi_quads.setter
    def hpi_quads(self, quads:list):
        if not isinstance(quads, list):
            raise TypeError("quadruples must be a list !")
        self._hpi_quads = quads
    
    

    def __PulseSchedule__(self, 
        new_degrees:dict,
        pi_druple_num:any,
        repetitions:int=1,
    )-> Schedule:
        qubits2read = list(new_degrees.keys())
        sched = Schedule("Pi half degree calibration", repetitions=repetitions)

        for acq_idx in range(array(new_degrees[qubits2read[0]]).shape[0]):
            align_pulse = sched.add(IdlePulse(4e-9))
            for q in qubits2read:
                new_theta = new_degrees[q][acq_idx]
                
                sched.add(Reset(q), ref_op=align_pulse)
                
                for pi_num in range(pi_druple_num):
                    for pi_idx in range(4):
                        half_pi = sched.add(Rxy(qubit=q, theta=new_theta, phi=0))

            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol="SSBIntegrationComplex", bin_mode=BinMode.AVERAGE))
                        
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(len(list(list(self._ro_elements.values())[0])))
        self.__theta_samples = {}
        
        for q in self._ro_elements:
            qubit_info:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
            self.__theta_samples[q] = 90 * self._ro_elements[q] 
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
           
        self.__theta = ManualParameter(name="theta", unit="degree", label="angle")
        self.__theta.batched = True
        self.__hpi_quadruples =  ManualParameter(name="half_pi_quadruples", unit="", label="number")
        self.__hpi_quadruples.batched = False


        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))
        self.__spec_sched_kwargs = dict(   
        new_degrees=self.__theta_samples,
        pi_druple_num=self.__hpi_quadruples,
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
            self.meas_ctrl.settables([self.__theta, self.__hpi_quadruples])
            self.meas_ctrl.setpoints_grid([self.__datapoint_idx, self._hpi_quads])
        
        else:
            preview_para = {}
            for q in self.__theta_samples:
                preview_para[q] = array([self.__theta_samples[q][1], self.__theta_samples[q][-2]])
            self.__spec_sched_kwargs['new_degrees']= preview_para
            self.__spec_sched_kwargs['pi_druple_num']= self._hpi_quads[0]
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rs_ds = self.meas_ctrl.run("half pi-amp calibration")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                coefs = 2*len(self._hpi_quads)*list(self._ro_elements[q])
                i_data = array(rs_ds[f'y{2*q_idx}']).reshape(len(self._hpi_quads),self.__datapoint_idx.shape[0])
                q_data = array(rs_ds[f'y{2*q_idx+1}']).reshape(len(self._hpi_quads),self.__datapoint_idx.shape[0])
                dict_[q] = (["mixer", "PiPairNum", "HalfPiCoef"],array([i_data,q_data]))
                dict_[f'{q}_HalfPIcoef'] = (["mixer", "PiPairNum", "HalfPiCoef"],array(coefs).reshape(2,len(self._hpi_quads),self.__datapoint_idx.shape[0]))

            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"PiPairNum":array(self._hpi_quads),"HalfPiCoef":self.__datapoint_idx})
            
            ds.attrs["execution_time"] = Data_manager().get_time_now()
            ds.attrs["method"] = "Average"
            ds.attrs["system"] = "qblox"
            self.dataset = ds
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)

