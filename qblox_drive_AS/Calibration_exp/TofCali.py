import os, time
from datetime import datetime
from numpy import  linspace, array, arange
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Measure
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from quantify_scheduler.backends.qblox import constants
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support import Data_manager, dh, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from quantify_scheduler.math import closest_number_ceil
from qblox_drive_AS.support.ExpFrames import ExpGovernment

class TOFcaliPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._target_q:str
        

    @property
    def target_q( self ):
        return self._target_q
    @target_q.setter
    def target_q(self, q:str):
        if not isinstance(q, str):
            raise TypeError("q must be a str !")
        self._target_q = q
    
    

    def __PulseSchedule__(self,
        qubit_name: str,
        repetitions: int = 1,
    ) -> Schedule:
        schedule = Schedule("Trace measurement schedule", repetitions=repetitions)
        schedule.add(Measure(qubit_name, acq_protocol="Trace"))
        return schedule
        
        

    def __SetParameters__(self, *args, **kwargs):
        
        qubit:BasicTransmonElement = self.QD_agent.quantum_device.get_element(self._target_q)
        self.__tof_t = ManualParameter(name="tof_t", unit="ns", label="Trace acquisition sample")
        self.__tof_t.batched = True
        self.__tof_t.batch_size = round(qubit.measure.integration_time() * constants.SAMPLING_RATE)

        self.__sched_kwargs = dict(
            qubit_name=self._target_q,
        )

    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=False,
                batched=True,
                num_channels=1,
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            tof_t_setpoints = arange(self.__tof_t.batch_size)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables(self.__tof_t)
            self.meas_ctrl.setpoints(tof_t_setpoints)
            
        else:
            self.__sched_kwargs['qubit_name']= self._target_q
            
    
    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run('TOF_calibration')
            self.tuid = dh.get_latest_tuid()
            
            dataset = dh.to_gridded_dataset(ds)
                
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)

class TofCalirator(ExpGovernment):
    
    def __init__(self,data_folder:str=None,JOBID:str=None):
        
        super().__init__()
        self.QD_path:str = ""
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID
        self.execution:bool = True
        self.q:str = "q0"
        

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self)->None:
        pass
        
        
        
    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
    
        # offset bias, LO and atte
        for q in [self.q]:
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=0)
        

    def RunMeasurement(self):
        meas = TOFcaliPS()
        meas.target_q = self.q
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset
        self.tuid = meas.tuid
        # dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"TOFcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                try:
                    dataset.to_netcdf(self.__raw_data_location)
                except:
                    print(f"Fail to save dataset, TUID: {self.tuid}")
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_QDagent:QDmanager=None,playback_delay:float=149e-9):
        """ User callable analysis function pack """
        from quantify_core.analysis.time_of_flight_analysis import TimeOfFlightAnalysis
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_QDagent is None:
                QD_savior = QDmanager(QD_file)
                QD_savior.QD_loader()
            else:
                QD_savior = new_QDagent
            
            tof_analysis = TimeOfFlightAnalysis(tuid=self.tuid)
            tof_analysis.run(playback_delay=playback_delay)
            fit_results = tof_analysis.quantities_of_interest
            nco_prop_delay = fit_results["nco_prop_delay"]
            measured_tof = fit_results["tof"]
            acq_delay = closest_number_ceil(measured_tof * constants.SAMPLING_RATE, constants.MIN_TIME_BETWEEN_OPERATIONS) / constants.SAMPLING_RATE
            for q in QD_savior.quantum_device.elements():
                qubit:BasicTransmonElement = QD_savior.quantum_device.get_element(q)
                qubit.measure.acq_delay(acq_delay)
            

            print(f"nco_prop_delay: {nco_prop_delay},\n acq_delay: {acq_delay}")
            QD_savior.Notewriter.save_NcoPropDelay(nco_prop_delay)

            QD_savior.QD_keeper()

           
            

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()