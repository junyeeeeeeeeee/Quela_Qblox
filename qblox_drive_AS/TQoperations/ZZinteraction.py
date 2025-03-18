import time, os
from datetime import datetime
from numpy import array, arange, ndarray
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager
from qblox_drive_AS.support import compose_para_for_multiplexing, Data_manager, check_OS_model_ready, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Readout, Multi_Readout, Integration, electrical_delay, Z
from quantify_scheduler.operations.gate_library import Reset
from qblox_drive_AS.support.ExpFrames import ExpGovernment
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer


class ZZinteractionPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._xy_elements:list = []
        self._zc_elements:list = []
        self._time_samples:dict = {}
        self._zamp_samples:ndarray = []
        self._os_mode:bool = False
        self._repeat:int = 1
        self._avg_n:int = 300

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
    def os_mode( self ):
        return self._os_mode
    @os_mode.setter
    def set_os_mode(self, os_mode:bool):
        if not isinstance(os_mode,bool):
            if os_mode in [0,1]:
                pass
            else:
                raise TypeError("Arg 'os_mode' must be a bool, 0 or 1 !")
        
        self._os_mode = os_mode
    @property
    def n_avg( self ):
        return self._avg_n
    @n_avg.setter
    def set_n_avg(self, avg_num:int):
        if not isinstance(avg_num,(int,float)):
            raise TypeError("Ard 'avg_num' must be a int or float !")
        self._avg_n = int(avg_num)
    
    @property
    def ro_elements( self ):
        return self._ro_elements
    @ro_elements.setter
    def set_ro_elements(self, ro_elements:list):
        if not isinstance(ro_elements, list):
            raise TypeError("ro_elements must be a list !")
        self._ro_elements = ro_elements
    @property
    def xy_elements( self ):
        return self._xy_elements
    @xy_elements.setter
    def set_xy_elements(self, xy_elements:list):
        if not isinstance(xy_elements, list):
            raise TypeError("xy_elements must be a list !")
        self._xy_elements = xy_elements
    @property
    def zc_elements( self ):
        return self._zc_elements
    @zc_elements.setter
    def set_zc_elements(self, zc_elements:list):
        if not isinstance(zc_elements, list):
            raise TypeError("zc_elements must be a list !")
        self._zc_elements = zc_elements
    @property
    def zamp_samples( self ):
        return self._zamp_samples
    @zamp_samples.setter
    def set_zamp_samples(self, zamp_samples:ndarray):
        if not isinstance(zamp_samples, ndarray):
            raise TypeError("zamp_samples must be a ndarray !")
        
        self._zamp_samples = zamp_samples

    def __PulseSchedule__(self,
        freeduration:dict,
        xy_qs:dict,
        z_cs:list,
        Z_amp:any,
        pi_amp: dict,
        pi_dura:dict,
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        repetitions:int=1,
        singleshot:bool=False
    ) -> Schedule:
        qubits2read = list(freeduration.keys())
        sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
        sched = Schedule("Ramsey", repetitions=repetitions)
        
        for acq_idx in range(sameple_idx):  
            freeDu = freeduration[qubits2read[0]][acq_idx] 
            for qubit_idx, q in enumerate(qubits2read):
                
                sched.add(Reset(q))

                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)

                
                first_half_pi = self.QD_agent.Waveformer.X_pi_2_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu=electrical_delay)
                  
                pi = self.QD_agent.Waveformer.X_pi_p(sched,pi_amp,q,pi_dura[q],first_half_pi,0.5*freeDu)
                
                self.QD_agent.Waveformer.X_pi_2_p(sched,pi_amp,q,pi_dura[q],pi,0.5*freeDu)
                

                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=singleshot,get_trace=False,trace_recordlength=0)

            for q in xy_qs:
                self.QD_agent.Waveformer.X_pi_p(sched,{q:xy_qs[q]["pi_amp"]},q,xy_qs[q]["pi_du"],pi,freeDu=0,ref_point='end')
                
            if freeDu != 0:
                for c in z_cs:
                    Z(sched, Z_amp, 0.5*freeDu, c, first_half_pi, freeDu=0)
                    Z(sched, Z_amp, 0.5*freeDu, c, pi, freeDu=0)

            
        return sched
        

    def __SetParameters__(self, *args, **kwargs):

        self.__time_data_idx = arange(len(list(self._time_samples.values())[0]))
        for q in self._time_samples:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            qubit_info.reset.duration(qubit_info.reset.duration()+max(self._time_samples[q]))
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__Para_zamp = ManualParameter(name="z_amp", unit="V", label="Amplitude")
        self.__Para_zamp.batched = False
        self.__repeat_data_idx = arange(self._repeat)


        drive_elements = {}
        for q in self._xy_elements:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            drive_elements[q] = {}
            drive_elements[q]['pi_du'] = qubit_info.rxy.duration()
            drive_elements[q]['pi_amp']= qubit_info.rxy.amp180()
        
        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")

        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        xy_qs=drive_elements,
        z_cs=self._zc_elements,
        Z_amp=self.__Para_zamp,
        singleshot=self._os_mode,
        pi_amp=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'d1'),
        pi_dura=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'d3'),
        R_amp=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r1'),
        R_duration=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r3'),
        R_integration=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r4'),
        R_inte_delay=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r2'),
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(self._ro_elements),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            if not self._os_mode:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__Para_zamp,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,self._zamp_samples,self.__repeat_data_idx))
            else:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__one_shot_para,self.__Para_zamp,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,arange(self._avg_n),self._zamp_samples,self.__repeat_data_idx))
        else:
            time_preview_para = {}
            for q in self._ro_elements:
                time_preview_para[q] = array([self._time_samples[q][0],self._time_samples[q][-1]])
            Z_amps = self._zamp_samples[1]
            self.__sched_kwargs['freeduration']= time_preview_para
            self.__sched_kwargs['Z_amp']= Z_amps
    
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
                
                
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"
            dataset.attrs["biased_cs"] = self._zc_elements
            dataset.attrs["probe_qs"] = self._xy_elements
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)

class ZZinteraction(ExpGovernment):
    def __init__(self,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path:str = ""
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID
        self.z_amp_samples:ndarray = []
        self.time_samples:ndarray = []
        self.read_qs:list = []
        self.probe_qs:list = []
        self.bias_cs:list = []
        self.OSmode:bool = False
        self.avg_n:int = 300
        self.execution:bool = True
        self.histo_count:int = 1

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self)->None:
        
        self.time_samples_dict = {}
        for q in self.read_qs:
            self.time_samples_dict[q] = array(self.time_samples)
        

        self.avg_n = int(self.avg_n)

        if self.histo_count <= 100:
            self.want_while = False
            self.histos = self.histo_count
        else:
            self.want_while = True
            self.histos = 1
        
        # FPGA memory limit guard
        if self.OSmode:
            for q in self.time_samples:
                if self.avg_n * self.time_samples[q].shape[0] > 131000:
                    raise MemoryError(f"Due to Qblox FPGA memory limit, time_pts * shots must be less than 131000. And you are trying to set {self.avg_n * self.time_samples[q].shape[0]} for {q}")
        
    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        if self.OSmode:
            check_OS_model_ready(self.QD_agent, self.read_qs)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.read_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        for q in self.probe_qs:
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q], xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))

    def RunMeasurement(self):
        meas = ZZinteractionPS()
        meas.set_zamp_samples = self.z_amp_samples
        meas.set_ro_elements = self.read_qs
        meas.set_xy_elements = self.probe_qs
        meas.set_zc_elements = self.bias_cs
        meas.set_time_samples = self.time_samples_dict
        meas.set_os_mode = self.OSmode
        meas.set_n_avg = self.avg_n
        meas.set_repeat = self.histos
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset
        # dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"ZZinteraction_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None,new_QDagent:QDmanager=None,new_pic_save_place:str=None):
        """ User callable analysis function pack """
        
        if self.execution:
            if new_QD_path is None:
                QD_file = self.QD_path
            else:
                QD_file = new_QD_path

            if new_file_path is None:
                file_path = self.__raw_data_location
                fig_path = self.save_dir
            else:
                file_path = new_file_path
                fig_path = os.path.split(new_file_path)[0]
            
            if new_QDagent is None:
                QD_savior = QDmanager(QD_file)
                QD_savior.QD_loader()
            else:
                QD_savior = new_QDagent
            
            if new_pic_save_place is not None:
                fig_path = new_pic_save_place

            ds = open_dataset(file_path)
            md = None
            
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    self.ANA = Multiplex_analyzer("m12")
                    if ds.attrs['method'].lower() == "shot":
                        md = QD_savior.StateDiscriminator.summon_discriminator(var)
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    self.ANA._import_data(ds,var_dimension=2,refIQ=ref)
                    self.ANA._start_analysis(var_name=var, OSmodel=md)
                    if self.save_pics:
                        self.ANA._export_result(fig_path)

                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_echoT2_for(self.ANA.fit_packs["median_T2"],var)
                   
            ds.close()
            if self.keep_QD:
                QD_savior.QD_keeper()
            

    def WorkFlow(self):
        while True:
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()   
            if not self.want_while:
                break