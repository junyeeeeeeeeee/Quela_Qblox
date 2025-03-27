import time, os
from datetime import datetime
from numpy import array, arange, ndarray, linspace
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support import Data_manager, check_OS_model_ready, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down, check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, IdlePulse, Measure, X, X90, SquarePulse, electrical_delay, ConditionalReset, BinMode
from quantify_scheduler.operations.gate_library import Reset
from qblox_drive_AS.support.ExpFrames import ExpGovernment
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
from qblox_drive_AS.SOP.FluxQubit import z_pulse_amp_OVER_const_z

class ZZinteractionPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._xy_elements:list = []
        self._zc_elements:list = []
        self._time_samples:dict = {}
        self._zamp_samples:ndarray = []


    @property
    def time_samples( self ):
        return self._time_samples
    @time_samples.setter
    def set_time_samples(self, time_samples:dict):
        if not isinstance(time_samples,dict):
            raise TypeError("Arg 'time_samples' must be a dict !")
        self._time_samples = time_samples
    
    
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
        xy_qs:list,
        z_cs:list,
        Z_amp:any,
        repetitions:int=1,
        activeReset:bool=False,
        singleshot:bool=False
    ) -> Schedule:
        qubits2read = list(freeduration.keys())
        freeDu_sameples = array(freeduration[qubits2read[0]])
        sched = Schedule("ZZ interaction", repetitions=repetitions)
        
        for acq_idx, freeDu in enumerate(freeDu_sameples): 
            align_pulse = sched.add(IdlePulse(4e-9))
            for q in qubits2read:
            
                if not activeReset:
                    reset = sched.add(Reset(q), ref_op=align_pulse)
                else:
                    reset = sched.add(
                        ConditionalReset(q, acq_index=acq_idx),
                        label=f"Reset {acq_idx}",
                        ref_op=align_pulse
                    )
                
                # first pi/2
                first_pulse = sched.add(X90(q), ref_op=reset)
                pi = sched.add(X(q), ref_op=first_pulse, rel_time=0.5*freeDu)
                final_pulse = sched.add(X90(q), ref_op=pi, rel_time=0.5*freeDu)
            
            # probe_pi
            for q in xy_qs:
                sched.add(X(q), ref_op=pi, ref_pt="start")
            
            # z 
            for q in z_cs:
                sched.add(SquarePulse(amp=Z_amp, duration=0.5*freeDu, port=f"{q}:fl", clock="cl0.baseband"), ref_op=first_pulse)
                sched.add(SquarePulse(amp=Z_amp, duration=0.5*freeDu, port=f"{q}:fl", clock="cl0.baseband"), ref_op=pi)
                
            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex' if not activeReset else 'ThresholdedAcquisition', bin_mode=BinMode.APPEND if singleshot else BinMode.AVERAGE), rel_time=electrical_delay, ref_op=final_pulse)
        
            
        return sched
        

    def __SetParameters__(self, *args, **kwargs):

        self.__time_data_idx = arange(len(list(self._time_samples.values())[0]))
        for q in self._time_samples:
            qubit_info:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
            qubit_info.reset.duration(qubit_info.reset.duration()+max(self._time_samples[q]))
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_zamp = ManualParameter(name="z_amp", unit="V", label="Amplitude")
        self.__Para_zamp.batched = False

        
        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")
        
        self.QD_agent = check_acq_channels(self.QD_agent, list(self._time_samples.keys()))
        
        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        xy_qs=self._xy_elements,
        z_cs=self._zc_elements,
        Z_amp=self.__Para_zamp,
        singleshot=self._os_mode,
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(self._time_samples),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            if not self._os_mode:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__Para_zamp])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,self._zamp_samples*z_pulse_amp_OVER_const_z))
            else:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__one_shot_para,self.__Para_zamp])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,arange(self._avg_n),self._zamp_samples*z_pulse_amp_OVER_const_z))
        else:
            time_preview_para = {}
            for q in self._time_samples:
                time_preview_para[q] = array([self._time_samples[q][1],self._time_samples[q][-1]])
            
            
            self.__sched_kwargs['freeduration']= time_preview_para
            self.__sched_kwargs['Z_amp']= self._zamp_samples[1]
    
    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run('RamseyT2')
            dict_ = {}
            if not self._os_mode:
                for q_idx, q in enumerate(self._time_samples):
                    i_data = array(ds[f'y{2*q_idx}']).reshape(self._zamp_samples.shape[0],self.__time_data_idx.shape[0])
                    q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._zamp_samples.shape[0],self.__time_data_idx.shape[0])
                    dict_[q] = (["mixer","bias","idx"],array([i_data,q_data]))
                    time_values = list(self._time_samples[q])*2*self._zamp_samples.shape[0]
                    dict_[f"{q}_x"] = (["mixer","bias","idx"],array(time_values).reshape(2,self._zamp_samples.shape[0],self.__time_data_idx.shape[0]))
                
                    dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":self._zamp_samples,"idx":self.__time_data_idx})
        
            else:
                dict_ = {}
                for q_idx, q in enumerate(self._time_samples):
                    i_data = array(ds[f'y{2*q_idx}']).reshape(self._zamp_samples.shape[0],self._avg_n,self.__time_data_idx.shape[0])
                    q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._zamp_samples.shape[0],self._avg_n,self.__time_data_idx.shape[0])
                    dict_[q] = (["mixer","prepared_state","bias","index","time_idx"],array([[i_data],[q_data]]))
                    time_values = list(self._time_samples[q])*2*self._zamp_samples.shape[0]*self._avg_n
                    dict_[f"{q}_x"] = (["mixer","prepared_state","bias","index","time_idx"],array(time_values).reshape(2,1,self._zamp_samples.shape[0],self._avg_n,self.__time_data_idx.shape[0]))

                dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":self._zamp_samples,"prepared_state":array([1]),"index":arange(self._avg_n),"time_idx":self.__time_data_idx})
                
                
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
        self.histos:int = 1
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

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self)->None:

        self.time_samples = sort_elements_2_multiples_of(self.time_samples*1e9, 2)*1e-9
        
        self.time_samples_dict = {}
        for q in self.read_qs:
            self.time_samples_dict[q] = self.time_samples
        

        self.avg_n = int(self.avg_n)
        
        # FPGA memory limit guard
        if self.OSmode:
            if self.avg_n * self.time_samples.shape[0] > 131000:
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
        meas.set_xy_elements = self.probe_qs
        meas.set_zc_elements = self.bias_cs
        meas.set_time_samples = self.time_samples_dict
        meas.set_os_mode = self.OSmode
        meas.n_avg = self.avg_n
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
            
            pass 
            # ds = open_dataset(file_path)
            # md = None
            
            # for var in ds.data_vars:
            #     if var.split("_")[-1] != 'x':
            #         self.ANA = Multiplex_analyzer("m12")
            #         if ds.attrs['method'].lower() == "shot":
            #             md = QD_savior.StateDiscriminator.summon_discriminator(var)
            #         if QD_savior.rotate_angle[var][0] != 0:
            #             ref = QD_savior.rotate_angle[var]
            #         else:
            #             eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
            #             ref = QD_savior.refIQ[var]
            #         self.ANA._import_data(ds,var_dimension=2,refIQ=ref)
            #         self.ANA._start_analysis(var_name=var, OSmodel=md)
            #         if self.save_pics:
            #             self.ANA._export_result(fig_path)

            #         """ Storing """
            #         if self.histos >= 50:
            #             QD_savior.Notewriter.save_echoT2_for(self.ANA.fit_packs["median_T2"],var)
                   
            # ds.close()
            # if self.keep_QD:
            #     QD_savior.QD_keeper()
            

    def WorkFlow(self):
        for i in range(self.histos):
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()   
            