import time, os
from datetime import datetime
from numpy import array, arange, ndarray, linspace, rad2deg, logspace
from numpy import pi as PI
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support.Notebook import Notebook
from qblox_drive_AS.support import Data_manager, check_OS_model_ready, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down, check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, IdlePulse, Measure, X90, Y90, DRAGPulse, Rxy, X, electrical_delay, ConditionalReset, BinMode
from quantify_scheduler.operations.gate_library import Reset
from qblox_drive_AS.support.ExpFrames import ExpGovernment
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of


class Ramsey12PS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._time_samples:dict = {}
        self._use_arti_detune:bool = False
        self._repeat:int = 1
        self._second_phase:str = 'x'

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
    def arti_detune_status(self):
        return self._use_arti_detune
    @arti_detune_status.setter
    def enable_arti_detune( self, switch:bool):
        if not isinstance(switch,bool):
            if switch in [0,1]:
                pass
            else:
                raise TypeError("Arg 'os_mode' must be a bool, 0 or 1 !")
        
        self._use_arti_detune = switch

    @property
    def second_phase( self ):
        return self._second_phase
    @second_phase.setter
    def set_second_phase(self, phase:str):
        if phase.lower() in ["x", "y"]:
            self._second_phase = phase.lower()
        else:
            print(f"Second phase should be 'x' or 'y', but '{phase}' was given. We choose to use 'x' instead.")
            self._second_phase = 'x'


    def __PulseSchedule__(self,
        freeduration:dict,
        Note:Notebook,
        repetitions:int=1,
        second_pulse_phase:str='x',
        arti_detune:dict={},
        singleshot:bool=False,
        
    ) -> Schedule:
        qubits2read = list(freeduration.keys())
        sched = Schedule("Ramsey", repetitions=repetitions)
        
        for acq_idx in range(array(freeduration[qubits2read[0]]).shape[0]): 
            align_pulse = sched.add(IdlePulse(4e-9))
            for qubit_idx, q in enumerate(qubits2read):
                
                freeDu = freeduration[q][acq_idx]
                
                sched.add(Reset(q), ref_op=align_pulse)
                # prepare |1>
                sched.add(X(q))
                # first pi/2
                first_pulse = sched.add(DRAGPulse(G_amp=Note.get_12ampFor(q)/2,D_amp=0,phase=0,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"))

                if q in arti_detune:
                    recovery_phase = rad2deg(2 * PI * arti_detune[q] * freeDu)
                    sched.add(
                        DRAGPulse(G_amp=Note.get_12ampFor(q)/2,D_amp=0,phase=recovery_phase,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"), ref_op=first_pulse, rel_time=freeDu
                    )
                else:
                    if second_pulse_phase.lower()=='x':
                        sched.add(DRAGPulse(G_amp=Note.get_12ampFor(q)/2,D_amp=0,phase=0,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"), ref_op=first_pulse, rel_time=freeDu)
                    else:
                        sched.add(DRAGPulse(G_amp=Note.get_12ampFor(q)/2,D_amp=0,phase=90,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"), ref_op=first_pulse, rel_time=freeDu)
            
                sched.add(X(q))
            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.APPEND if singleshot else BinMode.AVERAGE), rel_time=electrical_delay)
        
        return sched
        

    def __SetParameters__(self, *args, **kwargs):
        
        self.__time_data_idx = arange(len(list(self._time_samples.values())[0]))
        for q in self._time_samples:
            qubit_info:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
            
            qubit_info.reset.duration(qubit_info.reset.duration()+max(self._time_samples[q]))
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__repeat_data_idx = arange(self._repeat)

        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")

        self.QD_agent = check_acq_channels(self.QD_agent, list(self._time_samples.keys()))
        arti_detunes = {}
        if self._use_arti_detune:
            for q in list(self._time_samples.keys()):
                arti_detunes[q] = self.QD_agent.Notewriter.get_12artiDetuneFor(q)

        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        Note=self.QD_agent.Notewriter,
        singleshot=self._os_mode,
        second_pulse_phase=self._second_phase,
        arti_detune=arti_detunes
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(self._time_samples.keys())),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            if not self._os_mode:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,self.__repeat_data_idx))
            else:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__one_shot_para,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,arange(self._avg_n),self.__repeat_data_idx))
        else:
            preview_para = {}
            for q in self._time_samples:
                preview_para[q] = array([self._time_samples[q][0],self._time_samples[q][-1]])
            self.__sched_kwargs['freeduration']= preview_para
    
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
                
                
            for var in [i for i in dataset.data_vars if "_" not in i]:
                dataset[var].attrs["spin_num"] = 0
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["states"] = "EF"
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)


class EF_Ramsey(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID
        self.histos:int = 0
        self.sec_phase = 'x'

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, max_evo_time:float, target_qs:list, time_sampling_func:str, time_pts_or_step:int|float=100,histo_counts:int=1, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * max_evo_time: 100e-6\n
            * target_qs: ["q0", "q1", ...]
            * histo_counts: int, larger than 100 use while loop.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in target_qs:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(0, max_evo_time,time_pts_or_step)*1e9,1)*1e-9
        else:
            for q in target_qs: 
                self.time_samples[q] = sampling_func(0, max_evo_time,time_pts_or_step)

        self.avg_n = avg_n

        if histo_counts <= 100:
            self.want_while = False
            self.histos = histo_counts
        else:
            self.want_while = True
            self.histos = 1
        
        self.execution = execution
        self.OSmode = OSmode
        self.spin_num = {}
        self.target_qs = target_qs
        
        # FPGA memory limit guard
        if self.OSmode:
            for q in self.time_samples:
                if self.avg_n * self.time_samples[q].shape[0] > 131000:
                    raise MemoryError(f"Due to Qblox FPGA memory limit, time_pts * shots must be less than 131000. And you are trying to set {self.avg_n * self.time_samples[q].shape[0]} for {q}")
        
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        if self.OSmode:
            check_OS_model_ready(self.QD_agent, self.target_qs)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            slightly_print(f"{q} f12 arti-detune = {round(self.QD_agent.Notewriter.get_12artiDetuneFor(q)*1e-6,2)} MHz")
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        
        meas = Ramsey12PS()
        meas._use_arti_detune= True
        meas.set_time_samples = self.time_samples
        meas.set_os_mode = self.OSmode
        meas.n_avg = self.avg_n
        meas.set_repeat = self.histos
        meas.set_second_phase = self.sec_phase
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset

        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"Ramsey_EF_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
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
            self.corrected_detune = {}
            
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

                    if self.sec_phase.lower() == 'y':
                        if (self.ANA.fit_packs['phase'] % 360 + 360) % 360 > 180:
                            sign = -1
                        else:
                            sign = 1
                        
                        self.corrected_detune[var] = sign*self.ANA.fit_packs['freq']
                    self.ANA.fit_packs.update({"plot_item":self.ANA.plot_item})
                    self.sum_dict[var] = self.ANA.fit_packs
                    """ Storing """
                    if self.histos >= 50:
                        QD_savior.Notewriter.save_T2_for(self.ANA.fit_packs["median_T2"],var)
                   
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

class f12_calibration(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID
        self.histos:int = 0
        self.sec_phase = 'y'

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, max_evo_time:float, target_qs:list, time_sampling_func:str, time_pts_or_step:int|float=100,histo_counts:int=1, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * max_evo_time: 100e-6\n
            * target_qs: ["q0", "q1", ...]
            * histo_counts: int, larger than 100 use while loop.\n
            * time_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * time_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if time_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(time_sampling_func)
        else:
            raise ValueError(f"Can't recognize the given sampling function name = {time_sampling_func}")
        
        self.time_samples = {}
        if sampling_func in [linspace, logspace]:
            for q in target_qs:
                self.time_samples[q] = sort_elements_2_multiples_of(sampling_func(0, max_evo_time,time_pts_or_step)*1e9,1)*1e-9
        else:
            for q in target_qs: 
                self.time_samples[q] = sampling_func(0, max_evo_time,time_pts_or_step)

        self.avg_n = avg_n

        if histo_counts <= 100:
            self.want_while = False
            self.histos = histo_counts
        else:
            self.want_while = True
            self.histos = 1
        
        self.execution = execution
        self.OSmode = OSmode
        self.spin_num = {}
        self.target_qs = target_qs
        
        # FPGA memory limit guard
        if self.OSmode:
            for q in self.time_samples:
                if self.avg_n * self.time_samples[q].shape[0] > 131000:
                    raise MemoryError(f"Due to Qblox FPGA memory limit, time_pts * shots must be less than 131000. And you are trying to set {self.avg_n * self.time_samples[q].shape[0]} for {q}")
        
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        if self.OSmode:
            check_OS_model_ready(self.QD_agent, self.target_qs)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            slightly_print(f"{q} arti-detune = {round(self.QD_agent.Notewriter.get_artiT2DetuneFor(q)*1e-6,2)} MHz")
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
        
        meas = Ramsey12PS()
        meas.set_time_samples = self.time_samples
        meas.set_os_mode = self.OSmode
        meas.n_avg = self.avg_n
        meas.set_repeat = self.histos
        meas.set_second_phase = self.sec_phase
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset

        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"Ramsey_EF_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
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

            QD_savior = QDmanager(QD_file)
            QD_savior.QD_loader()

            ds = open_dataset(file_path)
            answer = {}
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    ANA = Multiplex_analyzer("c2")
                    if QD_savior.rotate_angle[var][0] != 0:
                        ref = QD_savior.rotate_angle[var]
                    else:
                        eyeson_print(f"{var} rotation angle is 0, use contrast to analyze.")
                        ref = QD_savior.refIQ[var]
                    
                    ANA._import_data(ds,var_dimension=2,refIQ=ref)
                    ANA._start_analysis(var_name=var)
                    ANA._export_result(fig_path)

                    if (ANA.fit_packs['phase'] % 360 + 360) % 360 > 180:
                        sign = -1
                    else:
                        sign = 1
                    
                    answer[var] = sign*ANA.fit_packs['freq']
                    highlight_print(f"{var}: actual detune = {round(answer[var]*1e-6,4)} MHz")
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).clock_freqs.f12(QD_savior.quantum_device.get_element(permi).clock_freqs.f12()-answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).clock_freqs.f12(QD_savior.quantum_device.get_element(q).clock_freqs.f12()-answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")
            

    def WorkFlow(self):
        while True:
            self.PrepareHardware()

            self.RunMeasurement()

            self.CloseMeasurement()   
            if not self.want_while:
                break

