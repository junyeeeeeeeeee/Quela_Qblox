import time, os
from datetime import datetime
from numpy import array, arange, ndarray, linspace
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support.Notebook import Notebook
from qblox_drive_AS.support import Data_manager, check_OS_model_ready, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down, check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, IdlePulse, Measure, X, DRAGPulse, ConditionalReset, BinMode
from quantify_scheduler.operations.gate_library import Reset
from qblox_drive_AS.support.ExpFrames import ExpGovernment
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of


class Rabi12PS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._RabiType:str = "power"
        self._variables:dict = {}
    

    @property
    def RabiType( self ):
        return self._RabiType
    @RabiType.setter
    def set_RabiType( self, type:str):
        if type.lower() in ["power", 'time']:
            self._RabiType = type.lower()
        else:
            raise ValueError("Arg 'type' must be given as 'time' or 'power' !")

    @property
    def samples( self ):
        return self._variables
    @samples.setter
    def set_samples( self, samples:dict):
        if not isinstance(samples, dict):
            raise TypeError("Arg 'samples' must be a dict !")
        self._variables = samples

    def __PulseSchedule__(self,
        Note:Notebook,
        pi_amp:dict,
        pi_dura:dict,
        XY_theta:str,
        repetitions:int=1,
        OS_or_not:bool=False
    ) -> Schedule:

        qubits2read = list(pi_amp.keys()) if self._RabiType.lower() == 'power' else list(pi_dura.keys())
        sample_len = pi_amp[qubits2read[0]].shape[0] if self._RabiType.lower() == 'power' else pi_dura[qubits2read[0]].shape[0]
        sched = Schedule("RabiOscillation",repetitions=repetitions)

        match XY_theta:
            case 'Y_theta':
                phi = 90.0
            case _:
                phi = 0

        
        for acq_idx in range(sample_len):
            align_pulse = sched.add(IdlePulse(4e-9))    
            for qubit_idx, q in enumerate(qubits2read):
                sched.add(Reset(q), ref_op=align_pulse)
                sched.add(X(q))

                if self._RabiType.lower() == 'power':
                    sched.add(DRAGPulse(G_amp=pi_amp[q][acq_idx],D_amp=0,phase=phi,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"))
                else:
                    sched.add(DRAGPulse(G_amp=Note.get_12ampFor(q),D_amp=0,phase=phi,duration=pi_dura[q][acq_idx],port=f"{q}:mw",clock=f"{q}.12"))

                sched.add(X(q))
            
            sched.add(Measure(*qubits2read, acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE))
        
        self.schedule = sched
        return sched

    def __SetParameters__(self, *args, **kwargs):
        self.amps, self.duras = {}, {}
        for q in self._variables:
            qubit_info:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
            eyeson_print(f"f12 = {round(qubit_info.clock_freqs.f12()*1e-9,3)} GHz")
            self.__variable_idx = arange(0,self._variables[q].shape[0])

        
        if self._RabiType == 'power':
            self.__Sweep_para = ManualParameter(name="PowerRabi", unit="V", label="amplitude")
            self.amps = self._variables
            for q in self._variables:
                self.duras[q] = self.QD_agent.quantum_device.get_element(q).rxy.duration()
        else:
            self.__Sweep_para = ManualParameter(name="TimeRabi", unit="sec", label="time")
            self.duras = self._variables
            for q in self._variables:
                self.amps[q] = self.QD_agent.quantum_device.get_element(q).rxy.amp180()
       
        self.__Sweep_para.batched = True
        
        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")
        
        self.QD_agent = check_acq_channels(self.QD_agent, list(self._variables.keys()))
        self.__sched_kwargs = dict(
            Note = self.QD_agent.Notewriter,
            pi_amp = self.amps,
            pi_dura = self.duras,
            XY_theta='X_theta',
            OS_or_not=self._os_mode
            )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(self._variables.keys())),
            )
        
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)

            if not self._os_mode:
                self.meas_ctrl.settables(self.__Sweep_para)
                self.meas_ctrl.setpoints(self.__variable_idx)
            else:
                self.meas_ctrl.settables([self.__Sweep_para, self.__one_shot_para])
                self.meas_ctrl.setpoints_grid((self.__variable_idx, arange(self._avg_n)))
        
        else:
            preview_para = {}
            for q in self._variables:
                preview_para[q] = array([self._variables[q][0],self._variables[q][-1]])
            
            if self._RabiType == "power":
                self.__sched_kwargs['pi_amp'] = preview_para
            else:
                self.__sched_kwargs['pi_dura'] = preview_para

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run("RabiOscillation")
            dict_ = {}

            if not self._os_mode:
                for idx, q in enumerate(list(self._variables.keys())):
                    I = array(ds[f'y{2*idx}'])
                    Q = array(ds[f'y{2*idx+1}'])
                    dict_[q] = (['mixer','var_idx'],array([I,Q]))
                    dict_[f"{q}_variable"] = (['mixer','var_idx'],array(2*list(self._variables[q])).reshape(2,self._variables[q].shape[0]))
                
                dataset = Dataset(dict_, coords={"mixer":array(["I","Q"]),"var_idx":self.__variable_idx})
                
            else:
                for idx, q in enumerate(list(self._variables.keys())):
                    I = array(ds[f'y{2*idx}']).reshape(self._avg_n,self._variables[q].shape[0])
                    Q = array(ds[f'y{2*idx+1}']).reshape(self._avg_n,self._variables[q].shape[0])
                    dict_[q] = (["mixer","prepared_state","index","var_idx"],array([[I],[Q]]))
                    variables =  list(self._variables[q])*2*self._avg_n
                    dict_[f"{q}_variable"] = (["mixer","prepared_state","index","var_idx"],array(variables).reshape(2,1,self._avg_n,self._variables[q].shape[0]))
                
                dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"prepared_state":array([1]),"index":arange(self._avg_n),"var_idx":self.__variable_idx})
                    
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["rabi_type"] = self._RabiType
            dataset.attrs["states"] = "EF"
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"            
            if self._RabiType == 'time':
                for q in self.amps:
                    dataset.attrs[f"{q}_piamp"] = self.amps[q]
            else:
                for q in self.duras:
                    dataset.attrs[f"{q}_pidura"] = self.duras[q]
            
            self.dataset = dataset
    
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)


class PowerRabi_12_Osci(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, pi_amp:dict, pi_amp_sampling_func:str, pi_amp_pts_or_step:float=100, avg_n:int=100, execution:bool=True, OSmode:bool=False):
        """ ### Args:
            * pi_amp: {"q0":[pi_amp_start, pi_amp_end], ...}\n
            * pi_amp_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * pi_amp_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        if pi_amp_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(pi_amp_sampling_func)
        else:
            sampling_func:callable = linspace
        
        self.pi_amp_samples = {}
        for q in pi_amp:
            self.pi_amp_samples[q] = sampling_func(*pi_amp[q],pi_amp_pts_or_step)
            
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(pi_amp.keys())

        # FPGA memory limit guard
        if self.OSmode:
            for q in self.pi_amp_samples:
                if self.avg_n * self.pi_amp_samples[q].shape[0] > 131000:
                    raise MemoryError(f"Due to Qblox FPGA memory limit, amp_pts * shots must be less than 131000. And you are trying to set {self.avg_n * self.pi_amp_samples[q].shape[0]} for {q}")

        
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and driving atte
        self.pi_dura = {}
        if self.OSmode:
            check_OS_model_ready(self.QD_agent, self.target_qs)
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    def RunMeasurement(self):
    
        meas = Rabi12PS()
        meas.execution = self.execution
        meas.set_samples = self.pi_amp_samples
        meas.set_RabiType = "power"
        meas.set_os_mode = self.OSmode
        meas.n_avg = self.avg_n
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        
        meas.run()
        dataset = meas.dataset
    
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"PowerRabi_12_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QDagent:QDmanager=None,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RabiOsci import conditional_update_qubitInfo
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

            ds = open_dataset(file_path)
            md = None

            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'variable':
                    ANA = Multiplex_analyzer("m11")
                    if ds.attrs['method'].lower() == "shot":
                        md = QD_savior.StateDiscriminator.summon_discriminator(var)   
                    
                    ANA._import_data(ds,1,QD_savior.refIQ[var] if QD_savior.rotate_angle[var][0] == 0 else QD_savior.rotate_angle[var])
                    ANA._start_analysis(var_name=var, OSmodel=md)
                    ANA._export_result(fig_path)
                    conditional_update_qubitInfo(QD_savior,ANA.fit_packs,var,fit_12=True)  

            ds.close()
            QD_savior.QD_keeper()
            



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   



class TimeRabi_12_Osci(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, pi_dura:dict, pi_dura_sampling_func:str, pi_dura_pts_or_step:float=100, avg_n:int=100, execution:bool=True, OSmode:bool=False):
        """ ### Args:
            * pi_amp: {"q0": pi-amp in V, ...}\n
            * pi_dura: {"q0":[pi_duration_start, pi_duration_end], ...}\n
            * pi_dura_sampling_func (str): 'linspace', 'arange', 'logspace'\n
            * pi_dura_pts_or_step: Depends on what sampling func you use, `linspace` or `logspace` set pts, `arange` set step. 
        """
        from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
        if pi_dura_sampling_func in ['linspace','logspace','arange']:
            sampling_func:callable = eval(pi_dura_sampling_func)
        else:
            sampling_func:callable = linspace
        
        self.pi_dura_samples = {}
        for q in pi_dura:
            if min(pi_dura[q]) == 0: pi_dura[q] = [4e-9, max(pi_dura[q])]
            self.pi_dura_samples[q] = sort_elements_2_multiples_of(sampling_func(*pi_dura[q],pi_dura_pts_or_step)*1e9,1)*1e-9
            
        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = list(pi_dura.keys())
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        if self.OSmode:
            check_OS_model_ready(self.QD_agent, self.target_qs)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        self.pi_amp = {}
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
            # self.pi_amp[q] = self.QD_agent.quantum_device.get_element(q).rxy.amp180()
          
    def RunMeasurement(self):
        meas = Rabi12PS()
        meas.set_samples = self.pi_dura_samples
        meas.set_RabiType = "time"
        meas.set_os_mode = self.OSmode
        meas.n_avg = self.avg_n
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset

        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"TimeRabi_12_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QDagent:QDmanager=None,new_QD_path:str=None,new_file_path:str=None):
        """ User callable analysis function pack """
        from qblox_drive_AS.SOP.RabiOsci import conditional_update_qubitInfo
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

            ds = open_dataset(file_path)
            md = None
        
            for var in ds.data_vars:
                if str(var).split("_")[-1] != 'variable':
                    ANA = Multiplex_analyzer("m11")  
                    if ds.attrs['method'].lower() == "shot":
                        md = QD_savior.StateDiscriminator.summon_discriminator(var)    
                    ANA._import_data(ds,1,QD_savior.refIQ[var] if QD_savior.rotate_angle[var][0] == 0 else QD_savior.rotate_angle[var])
                    ANA._start_analysis(var_name=var, OSmodel=md)
                    ANA._export_result(fig_path)
                    conditional_update_qubitInfo(QD_savior,ANA.fit_packs,var,fit_12=True)  
                    
            ds.close()
            # QD_savior.QD_keeper(new_QD_dir)
            



    def WorkFlow(self):
    
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()   


