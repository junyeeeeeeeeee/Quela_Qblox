import time, os
from datetime import datetime
from numpy import array, arange, ndarray, linspace, rad2deg, logspace
from numpy import pi as PI
from numpy import nan as NaN
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support.Notebook import Notebook
from qblox_drive_AS.support import Data_manager, check_OS_model_ready, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down, check_acq_channels, sort_dict_with_qidx
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, IdlePulse, Measure, X90, Y90, DRAGPulse, Rxy, X, ConditionalReset, BinMode, ClockResource, SetClockFrequency
from quantify_scheduler.operations.gate_library import Reset
from qblox_drive_AS.support.ExpFrames import ExpGovernment
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer


class ROfreqGEFPS(ScheduleConductor):
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
        frequencies: dict,
        Note:Notebook,
        state:int,
        repetitions:int=1,    
    ) -> Schedule:
        
        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
        sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)

        for acq_idx in range(sameple_idx):    
            align_pulse = sched.add(IdlePulse(4e-9))
            for qubit_idx, q in enumerate(qubits2read):
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
                
                sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
                sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")
                reset = sched.add(Reset(q), ref_op=align_pulse)
                if state == 1:
                    pi_pulse = sched.add(X(q), ref_op=reset)
                elif state == 2:
                    sched.add(X(q))
                    pi_pulse = sched.add(DRAGPulse(G_amp=Note.get_12ampFor(q),D_amp=0,phase=0,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"))
                else:
                    pi_pulse = sched.add(X(q, amp180=0), ref_op=reset)

            sched.add(Measure(*qubits2read,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.AVERAGE), ref_op=pi_pulse if state else reset)
                
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(0,len(list(list(self._ro_elements.values())[0])))
        self.__prepared_states = array([0, 1, 2])
        
        self.__ro_f_origin = {}
        for q in self._ro_elements:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            self.__ro_f_origin[q] = qubit_info.clock_freqs.readout()
            qubit_info.clock_freqs.readout(NaN)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
           
        self.__freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
        self.__freq.batched = True
        self.__state =  ManualParameter(name="prepared_state", unit="", label="state")
        self.__state.batched = False


        self.QD_agent = check_acq_channels(self.QD_agent, list(self._ro_elements.keys()))
        self.__spec_sched_kwargs = dict(
        Note=self.QD_agent.Notewriter,   
        frequencies=self._ro_elements,
        state=self.__state
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
            self.meas_ctrl.settables([self.__freq, self.__state])
            self.meas_ctrl.setpoints_grid([self.__datapoint_idx, self.__prepared_states])
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
        
            self.__spec_sched_kwargs['frequencies']= preview_para
            self.__spec_sched_kwargs['state']= array([1])
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rs_ds = self.meas_ctrl.run("ROF calibration")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                freq_values = 2*self.__prepared_states.shape[0]*list(self._ro_elements[q])
                i_data = array(rs_ds[f'y{2*q_idx}']).reshape(self.__prepared_states.shape[0],array(self._ro_elements[q]).shape[0])
                q_data = array(rs_ds[f'y{2*q_idx+1}']).reshape(self.__prepared_states.shape[0],array(self._ro_elements[q]).shape[0])
                dict_[q] = (["mixer","state","rof"],array([i_data,q_data]))
                dict_[f'{q}_rof'] = (["mixer","state","rof"],array(freq_values).reshape(2,self.__prepared_states.shape[0],array(self._ro_elements[q]).shape[0]))

            ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"state":self.__prepared_states,"rof":self.__datapoint_idx})
            
            ds.attrs["execution_time"] = Data_manager().get_time_now()
            ds.attrs["method"] = "Average"
            ds.attrs["system"] = "qblox"
            for q in self.__ro_f_origin:
                ds.attrs[f"{q}_ori_rof"] = self.__ro_f_origin[q]
            self.dataset = ds
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)


class GEF_ROFcali(ExpGovernment):
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self, freq_span_range:dict, freq_pts:int=100, avg_n:int=100, execution:bool=True, OSmode:bool=False)->None:
        """ ### Args:
            * freq_span_range: {"q0":[freq_span_start, freq_span_end], ....}\n
        """
        self.freq_samples = {}
        self.tempor_freq = [freq_span_range, freq_pts]

        self.avg_n = avg_n
        self.execution = execution
        self.OSmode = OSmode
        self.target_qs = sort_dict_with_qidx(list(freq_span_range.keys()))
        

    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
            rof = self.QD_agent.quantum_device.get_element(q).clock_freqs.readout()
            self.freq_samples[q] = linspace(rof+self.tempor_freq[0][q][0],rof+self.tempor_freq[0][q][1],self.tempor_freq[1])
        
    def RunMeasurement(self):
        
        meas = ROfreqGEFPS()
        meas.ro_elements = self.freq_samples
        meas.execution = self.execution
        meas.n_avg = self.avg_n
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.run()
        dataset = meas.dataset
        
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"GEF_ROFcali_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
                
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None):
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
                if var.split("_")[-1] != 'rof':
                    ANA = Multiplex_analyzer("s5")
                    ANA._import_data(ds,var_dimension=1)
                    ANA._start_analysis(var_name = var)
                    ANA._export_result(fig_path)
                    answer[var] = ANA.fit_packs[var]["optimal_rof"]
            ds.close()

            permi = mark_input(f"What qubit can be updated ? {list(answer.keys())}/ all/ no ").lower()
            if permi in list(answer.keys()):
                QD_savior.quantum_device.get_element(permi).clock_freqs.readout(answer[permi])
                QD_savior.QD_keeper()
            elif permi in ["all",'y','yes']:
                for q in answer:
                    QD_savior.quantum_device.get_element(q).clock_freqs.readout(answer[q])
                QD_savior.QD_keeper()
            else:
                print("Updating got denied ~")


    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

if __name__ == "__main__":
    
    EXP = GEF_ROFcali("")
    EXP.execution = True
    EXP.RunAnalysis("qblox_drive_AS/QD_backup/20250418/DR1#11_SumInfo.pkl","qblox_drive_AS/Meas_raw/20250418/H16M26S32/GEF_ROFcali_20250418162924.nc")