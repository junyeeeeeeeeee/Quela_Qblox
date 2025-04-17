import time, os
from datetime import datetime
from numpy import array, arange, ndarray, linspace, moveaxis, median, std
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset, DataArray
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
import traceback


class GEF_ROFidelityPS(ScheduleConductor):
    def __init__( self ):
        super().__init__()
        self._target_qs:list = []
        self._shots:int = 5001
        self._repeat:int = 1

    @property
    def target_qs( self ):
        return self._target_qs
    @target_qs.setter
    def set_target_qs( self, target_qs:list):
        if not isinstance(target_qs, list):
            raise TypeError("Target_qs must be given as a list !")
        else:
            self._target_qs = target_qs
    @property
    def shots( self ):
        return self._shots
    @shots.setter
    def set_shots( self, shot_num:int):
        self._shots = shot_num
    @property
    def repeat( self ):
        return self._repeat
    @repeat.setter
    def set_repeat(self, repeat_num:int):
        if not isinstance(repeat_num,(int,float)):
            raise TypeError("Arg 'repeat_num' must be a int or float !")
        self._repeat = int(repeat_num)

    def __PulseSchedule__(self,
        meas_qs: list,
        Note:Notebook,
        states:ndarray = array([0,1,2]),
        repetitions:int=1,
        activeReset:bool=False
        ) -> Schedule:
        
        sameple_idx = states.shape[0]
        sched = Schedule("Single shot", repetitions=repetitions)
        for acq_idx in range(sameple_idx):
            align_pulse = sched.add(IdlePulse(4e-9))
            for qubit_idx, q in enumerate(meas_qs):

                reset = sched.add(Reset(q), ref_op=align_pulse)
                
                if acq_idx == 1:
                    pu = sched.add(X(q))
                elif acq_idx == 2:
                    sched.add(X(q))
                    pu = sched.add(DRAGPulse(G_amp=Note.get_12ampFor(q),D_amp=0,phase=0,duration=Note.get_12durationFor(q),port=f"{q}:mw",clock=f"{q}.12"))
                
            sched.add(Measure(*meas_qs,  acq_index=acq_idx, acq_protocol='SSBIntegrationComplex', bin_mode=BinMode.APPEND), ref_op=pu if acq_idx in [1,2] else reset)
            
        self.schedule = sched
        return sched
    
    def __SetParameters__(self, *args, **kwargs):
        for q in self._target_qs:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
            eyeson_print(f"{q} Integration time: {round(qubit_info.measure.integration_time()*1e6,1)} µs")
        
        self.__Para_state = ManualParameter(name="state")
        self.__Para_state.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__one_shot_para =  ManualParameter(name="Shot")

        self.__repeat_data_idx = arange(self._repeat)
        self.__state_idx = array([0, 1, 2])
        
        self.QD_agent = check_acq_channels(self.QD_agent, self._target_qs)
        
        self.__sched_kwargs = dict(   
            Note=self.QD_agent.Notewriter,
            meas_qs=self._target_qs,
            states=self.__state_idx
        )
    

    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(self._target_qs),
                )
            
            
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._shots)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables([self.__Para_state,self.__one_shot_para,self.__Para_repeat])
            self.meas_ctrl.setpoints_grid((self.__state_idx,arange(self._shots),self.__repeat_data_idx))
    
    def __RunAndGet__(self, *args, **kwargs):
        if self._execution:
            ds = self.meas_ctrl.run('Readout_Fidelity')
            dict_ = {}
            for q_idx, q in enumerate(self._target_qs):
                i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._shots,3)
                q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._shots,3)
                raw_data = moveaxis(moveaxis(array([i_data,q_data]),0,1),2,-1)  # (mixer, repeat, index, prepared_state) -> (repeat, mixer, prepared_state, index)
                dict_[q] = (["repeat","mixer","prepared_state","index"],raw_data)
        
            dataset = Dataset(dict_,coords={"repeat":self.__repeat_data_idx,"mixer":array(["I","Q"]),"prepared_state":array([0, 1, 2]),"index":arange(self._shots)})
                
                
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["method"] = "Shot"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)



class GEF_ROFidelity(ExpGovernment):
    """ Helps you get the **Dressed** cavities. """
    def __init__(self,QD_path:str,data_folder:str=None,JOBID:str=None):
        super().__init__()
        self.QD_path = QD_path
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID

    @property
    def RawDataPath(self):
        return self.__raw_data_location
    
    @RawDataPath.setter
    def RawDataPath(self,path:str):
        self.__raw_data_location = path

    def SetParameters(self, target_qs:list, histo_counts:int=1, shots:int=10000, execution:bool=True):
        """ 
        ### Args:\n
        * target_qs: list, like ["q0", "q1", ...]
        """
        self.use_time_label:bool = False
        self.avg_n = shots
        self.execution = execution
        self.target_qs = target_qs
        self.histos = histo_counts
        self.counter:int = 0
        if self.histos > 1:
            self.use_time_label = True



    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and driving atte
        for q in self.target_qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
    
    def RunMeasurement(self):
        
        meas = GEF_ROFidelityPS()
        meas.set_target_qs = self.target_qs
        meas.set_shots = self.avg_n
        meas.set_repeat = self.histos
        
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset

        
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"GEF_SingleShot_{datetime.now().strftime('%Y%m%d%H%M%S') if (self.JOBID is None or self.use_time_label) else self.JOBID}")
                self.__raw_data_location = self.save_path + ".nc"
                dataset.to_netcdf(self.__raw_data_location)
            else:
                self.save_fig_path = None
        
    def CloseMeasurement(self):
        shut_down(self.cluster,self.Fctrl)


    def RunAnalysis(self,new_QD_path:str=None,new_file_path:str=None,new_QDagent:QDmanager=None,new_pic_save_place:str=None):
        """ if histo_ana, it will check all the data in the same folder with the given new_file_path """

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

            if new_pic_save_place is not None:
                fig_path = new_pic_save_place

            if new_QDagent is None:
                QD_savior = QDmanager(QD_file)
                QD_savior.QD_loader()
            else:
                QD_savior = new_QDagent

            parent_dir = os.path.dirname(file_path)  # Get the parent directory
            date_part = os.path.basename(os.path.dirname(parent_dir))  # "20250122"
            time_part = os.path.basename(parent_dir) 
            ds = open_dataset(file_path)
            
            if self.histos == 1:
                for var in ds.data_vars:
                    try:
                        self.ANA = Multiplex_analyzer("m14")
                        pic_path = os.path.join(fig_path,f"{var}_GEF_SingleShot_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
                        self.ANA._import_data(ds[var]*1000,var_dimension=0,fq_Hz=QD_savior.quantum_device.get_element(var).clock_freqs.f01())
                        self.ANA._start_analysis(var_name=var)
                        if self.save_pics:
                            self.ANA._export_result(pic_path)

                        
                        self.sum_dict[var] = self.ANA.fit_packs

                    
                    except BaseException as err:
                        print(f"Get error while analyze your one-shot data: {err}")
                        traceback.print_exc()
                        eyeson_print("Trying to plot the raw data now... ")
                        self.ANA = Multiplex_analyzer("m14b")
                        self.ANA._import_data(ds[var]*1000,var_dimension=0)
                        pic_path = os.path.join(fig_path,f"{var}_GEF_SingleShot_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}.png")
                        if self.save_pics:
                            self.ANA._export_result(pic_path)
                        
            

                ds.close()
                if self.keep_QD:
                    QD_savior.QD_keeper()

            else:
                for var in ds.data_vars:
                    self.ANA = Multiplex_analyzer("m14")
                    self.ANA._import_data(ds[var]*1000,var_dimension=0,fq_Hz=QD_savior.quantum_device.get_element(var).clock_freqs.f01())
                    self.ANA._start_analysis(var_name=var)
                
                    highlight_print(f"{var}: {round(median(array(self.ANA.fit_packs['effT_mK'])),1)} +/- {round(std(array(self.ANA.fit_packs['effT_mK'])),1)} mK")
                    Data_manager().save_histo_pic(QD_savior,array(self.ANA.fit_packs["effT_mK"]),var,mode="ss",pic_folder=fig_path)
                    Data_manager().save_histo_pic(QD_savior,array(self.ANA.fit_packs["thermal_population"])*100,var,mode="pop",pic_folder=fig_path)

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement() 

if __name__ == "__main__":
    EXP = GEF_ROFidelity("")
    EXP.execution = True
    EXP.histos = 1
    EXP.RunAnalysis(new_QD_path="qblox_drive_AS/QD_backup/20250411/DR4#81_SumInfo.pkl", new_file_path="qblox_drive_AS/Meas_raw/20250411/H17M21S09/GEF_SingleShot_20250411172138.nc")

    