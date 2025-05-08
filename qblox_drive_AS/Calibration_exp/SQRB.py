import os, time
from datetime import datetime
from typing import Iterable
from numpy import asarray, ndarray, linspace, array, arange
from quantify_scheduler.backends.qblox.constants import MIN_TIME_BETWEEN_OPERATIONS
from qblox_drive_AS.SQRB_utils.pycqed_randomized_benchmarking.randomized_benchmarking import randomized_benchmarking_sequence
from qblox_drive_AS.SQRB_utils.pycqed_randomized_benchmarking.two_qubit_clifford_group import SingleQubitClifford, common_cliffords
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Measure, Reset, X, X90, Y, Y90, Rxy, CZ, IdlePulse
from quantify_scheduler.gettables import ScheduleGettable
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support import Data_manager, check_OS_model_ready, init_meas, coupler_zctrl, set_LO_frequency, init_system_atte, shut_down, check_acq_channels
from qblox_drive_AS.support.Pulse_schedule_library import  pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.ExpFrames import ExpGovernment
from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
from qblox_drive_AS.analysis.SQRBana import SQRB_ana

class SQRBPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._target_qs:Iterable = ["q0"]
        self._GateNum_samples:ndarray = []
        self._circuits_num:int = 10


    @property
    def target_qs( self ):
        return self._target_qs
    @target_qs.setter
    def target_qs(self, qs:Iterable):
        if not isinstance(qs, Iterable):
            raise TypeError("q must be an iterable type !")
        self._target_qs = qs
    
    @property
    def GateeNum_samples( self ):
        return self._GateNum_samples
    @GateeNum_samples.setter
    def set_gate_samples(self, gates_samples:ndarray):
        if not isinstance(gates_samples, ndarray):
            raise TypeError("gates_samples must be a ndarray !")
        
        self._GateNum_samples = gates_samples
    @property
    def randoms(self):
        return self._circuits_num
    @randoms.setter
    def set_random_times(self, num:int):
        if not isinstance(num,int):
            raise TypeError("Random number must be a int !")

    def __PulseSchedule__(self,
        qubit_names: Iterable[str],
        lengths: Iterable[int],
        desired_net_clifford_index: int | None = common_cliffords["I"],
        seed: int | None = None,
        repetitions: int = 1,
    ) -> Schedule:
        """
        Generate a randomized benchmarking schedule.

        All Clifford gates in the schedule are decomposed into products
        of the following unitary operations:

            {'CZ', 'I', 'Rx(pi)', 'Rx(pi/2)', 'Ry(pi)', 'Ry(pi/2)', 'Rx(-pi/2)', 'Ry(-pi/2)'}

        Parameters
        ----------
        qubit_specifier
            String or iterable of strings specifying which qubits to conduct the
            experiment on. If one name is specified, then single qubit randomized
            benchmarking is performed. If two names are specified, then two-qubit
            randomized benchmarking is performed.
        lengths
            Array of non-negative integers specifying how many Cliffords
            to apply before each recovery and measurement. If lengths is of size M
            then there will be M recoveries and M measurements in the schedule.
        desired_net_clifford_index
            Optional index specifying what the net Clifford gate should be. If None
            is specified, then no recovery Clifford is calculated. The default index
            is 0, which corresponds to the identity gate. For a map of common Clifford
            gates to Clifford indices, please see: two_qubit_clifford_group.common_cliffords
        seed
            Optional random seed to use for all lengths m. If the seed is None,
            then a new seed will be used for each length m.
        repetitions
            Optional positive integer specifying the amount of times the
            Schedule will be repeated. This corresponds to the number of averages
            for each measurement.

        """

        # ---- PycQED mappings ----#
        pycqed_qubit_map = {name: name for name in qubit_names}
        
        pycqed_operation_map = {
            "X180": lambda q: X(pycqed_qubit_map[q]),
            "X90": lambda q: X90(pycqed_qubit_map[q]),
            "Y180": lambda q: Y(pycqed_qubit_map[q]),
            "Y90": lambda q: Y90(pycqed_qubit_map[q]),
            "mX90": lambda q: Rxy(qubit=pycqed_qubit_map[q], phi=0.0, theta=-90.0),
            "mY90": lambda q: Rxy(qubit=pycqed_qubit_map[q], phi=90.0, theta=-90.0),
        }

        # ---- Build RB schedule ----#
        sched = Schedule(
            "Randomized benchmarking on " + " and ".join(qubit_names), repetitions=repetitions
        )
        

        # two-qubit RB needs buffer time for phase corrections on drive lines
        operation_buffer_time = 0.0

        for idx, m in enumerate(lengths):
            reset = sched.add(Reset(*qubit_names))
            if m > 0:
                # m-sized random sample of representatives in the quotient group C_n / U(1)
                # where C_n is the n-qubit Clifford group and U(1) is the circle group
                rb_sequence_m: list[int] = randomized_benchmarking_sequence(
                    m, number_of_qubits=1, seed=seed, desired_net_cl=desired_net_clifford_index
                )
                for q in qubit_names:   
                    align_reset = 1
                    for clifford_gate_idx in rb_sequence_m:
                        cl_decomp = SingleQubitClifford(clifford_gate_idx).gate_decomposition
                        for gate, _ in  cl_decomp:
                            if gate != "I":
                                sched.add(pycqed_operation_map[gate](q), rel_time=operation_buffer_time, ref_op=reset if align_reset else None)
                        align_reset = 0
                                

                    sched.add(Measure(q, acq_index=idx))
            else:
                sched.add(Measure(*qubit_names, acq_index=idx))

        return sched
        

    def __SetParameters__(self, *args, **kwargs):
        
        for q in self._target_qs:
            qubit_info:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
            eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        self.__circuits = arange(self._circuits_num)
        self.__length = ManualParameter(name="length", unit="#", label="Number of Clifford gates")
        self.__length.batched = True
        self.__length.batch_size = 10
        self.__random_circuit = ManualParameter(name="Circuit", unit="#", label="randoms")
        self.__random_circuit.batched = False
        self.QD_agent = check_acq_channels(self.QD_agent, self._target_qs)
        self.__sched_kwargs = {"qubit_names": self._target_qs, "lengths": self.__length, "seed":self.__random_circuit}
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(self._target_qs)),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            self.meas_ctrl.settables([self.__length, self.__random_circuit])
            self.meas_ctrl.setpoints_grid([self._GateNum_samples, self.__circuits])
            
        else:
            self.__sched_kwargs['lengths']= array([self._GateNum_samples[0], self._GateNum_samples[1], self._GateNum_samples[-1]])
            self.__sched_kwargs['seed'] = self._circuits_num
    
    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run('SQRB')
            dict_ = {}
            
            for q_idx, q in enumerate(self._target_qs):
                i_data = ds[f'y{2*q_idx}'].values.reshape(self._circuits_num, self._GateNum_samples.shape[0])
                q_data = ds[f'y{2*q_idx+1}'].values.reshape(self._circuits_num, self._GateNum_samples.shape[0])
                dict_[q] = (["mixer","random_circuits","gate_length"],array([i_data,q_data]))
            dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"random_circuits":self.__circuits,"gate_length":self._GateNum_samples})
                
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)

class SQRB(ExpGovernment):
    
    def __init__(self,data_folder:str=None,JOBID:str=None):
        
        super().__init__()
        self.QD_path:str = ""
        self.save_dir = data_folder
        self.__raw_data_location:str = ""
        self.JOBID = JOBID
        self.avg_n:int = 300
        self.execution:bool = True
        self.qs:list = ["q0"]
        self.max_gate_num:int = 100
        self.circuits_num:int=10
        self.gate_pts:int = 50

    @property
    def RawDataPath(self):
        return self.__raw_data_location

    def SetParameters(self)->None:
        self.gate_length = linspace(0, self.max_gate_num, self.gate_pts)
        
        
        
    def PrepareHardware(self):
        self.QD_agent, self.cluster, self.meas_ctrl, self.ic, self.Fctrl = init_meas(QuantumDevice_path=self.QD_path)
        
        # bias coupler
        self.Fctrl = coupler_zctrl(self.Fctrl,self.QD_agent.Fluxmanager.build_Cctrl_instructions([cp for cp in self.Fctrl if cp[0]=='c' or cp[:2]=='qc'],'i'))
        # offset bias, LO and atte
        for q in self.qs:
            self.Fctrl[q](self.QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q))
            IF_minus = self.QD_agent.Notewriter.get_xyIFFor(q)
            xyf = self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()
            set_LO_frequency(self.QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf-IF_minus)
            init_system_atte(self.QD_agent.quantum_device,[q],ro_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q, 'ro'), xy_out_att=self.QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        

    def RunMeasurement(self):
        meas = SQRBPS()
        meas.target_qs = self.qs
        meas.set_random_times = self.circuits_num
        meas.set_gate_samples = self.gate_length
        meas.n_avg = self.avg_n
        meas.meas_ctrl = self.meas_ctrl
        meas.QD_agent = self.QD_agent
        meas.execution = self.execution
        
        meas.run()
        dataset = meas.dataset
        # dataset = Ramsey(self.QD_agent,self.meas_ctrl,self.time_samples,self.spin_num,self.histos,self.avg_n,self.execution)
        if self.execution:
            if self.save_dir is not None:
                self.save_path = os.path.join(self.save_dir,f"SQRB_{datetime.now().strftime('%Y%m%d%H%M%S') if self.JOBID is None else self.JOBID}")
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
            SQRB_ana(ds, QD_savior.rotate_angle, fig_path, QD_savior.quantum_device)
            

    def WorkFlow(self):
        
        self.PrepareHardware()

        self.RunMeasurement()

        self.CloseMeasurement()