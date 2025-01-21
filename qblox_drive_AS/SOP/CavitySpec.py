"""
Use the results from m1 and a light attenuation (10 ~ 16 is recommended) to find the BARE cavity frequency.\n
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from xarray import Dataset
import matplotlib.pyplot as plt
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from numpy import array, arange, real, imag, arctan2
from qcat.analysis.resonator.photon_dep.res_data import ResonatorData
from qblox_drive_AS.support import Data_manager, QDmanager, compose_para_for_multiplexing 
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Readout, Multi_Readout, Integration, pulse_preview
from quantify_scheduler.operations.gate_library import Reset
from quantify_scheduler.operations.pulse_library import IdlePulse,SetClockFrequency
from quantify_scheduler.resources import ClockResource


def QD_RO_init(QD_agent:QDmanager, ro_elements:dict):
    for target_q in list(ro_elements.keys()):
        qubit = QD_agent.quantum_device.get_element(target_q)
        qubit.measure.pulse_duration(100e-6)
        qubit.measure.integration_time(100e-6)


def multiplexing_CS_ana(QD_agent:QDmanager, ds:Dataset, save_pic_folder:str=None)->dict:
    """
    # Return\n
    A dict sorted by q_name with its fit results.\n
    Ex. {'q0':{..}, ...}\n
    ----------------------------
    # fit results key names: \n
    ['Qi_dia_corr', 'Qi_no_corr', 'absQc', 'Qc_dia_corr', 'Ql', 'fr', 'theta0', 'phi0', 'phi0_err', 'Ql_err', 'absQc_err', 'fr_err', 'chi_square', 'Qi_no_corr_err', 'Qi_dia_corr_err', 'A', 'alpha', 'delay', 'input_power']
    """
    fit_results = {}
    for idx, q in enumerate(ds.data_vars):
        if str(q).split("_")[-1] != "freq":
            S21 = array(ds[q])[0] + array(ds[q])[1]*1j
            freq = array(ds[f"{q}_freq"])[0][5:]
            res_er = ResonatorData(freq=freq,zdata=array(S21)[5:])
            result, data2plot, fit2plot = res_er.fit()
            fig, ax = plt.subplots(2,2,figsize=(12,12))
            ax0:plt.Axes = ax[0][0] 
            ax0.grid()       
            ax0.plot(freq,result['A']*abs(data2plot))
            ax0.plot(freq,result['A']*abs(fit2plot),c="red",label='fitting')
            ax0.vlines(result['fr'],result['A']*min(data2plot),result['A']*max(data2plot),linestyles="--")
            ax0.set_title(f"{q} cavity @ {round(float(result['fr'])*1e-9,5)} GHz")
            ax0.legend()
            ax1:plt.Axes = ax[0][1] 
            ax1.grid()       
            ax1.plot(freq,arctan2(imag(data2plot),real(data2plot)))
            ax1.plot(freq,arctan2(imag(fit2plot),real(fit2plot)),c="red",label='fitting')
            ax1.set_title("Phase")
            ax1.legend()
            ax2:plt.Axes = ax[1][0]  
            ax2.grid()      
            ax2.scatter(real(array(S21)[1:]),imag(array(S21)[1:]),label='data')
            ax2.set_title("S21 raw data")
            ax2.legend()
            ax3:plt.Axes = ax[1][1]  
            ax3.grid()      
            ax3.scatter(result['A']*real(data2plot),result['A']*imag(data2plot),label='data')
            ax3.scatter(result['A']*real(fit2plot),result['A']*imag(fit2plot),label='fit',c='red',s=10)
            ax3.set_title("S21 after fit")
            ax3.legend()
            plt.tight_layout()
            if save_pic_folder is not None:
                Data_manager().save_multiplex_pics(QD_agent,q,'cs',fig,save_pic_folder)
                plt.close()
            else:
                plt.show()
                
            fit_results[q] = result

    return fit_results


def CS_ana(QD_agent:QDmanager, cs_ds:Dataset, pic_save_folder:str=None, keep_bare:bool=True):
    Quality_values = ["Qi_dia_corr", "Qc_dia_corr", "Ql"]
    Quality_errors = ["Qi_dia_corr_err", "absQc_err", "Ql_err"]
    CS_results = multiplexing_CS_ana(QD_agent,cs_ds, Data_manager().get_today_picFolder() if pic_save_folder is None else pic_save_folder)
    for qubit in CS_results:
        qu = QD_agent.quantum_device.get_element(qubit)
        qu.clock_freqs.readout(float(CS_results[qubit]['fr']))
        if keep_bare:
            QD_agent.Notewriter.save_bareFreq_for(target_q=qubit,bare_freq=CS_results[qubit]['fr'])
        print(f"{qubit}:")
        print("Res @ ",round(qu.clock_freqs.readout()*1e-9,4)," GHz")
        for Qua_idx, Qua in enumerate(Quality_values):
            print(f"{Qua[:2]} = {round(float(CS_results[qubit][Qua])/1000,2)} 土 {round(float(CS_results[qubit][Quality_errors[Qua_idx]])/1000,2)} k")


class CavitySearch(ScheduleConductor):
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
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        repetitions:int=1,    
    ) -> Schedule:
        
        qubits2read = list(frequencies.keys())
        sameple_idx = array(frequencies[qubits2read[0]]).shape[0]
        sched = Schedule("One tone multi-spectroscopy (NCO sweep)",repetitions=repetitions)


        for acq_idx in range(sameple_idx):    

            for qubit_idx, q in enumerate(qubits2read):
                freq = frequencies[q][acq_idx]
                if acq_idx == 0:
                    sched.add_resource(ClockResource(name=q+ ".ro", freq=array(frequencies[q]).flat[0]))
                
                sched.add(Reset(q))
                sched.add(SetClockFrequency(clock=q+ ".ro", clock_freq_new=freq))
                sched.add(IdlePulse(duration=4e-9), label=f"buffer {qubit_idx} {acq_idx}")

                
                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
                
                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_index=acq_idx,acq_channel=qubit_idx,single_shot=False,get_trace=False,trace_recordlength=0)
        
        self.schedule =  sched  
        return sched
        
    def __SetParameters__(self, *args, **kwargs):
         
        self.__datapoint_idx = arange(0,len(list(list(self._ro_elements.values())[0])))

        quantum_device = self.QD_agent.quantum_device
        for q in self._ro_elements:
            quantum_device.get_element(q).clock_freqs.readout(NaN) # avoid cluster clock warning

        self.__freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
        self.__freq.batched = True

        self.__spec_sched_kwargs = dict(   
        frequencies=self._ro_elements,
        R_amp=compose_para_for_multiplexing(self.QD_agent,self._ro_elements,'r1'),
        R_duration=compose_para_for_multiplexing(self.QD_agent,self._ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(self.QD_agent,self._ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(self.QD_agent,self._ro_elements,'r2')
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
            self.meas_ctrl.settables(self.__freq)
            self.meas_ctrl.setpoints(self.__datapoint_idx)
        
        else:
            n_s = 2
            preview_para = {}
            for q in self._ro_elements:
                preview_para[q] = self._ro_elements[q][:n_s]
        
            self.__spec_sched_kwargs['frequencies']= preview_para
        

    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            rs_ds = self.meas_ctrl.run("One-tone")
            dict_ = {}
            for q_idx, q in enumerate(list(self._ro_elements.keys())):
                i_data = array(rs_ds[f'y{2*q_idx}'])
                q_data = array(rs_ds[f'y{2*q_idx+1}'])
                dict_[q] = (["mixer","freq"],array([i_data,q_data]))
                dict_[f'{q}_freq'] = (["mixer","freq"],array([self._ro_elements[q],self._ro_elements[q]]))
            
            self.dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"freq":self.__datapoint_idx})
        
        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__spec_sched_kwargs)


if __name__ == "__main__":
    # meas = CavitySearch()
    # ps = meas.get_adjsutable_paras(display=True)
    QD_agent = QDmanager("qblox_drive_AS/QD_backup/20250120/DR1#11_SumInfo.pkl")
    QD_agent.QD_loader()
    from xarray import open_dataset
    ds = open_dataset("qblox_drive_AS/Meas_raw/20250120/H14M32S19/zoomCS_20250120143243.nc")

    multiplexing_CS_ana(QD_agent,ds)