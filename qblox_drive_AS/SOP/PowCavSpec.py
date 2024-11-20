""" 
Base on a BARE cavity observe a dispersive shift in RO-freq with the variable RO-amp.  
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, linspace, sqrt, mean, arange, ndarray
from xarray import open_dataset, Dataset
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager, compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
import matplotlib.pyplot as plt
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down
from qblox_drive_AS.support.Pulse_schedule_library import One_tone_multi_sche, pulse_preview


def PowerDep_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,power_samples:ndarray,n_avg:int=100,run:bool=True)->Dataset:

    sche_func = One_tone_multi_sche
    freq_datapoint_idx = arange(0,len(list(list(ro_elements.values())[0])))
    original_rof = {}
     
    from numpy import NaN
    for q in ro_elements:
        
        qubit_info = QD_agent.quantum_device.get_element(q)
        original_rof[q] = qubit_info.clock_freqs.readout()
        # avoid frequency conflicts
        qubit_info.clock_freqs.readout(NaN)

    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
    ro_pulse_amp.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        R_amp=ro_pulse_amp,
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
        powerDep=True,
    )

    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((freq_datapoint_idx,power_samples)) # -> (x0, x1) if do dh.to_gridded_dataset(ds)
        
        
        
        rp_ds = meas_ctrl.run("One-tone-powerDep")
        dict_ = {}
        for q_idx, q in enumerate(ro_elements):
            i_data = array(rp_ds[f'y{2*q_idx}']).reshape(power_samples.shape[0],ro_elements[q].shape[0])
            q_data = array(rp_ds[f'y{2*q_idx+1}']).reshape(power_samples.shape[0],ro_elements[q].shape[0])
            dict_[q] = (["mixer","ro_amp","freq"],array([i_data,q_data]))
            dict_[f'{q}_freq'] = (["mixer","ro_amp","freq"],array(list(ro_elements[q])*2*power_samples.shape[0]).reshape(2,power_samples.shape[0],ro_elements[q].shape[0]))
            
        DS = Dataset(dict_,coords={"mixer":array(["I","Q"]),"ro_amp":power_samples,"freq":freq_datapoint_idx})

        

        
    else:
        n_s = 2
        preview_para = {}
        for q in ro_elements:
            preview_para[q] = ro_elements[q][:n_s]
        sweep_para2 = array(power_samples[:2])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['R_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        DS = {}
    
    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        original_rof[q] = qubit_info.clock_freqs.readout(original_rof[q])

    return DS

def plot_powerCavity_S21(ds:Dataset, QD_agent:QDmanager=None, save_fig_folder:str=None):
    """
    Plot |S21| from a given power cavity nc file and save it in the pic folder within the same day.
    """
    power = array(ds.coords['ro_amp'])
    # ds.x0 = freq. ; ds.x1 = power
    for idx, q in enumerate(ds.data_vars):
        if str(q).split("_")[-1] != 'freq':
            amp:ndarray = sqrt(array(ds[q])[0]**2+array(ds[q])[1]**2)
            s21 = []
            for i in range(amp.shape[0]):
                if power[i] != 0:
                    s21.append(list(array(amp[i])/power[i]))
                else:
                    s21.append(list(amp[i]))
            s21 = array(s21)
            freq = array(ds[f"{q}_freq"])[0][0]

            fig, ax = plt.subplots()
            ax:plt.Axes
            d = ax.pcolormesh(freq*1e-9, power, s21, shading='gouraud',cmap='RdBu')
           
            ax.vlines(mean(freq)*1e-9 if QD_agent is None else QD_agent.quantum_device.get_element(q).clock_freqs.readout()*1e-9,ymin=min(power),ymax=max(power),linestyles='--',colors='#FF00FF',label='window_shift_baseline')
            fig.colorbar(d, ax=ax)
            plt.xlabel("frequency (GHz)")
            plt.ylabel("Power (V)")
            plt.minorticks_on()
            if QD_agent is not None:
                ro_atte:str = f"_{int(QD_agent.Notewriter.get_DigiAtteFor(q,'ro'))}dB"
            else:
                ro_atte:str = "" 

            plt.title(f"PowerCavity_{q}{ro_atte}")
            plt.grid()
            plt.tight_layout()
            plt.legend()

            if save_fig_folder is not None:
                plt.savefig(os.path.join(save_fig_folder,f"PowerCavity_{q}{ro_atte}.png"))
                plt.close()
            else:
                plt.show()


    
