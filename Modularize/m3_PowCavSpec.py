""" 
Base on a BARE cavity observe a dispersive shift in RO-freq with the variable RO-amp.  
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, linspace, sqrt, moveaxis, arange, cos, sin, deg2rad, real, imag
from xarray import open_dataset, Dataset
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import Data_manager, QDmanager, compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from Modularize.support.Pulse_schedule_library import One_tone_multi_sche, pulse_preview

def PowerDep_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,ro_p_min:float=0.01,ro_p_max:float=0.5,p_points:int=20,n_avg:int=100,run:bool=True,Experi_info:dict={}, rof_marks:dict={})->Dataset:

    sche_func = One_tone_multi_sche
    freq_datapoint_idx = arange(0,len(list(list(ro_elements.values())[0])))
    original_rof = {}
     
    from numpy import NaN
    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        original_rof[q] = qubit_info.clock_freqs.readout()
        # avoid frequency conflicts
        qubit_info.clock_freqs.readout(NaN)

    ro_p_samples = linspace(ro_p_min,ro_p_max,p_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
    ro_pulse_amp.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        R_amp=ro_pulse_amp,
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,3),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,4),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,2),
        powerDep=True,
    )
    # exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
    #                  Power=['start '+'%E' %ro_p_samples[0],'end '+'%E' %ro_p_samples[-1]])
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((freq_datapoint_idx,ro_p_samples)) # -> (x0, x1) if do dh.to_gridded_dataset(ds)
        
        
        
        rp_ds = meas_ctrl.run("One-tone-powerDep")
        # Save the raw data into netCDF
        nc_file = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rp_ds,qb="MultiQ",exp_type='PD',get_data_loc=True)

        # analysis_result[q] = Basic2DAnalysis(tuid=rp_ds.attrs["tuid"], dataset=rp_ds).run()
        # show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        # try: 
        plot_powerCavity_S21(nc_file, ro_elements, rof_marks)
        # except:
        #     from Modularize.support.UserFriend import warning_print
        #     warning_print("Beware! plot_powerCavity_S21 didn't work~")
        
    else:
        n_s = 2
        preview_para = {}
        for q in ro_elements:
            preview_para[q] = ro_elements[q][:n_s]
        sweep_para2 = array(ro_p_samples[:2])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['R_amp']= {q:sweep_para2.reshape(sweep_para2.shape or (1,))[0]}
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)

        # show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    
    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        original_rof[q] = qubit_info.clock_freqs.readout(original_rof[q])

    return rp_ds

def plot_powerCavity_S21(PC_nc_file:str, ro_elements:dict, rof_pos:dict={}):
    """
    Plot |S21| from a given power cavity nc file and save it in the pic folder within the same day.
    """
    import quantify_core.data.handling as dh
    import matplotlib.pyplot as plt
    title = f"{os.path.split(PC_nc_file)[-1].split('.')[0]}"
    pic_dir = os.path.join(os.path.split(PC_nc_file)[0],"PowerCav_pic")
    if not os.path.exists(pic_dir):         
        os.mkdir(pic_dir)
    ds = dh.to_gridded_dataset(open_dataset(PC_nc_file))
    # ds.x0 = freq. ; ds.x1 = power
    power = array(ds.x1)
    for idx, q in enumerate(ro_elements):
        S21 = ds[f"y{2*idx}"] * cos(
                deg2rad(ds[f"y{2*idx+1}"])
            ) + 1j * ds[f"y{2*idx}"] * sin(deg2rad(ds[f"y{2*idx+1}"]))
        I, Q = real(S21), imag(S21)
        amp = moveaxis(array(sqrt(I**2+Q**2)),0,-1) # make shape in (power, freq)
        s21 = []
        for i in range(amp.shape[0]):
            s21.append(list(array(amp[i])/power[i]))
        s21 = array(s21)
        freq = ro_elements[q]
        fig, ax = plt.subplots()
        ax:plt.Axes
        d = ax.pcolormesh(freq*1e-9, power, s21, shading='gouraud',cmap='RdBu')
        ax.vlines(rof_pos[q]*1e-9,ymin=min(power),ymax=max(power),linestyles='--',colors='#FF00FF',label='window_shift_baseline')
        fig.colorbar(d, ax=ax)
        plt.xlabel("frequency (GHz)")
        plt.ylabel("Power (V)")
        plt.minorticks_on()
        plt.title(title+f"_{q}")
        plt.grid()
        plt.tight_layout()
        plt.legend()
        plt.savefig(os.path.join(pic_dir,f"{title}_{q}.png"))
        plt.close()


def powerCavity_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Fctrl:dict,qubits:list,ro_atte:int,ro_span_Hz:float=3e6,max_power:float=0.7,run:bool=True,sweet_spot:bool=False,fpts:int=30,ppts:int=30,avg_n:int=10,dressHbare:bool=True):
    ro_elements = {}
    rof_baselines = {}
    for qb in qubits:
        QD_agent.Notewriter.save_DigiAtte_For(ro_atte,qb,'ro')
        rof_baselines[qb] = QD_agent.quantum_device.get_element(qb).clock_freqs.readout()
        if dressHbare:
            rof = rof_baselines[qb]-1e6
            ro_elements[qb] = linspace(rof, rof+2*ro_span_Hz, fpts)
        else:
            rof = rof_baselines[qb]
            ro_elements[qb] = linspace(rof-ro_span_Hz, rof+ro_span_Hz, fpts)
    init_system_atte(QD_agent.quantum_device,qubits,ro_out_att=ro_atte)
    
    if run:
        if sweet_spot:
            Fctrl[list(ro_elements.keys())[0]](QD_agent.Fluxmanager.get_sweetBiasFor(target_q=list(ro_elements.keys())[0]))
        PD_ds = PowerDep_spec(QD_agent,meas_ctrl,ro_elements, ro_p_max=max_power, p_points=ppts, n_avg=avg_n, rof_marks=rof_baselines)
        if sweet_spot:
            Fctrl[list(ro_elements.keys())[0]](0.0)
    else:
        PD_ds = PowerDep_spec(QD_agent,meas_ctrl,ro_elements,run=False,ro_p_max=max_power)

    

if __name__ == "__main__":
    
    """ fill in """
    execution:bool = True
    sweetSpot_dispersive:bool = 0 # if true, only one qubit should be in the ro_elements 
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements =["q2"]     # measurement target q from this dict # q1, q2 44dB 0.2
    ro_atte_for_all:int=20 

    """ Optional paras"""
    maxima_power = 0.6
    half_ro_freq_window_Hz = 2e6
    freq_data_points = 100
    power_data_points = 30

    """ in Case paras """
    dress_higher_bare:bool = 0


    """ preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')


    """ Running """
    powerCavity_executor(QD_agent,meas_ctrl,Fctrl,qubits=ro_elements,ro_atte=ro_atte_for_all,run=execution,sweet_spot=sweetSpot_dispersive,max_power=maxima_power,ro_span_Hz=half_ro_freq_window_Hz, fpts=freq_data_points, ppts=power_data_points,dressHbare=dress_higher_bare,avg_n=10)
    cluster.reset()
    

    """ Storing """
    if execution: 
        QD_agent.refresh_log('after PowerDep')
        QD_agent.QD_keeper()
    

    """ Close """
    print('Power dependence done!')
    shut_down(cluster,Fctrl)

    