"""
Use the results from m1 and a light attenuation (10 ~ 16 is recommended) to find the BARE cavity frequency.\n
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from xarray import Dataset
import matplotlib.pyplot as plt
from Modularize.support import cds
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from numpy import array, linspace, arange, cos, sin, deg2rad, real, imag, sqrt
from quantify_core.measurement.control import MeasurementControl
from qcat.analysis.resonator.photon_dep.res_data import ResonatorData
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support import Data_manager, QDmanager, compose_para_for_multiplexing
from Modularize.support.Pulse_schedule_library import One_tone_multi_sche, pulse_preview


def Cavity_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,n_avg:int=10,run:bool=True,Experi_info:dict={},particular_folder:str="",ro_amps:dict={})->Dataset:
    """
        Doing the multiplexing cavity search according to the arg `ro_elements`\n
        Please fill up the initial value about measure for qubit in QuantumDevice first, like: amp, duration, integration_time and acqusition_delay!\n
        ----
        ### Args:\n
        * ro_elements: {'q0': data point array, 'q1':[], ...}\n
        ----
        ## Warning:\n
        The sweep frequency data-point for each cavity in ro_elements is better to set equally.
    """
    datapoint_idx = arange(0,len(list(list(ro_elements.values())[0])))

    quantum_device = QD_agent.quantum_device
    sche_func = One_tone_multi_sche
    for q in ro_elements:
        quantum_device.get_element(q).clock_freqs.readout(NaN) # avoid cluster clock warning
    

    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True

    amps = {}
    if ro_amps != {}:
        init = compose_para_for_multiplexing(QD_agent,ro_elements,1)
        for q in ro_elements:
            if q in list(ro_amps.keys()):
                amps[q] = ro_amps[q]
            else:
                amps[q] = init[q]
    else:
        amps = compose_para_for_multiplexing(QD_agent,ro_elements,1)


    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        R_amp=amps,
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,3),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,4),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,2),
        powerDep=False,
    )
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(datapoint_idx)
        
        rs_ds = meas_ctrl.run("One-tone")
        for idx, q in enumerate(list(ro_elements.keys())):
            rs_ds[f'x{idx}'] = ro_elements[q]

        Data_manager().save_raw_data(QD_agent=QD_agent,ds=rs_ds,qb=q,exp_type='CS',specific_dataFolder=particular_folder)
        
        print(f"{q} Cavity:")
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s = 2
        preview_para = {}
        for q in ro_elements:
            preview_para[q] = ro_elements[q][:n_s]
        
        spec_sched_kwargs['frequencies']= preview_para
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
    
        if Experi_info != {}:
            show_args(Experi_info(q))

    return rs_ds

def QD_RO_init(QD_agent:QDmanager, ro_elements:dict):
    for target_q in list(ro_elements.keys()):
        qubit = QD_agent.quantum_device.get_element(target_q)
        qubit.reset.duration(250e-6)
        qubit.measure.acq_delay(0)
        qubit.measure.pulse_amp(0.15)
        qubit.measure.pulse_duration(100e-6)
        qubit.measure.integration_time(100e-6-4e-9)


def multiplexing_CS_ana(QD_agent:QDmanager, q_label_inorder:list, ds:Dataset, save_pic:bool=True)->dict:
    """
    # Return\n
    A dict sorted by q_name with its fit results.\n
    Ex. {'q0':{..}, ...}\n
    ----------------------------
    # fit results key names: \n
    ['Qi_dia_corr', 'Qi_no_corr', 'absQc', 'Qc_dia_corr', 'Ql', 'fr', 'theta0', 'phi0', 'phi0_err', 'Ql_err', 'absQc_err', 'fr_err', 'chi_square', 'Qi_no_corr_err', 'Qi_dia_corr_err', 'A', 'alpha', 'delay', 'input_power']
    """
    fit_results = {}
    for idx, q in enumerate(q_label_inorder):
        S21 = ds[f"y{2*idx}"] * cos(
                deg2rad(ds[f"y{2*idx+1}"])
            ) + 1j * ds[f"y{2*idx}"] * sin(deg2rad(ds[f"y{2*idx+1}"]))
        coord = f'x{idx}'
        freq = array(ds[coord])[5:]
        res_er = ResonatorData(freq=freq,zdata=array(S21)[5:])
        result, data2plot, fit2plot = res_er.fit()
        fig, ax = plt.subplots(1,2,figsize=(9,6))
        ax0:plt.Axes = ax[0]        
        ax0.plot(freq,abs(data2plot))
        ax0.plot(freq,abs(fit2plot),c="red",label='fitting')
        # ax.vlines(float(result['fr']),min(data2plot),max(data2plot),linestyles="--")
        ax0.set_title(f"{q} cavity @ {round(float(result['fr'])*1e-9,5)} GHz")
        ax0.legend()
        ax1:plt.Axes = ax[1]        
        ax1.scatter(real(data2plot),imag(data2plot),label='data')
        ax1.scatter(real(fit2plot),imag(fit2plot),label='fit',c='red',s=10)
        ax1.set_title("S21")
        ax1.legend()
        plt.tight_layout()
        if save_pic:
            Data_manager().save_multiplex_pics(QD_agent, q, 'CS', fig)
        else:
            plt.show()
            plt.close()
        fit_results[q] = result

    return fit_results

# execution pack
def cavitySpectro_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_bare_guess:dict,ro_span_Hz:float=10e6,run:bool=True,fpts:int=101,avg_times:int=10)->dict:
    ro_elements = {}
    Quality_values = ["Qi_dia_corr", "Qc_dia_corr", "Ql"]
    Quality_errors = ["Qi_dia_corr_err", "absQc_err", "Ql_err"]
    for qb in list(ro_bare_guess.keys()):
        ro_elements[qb] = linspace(ro_bare_guess[qb]-ro_span_Hz, ro_bare_guess[qb]+ro_span_Hz, fpts)
    if run:
        cs_ds = Cavity_spec(QD_agent,meas_ctrl,ro_elements,n_avg=avg_times)
        CS_results = multiplexing_CS_ana(QD_agent, list(ro_elements.keys()),cs_ds, ro_elements)
        for qubit in CS_results:
            qu = QD_agent.quantum_device.get_element(qubit)
            if qubit != 'q0':
                qu.clock_freqs.readout(float(CS_results[qubit]['fr']))
            else:
                qu.clock_freqs.readout(5.434e9)
            print(f"{qubit}:")
            print("Res @ ",round(qu.clock_freqs.readout()*1e-9,4)," GHz")
            for Qua_idx, Qua in enumerate(Quality_values):
                print(f"{Qua[:2]} = {round(float(CS_results[qubit][Qua])/1000,2)} 土 {round(float(CS_results[qubit][Quality_errors[Qua_idx]])/1000,2)} k")
            
            

    else:
        # For pulse preview
        cs_ds = Cavity_spec(QD_agent,meas_ctrl,ro_elements,run=False)
        CS_results = {}
        
    return CS_results


if __name__ == "__main__":
    
    """ fill in part """
    # Basic info of measuring instrument, chip
    # e.g. QD_path, dr, ip, mode, chip_name, chip_type = '', 'dr3', '13', 'n','20240430_8_5Q4C', '5Q4C'
    QD_path, dr, mode, chip_name, chip_type = '', 'drke', 'n','WJ3CQ3FQ_aivon', '5Q4C'
    execution:bool = 1
    chip_info_restore:bool = 0
    # RO attenuation
    init_RO_DigiAtte = 20 # multiple of 2, 10 ~ 16 recommended

    ro_bare=dict(
        q0=5.434e9,
        q1=5.9757e9,
        q2=6.0299e9,
        q3=5.8765e9,
        q4=5.9268e9,
    )

    """ Optional paras """ 
    coupler_number:int = 0
    qubit_num:int = 5
    freq_data_points = 401
    half_freq_window_Hz = 20e6
    n_avg: int = 100



    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,
                                                        dr_loc=dr,
                                                        mode=mode,
                                                        chip_name=chip_name,
                                                        chip_type=chip_type,
                                                        qubit_number=qubit_num,
                                                        coupler_number=coupler_number)
    # Create or Load chip information
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    # Set the system attenuations
    if QD_path == '': QD_RO_init(QD_agent,ro_bare)
    init_system_atte(QD_agent.quantum_device,list(ro_bare.keys()),ro_out_att=init_RO_DigiAtte)
    
    """ Measurements """
    CS_results = cavitySpectro_executor(QD_agent=QD_agent,meas_ctrl=meas_ctrl,ro_bare_guess=ro_bare,run = execution,ro_span_Hz=half_freq_window_Hz,fpts=freq_data_points,avg_times=n_avg)
    
    
    """ Storing """
    if execution:
        QD_agent.refresh_log("After cavity search")
        QD_agent.QD_keeper()
        print('CavitySpectro done!')
        
        # Chip info!
        if chip_info_restore:
            chip_info.update_CavitySpec(result=CS_results)


    """ Close """
    shut_down(cluster,Fctrl)
    

