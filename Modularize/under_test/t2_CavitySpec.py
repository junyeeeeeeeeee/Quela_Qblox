"""
Use the results from m1 and a light attenuation (10 ~ 16 is recommended) to find the BARE cavity frequency.\n
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from Modularize.support import uw
from numpy import array, linspace, arange, cos, sin, deg2rad, real, imag, sqrt
from Modularize.support import cds
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import Data_manager, QDmanager, compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support import init_meas, init_system_atte, shut_down, QRM_nco_init
from Modularize.support.Pulse_schedule_library import One_tone_multi_sche, pulse_preview
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
import matplotlib.pyplot as plt
from xarray import Dataset

def Cavity_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,n_avg:int=10,run:bool=True,Experi_info:dict={},particular_folder:str="")->Dataset:
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

    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,1),
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
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=rs_ds,qb=q,exp_type='CS',specific_dataFolder=particular_folder)
        # analysis_result = ResonatorSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
        # save the xarrry into netCDF
        
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
        qubit.reset.duration(150e-6)
        qubit.measure.acq_delay(280e-9)
        qubit.measure.pulse_amp(0.15)
        qubit.measure.pulse_duration(10e-6)
        qubit.measure.integration_time(10e-6-4e-9)


def multiplexing_CS_ana(QD_agent:QDmanager, ds:Dataset, ro_elements:dict):
    for idx, q in enumerate(ro_elements):
        S21 = ds[f"y{2*idx}"] * cos(
                deg2rad(ds[f"y{2*idx+1}"])
            ) + 1j * ds[f"y{2*idx}"] * sin(deg2rad(ds[f"y{2*idx+1}"]))
        I, Q = real(S21), imag(S21)
        amp = sqrt(I**2+Q**2)
        fig, ax = plt.subplots()
        ax:plt.Axes        
        ax.plot(ro_elements[q],amp)
        Data_manager().save_multiplex_pics(QD_agent, q, 'CS', fig)

# execution pack
def cavitySpectro_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_bare_guess:dict,ro_span_Hz:float=10e6,run:bool=True,fpts:int=101):
    ro_elements = {}
    for qb in list(ro_bare_guess.keys()):
        ro_elements[qb] = linspace(ro_bare_guess[qb]-ro_span_Hz, ro_bare_guess[qb]+ro_span_Hz, fpts)
    if run:
        cs_ds = Cavity_spec(QD_agent,meas_ctrl,ro_elements)
        multiplexing_CS_ana(QD_agent, cs_ds, ro_elements)

    else:
        # For pulse preview
        cs_ds = Cavity_spec(QD_agent,meas_ctrl,ro_elements,run=False)
        
    return cs_ds


if __name__ == "__main__":
    
    """ fill in part """
    # Basic info of measuring instrument, chip
    # e.g. QD_path, dr, ip, mode, chip_name, chip_type = '', 'dr3', '13', 'n','20240430_8_5Q4C', '5Q4C'
    QD_path, dr, ip, mode, chip_name, chip_type = '', 'dr3', '', 'n','20240702_AS1606_multi', '5Q4C'
    execution:bool = 1
    chip_info_restore:bool = 0
    # RO attenuation
    init_RO_DigiAtte = 18 # multiple of 2, 10 ~ 16 recommended

    ro_bare=dict(
        q0=5973306404,
        # q1=6083525235,
        # q2=5920038579,
        # q3=6099702007,
        # q4=6010586222
    )

    """ Optional paras """
    freq_data_points = 101
    half_freq_window_Hz = 10e6


    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,
                                                        dr_loc=dr,
                                                        cluster_ip=ip,
                                                        mode=mode,
                                                        chip_name=chip_name,
                                                        chip_type=chip_type,
                                                        qubit_number=5,
                                                        coupler_number=4)
    # Create or Load chip information
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    # Set the system attenuations
    if QD_path == '': QD_RO_init(QD_agent,ro_bare)
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=init_RO_DigiAtte)
    
    """ Measurements """
    CS_results = cavitySpectro_executor(QD_agent=QD_agent,meas_ctrl=meas_ctrl,ro_bare_guess=ro_bare,run = execution,ro_span_Hz=half_freq_window_Hz,fpts=freq_data_points)
    
    
    """ Storing """
    if execution:
        QD_agent.refresh_log("After cavity search")
        # QD_agent.QD_keeper()
        print('CavitySpectro done!')
        
        # Chip info!
        if chip_info_restore:
            chip_info.update_CavitySpec(result=CS_results)

    """ Close """
    shut_down(cluster,Fctrl)
    

