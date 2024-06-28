"""
Use the results from m1 and a light attenuation (10 ~ 16 is recommended) to find the BARE cavity frequency.
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from Modularize.support import uw
from numpy import array, linspace
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

def Cavity_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,n_avg:int=300,run:bool=True,Experi_info:dict={},particular_folder:str="")->dict:
    """
        Do the cavity search by the given QuantumDevice with a given target qubit q. \n
        Please fill up the initial value about measure for qubit in QuantumDevice first, like: amp, duration, integration_time and acqusition_delay! 
    """
    quantum_device = QD_agent.quantum_device
    sche_func = One_tone_multi_sche
    for q in ro_elements:
        quantum_device.get_element(q).clock_freqs.readout(NaN) # avoid cluster clock warning
    analysis_result = {}

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
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(ro_f_samples)
        
        rs_ds = meas_ctrl.run("One-tone")
        analysis_result[q] = ResonatorSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
        # save the xarrry into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=rs_ds,qb=q,exp_type='CS',specific_dataFolder=particular_folder)

        print(f"{q} Cavity:")
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s=2 
        sweep_para= array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result

def QD_RO_init(QD_agent:QDmanager, target_q:str):
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.reset.duration(150e-6)
    qubit.measure.acq_delay(280e-9)
    qubit.measure.pulse_amp(0.15)
    qubit.measure.pulse_duration(1e-6)
    qubit.measure.integration_time(0.5e-6-4e-9)


# execution pack
def cavitySpectro_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_bare_guess:dict,qb:str,ro_span_Hz:float=10e6,run:bool=True):
    if run:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess,q=qb,ro_span_Hz=ro_span_Hz,ro_amp=1)[qb]
        if qb_CSresults != {}:
            print(f'Cavity {qb} @ {qb_CSresults.quantities_of_interest["fr"].nominal_value} Hz')
            QD_agent.quantum_device.get_element(qb).clock_freqs.readout(qb_CSresults.quantities_of_interest["fr"].nominal_value)
        else:
            print(f"Cavity Spectroscopy error qubit: {qb}")

    else:
        # For pulse preview
        for qb in [list(ro_bare_guess.keys())[0]]:
            qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess,q=qb,ro_span_Hz=ro_span_Hz,run=False)
        
    return qb_CSresults


if __name__ == "__main__":
    
    """ fill in part """
    # Basic info of measuring instrument, chip
    # e.g. QD_path, dr, ip, mode, chip_name, chip_type = '', 'dr3', '13', 'n','20240430_8_5Q4C', '5Q4C'
    QD_path, dr, ip, mode, chip_name, chip_type = '', 'dr1sca', '', 'n','20240611_ScalinQ', '5Q4C'
    # 1 = Run the measurement
    # 0 = plot the output signal
    execution:bool = True
    chip_info_restore:bool = 1
    # RO attenuation
    # 10 ~ 16
    init_RO_DigiAtte = 10 # multiple of 2 

    ro_bare=dict(
        q0=5.95388e9   
    )


    """ Preparations """
    
    # Reload the QuantumDevice or build up a new one
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,
                                                        dr_loc=dr,
                                                        cluster_ip=ip,
                                                        mode=mode,
                                                        chip_name=chip_name,
                                                        chip_type=chip_type,
                                                        qubit_number=1,
                                                        coupler_number=1)
    # Create or Load chip information
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    # Set the system attenuations
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=init_RO_DigiAtte)
    
    """ Measurements """
    CS_results = {}
    for qubit in ro_bare:
        if QD_path == '': QD_RO_init(QD_agent,qubit)
        CS_results[qubit] = cavitySpectro_executor(QD_agent=QD_agent,meas_ctrl=meas_ctrl,ro_bare_guess=ro_bare,qb=qubit,run = execution,ro_span_Hz=10e6)
        if not execution:
            break
    
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
    

