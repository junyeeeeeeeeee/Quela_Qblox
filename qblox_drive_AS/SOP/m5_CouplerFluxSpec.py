import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, linspace, pi, arange, sqrt, ndarray
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds, compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from qblox_drive_AS.Configs.ClusterAddress_rec import ip_register
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from qblox_drive_AS.support.Pulse_schedule_library import RabiSplitting_multi_sche, pulse_preview
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.analysis.raw_data_demolisher import fluxCoupler_dataReducer
from xarray import Dataset

z_pulse_amp_OVER_const_z = sqrt(2)/2.5


def fluxCoupler_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,bias_elements:dict,n_avg:int=300,run:bool=True,Experi_info:dict={}):
    from numpy import NaN
    sche_func = RabiSplitting_multi_sche
    freq_datapoint_idx = arange(0,len(list(list(ro_elements.values())[0])))
    original_rof = {}
    flux_dura = 0
    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.measure.pulse_duration(100e-6)
        qubit_info.measure.integration_time(100e-6)
        qubit_info.reset.duration(250e-6)
        flux_dura = qubit_info.reset.duration()+qubit_info.measure.integration_time()
        original_rof[q] = qubit_info.clock_freqs.readout()
        qubit_info.clock_freqs.readout(NaN)

    flux_samples:ndarray = bias_elements["bias_samples"] 
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    bias = ManualParameter(name="bias", unit="V", label="Flux voltage")
    bias.batched = False
    
    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        bias_couplers=bias_elements["target_couplers"],
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
        powerDep=False,
        bias = bias,
        bias_dura = flux_dura
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
        meas_ctrl.settables((freq,bias))
        meas_ctrl.setpoints_grid((freq_datapoint_idx,flux_samples)) # x0, x1
        
        ds = meas_ctrl.run("One-tone-Flux")
        dict_ = {}
        for idx, q in enumerate(ro_elements):   
            freq_values = 2*flux_samples.shape[0]*list(ro_elements[q])
        
            i_data = array(ds[f'y{2*idx}']).reshape(flux_samples.shape[0],array(ro_elements[q]).shape[0])
            q_data = array(ds[f'y{2*idx+1}']).reshape(flux_samples.shape[0],array(ro_elements[q]).shape[0])
            dict_[q] = (["mixer","bias","freq"],array([i_data,q_data]))
            
            dict_[f"{q}_freq"] = (["mixer","bias","freq"],array(freq_values).reshape(2,flux_samples.shape[0],array(ro_elements[q]).shape[0]))

        
        rfs_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"bias":flux_samples/z_pulse_amp_OVER_const_z,"freq":arange(flux_samples.shape[0])})
        rfs_ds.attrs["execution_time"] = Data_manager().get_time_now()
        rfs_ds.attrs["cntrl_couplers"] =  "_".join(bias_elements["target_couplers"])
        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rfs_ds,qb="multiQ",exp_type='FD',get_data_loc=True)
        
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        
    else:
        preview_para = {}
        for q in ro_elements:
            preview_para[q] = ro_elements[q][:2]
        sweep_para2= array([flux_samples[0],flux_samples[-1]])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['bias']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)

        if Experi_info != {}:
            show_args(Experi_info(q))
        nc_path = ''

    for q in ro_elements:
        QD_agent.quantum_device.get_element(q).clock_freqs.readout(original_rof[q])
    
    return nc_path

def update_flux_info_in_results_for(QD_agent:QDmanager,qb:str,FD_results:dict):
    qubit = QD_agent.quantum_device.get_element(qb)
    qubit.clock_freqs.readout(FD_results[qb].quantities_of_interest["freq_0"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=FD_results[qb].quantities_of_interest["offset_0"].nominal_value)
    QD_agent.Fluxmanager.save_period_for(target_q=qb, period=2*pi/FD_results[qb].quantities_of_interest["frequency"].nominal_value)
    QD_agent.Fluxmanager.save_tuneawayBias_for(target_q=qb,mode='auto')
    QD_agent.Fluxmanager.save_cavFittingParas_for(target_q=qb,
        f=FD_results[qb].quantities_of_interest["frequency"].nominal_value,
        amp=FD_results[qb].quantities_of_interest["amplitude"].nominal_value,
        phi=FD_results[qb].quantities_of_interest["shift"].nominal_value,
        offset=FD_results[qb].quantities_of_interest["offset"].nominal_value
    )

def update_coupler_bias(QD_agent:QDmanager,cp_elements:dict):
    """
    Update the idle bias in Fluxmanager for couplers.\n
    --------------------------
    ### Args:\n
    cp_elements = {"c0":0.2}
    """
    for cp in cp_elements:
        QD_agent.Fluxmanager.save_idleBias_for(cp, cp_elements[cp])

def fluxCoupler_waiter(QD_agent:QDmanager,freq_halfspan:dict, freq_shift:dict, fpts:int, target_couplers:list, flux_span:float=0.3, flux_pts:int=30)->tuple[dict,dict]:
    """
    Generate the frequencies samples array with the rule: array = linspace(ROfreq+freq_shift-freq_halfspan, ROfreq+freq_shift+freq_halfspan, fpts)
    * frequency unit in Hz
    """
    # RO freq settings
    ro_elements = {}
    for q in freq_halfspan:
        rof = QD_agent.quantum_device.get_element(q).clock_freqs.readout()
        ro_elements[q] = linspace(rof+freq_shift[q]-freq_halfspan[q], rof+freq_shift[q]+freq_halfspan[q], fpts)
    
    # bias settings
    if flux_span > 0.4: mark_input("Attempting to set flux voltage higher than 0.4 V, press any key to continue...")
    bias_elements = {"bias_samples":linspace(-flux_span,flux_span,flux_pts)*z_pulse_amp_OVER_const_z,"target_couplers":target_couplers}
    return ro_elements, bias_elements




def fluxCoupler_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,bias_elements:dict,run:bool=True,avg_n=20)->str:
    
    nc_path = fluxCoupler_spec(QD_agent,meas_ctrl,ro_elements,bias_elements=bias_elements,run=run,n_avg=avg_n)

    return nc_path

# accident: q2, q3, q4

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = True
    chip_info_restore:bool = 0 # <- haven't been modified with the multiplexing version
    DRandIP = {"dr":"dr2","last_ip":"10"}
    cp_ctrl = ["c0", "c1"]
    freq_half_span = {
        "q0":5e6,
        
    }
    freq_shift = {
        "q0":0e6,
        
    }

    """ Optional paras """
    freq_data_points = 40
    flux_half_window_V  = 0.3
    flux_data_points = 40
    avg_n = 50
    

    
    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    ro_elements, bias_elements = fluxCoupler_waiter(QD_agent, freq_half_span, freq_shift, freq_data_points, cp_ctrl, flux_half_window_V, flux_data_points)
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))


    """ Running """
    update = False
    FD_results = {}
    nc_path = fluxCoupler_executor(QD_agent,meas_ctrl,ro_elements,bias_elements,run=execution, avg_n=avg_n)
    slightly_print(f" Nc saved located:\n{nc_path}")
    cluster.reset()


    """ Analysis """
    if execution:
        ds = fluxCoupler_dataReducer(nc_path)
        for var in ds.data_vars:
            ANA = Multiplex_analyzer("m5")
            if var.split("_")[-1] != 'freq':
                ANA._import_data(ds,2)
                ANA._start_analysis(var_name=var)
                pic_folder = Data_manager().get_today_picFolder()
                ANA._export_result(pic_folder)


    """ Close """
    print('Flux dependence done!')
    shut_down(cluster,Fctrl)