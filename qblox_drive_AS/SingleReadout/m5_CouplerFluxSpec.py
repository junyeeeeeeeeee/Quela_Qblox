"""
### This program focus on a readout freq. and tune the flux on the given coupler.\n
We expects to see 2 rabi-splitting poles.\n
Once you get the bias about these poles you calculate the max coupler freq bias and use it in m6 program.
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from numpy import array, linspace, pi
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from qblox_drive_AS.support.Pulse_schedule_library import One_tone_sche, pulse_preview
from qblox_drive_AS.support.QuFluxFit import plot_QbFlux



def CpCav_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,flux_ctrl:dict,ro_span_Hz:int=3e6,flux_span:float=0.3,n_avg:int=10,f_points:int=20,flux_points:int=20,run:bool=True,q:str='q1',c:str="c0",Experi_info:dict={}):

    sche_func = One_tone_sche
        
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    ro_f_center = qubit_info.clock_freqs.readout()
    # avoid frequency conflicts 
    from numpy import NaN
    qubit_info.clock_freqs.readout(NaN)
    qubit_info.measure.pulse_duration(100e-6)
    qubit_info.measure.integration_time(100e-6-1e-6)
    qubit_info.reset.duration(250e-6)
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,f_points)
    flux_samples = linspace(-flux_span,flux_span,flux_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        powerDep=False,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     Flux=['start '+'%E' %flux_samples[0],'end '+'%E' %flux_samples[-1]])
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,flux_ctrl[c]])
        meas_ctrl.setpoints_grid((ro_f_samples,flux_samples))
        
        
        
        rfs_ds = meas_ctrl.run("One-tone-Flux")
        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rfs_ds,qb=q,exp_type='FD',get_data_loc=True)
        analysis_result[q] = ResonatorFluxSpectroscopyAnalysis(tuid=rfs_ds.attrs["tuid"], dataset=rfs_ds).run(sweetspot_index=1)
        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        plot_QbFlux(QD_agent,nc_path,q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        # reset flux bias
        flux_ctrl[c](0.0)
        qubit_info.clock_freqs.readout(ro_f_center)
        
    else:
        n_s = 2
        sweep_para= array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)

        show_args(exp_kwargs, title="One_tone_FluxDep_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result



def update_coupler_bias(QD_agent:QDmanager,cp_elements:dict):
    """
    Update the idle bias in Fluxmanager for couplers.\n
    --------------------------
    ### Args:\n
    cp_elements = {"c0":0.2}
    """
    for cp in cp_elements:
        QD_agent.Fluxmanager.save_idleBias_for(cp, cp_elements[cp])


def fluxCavity_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,specific_cp:str,run:bool=True,flux_span:float=0.4,ro_span_Hz=3e6,zpts=20,fpts=30):
    
    if run:
        print(f"{specific_qubits} are under the measurement ...")
        FD_results = CpCav_spec(QD_agent,meas_ctrl,Fctrl,ro_span_Hz=ro_span_Hz,q=specific_qubits,c=specific_cp,flux_span=flux_span,flux_points=zpts,f_points=fpts)[specific_qubits]
        if FD_results == {}:
            print(f"Flux dependence error qubit: {specific_qubits}")
        
    else:
        FD_results = CpCav_spec(QD_agent,meas_ctrl,Fctrl,ro_span_Hz=ro_span_Hz,q=specific_qubits,c=specific_cp,flux_span=flux_span,run=False,flux_points=zpts,f_points=fpts)

    return FD_results

# accident: q2, q3, q4

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 1
    sweet_spot:bool = 0
    DRandIP = {"dr":"drke","last_ip":"242"}
    ro_elements = {"q1":["c1"]}

    
    """ Optional paras """
    freq_half_window_Hz = 7e6
    flux_half_window_V  = 0.4
    freq_data_points = 80
    flux_data_points = 80
    freq_center_shift = 0e6 # freq axis shift

    
    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    
    """ Running """
    from qblox_drive_AS.Configs.ClusterAddress_rec import get_CouplerController, ip_register
    update = False
    Cctrl = get_CouplerController(cluster, ip=ip_register[DRandIP["dr"]])
    FD_results = {}
    for qubit in ro_elements:
        qu = QD_agent.quantum_device.get_element(qubit)
        qu.clock_freqs.readout(qu.clock_freqs.readout()+freq_center_shift)
        for coupler in ro_elements[qubit]:
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
            if sweet_spot:
                Fctrl[qubit](QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qubit))
            FD_results[qubit] = fluxCavity_executor(QD_agent,meas_ctrl,Cctrl,qubit,coupler,run=execution,flux_span=flux_half_window_V,ro_span_Hz=freq_half_window_Hz, zpts=flux_data_points,fpts=freq_data_points)
            if sweet_spot:
                Fctrl[qubit](0)
            cluster.reset()
        
        """ Storing """
        Fctrl.update(Cctrl)
    """ Close """
    print('Rabi-spliting done!')

    shut_down(cluster,Fctrl)