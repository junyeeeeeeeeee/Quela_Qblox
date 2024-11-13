import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from numpy import NaN
from numpy import array, linspace
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, reset_offset, coupler_zctrl
from qblox_drive_AS.support.QuFluxFit import plot_QbFlux
from qblox_drive_AS.support.Pulse_schedule_library import Z_gate_two_tone_sche, set_LO_frequency, pulse_preview




def Zgate_two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Z_amp_start:float,Z_amp_end:float,IF:int=200e6,xyf:float=0e9,xyf_span_Hz:float=400e6,n_avg:int=1000,RO_z_amp:float=0,Z_points:int=40,f_points:int=60,run:bool=True,q:str='q1',Experi_info={},get_data_path:bool=False,analysis:bool=True):
    print("Zgate 2tone start")
    trustable = True
    sche_func = Z_gate_two_tone_sche

    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    original_f01 = qubit_info.clock_freqs.f01()
    print(original_f01)

    if xyf == 0:
        xyf_highest = original_f01+IF
    else:
        xyf_highest = xyf + IF
    qubit_info.clock_freqs.f01(NaN)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf_highest)
    f01_samples = linspace(xyf_highest-xyf_span_Hz,xyf_highest,f_points)
    
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    
    # temperature quard
    if Z_amp_end > 0.4:
        Z_amp_end = 0.4
    elif Z_amp_end < -0.4:
        Z_amp_end = -0.4
    else:
        pass

    if Z_amp_start > 0.4:
        Z_amp_start = 0.4
    elif Z_amp_start < -0.4:
        Z_amp_start = -0.4
    else:
        pass 

    Z_samples = linspace(Z_amp_start,Z_amp_end,Z_points)
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        Z_amp=Z_bias,
        spec_amp=QD_agent.Notewriter.get_2tone_piampFor(q),
        spec_Du=50*1e-6,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        Z_ro_amp=RO_z_amp,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
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
        meas_ctrl.settables([freq,Z_bias])
        meas_ctrl.setpoints_grid((f01_samples,Z_samples))
        qs_ds = meas_ctrl.run("Zgate-two-tone")
        
        # Save the raw data into netCDF
        if get_data_path:
            path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='F2tone',get_data_loc=get_data_path)
        else:
            path = ''
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='F2tone',get_data_loc=get_data_path)

        if analysis:
            try:
                analysis_result[q] = QubitFluxSpectroscopyAnalysis(tuid=qs_ds.attrs["tuid"], dataset=qs_ds).run()
            except:
                analysis_result[q] = {}
                print("Qb vs Flux fitting failed! Raw data had been saved.")
                trustable = False
        
        
        show_args(exp_kwargs, title="Zgate_two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s = 2
        sweep_para1= array(f01_samples[:n_s])
        sweep_para2= array(Z_samples[:2])
        spec_sched_kwargs['frequencies']= sweep_para1.reshape(sweep_para1.shape or (1,))
        spec_sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        
        
        show_args(exp_kwargs, title="Zgate_two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        path = ''

    qubit_info.clock_freqs.f01(original_f01)

    return analysis_result, path, trustable


def update_by_fluxQubit(QD_agent:QDmanager,correct_results:dict,target_q:str):
    """
    correct_results dict in the form: {"xyf":float,"sweet_bias":float}
    """
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.clock_freqs.f01(correct_results["xyf"])
    QD_agent.Fluxmanager.check_offset_and_correctFor(target_q=target_q,new_offset=correct_results["sweet_bias"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=target_q,bias=correct_results["sweet_bias"])



def fluxQubit_executor(QD_agent:QDmanager,Fctrl:dict,meas_ctrl:MeasurementControl,specific_qubits:str,run:bool=True,z_shifter:float=0,zpts:int=5,fpts:int=40,span_priod_factor:int=12,f_sapn_Hz=400e6,avg_times:int=1000,xy_IF:float=200e6):
    center = QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=specific_qubits)
    partial_period = QD_agent.Fluxmanager.get_PeriodFor(target_q=specific_qubits)/span_priod_factor

    if run:
        
        Fctrl[specific_qubits](center)
        results, nc_path, trustable= Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=-partial_period+z_shifter,Z_points=zpts,f_points=fpts,Z_amp_end=partial_period+z_shifter,q=specific_qubits,run=True,get_data_path=True,xyf_span_Hz=f_sapn_Hz,IF=xy_IF,n_avg=avg_times)
        reset_offset(Fctrl)
        if trustable:
            plot_QbFlux(QD_agent,nc_path,specific_qubits)
            permission = mark_input("Update the QD with this result ? [y/n]") 
            if permission.lower() in ['y','yes']:
                return trustable, {"xyf":results[specific_qubits].quantities_of_interest["freq_0"].nominal_value,"sweet_bias":results[specific_qubits].quantities_of_interest["offset_0"].nominal_value+center}
            else:
                return False, {}
        else:
            plot_QbFlux(QD_agent,nc_path,specific_qubits)
            trustable = False
            return False, {}

    else:
        results, _, trustable= Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=center-partial_period+z_shifter,Z_points=10,Z_amp_end=center+partial_period+z_shifter,q=specific_qubits,run=False)
        return False, {}


if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = False
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = ['q1']
    couplers = []
    z_shifter = 0.0 # V

    
    """ Optional paras """
    span_period_factor:float = 12 # range in [sweet - period/span_period_factor, sweet + period/span_period_factor]
    flux_pts:int = 40
    freq_pts:int = 40
    freq_span_Hz:float = 500e6
    sweet_flux_shifter:float = 0
    xy_IF = 100e6
    avg_n:int = 500



    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    chip_info = cds.Chip_file(QD_agent=QD_agent)


    """ Running """
    FQ_results = {}
    check_again =[]
    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        if not QD_agent.Fluxmanager.get_offsweetspot_button(qubit):
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            trustable, new_ans = fluxQubit_executor(QD_agent,Fctrl,meas_ctrl,qubit,run=execution,z_shifter=z_shifter,zpts=flux_pts,fpts=freq_pts,span_priod_factor=span_period_factor,f_sapn_Hz=freq_span_Hz,avg_times=avg_n,xy_IF=xy_IF)
            cluster.reset()

            """ Storing """
            if  trustable:
                update_by_fluxQubit(QD_agent,new_ans,qubit)
                QD_agent.QD_keeper()
                if chip_info_restore:
                    chip_info.update_FluxQubit(qb=qubit, result=new_ans)
            else:
                check_again.append(qubit)    

    """ Close """
    print('Flux qubit done!')
    if len(check_again) != 0:
        warning_print(f"qubits to check again: {check_again}")
    shut_down(cluster,Fctrl)

  
    


