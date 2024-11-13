import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from numpy import array, linspace, arange, sqrt
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, reset_offset, coupler_zctrl, compose_para_for_multiplexing
from qblox_drive_AS.analysis.raw_data_demolisher import fluxQub_dataReductor
from qblox_drive_AS.support.Pulse_schedule_library import multi_Z_gate_two_tone_sche, set_LO_frequency, pulse_preview, QS_fit_analysis
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer

z_pulse_amp_OVER_const_z = sqrt(2)/2.5


def Zgate_two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,XYFs:dict,Bias:dict,ref_z:dict,n_avg:int=1000,run:bool=True):
    print("Zgate 2tone start")
    
    sche_func = multi_Z_gate_two_tone_sche

    original_xyfs = {}
    for q in XYFs:
        freq_datapoint_idx = arange(0,XYFs[q].shape[0])
        qubit_info = QD_agent.quantum_device.get_element(q)  
        qubit_info.reset.duration(250e-6) 
        original_xyfs[q] = qubit_info.clock_freqs.f01()
        qubit_info.clock_freqs.f01(NaN)
        spec_pulse_amp = QD_agent.Notewriter.get_2tone_piampFor(q)
        qubit_info.measure.pulse_duration(1.5e-6)
        qubit_info.measure.integration_time(1.5e-6)
    
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    Z_bias = ManualParameter(name="Flux", unit="V", label="Z bias")
    Z_bias.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=XYFs,
        bias_qs=Bias["bias_qs"],
        Z_amp=Z_bias,
        spec_amp=spec_pulse_amp,
        spec_Du=10e-6,
        R_amp=compose_para_for_multiplexing(QD_agent,XYFs,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,XYFs,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,XYFs,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,XYFs,'r2'),
    )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
            num_channels=len(list(XYFs.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,Z_bias])
        meas_ctrl.setpoints_grid((freq_datapoint_idx,Bias["bias_span"]))
        qs_ds = meas_ctrl.run("Zgate-two-tone")
        qs_ds.attrs["RO_qs"] = ""
        qs_ds.attrs["ref_z"] = ""
        qs_ds.attrs["execution_time"] = Data_manager().get_time_now()
        
        for idx, q in enumerate(XYFs):
            attr_0 = qs_ds['x0'].attrs
            qs_ds[f'x{2*idx}'] = array(list(XYFs[q])*array(Bias["bias_span"]).shape[0])
            qs_ds[f'x{2*idx}'].attrs = attr_0
            
            attr_1 = qs_ds['x1'].attrs
            qs_ds[f'x{2*idx+1}'] = array([[i]*XYFs[q].shape[0] for i in array(Bias["bias_span"])/z_pulse_amp_OVER_const_z]).reshape(-1)
            qs_ds[f'x{2*idx+1}'].attrs = attr_1

            qs_ds.attrs["RO_qs"] += f" {q}"
        for idx, q in enumerate(Bias["bias_qs"]):
            qs_ds.attrs["ref_z"] += f"{q}_{round(ref_z[q],4)}_"
        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb="multiQ",exp_type='F2tone',get_data_loc=True)
        


    else:
        n_s = 2
        preview_para = {}
        for q in XYFs:
            preview_para[q] = XYFs[q][:n_s]
        sweep_para2 = array(Bias["bias_span"][:2])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        nc_path = ""


        
    for q in XYFs:
        qubit_info = QD_agent.quantum_device.get_element(q)   
        qubit_info.clock_freqs.f01(original_xyfs[q])

    return nc_path


def update_by_fluxQubit(QD_agent:QDmanager,correct_results:dict,target_q:str):
    """
    correct_results dict in the form: {"xyf":float,"sweet_bias":float}
    """
    qubit = QD_agent.quantum_device.get_element(target_q)
    qubit.clock_freqs.f01(correct_results["xyf"])
    QD_agent.Fluxmanager.check_offset_and_correctFor(target_q=target_q,new_offset=correct_results["sweet_bias"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=target_q,bias=correct_results["sweet_bias"])


def flux_qubitspectro_waiter(QD_agent:QDmanager,ro_elements:dict,bias_elements:dict)->tuple[dict,dict]:
    # frequencies
    freqs = {}
    xyf_pts = ro_elements["xyf_pts"]
    for q in ro_elements:
        if q[0] == "q":
            if len(ro_elements[q]["assigned_xyf_range"]) != 0:
                freqs[q] = linspace(min(ro_elements[q]["assigned_xyf_range"]),max(ro_elements[q]["assigned_xyf_range"]),xyf_pts)
                set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=max(ro_elements[q]["assigned_xyf_range"]))
            else:
                qubit_info = QD_agent.quantum_device.get_element(q)   
                original_f01 = qubit_info.clock_freqs.f01()
                set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=original_f01+ro_elements[q]["xy_IF"])
                freqs[q] = linspace(original_f01+ro_elements[q]["xy_IF"]-ro_elements[q]["freq_span"],original_f01+(ro_elements[q]["xy_IF"]),xyf_pts)
    # flux
    bias_elements["bias_qs"] = bias_elements["bias_qs"] if len(bias_elements["bias_qs"]) != 0 else list(freqs.keys())
    bias_elements["bias_span"] = linspace(-bias_elements["bias_span"]*z_pulse_amp_OVER_const_z,bias_elements["bias_span"]*z_pulse_amp_OVER_const_z,bias_elements["bias_pts"])
    
    return freqs, bias_elements



def fluxQubit_executor(QD_agent:QDmanager,Fctrl:dict,meas_ctrl:MeasurementControl,XYFs:dict,bias:dict,run:bool=True,avg_times:int=1000):
    nc_path = ""
    if run:
        ref_z = {}
        for q in bias["bias_qs"]:
            center = QD_agent.Fluxmanager.get_proper_zbiasFor(target_q=q)
            Fctrl[q](center)
            ref_z[q] = center
        
        nc_path = Zgate_two_tone_spec(QD_agent,meas_ctrl,XYFs,bias,ref_z,run=True,n_avg=avg_times)
        reset_offset(Fctrl)
    else:
        _ = Zgate_two_tone_spec(QD_agent,meas_ctrl,XYFs,bias,ref_z,run=False,n_avg=avg_times)
    return nc_path

if __name__ == "__main__":
    #?  NOTE: Readout the qubit in the ro_elements at the same time, bias the qubit in the bias_element at the same time 
    #// NOTE: If the "bias_qs" in the bias_element is empty, bias the qubit in the ro_elements at the same time
    
    """ Fill in """
    execution:bool = True
    chip_info_restore:bool = 0
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {'q0':{"freq_span":500e6,"xy_IF":200e6,"assigned_xyf_range":[]},'q1':{"freq_span":500e6,"xy_IF":200e6,"assigned_xyf_range":[]},"xyf_pts":50}
    bias_element = {"bias_qs":['q0',"q1"],"bias_span":0.05,"bias_pts":20}
    couplers = []
    

    """ Optional paras """
    avg_n:int = 500


    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    if ro_elements == 'all':
        ro_elements = list(Fctrl.keys())
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    freq_elements, bias_elements = flux_qubitspectro_waiter(QD_agent,ro_elements,bias_element)


    """ Running """
    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        if qubit[0] == "q":
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
    
    nc_path = fluxQubit_executor(QD_agent,Fctrl,meas_ctrl,freq_elements,bias_elements,run=execution,avg_times=avg_n)
    slightly_print(f"Raw data loc:\n{nc_path}")
    cluster.reset()

    
    """ Analysis """
    ana_results = {}
    dss = fluxQub_dataReductor(nc_path)
    ANA = Multiplex_analyzer("m9") 
    for q in dss:
        ANA.target_q = q
        ANA._import_data(dss[q],2,QD_agent.refIQ[q],QS_fit_analysis)
        ANA._start_analysis()
        ANA._export_result(Data_manager().get_today_picFolder())
        ana_results[q] = ANA.fit_packs
        if len(list(ana_results[q].keys())) != 0:
            update_by_fluxQubit(QD_agent,ana_results[q],q)
    
    
    """ Storing """
    trustable = mark_input("Keep these results ? [y/n]")
    if  trustable.lower() in ['y', 'yes']:
        QD_agent.QD_keeper()
        if chip_info_restore:
            chip_info.update_FluxQubit(qb=qubit, result=ana_results[q])
    else:
        warning_print(f"The results were deleted !")


    """ Close """
    print('Flux qubit done!')
    shut_down(cluster,Fctrl)

  
    


