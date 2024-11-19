"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from qblox_drive_AS.support.UserFriend import *
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import linspace, array, arange, NaN, ndarray
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing, reset_offset
from qblox_drive_AS.support.Pulse_schedule_library import multi_PowerRabi_sche, multi_TimeRabi_sche, set_LO_frequency, pulse_preview
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.analysis.raw_data_demolisher import Rabi_dataReducer
from xarray import Dataset

#? The way to merge two dict a and b : c = {**a,**b}

def round_to_nearest_multiple_of_multipleNS(x, specific_multiple:int=None):
    if specific_multiple is None:
        specific_multiple = 4
    return specific_multiple * round(x / specific_multiple)


def PowerRabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,pi_amp:dict,pi_dura:dict,n_avg:int=300,run:bool=True,OSmode:bool=False):
    
    sche_func= multi_PowerRabi_sche

    for q in pi_amp:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
        pi_amp_idxes = arange(0,pi_amp[q].shape[0])
        
    
    Sweep_para = ManualParameter(name="PowerRabi", unit="Volt", label="Amp")
    Sweep_para.batched = True
    
    
    sched_kwargs = dict(
        pi_amp = pi_amp,
        pi_dura=pi_dura,
        R_amp=compose_para_for_multiplexing(QD_agent,pi_amp,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,pi_amp,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,pi_amp,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,pi_amp,'r2'),
        XY_theta='X_theta',
        OS_or_not=OSmode
        )
    
    
    if run:
        gettable = ScheduleGettable(
        QD_agent.quantum_device,
        schedule_function=sche_func,
        schedule_kwargs=sched_kwargs,
        real_imag=True,
        batched=True,
        num_channels=len(list(pi_amp.keys())),
        )
        
   
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(pi_amp_idxes)
    
       
        rabi_ds = meas_ctrl.run("RabiOsci")

        dict_ = {}
        for idx, q in enumerate(list(pi_amp.keys())):
            I = array(rabi_ds[f'y{2*idx}'])
            Q = array(rabi_ds[f'y{2*idx+1}'])
            dict_[q] = (['mixer','pi_amp'],array([I,Q]))
            dict_[f"{q}_piamp"] = (['mixer','pi_amp'],array(2*list(pi_amp[q])).reshape(2,pi_amp[q].shape[0]))
        ds = Dataset(dict_, coords={"mixer":array(["I","Q"]),"pi_amp":pi_amp_idxes})
        ds.attrs["execution_time"] = Data_manager().get_time_now()
        ds.attrs["rabi_type"] = "PowerRabi"
        ds.attrs["OS_mode"] = OSmode
        for q in pi_dura:
            ds.attrs[f"{q}_pidura"] = pi_dura[q]

     
    else:

        preview_para_pi_amp = {}
        preview_para_pi_dura = {}
        for q in pi_amp:
            preview_para_pi_amp[q] = array([pi_amp[q][0],pi_amp[q][-1]])
            preview_para_pi_dura[q] = pi_dura[q]
        
        sched_kwargs['pi_amp']= preview_para_pi_amp
        sched_kwargs['pi_dura']= preview_para_pi_dura
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        ds = ""

    return ds
    
def TimeRabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,pi_amp:dict,pi_dura:dict,n_avg:int=300,run:bool=True,OSmode:bool=False):
    
    sche_func= multi_TimeRabi_sche

    for q in pi_dura:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
        pi_dura_idxes = arange(0,pi_dura[q].shape[0])
        
    
    Sweep_para = ManualParameter(name="TimeRabi", unit="Sec", label="time")
    Sweep_para.batched = True
    
    
    sched_kwargs = dict(
        pi_amp = pi_amp,
        pi_dura=pi_dura,
        R_amp=compose_para_for_multiplexing(QD_agent,pi_dura,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,pi_dura,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,pi_dura,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,pi_dura,'r2'),
        XY_theta='X_theta',
        OS_or_not=OSmode
        )
    
    
    if run:
        gettable = ScheduleGettable(
        QD_agent.quantum_device,
        schedule_function=sche_func,
        schedule_kwargs=sched_kwargs,
        real_imag=True,
        batched=True,
        num_channels=len(list(pi_dura.keys())),
        )
        
   
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(pi_dura_idxes)
    
       
        rabi_ds = meas_ctrl.run("RabiOsci")

        dict_ = {}
        for idx, q in enumerate(list(pi_dura.keys())):
            I = array(rabi_ds[f'y{2*idx}'])
            Q = array(rabi_ds[f'y{2*idx+1}'])
            dict_[q] = (['mixer','pi_dura'],array([I,Q]))
            dict_[f"{q}_pidura"] = (['mixer','pi_dura'],array(2*list(pi_dura[q])).reshape(2,pi_dura[q].shape[0]))
        ds = Dataset(dict_, coords={"mixer":array(["I","Q"]),"pi_dura":pi_dura_idxes})
        ds.attrs["execution_time"] = Data_manager().get_time_now()
        ds.attrs["rabi_type"] = "TimeRabi"
        ds.attrs["OS_mode"] = OSmode
        for q in pi_amp:
            ds.attrs[f"{q}_piamp"] = pi_amp[q]

     
    else:

        preview_para_pi_amp = {}
        preview_para_pi_dura = {}
        for q in pi_amp:
            preview_para_pi_dura[q] = array([pi_dura[q][0],pi_dura[q][-1]])
            preview_para_pi_amp[q] = pi_amp[q]
        
        sched_kwargs['pi_amp']= preview_para_pi_amp
        sched_kwargs['pi_dura']= preview_para_pi_dura
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        ds = ""

    return ds

def conditional_update_qubitInfo(QD_agent:QDmanager,fit_results:Dataset,target_q:str):
    if fit_results.attrs['pi_2'] >= min(fit_results.coords['samples']) and fit_results.attrs['pi_2'] <= max(fit_results.coords['samples']) :
        qubit = QD_agent.quantum_device.get_element(target_q)
        match str(fit_results.attrs['Rabi_type']).lower():
            case 'powerrabi':
                qubit.rxy.amp180(fit_results.attrs['pi_2'])
                qubit.rxy.duration(fit_results.attrs["fix_variable"])
            case 'timerabi':
                qubit.rxy.amp180(fit_results.attrs["fix_variable"])
                qubit.rxy.duration(fit_results.attrs['pi_2'])
    else:
        warning_print(f"Results for {target_q} didn't satisfy the update condition !")

    
  
def RabiOsci_waiter(QD_agent:QDmanager,ro_elements:dict):
    """
    There must be the keys in the ro_elements[qubit] accordingly named with 'XY_amp' and 'XY_duration'.\n
    For example: ro_elements = {"q0":{"XY_amp":[],"XY_duration":40e-9,"XY_IF":150e6}, ...}\n
    Which name contains an Arraylike value will be the varable to execute. And the other one must be a float number.
    """ 
    # Input paras check
    for qubit in ro_elements:
        if "XY_amp" in list(ro_elements[qubit].keys()) or "xy_amp" in list(ro_elements[qubit].keys()):
            if "XY_duration" in list(ro_elements[qubit].keys()) or "xy_duration" in list(ro_elements[qubit].keys()):
                pass
            else:
                raise KeyError(f"Check your varables in ro_elements for {qubit}, there must be 'XY_duration' ! ")
        else:
            raise KeyError(f"Check your varables in ro_elements for {qubit}, there must be 'XY_amp' ! ")

    new_ro_elements = {}
    for qubit in ro_elements:
        set_LO_mark:bool = False
        xy_freq = QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()
        new_ro_elements[qubit] = {}
        for item in ro_elements[qubit]:
            if type(ro_elements[qubit][item]) == float:
                match str(item).lower():
                    case "xy_duration":
                        new_ro_elements[qubit]["XY_duration"] = ro_elements[qubit][item]
                    case "xy_amp":
                        new_ro_elements[qubit]["XY_amp"] = ro_elements[qubit][item]
                    case "xy_if":
                        set_LO_mark = True
                        LO = xy_freq+ro_elements[qubit][item]
                        set_LO_frequency(QD_agent.quantum_device,q=qubit,module_type='drive',LO_frequency=LO)
                    case _:
                        pass
            elif type(ro_elements[qubit][item]) in [ndarray, list]:
                match str(item).lower():
                    case "xy_duration":
                        time_sample_sec = list(set([round_to_nearest_multiple_of_4ns(x) for x in array(ro_elements[qubit][item])*1e9]))
                        time_sample_sec.sort()
                        new_ro_elements[qubit]["XY_duration"] = array(time_sample_sec)*1e-9
                    case "xy_amp":
                        new_ro_elements[qubit]["XY_amp"] = ro_elements[qubit][item]
                    case _:
                        pass
            else:
                pass
        if not set_LO_mark :
            set_LO_frequency(QD_agent.quantum_device,q=qubit,module_type='drive',LO_frequency=xy_freq+150e6)

    return new_ro_elements


if __name__ == "__main__":
    #? please use `np.linspace()` to sample your varables, for example: {"q0":{"XY_amp":np.linspace(-0.6,0.6,100),...}}
    #! Warning: For all the varables, the sample number must be the same!
    #// Because time rabi resolution is limited as 4ns, it will makes sample number different between different qubits,
    #// please do the same type rabi oscillation for all the qubits till the time resolution limit got released.
    
    """ Fill in """
    execution:bool = True
    chip_info_restore:bool = 0
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {"q0":{"XY_amp":linspace(-0.6,0.6,100),"XY_duration":40e-9,"XY_IF":150e6},
                   "q1":{"XY_amp":linspace(-0.6,0.6,100),"XY_duration":40e-9,"XY_IF":150e6}}
    couplers = []

    
    """ Optional paras """
    avg_n:int = 300
    OneShot:bool = False
    xy_atte:int = 10

    
    """ Preparations """
    start_time = time.time()
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    ro_elements = RabiOsci_waiter(QD_agent,ro_elements)
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    
    for qubit in ro_elements:
        QD_agent.Notewriter.save_DigiAtte_For(xy_atte,qubit,'xy')
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
    
    
    """Running """
    nc_path = rabi_executor(QD_agent,cluster,meas_ctrl,Fctrl,ro_elements,run=execution,OneShot=OneShot,avg_times=avg_n)
    slightly_print(f"Raw data loc:\n{nc_path}")


    """ Analysis """
    dss = Rabi_dataReducer(nc_path)  
    ANA = Multiplex_analyzer("m11")      
    for q in dss:
        ANA._import_data(dss[q],1,QD_agent.refIQ)
        ANA._start_analysis()
        ANA._export_result(Data_manager().get_today_picFolder())
        conditional_update_qubitInfo(QD_agent,ANA.fit_packs,ro_elements,q)


    """ Storing """
    QD_agent.refresh_log("after Rabi")
    QD_agent.QD_keeper()
    if chip_info_restore:
        chip_info.update_RabiOsci(qb=qubit)


    """ Close """
    print('Rabi osci done!')
    shut_down(cluster,Fctrl)
    end_time = time.time()
    slightly_print(f"time cost: {round(end_time-start_time,1)} secs")