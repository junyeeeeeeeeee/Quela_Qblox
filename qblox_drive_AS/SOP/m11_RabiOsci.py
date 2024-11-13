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
from qblox_drive_AS.support.Pulse_schedule_library import multi_Rabi_sche, set_LO_frequency, pulse_preview
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.analysis.raw_data_demolisher import Rabi_dataReducer
from xarray import Dataset

#? The way to merge two dict a and b : c = {**a,**b}

def round_to_nearest_multiple_of_4ns(x):
    return 4 * round(x / 4)


def Rabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Paras:dict,n_avg:int=300,run:bool=True,OSmode:bool=False,specific_data_folder:str=''):
    nc_path = ''
    sche_func= multi_Rabi_sche
    pulse_review_paras = {}
    ds_restore = {}

    for q in Paras:
        pulse_review_paras[q],ds_restore[q] = {}, {}
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
        if isinstance(Paras[q]["XY_amp"],ndarray):
            dummy_sweep_varabel = arange(0,Paras[q]["XY_amp"].shape[0])
            pulse_review_paras[q]["XY_amp"] = Paras[q]["XY_amp"][-2:]
            pulse_review_paras[q]["XY_duration"] = Paras[q]["XY_duration"]
            ds_restore[q]["XY_amp"] = Paras[q]["XY_amp"]
        elif isinstance(Paras[q]["XY_duration"],ndarray):
            dummy_sweep_varabel = arange(0,Paras[q]["XY_duration"].shape[0])
            pulse_review_paras[q]["XY_amp"] = Paras[q]["XY_amp"]
            pulse_review_paras[q]["XY_duration"] = Paras[q]["XY_duration"][-2:]
            ds_restore[q]["XY_duration"] = Paras[q]["XY_duration"]
        else:
            pass
    
    Sweep_para = ManualParameter(name="RabiVarable", unit="SecOrVolt", label="TimeOrAmp")
    Sweep_para.batched = True
    
    
    sched_kwargs = dict(
        Paras = Paras if run else pulse_review_paras,
        R_amp=compose_para_for_multiplexing(QD_agent,Paras,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,Paras,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,Paras,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,Paras,'r2'),
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
        num_channels=len(list(Paras.keys())),
        )
        
   
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(dummy_sweep_varabel)
    
       
        rabi_ds = meas_ctrl.run("RabiOsci")
        rabi_ds.attrs["RO_qs"] = ""
        rabi_ds.attrs["OS_mode"] = OSmode
        rabi_ds.attrs["execution_time"] = Data_manager().get_time_now()
        name_dict = {"XY_duration":{"name":"Time","long_name":"TimeRabi","unit":"s"},"XY_amp":{"name":"Amplitude","long_name":"PowerRabi","unit":"V"}}

        for idx, q in enumerate(list(ro_elements.keys())):
            rabi_ds.attrs["RO_qs"] += f" {q}"
            attr = rabi_ds['x0'].attrs
            rabi_ds[f'x{idx}'] = array(list(ds_restore[q].values())[0])
            rabi_ds[f'x{idx}'].attrs = attr
            rabi_ds[f'x{idx}'].attrs['name'] = name_dict[list(ds_restore[q].keys())[0]]["name"]
            rabi_ds[f'x{idx}'].attrs['long_name'] = name_dict[list(ds_restore[q].keys())[0]]["long_name"]
            rabi_ds[f'x{idx}'].attrs['unit'] = name_dict[list(ds_restore[q].keys())[0]]["unit"]
        
                
        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rabi_ds,qb="multiQ",exp_type="Rabi",specific_dataFolder=specific_data_folder,get_data_loc=True)
        
        
    else:
        sched_kwargs = sched_kwargs
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)


    return nc_path
    

def rabi_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,ro_elements:dict,run:bool=True,OneShot:bool=False,avg_times:int=500,data_folder:str=''):
    
    if run:
        for q in ro_elements:
            Fctrl[q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(q)))
        
        nc_path = Rabi(QD_agent,meas_ctrl,ro_elements,run=True,n_avg=avg_times,specific_data_folder=data_folder,OSmode=OneShot)
        reset_offset(Fctrl)
        cluster.reset()
        
    else:
        nc_path = Rabi(QD_agent,meas_ctrl,ro_elements,run=False,n_avg=avg_times,specific_data_folder=data_folder,OSmode=OneShot)
    
    return nc_path


def conditional_update_qubitInfo(QD_agent:QDmanager,fit_results:Dataset,ro_elements:dict,target_q:str):
    if fit_results.attrs['pi_2'] >= min(fit_results.coords['samples']) and fit_results.attrs['pi_2'] <= max(fit_results.coords['samples']) :
        qubit = QD_agent.quantum_device.get_element(target_q)
        match str(fit_results.attrs['Rabi_type']).lower():
            case 'powerrabi':
                qubit.rxy.amp180(fit_results.attrs['pi_2'])
                qubit.rxy.duration(ro_elements[target_q]["XY_duration"])
            case 'timerabi':
                qubit.rxy.amp180(ro_elements[target_q]["XY_amp"])
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