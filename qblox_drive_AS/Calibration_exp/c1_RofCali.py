import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from numpy import linspace, array, where, max, ndarray, sqrt, arctan2,NaN, arange, moveaxis
from qblox_drive_AS.support import QDmanager, Data_manager, init_meas, shut_down, init_system_atte,coupler_zctrl, reset_offset, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_ROF_Cali_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.analysis.raw_data_demolisher import rofcali_dataReducer

def rofCali(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,n_avg:int=500,run:bool=True,data_folder:str=''):
    sche_func= multi_ROF_Cali_sche
    nc_path = ""
    ro_f_origin = {}
    for q in ro_elements["rof_samples"]:
        qubit_info = QD_agent.quantum_device.get_element(q)
        ro_f_origin[q] = qubit_info.clock_freqs.readout()
        qubit_info.clock_freqs.readout(NaN)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        rof_data_idx = arange(ro_elements["rof_samples"][q].shape[0])
    
    
    option_rof = ManualParameter(name="freq", unit="Hz", label="Frequency")
    option_rof.batched = True

    def state_dep_sched(ini_state:str):
        sched_kwargs = dict(
            ro_freq=ro_elements["rof_samples"],
            ini_state=ini_state,
            pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements["rof_samples"],'d1'),
            pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements["rof_samples"],'d3'),
            R_amp=compose_para_for_multiplexing(QD_agent,ro_elements["rof_samples"],'r1'),
            R_duration=compose_para_for_multiplexing(QD_agent,ro_elements["rof_samples"],'r3'),
            R_integration=compose_para_for_multiplexing(QD_agent,ro_elements["rof_samples"],'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements["rof_samples"],'r2'),
            )
        
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func,
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(ro_elements["rof_samples"].keys())),
            )
            
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            meas_ctrl.gettables(gettable)
            meas_ctrl.settables(option_rof)
            meas_ctrl.setpoints(rof_data_idx)
            
            ds = meas_ctrl.run("Rof-Calibrate")
            dict_ = {}
            for q_idx, q in enumerate(ro_elements["rof_samples"]):
                i_data = array(ds[f'y{2*q_idx}'])
                q_data = array(ds[f'y{2*q_idx+1}'])
                dict_[q] = array([i_data.tolist(),q_data.tolist()]) # shape (mixer, rof_idx)
                dict_[f"{q}_rof"] = array([array(ro_elements["rof_samples"][q])]*2) # shape (mixer, rof_idx)
            
            return dict_
    
        else:
            preview_para = {}
            for q in ro_elements:
                preview_para[q] = ro_elements[q][:2]
            sched_kwargs['ro_freq']= preview_para
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

            return {}
    
    slightly_print("Running |0>")
    data_dict_g = state_dep_sched('g')
    slightly_print("Running |1>")
    data_dict_e = state_dep_sched('e')

    if run:
        dict_ = {}
        for var in data_dict_g :
            if var.split("_")[-1] != "rof":
                state_data = []
                rof_data = [array(data_dict_g[f"{var}_rof"]).tolist()]*2
                for item in [array(data_dict_g[var]).tolist(), array(data_dict_e[var]).tolist()]: # shape (state, mixer, rof)
                    state_data.append(item)
                dict_[var] = (["mixer","state","rof"],moveaxis(array(state_data),0,1)) 
                dict_[f"{var}_rof"] = (["mixer","state","rof"],moveaxis(array(rof_data),0,1))
        
        rofcali_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"state":["g","e"],"rof":rof_data_idx})
        rofcali_ds.attrs["execution_time"] = Data_manager().get_time_now()
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rofcali_ds,qb="multiQ",exp_type='RofCali',specific_dataFolder=data_folder,get_data_loc=True)

        for q in ro_f_origin:
            qubit_info = QD_agent.quantum_device.get_element(q)
            qubit_info.clock_freqs.readout(ro_f_origin[q])
    
    return nc_path


def rofCali_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,ro_elements:dict,execution:bool=True,n_avg:int=300):
    if execution:
        for q in ro_elements["rof_samples"]:
            Fctrl[q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(q)))
        optimal_rof = rofCali(QD_agent,meas_ctrl,ro_elements,run=execution,n_avg=n_avg)
        reset_offset(Fctrl)
        cluster.reset()
    else:
        optimal_rof = rofCali(QD_agent,meas_ctrl,ro_elements,run=execution,n_avg=n_avg)

    return optimal_rof


def rofCali_waiter(QD_agent:QDmanager,ro_element:dict,fpts:int=100)->dict:
    ro_elements = {"rof_samples":{}}
    for q in ro_element:
        # rof samples settings
        ro_f_origin= QD_agent.quantum_device.get_element(q).clock_freqs.readout()
        ro_elements["rof_samples"][q] = linspace(ro_f_origin-ro_element[q]["span_Hz"],ro_f_origin+ro_element[q]["span_Hz"],fpts)
        # driving LO settings
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=QD_agent.quantum_device.get_element(q).clock_freqs.f01()+250e6)

    return ro_elements





if __name__ == '__main__':

    """ Fill in """
    execute:bool = True
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {'q0':{"span_Hz":8e6}}
    couplers = []

    """ Optional paras """
    freq_pts:int = 100
    avg_n:int = 300


    """ Preparation """
    keep = False
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    ro_elements = rofCali_waiter(QD_agent,ro_elements,freq_pts)
    for qubit in ro_elements["rof_samples"]:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
    


    """ Running """
    nc_path = rofCali_executor(QD_agent,cluster,meas_ctrl,Fctrl,ro_elements,execution=execute,n_avg=avg_n)
    slightly_print(f"data saved located:\n{nc_path}")
    if execute:
        ds = rofcali_dataReducer(nc_path)
        qubit = [x for x in list(ds.data_vars) if x.split("_")[-1] != "rof"]
        optimal_rof = {}
        for var in qubit:
            ANA = Multiplex_analyzer("c1")
            ANA._import_data(ds,var_dimension=1)
            ANA._start_analysis(var_name = var)
            fit_pic_folder = Data_manager().get_today_picFolder()
            ANA._export_result(fit_pic_folder)
            optimal_rof[var] = ANA.fit_packs[var]["optimal_rof"]
        
        permission = mark_input(f"Update the optimal ROF for {qubit}?[y/n] or a qubit name to update... ")
        match permission.lower()[0]:
            case 'y':
                for qubit in ro_elements["rof_samples"]:
                    QD_agent.quantum_device.get_element(qubit).clock_freqs.readout(optimal_rof[qubit])
                keep = True
            case 'q':
                QD_agent.quantum_device.get_element(permission.lower()).clock_freqs.readout(optimal_rof[permission.lower()])
                keep = True
            case _:
                warning_print("QD updates got denied !")


    """ Storing """ 
    if execute:
        if keep:
            QD_agent.QD_keeper()


    """ Close """    
    shut_down(cluster,Fctrl)
    