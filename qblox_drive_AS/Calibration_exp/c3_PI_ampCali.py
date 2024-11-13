"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from qblox_drive_AS.support.UserFriend import *
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from xarray import Dataset
from numpy import linspace, array, arange, NaN, ndarray, sin, mean, cos, moveaxis
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing, reset_offset
from qblox_drive_AS.support.Pulse_schedule_library import multi_PI_amp_cali_sche, set_LO_frequency, pulse_preview
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.analysis.raw_data_demolisher import piampcali_dataReducer


def pi_amp_cali(QD_agent:QDmanager,meas_ctrl:MeasurementControl, ro_elements:dict, pi_pair_num:list=[3,5],n_avg:int=300,run:bool=True,specific_data_folder:str=''):
    results = {}
    sche_func= multi_PI_amp_cali_sche
    nc_path = ""
    for q in ro_elements["amp_samples"]:
        results[q], results[f"{q}_PIcoef"] = [], []
        data_sample_idx = arange(ro_elements["amp_samples"][q].shape[0])
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")

    
    Sweep_para = ManualParameter(name="XY_Amp_coef")
    Sweep_para.batched = True
    
    def pi_pair_dep_exe(pi_pair_num:int)->dict:
        dict_ = {}
        sched_kwargs = dict(
            pi_amp_coefs=ro_elements["amp_samples"],
            pi_pair_num=pi_pair_num,
            XY_amp=compose_para_for_multiplexing(QD_agent,ro_elements["amp_samples"],'d1'),
            XY_duration=compose_para_for_multiplexing(QD_agent,ro_elements["amp_samples"],'d3'),
            R_amp=compose_para_for_multiplexing(QD_agent,ro_elements["amp_samples"],'r1'),
            R_duration=compose_para_for_multiplexing(QD_agent,ro_elements["amp_samples"],'r3'),
            R_integration=compose_para_for_multiplexing(QD_agent,ro_elements["amp_samples"],'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements["amp_samples"],'r2'),
            )
        
        
        if run:
            gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements["amp_samples"].keys())),
            )
            
    
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            meas_ctrl.gettables(gettable)
            meas_ctrl.settables(Sweep_para)
            meas_ctrl.setpoints(data_sample_idx)
        
        
            ds = meas_ctrl.run("Pi amp calibration")
            
            for q_idx, q in enumerate(ro_elements["amp_samples"]):
                I_data, Q_data = array(ds[f"y{2*q_idx}"]).tolist(), array(ds[f"y{2*q_idx+1}"]).tolist()
                dict_[q] = [I_data,Q_data] # shape in (mixer, pi_amp)
                dict_[f"{q}_PIcoef"] = [list(ro_elements["amp_samples"][q])]*2
        else:
            preview_para = {}
            for q in ro_elements["amp_samples"]:
                preview_para[q] = array([ro_elements["amp_samples"][q][0],ro_elements["amp_samples"][q][-1]])
            sched_kwargs['pi_amp_coefs']= preview_para
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

        return dict_
    
    for idx, pi_pair_number in enumerate(pi_pair_num):
        slightly_print(f"{pi_pair_number}*2 pi-pulses ...")
        dataDict = pi_pair_dep_exe(pi_pair_number) # {"q0":[],"q0_PIcoef":[], ...}
        for var in dataDict:
            results[var].append(dataDict[var])
            if idx == len(pi_pair_num) - 1: results[var] = (["mixer", "PiPairNum", "PiAmpCoef"],moveaxis(array(results[var]),0,1)) # shape (pi-pair_num, mixer, rof) -> (mixer, pi-pair_num, rof)

    
    if run:
        dataset_2_nc = Dataset(results,coords={"mixer":array(["I","Q"]),"PiPairNum":array(pi_pair_num),"PiAmpCoef":data_sample_idx})
        dataset_2_nc.attrs["execution_time"] = Data_manager().get_time_now()
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=dataset_2_nc,qb="multiQ",exp_type="xylcali",specific_dataFolder=specific_data_folder,get_data_loc=True)
   
    return nc_path


def pi_amp_calibrator(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,ro_elements:dict,pi_pair_num:list=[3,5],run:bool=True,avg_times:int=500,data_folder:str=''):
    if run:
        for q in ro_elements["amp_samples"]:
            Fctrl[q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(q)))
        nc_path = pi_amp_cali(QD_agent,meas_ctrl,ro_elements,run=True,pi_pair_num=pi_pair_num,n_avg=avg_times,specific_data_folder=data_folder)
        reset_offset(Fctrl)
        cluster.reset()
        
        return nc_path
    else:
        _ = pi_amp_cali(QD_agent,meas_ctrl,ro_elements,run=False,pi_pair_num=pi_pair_num,n_avg=avg_times,specific_data_folder=data_folder)
        return 0
  
def pi_cali_waiter(QD_agent:QDmanager, ro_element:dict, amp_pts:int=100)->dict:
    """
    ro_elements = {"q0":{"pi_amp_coef_span":0.15,"xy_IF":250e6}} for example.
    """
    ro_elements = {"amp_samples":{}}
    for q in ro_element :
        # LO setting
        IF_key = [name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'][0] if len([name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'])==1 else ""      
        IF = ro_element[q][IF_key] if IF_key != "" else 150e6
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=QD_agent.quantum_device.get_element(q).clock_freqs.f01()+IF)
        # pi-amp samples settings 
        ro_elements["amp_samples"][q] = linspace(1-ro_element[q]["pi_amp_coef_span"],1+ro_element[q]["pi_amp_coef_span"],amp_pts) 

    return ro_elements

if __name__ == "__main__":
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {"q0":{"pi_amp_coef_span":0.2,"xy_IF":250e6},"q1":{"pi_amp_coef_span":0.2,"xy_IF":250e6}}
    couplers = []


    """ Optional paras """
    pi_pair_num:list = [2,3]
    data_pts:int = 80
    avg_n:int = 500
    
    

    """ Preparations """
    start_time = time.time()
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    ro_elements = pi_cali_waiter(QD_agent,ro_elements,data_pts)
    for q in ro_elements["amp_samples"]:
        init_system_atte(QD_agent.quantum_device,list([q]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
    

    """Running """
    nc_path = pi_amp_calibrator(QD_agent,cluster,meas_ctrl,Fctrl,ro_elements,run=execution,avg_times=avg_n,pi_pair_num=pi_pair_num)
    slightly_print(f"File sved located:\n{nc_path}")

    """ Storing """
    ds = piampcali_dataReducer(nc_path)
    qubit = [x for x in list(ds.data_vars) if x.split("_")[-1] != "PIcoef"]
    for var in qubit:
        ANA = Multiplex_analyzer("c3")
        ANA._import_data(ds,var_dimension=1,refIQ=QD_agent.refIQ)
        ANA._start_analysis(var_name = var)
        fit_pic_folder = Data_manager().get_today_picFolder()
        ANA._export_result(fit_pic_folder)
        fit_packs = ANA.fit_packs

    """ Close """
    shut_down(cluster,Fctrl)
    end_time = time.time()
    slightly_print(f"time cost: {round(end_time-start_time,1)} secs")

    

