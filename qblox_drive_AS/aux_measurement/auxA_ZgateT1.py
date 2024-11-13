import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from numpy import mean, array, arange, std, linspace, ndarray
from utils.tutorial_utils import show_args
from qblox_drive_AS.support.UserFriend import *
from xarray import Dataset, open_dataset
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing, reset_offset
from qblox_drive_AS.support.Pulse_schedule_library import multi_Zgate_T1_sche, set_LO_frequency, pulse_preview
from qblox_drive_AS.SOP.m9_FluxQubit import z_pulse_amp_OVER_const_z
from qblox_drive_AS.analysis.raw_data_demolisher import ZgateT1_dataReducer
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer

def Zgate_T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,bias_elements:dict,n_avg:int=500,run:bool=True,exp_idx:int=0,no_pi_pulse:bool=False,data_folder:str=''):
    sche_func= multi_Zgate_T1_sche
    origin_pi_amp = {}
    for q in ro_elements["time_samples"]:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        origin_pi_amp[q] = qubit_info.rxy.amp180()
        if no_pi_pulse:
            qubit_info.rxy.amp180(0)
        time_data_idx = arange(ro_elements["time_samples"][q].shape[0])
    

    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples:ndarray = bias_elements["flux_samples"]
    
    sched_kwargs = dict(
        freeduration=ro_elements["time_samples"],
        Z_amp=Z_bias,
        pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements["time_samples"],'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements["time_samples"],'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements["time_samples"],'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements["time_samples"],'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements["time_samples"],'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements["time_samples"],'r2'),
        )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements["time_samples"].keys())),
        )
        
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([Para_free_Du,Z_bias])
        meas_ctrl.setpoints_grid((time_data_idx,Z_samples))


        ds = meas_ctrl.run('Zgate_T1')
        dict_ = {}
        ref_bias = ""
        for q_idx, q in enumerate(ro_elements["time_samples"]):   
            i_data = array(ds[f'y{2*q_idx}']).reshape(Z_samples.shape[0],ro_elements["time_samples"][q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(Z_samples.shape[0],ro_elements["time_samples"][q].shape[0])
            dict_[q] = (["mixer","z_voltage","idx"],array([i_data,q_data]))

            time_values = list(ro_elements["time_samples"][q])*2*Z_samples.shape[0]
            dict_[f"{q}_time"] = (["mixer","z_voltage","time"],array(time_values).reshape(2,Z_samples.shape[0],array(ro_elements["time_samples"][q]).shape[0]))
            ref_bias += f"_{q}_{round(QD_agent.Fluxmanager.get_proper_zbiasFor(q),3)}"
        
        zT1_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"z_voltage":Z_samples/z_pulse_amp_OVER_const_z,"time":arange(time_data_idx.shape[0])})
        zT1_ds.attrs["execution_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) # '2024-10-25 13:02:39'
        zT1_ds.attrs["ref_bias"] = ref_bias
        zT1_ds.attrs["prepare_excited"] = not no_pi_pulse

        Data_manager().save_raw_data(QD_agent=QD_agent,ds=zT1_ds,label=exp_idx,qb="multiQ",exp_type='zT1',specific_dataFolder=data_folder)

    else:
        preview_para = {}
        for q in ro_elements["time_samples"]:
            preview_para[q] = ro_elements["time_samples"][q][:2]
        sweep_para2 = array([bias_elements["flux_samples"][0],bias_elements["flux_samples"][-1]])
        sched_kwargs['freeduration']= preview_para
        sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

    for q in origin_pi_amp:
        QD_agent.quantum_device.get_element(q).rxy.amp180(origin_pi_amp[q])
    

def zgateT1_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,ro_elements:dict,bias_elements:dict,run:bool=True,specific_folder:str='',ith:int=0,pi_pulse:bool=True,n_avg:int=300):
    flux_guard_Volt = 0.4
    for q in ro_elements["time_samples"]:
        offset_flux = float(QD_agent.Fluxmanager.get_proper_zbiasFor(q))
        if abs(offset_flux+max(bias_elements["flux_samples"]))>flux_guard_Volt or abs(offset_flux+min(bias_elements["flux_samples"]))>flux_guard_Volt :
            raise ValueError(f"Z span from {max(abs(offset_flux+bias_elements['flux_samples']))} is not an appropriate voltage !")
            
    if run:
        for q in ro_elements["time_samples"]:
            offset_flux = float(QD_agent.Fluxmanager.get_proper_zbiasFor(q))
            Fctrl[q](offset_flux)

        every_start = time.time()
        slightly_print(f"The {ith}-th T1:")
        Zgate_T1(QD_agent,meas_ctrl,ro_elements,bias_elements, run=True,exp_idx=ith,data_folder=specific_folder,no_pi_pulse = not pi_pulse,n_avg=n_avg)
        reset_offset(Fctrl)
        cluster.reset()
        
        every_end = time.time()
        slightly_print(f"time cost: {round(every_end-every_start,1)} secs")
    else:
        Zgate_T1(QD_agent,meas_ctrl,ro_elements,bias_elements, run=False,exp_idx=ith,data_folder=specific_folder,no_pi_pulse = not pi_pulse,n_avg=n_avg)


def zgateT1_waiter(QD_agent:QDmanager,ro_element:dict,flux_span_factor:float,flux_pts:int=30, freetime_pts:int=50):
    """ 
    Set the flux sweep rasnge by arg `flux_span_factor`, that sweeps in [sweet_bias-period/(2*flux_span_factor), sweet_bias+period/(2*flux_span_factor)] 
    """
    ro_elements, bias_elements = {}, {}
    # z pulse amplitude settings
    flux_period:ndarray = array([QD_agent.Fluxmanager.get_PeriodFor(target_q=q) for q in ro_element])
    bias_elements["flux_samples"] = linspace(-max(flux_period)*z_pulse_amp_OVER_const_z/flux_span_factor,max(flux_period)*z_pulse_amp_OVER_const_z/flux_span_factor,flux_pts)
    # free evolution time settings
    ro_elements["time_samples"] = {}
    for q in ro_element:
        gap = (ro_element[q]["evoT"]*1e9 // freetime_pts) + ((ro_element[q]["evoT"]*1e9 // freetime_pts)%4)
        ro_elements["time_samples"][q] = arange(0,ro_element[q]["evoT"],gap*1e-9)
    # LO settings 
    for q in ro_element:
        IF_key = [name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'][0] if len([name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'])==1 else ""      
        IF = ro_element[q][IF_key] if IF_key != "" else 150e6
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=QD_agent.quantum_device.get_element(q).clock_freqs.f01()+IF)


    return ro_elements, bias_elements

    
   

if __name__ == "__main__":
    # (flux, evo, avg): (50, 50, 100) cost 120 secs, (100, 50, 100) cost 240 secs
    #                   (50, 50, 300) cost 280 secs

    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 0
    prepare_excited:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        "q0":{"evoT":100e-6,"xy_IF":250e6},
        "q1":{"evoT":60e-6,"xy_IF":250e6},
    }
    couplers = []


    """ Optional paras """
    repeat:int = 2
    flux_span_period_factor = 10 # flux axis = [sweet-period/(2*flux_span_period_factor), sweet+period/(2*flux_span_period_factor)]
    flux_data_points = 10
    evotime_data_points = 50
    avg_times:int = 500
    plot_time_vs_flux:bool = False # True for time_trace. Otherwise,  AVG on same bias
    

    """ Preparations """
    specific_folder = Data_manager().creat_datafolder_today(f"ZgateT1_MultiQ_{Data_manager().get_time_now()}")
    for ith_histo in range(repeat):
        
        
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
        chip_info = cds.Chip_file(QD_agent=QD_agent)
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        new_ro_elements, bias_elements = zgateT1_waiter(QD_agent,ro_elements,flux_span_period_factor,flux_data_points, evotime_data_points)
        for q in new_ro_elements["time_samples"]:
            init_system_atte(QD_agent.quantum_device,list([q]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        
        
        """ Running """
        zgateT1_executor(QD_agent,cluster,meas_ctrl,Fctrl,new_ro_elements,bias_elements, run=execution,ith=ith_histo,pi_pulse=prepare_excited,specific_folder=specific_folder,n_avg=avg_times)


        """ Analysis """
        if ith_histo == repeat-1:
            nc_paths = ZgateT1_dataReducer(specific_folder)
            for q in nc_paths:
                ds = open_dataset(nc_paths[q])
                ANA = Multiplex_analyzer("auxA")
                ANA._import_data(ds,var_dimension=2,refIQ=QD_agent.refIQ[q])
                ANA._start_analysis(time_sort=plot_time_vs_flux)
                ANA._export_result(nc_paths[q])


        """ Close """
        print('Zgate T1 done!')
        shut_down(cluster,Fctrl)
            

         

        

        

                

            
       

    
        