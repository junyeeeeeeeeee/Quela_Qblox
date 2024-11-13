import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from xarray import Dataset
from numpy import mean, array, arange, std
from utils.tutorial_utils import show_args
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import QDmanager, Data_manager, multiples_of_x, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing, reset_offset
from qblox_drive_AS.support.Pulse_schedule_library import multi_T1_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T1_fit_analysis, Fit_analysis_plot
from qblox_drive_AS.analysis.raw_data_demolisher import T1_dataReducer
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer

def T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,repeat:int=1,n_avg:int=300,run:bool=True,exp_idx:int=0,data_folder:str=''):
    sche_func= multi_T1_sche
    nc_path = ""

    for q in ro_elements["time_samples"]:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        time_data_idx = arange(ro_elements["time_samples"][q].shape[0])
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
    Para_repeat.batched = False
    repeat_data_idx = arange(repeat)

    sched_kwargs = dict(
        freeduration=ro_elements["time_samples"],
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
        meas_ctrl.settables([Para_free_Du,Para_repeat])
        meas_ctrl.setpoints_grid((time_data_idx,repeat_data_idx))

        ds = meas_ctrl.run('T1')
        dict_ = {}
        for q_idx, q in enumerate(ro_elements["time_samples"]):
            i_data = array(ds[f'y{2*q_idx}']).reshape(repeat,ro_elements["time_samples"][q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(repeat,ro_elements["time_samples"][q].shape[0])
            dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
            time_values = list(ro_elements["time_samples"][q])*2*repeat
            dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,repeat,ro_elements["time_samples"][q].shape[0]))
        
        T1_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":repeat_data_idx,"idx":time_data_idx})
        T1_ds.attrs["execution_time"] = Data_manager().get_time_now()
        T1_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=T1_ds,label=exp_idx,qb="multiQ",exp_type='T1',specific_dataFolder=data_folder,get_data_loc=True)
        

    else:
        preview_para = {}
        for q in ro_elements["time_samples"]:
            preview_para[q] = array([ro_elements["time_samples"][q][0],ro_elements["time_samples"][q][-1]])
        sched_kwargs['freeduration']= preview_para
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
    
    return nc_path


def T1_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,ro_elements:dict,run:bool=True,repeat:int=1,specific_folder:str='',ith:int=0,avg_times:int=500):
    if run:
        slightly_print(f"The {ith}-th T1:")
        for q in ro_elements["time_samples"]:
            Fctrl[q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(q)))
        
        nc_path = T1(QD_agent,meas_ctrl,ro_elements,repeat,run=True,exp_idx=ith,data_folder=specific_folder,n_avg=avg_times)
        reset_offset(Fctrl)
        cluster.reset()

    else:
        nc_path = T1(QD_agent,meas_ctrl,ro_elements,repeat=1,run=False,exp_idx=ith,data_folder=specific_folder,n_avg=avg_times)
        

    return nc_path

def T1_waiter(QD_agent:QDmanager,ro_element:dict,time_pts:int=100)->dict:
    """ `'evoT'` and `'xy_IF'` should be there in the ro_element """
    ro_elements = {"time_samples":{}}
    
    for q in ro_element:
        # LO settings
        IF_key = [name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'][0] if len([name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'])==1 else ""      
        IF = ro_element[q][IF_key] if IF_key != "" else 150e6
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=QD_agent.quantum_device.get_element(q).clock_freqs.f01()+IF)

        # evolution time settings 
        gap = (ro_element[q]["evoT"]*1e9 // time_pts) + ((ro_element[q]["evoT"]*1e9 // time_pts)%4)
        ro_elements["time_samples"][q] = arange(0,ro_element[q]["evoT"],gap*1e-9)

    return ro_elements


if __name__ == "__main__":
    

    """ Fill in """
    execution:bool = 0
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        "q0":{"evoT":120e-6,"xy_IF":250e6},
        "q1":{"evoT":100e-6,"xy_IF":250e6},
    }
    couplers = []

    """ Optional paras """
    histo_counts:int = 1 # > 100 use while loop outside the FPGA
    time_data_points:int = 100
    avg_n = 500
  

    """ Iterations """
    if histo_counts > 100: # use while loop
        repeat = 1
        DM = Data_manager()
        repeat_folder_path = DM.creat_datafolder_today(f"MultiQ_T1collections_{DM.get_time_now()}")
    
    ith_histo = 0
    dont_stop = True
    while dont_stop:
        if histo_counts <= 100 :
            repeat = histo_counts # repeat into FPGA
            dont_stop = False
            repeat_folder_path = ''
        start_time = time.time()


        """ Preparations """
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
        for q in ro_elements:
            init_system_atte(QD_agent.quantum_device,list([q]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        chip_info = cds.Chip_file(QD_agent=QD_agent)
        ro_elements = T1_waiter(QD_agent,ro_elements,time_data_points)
        

        """ Running """
        nc_path = T1_executor(QD_agent,cluster,meas_ctrl,Fctrl,ro_elements,run=execution,repeat=repeat,ith=ith_histo,avg_times=avg_n,specific_folder=repeat_folder_path)
        slightly_print(f"T1 data saved located:\n{nc_path}")

        
        """ Analysis """
        if not dont_stop:
            ds = T1_dataReducer(nc_path)
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    time_data = array(ds[f"{var}_x"])[0][0]
                    ANA = Multiplex_analyzer("m13")
                    ANA._import_data(ds[var],var_dimension=2,refIQ=QD_agent.refIQ[var])
                    ANA._import_2nddata(time_data)
                    ANA._start_analysis()

                    fit_pic_folder = Data_manager().get_today_picFolder()
                    ANA._export_result(fit_pic_folder)

                    """ Storing """
                    if histo_counts >= 50:
                        QD_agent.Notewriter.save_T1_for(ANA.fit_packs["median_T1"],var)
                        QD_agent.QD_keeper()
                        if chip_info_restore:
                            chip_info.update_T1(qb=var, T1=f'{ANA.fit_packs["median_T1"]} +- {ANA.fit_packs["std_T1"]}')

                
        """ Close """
        print('T1 done!')
        every_end = time.time()
        shut_down(cluster,Fctrl)
        slightly_print(f"time cost: {round(every_end-start_time,1)} secs")
        ith_histo += 1
        
        


            
       

    
        