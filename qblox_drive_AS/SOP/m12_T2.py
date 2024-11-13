import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
import matplotlib.pyplot as plt
from xarray import Dataset
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from numpy import std, arange, array, average, mean, ndarray, pi, arange, reshape
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl, compose_para_for_multiplexing, reset_offset
from qblox_drive_AS.support.Pulse_schedule_library import multi_ramsey_sche, set_LO_frequency, pulse_preview
from qblox_drive_AS.analysis.raw_data_demolisher import T2_dataReducer
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer


def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,repeat:int=1,n_avg:int=1000,run:bool=True,exp_idx:int=0,data_folder:str='', second_phase:str='x'):
    
    sche_func= multi_ramsey_sche
    nc_path = ''
    for q in ro_elements['time_samples']:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        time_data_idx = arange(ro_elements['time_samples'][q].shape[0])

    repeat_data_idx = arange(repeat)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
    Para_repeat.batched = False
    

    
   
    sched_kwargs = dict(
        New_fxy=ro_elements["xyf"],
        freeduration=ro_elements["time_samples"],
        pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements["xyf"],'d1'),
        pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements["xyf"],'d3'),
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements["xyf"],'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements["xyf"],'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements["xyf"],'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements["xyf"],'r2'),
        echo_pi_num=ro_elements["spin_num"],
        second_pulse_phase=second_phase
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
        
        
        ds = meas_ctrl.run('Ramsey')
        # Save the raw data into netCDF
        dict_ = {}
        for q_idx, q in enumerate(ro_elements["time_samples"]):
            i_data = array(ds[f'y{2*q_idx}']).reshape(repeat,ro_elements["time_samples"][q].shape[0])
            q_data = array(ds[f'y{2*q_idx+1}']).reshape(repeat,ro_elements["time_samples"][q].shape[0])
            dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
            time_values = list(ro_elements["time_samples"][q])*2*repeat
            dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,repeat,ro_elements["time_samples"][q].shape[0]))
        
        ramsey_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":repeat_data_idx,"idx":time_data_idx})
        for var in ramsey_ds.data_vars:
            ramsey_ds[var].attrs["spin_num"] = ro_elements["spin_num"][var[:2]]
        ramsey_ds.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        ramsey_ds.attrs["execution_time"] = Data_manager().get_time_now()
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=ramsey_ds,label=exp_idx,qb="multiQ",exp_type='T2',specific_dataFolder=data_folder,get_data_loc=True)
        
        # I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
        
        # data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        
    else:
        preview_para = {}
        for q in ro_elements["time_samples"]:
            preview_para[q] = array([ro_elements["time_samples"][q][0],ro_elements["time_samples"][q][-1]])
        sched_kwargs['freeduration']= preview_para
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        
    return nc_path


def modify_time_point(ary:ndarray,factor1:int, factor2:int=0):
    x = []
    for i in ary:
        if i % factor1 == 0:
            ii = i

        else:
            multiples = i // factor1
            ii = factor1*(multiples+1)
        
        if factor2 != 0:
            if ii % factor2 == 0: 
                    if ii not in x :
                        x.append(ii)
                    else:
                        pass
            else:
                multiples = ii // factor2
                multiples_of_factor = factor2*(multiples+1)

                if multiples_of_factor % factor1 == 0:
                    if multiples_of_factor not in x :
                        x.append(multiples_of_factor)
                    else:
                        pass
        else:
            x.append(ii)


    return array(x)


def ramsey_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,ro_elements:dict,repeat:int=1,ith:int=1,run:bool=True,specific_folder:str='', avg_n:int=800, second_phase:str='x'):
    if run:
        slightly_print(f"The {ith}-th T2:")
        for q in ro_elements["time_samples"]:
            Fctrl[q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(q)))
        
        nc_path = Ramsey(QD_agent,meas_ctrl,ro_elements,repeat,n_avg=avg_n,run=True,exp_idx=ith,data_folder=specific_folder,second_phase=second_phase)
        reset_offset(Fctrl)
        cluster.reset()

    else:
        nc_path = Ramsey(QD_agent,meas_ctrl,ro_elements,repeat,n_avg=1000,run=False)

    return nc_path


def T2_waiter(QD_agent:QDmanager,ro_element:dict,time_pts:int)->dict:
    """
    The ro_element must satisfy these key_name inside the item for each qubit: `'detune'`,`'evo_T'`,`'spin_num'`. And alternatively: `'xy_IF'`.\n
    For example: `ro_element = {"q0":{'detune':5e6,'evo_T':10e-6, 'spin_num':0,'xy_IF':150e6},...}`.
    """
    ro_elements = {"xyf":{},"time_samples":{},"spin_num":{}}
    
    for q in ro_element:
        
        
        # XY-freq setting
        IF_key = [name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'][0] if len([name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'])==1 else ""      
        IF = ro_element[q][IF_key] if IF_key != "" else 150e6
        new_freq = QD_agent.quantum_device.get_element(q).clock_freqs.f01() + (ro_element[q]['detune'] if 'detune' in list(ro_element[q].keys()) else 0)
        ro_elements["xyf"][q] = new_freq
        # LO settings
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=new_freq+IF)
        # Spin echo pi number settings
        ro_elements["spin_num"][q] = ro_element[q]['spin_num'] if 'spin_num' in list(ro_element[q].keys()) else 0
        # Time sample settings
        gap = (ro_element[q]['evo_T'])*1e9 // time_pts + (((ro_element[q]['evo_T'])*1e9 // time_pts) %(4*(ro_elements["spin_num"][q]+1))) # multiple by 8 ns
        time_samples = modify_time_point(arange(0,ro_element[q]['evo_T'],gap*1e-9), ro_elements["spin_num"][q]*8e-9 if ro_elements["spin_num"][q]>=1 else 4e-9)
        ro_elements["time_samples"][q] = time_samples
        

    return ro_elements




if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 0
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        "q0":{"detune":0.5e6, 'evo_T':60e-6, 'spin_num':1, "xy_IF":250e6},
        "q1":{"detune":0.5e6, 'evo_T':60e-6, 'spin_num':2, "xy_IF":250e6},
    }
    couplers = []

    """ Optional paras """
    histo_counts:int = 1       # if > 100, use while loop outside the FPGA
    time_data_points:int = 100
    avg_n:int = 500


    """ Iteration """
    if histo_counts > 100: # use while loop
        repeat = 1
        DM = Data_manager()
        repeat_folder_path = DM.creat_datafolder_today(f"MultiQ_T2collections_{DM.get_time_now()}")
    
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
        chip_info = cds.Chip_file(QD_agent=QD_agent)
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        new_ro_elements = T2_waiter(QD_agent,ro_elements,time_data_points)
    
        
        """ Running """
        for qubit in new_ro_elements["xyf"]:
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        
        nc_path = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,new_ro_elements,repeat,ith=ith_histo,run=execution,avg_n=avg_n,specific_folder=repeat_folder_path)
        slightly_print(f"Data saved located:\n{nc_path}")

        
        """ Analysis """
        if not dont_stop:
            ds = T2_dataReducer(nc_path)
            for var in ds.data_vars:
                if var.split("_")[-1] != 'x':
                    time_data = array(ds[f"{var}_x"])[0][0]
                    ANA = Multiplex_analyzer("m12")
                    ANA._import_data(ds[var],var_dimension=2,refIQ=QD_agent.refIQ[var])
                    ANA._import_2nddata(time_data)
                    ANA._start_analysis()

                    fit_pic_folder = Data_manager().get_today_picFolder()
                    ANA._export_result(fit_pic_folder)

                    """ Storing """
                    if histo_counts >= 50:
                        QD_agent.Notewriter.save_T2_for(ANA.fit_packs["median_T2"],var)
                        QD_agent.QD_keeper()
                        if chip_info_restore:
                            chip_info.update_T2(qb=qubit, T2=f'{ANA.fit_packs["median_T2"]} +- {ANA.fit_packs["std_T2"]}')


        """ Close """
        print('T2 done!')
        shut_down(cluster,Fctrl)
        end_time = time.time()
        slightly_print(f"time cost: {round(end_time-start_time,1)} secs")
        ith_histo += 1
        
    
    
        
            
        


    