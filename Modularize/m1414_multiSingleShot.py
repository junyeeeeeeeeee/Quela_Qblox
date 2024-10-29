import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from qblox_instruments import Cluster
from Modularize.support.UserFriend import *
from numpy import array, linspace, median, std, moveaxis, arange
from quantify_scheduler.gettables import ScheduleGettable
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import QDmanager, Data_manager,init_system_atte, init_meas, shut_down, coupler_zctrl, compose_para_for_multiplexing,reset_offset
from Modularize.support.Pulse_schedule_library import multi_Qubit_SS_sche, set_LO_frequency, pulse_preview
from Modularize.analysis.Multiplexing_analysis import Multiplex_analyzer
from Modularize.analysis.raw_data_demolisher import OneShot_dataReducer


def Qubit_state_single_shot(QD_agent:QDmanager,ro_elements:dict,shots:int=1000,run:bool=True,exp_idx:int=0,parent_datafolder:str=''):
    sche_func = multi_Qubit_SS_sche 

    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.measure.integration_time(1e-6)
        print("Integration time ",qubit_info.measure.integration_time()*1e6, "µs")
        print("Reset time ", qubit_info.reset.duration()*1e6, "µs")

        qubit_info.measure.pulse_amp(float(ro_elements[q]["ro_amp_factor"])*qubit_info.measure.pulse_amp())
        if float(ro_elements[q]["ro_amp_factor"]) != 1:
            eyeson_print(f"The new RO amp = {round(qubit_info.measure.pulse_amp(),2)}")
    
    folder = []

    def state_dep_sched(ini_state:str):
        slightly_print(f"Shotting for |{ini_state}>")
        sched_kwargs = dict(   
            ini_state=ini_state,
            pi_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'d1'),
            pi_dura=compose_para_for_multiplexing(QD_agent,ro_elements,'d3'),
            R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
            R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
            R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
        )
        
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(ro_elements.keys())),
            )
            QD_agent.quantum_device.cfg_sched_repetitions(shots)
            ss_da= gettable.get() # DataArray (2*ro_q, shots)
            reshaped_data = list(array(ss_da).reshape(len(list(ro_elements.keys())),2,shots)) # (ro_q, IQ, shots)
            folder.append(reshaped_data) # (state, ro_q, IQ, shots)

        else:
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

            
    state_dep_sched('g')
    state_dep_sched('e')
    folder = moveaxis(array(folder),0,2) # (state,ro_q,IQ,shots) -> (ro_q, IQ, state, shots)
    output_dict = {}
    for q_idx, q_name in enumerate(ro_elements):
        output_dict[q_name] = (["mixer","prepared_state","index"],folder[q_idx])

    SS_ds = Dataset(output_dict, coords= {"mixer":array(["I","Q"]), "prepared_state":array([0,1]),"index":arange(shots)})
    SS_ds.attrs["execution_time"] = Data_manager().get_time_now()

    if run:
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=SS_ds,qb="multiQ",exp_type='ss',label=exp_idx,specific_dataFolder=parent_datafolder,get_data_loc=True)
    else:
        nc_path = ""
    
    return nc_path


def SS_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict, ro_elements:dict, shots:int=10000,execution:bool=True,data_folder='',exp_label:int=0):
    
    for q in ro_elements:
        Fctrl[q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(q)))

    nc_path = Qubit_state_single_shot(QD_agent,ro_elements,shots=shots,run=execution,parent_datafolder=data_folder,exp_idx=exp_label)
    reset_offset(Fctrl)
    cluster.reset()
    
    # thermal_p, effT_mk, ro_fidelity, rotate_angle = a_OSdata_analPlot(QD_agent,target_q,nc,plot,save_pic=save_every_pic)
    return nc_path


def SS_waiter(QD_agent:QDmanager,ro_element:dict, repeat:int=1)->tuple[dict,bool]:
    auto_save:bool = True
    ro_elements = {}
    for q in ro_element:
        # RO-amp
        if 'ro_amp_factor' in list(ro_element[q].keys()):
            ro_elements[q] = ro_element[q]
            if ro_element[q]['ro_amp_factor'] != 1 :
                auto_save = False
                if repeat != 1:
                    raise ValueError(f"ro_amp_amplifier for {q} is not 1 so you can not repeat the exp.")
        else:
            ro_elements[q] = {"ro_amp_factor":1}

        # LO setting
        IF_key = [name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'][0] if len([name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'])==1 else ""      
        IF = ro_element[q][IF_key] if IF_key != "" else 150e6
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=QD_agent.quantum_device.get_element(q).clock_freqs.f01()+IF)

    return ro_elements, auto_save


if __name__ == '__main__':
    
    """ Fill in """
    execute:bool = True
    repeat:int = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        'q0':{"ro_amp_factor":1,"xy_IF":150e6},
        'q1':{"ro_amp_factor":1,"xy_IF":150e6}
    } 
    couplers = []


    """ Optional paras (don't use is better) """
    shot_num:int = 10000
    


    """ Iteration """
    snr_rec, effT_rec, thermal_pop = {}, {}, {}
    
    for i in range(repeat):
        start_time = time.time()

        """ Preparation """
        slightly_print(f"The {i}th OS:")
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
        Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        ro_elements, auto_save = SS_waiter(QD_agent,ro_elements,repeat)
        for qubit in ro_elements:
            init_system_atte(QD_agent.quantum_device,list([qubit]),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))


        """ Running """
        nc_path = SS_executor(QD_agent,cluster,Fctrl,ro_elements,execution=execute,shots=shot_num,exp_label=i)
        slightly_print(nc_path)
        
        """ Analysis """
        ds = OneShot_dataReducer(nc_path)
        ANA = Multiplex_analyzer("m1414")
        for var in ds.data_vars:
            ANA._import_data(ds[var],var_dimension=0,fq_Hz=QD_agent.quantum_device.get_element(var).clock_freqs.f01())
            ANA._start_analysis()
            pic_path = os.path.join(Data_manager().get_today_picFolder(),f"{var}_SingleShot_{ds.attrs['execution_time']}")
            ANA._export_result(pic_path)


        """ Storing """ 
        if execute and repeat == 1:
            if auto_save:
                QD_agent.QD_keeper() 
                
                
        """ Close """    
        shut_down(cluster,Fctrl,Cctrl)
        end_time = time.time()
        slightly_print(f"time cose: {round(end_time-start_time,1)} secs")

    # for qubit in effT_rec:
    #     highlight_print(f"{qubit}: {round(median(array(effT_rec[qubit])),1)} +/- {round(std(array(effT_rec[qubit])),1)} mK")
    #     Data_manager().save_histo_pic(QD_agent,effT_rec,qubit,mode="ss")
    #     Data_manager().save_histo_pic(QD_agent,thermal_pop,qubit,mode="pop")
        
        

        
    
