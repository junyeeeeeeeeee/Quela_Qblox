import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
import matplotlib.pyplot as plt
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from numpy import std, arange, array, average, mean, ndarray, pi
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot, Fit_T2_cali_analysis_plot, T1_fit_analysis


def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,arti_detune:int=0,IF:int=250e6,n_avg:int=1000,points:int=101,run:bool=True,q='q1', ref_IQ:list=[0,0],Experi_info:dict={},exp_idx:int=0,data_folder:str='',spin:int=0, second_phase:str='x'):
    
    T2_us = {}
    analysis_result = {}
    Real_detune= {}
    
    qubit_info = QD_agent.quantum_device.get_element(q)

    
    # qubit_info.reset.duration(qubit_info.reset.duration()*2)
    print("Integration time ",qubit_info.measure.integration_time()*1e6, "µs")
    print("Reset time ", qubit_info.reset.duration()*1e6, "µs")
    # Manually change f01
    # f01 = qubit.clock_freqs.f01()
    # qubit.clock_freqs.f01(f01-2.47e6)
    
    New_fxy= qubit_info.clock_freqs.f01()+arti_detune
    
    LO= New_fxy+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    
    gap = (freeduration)*1e9 // points + (((freeduration)*1e9 // points) %(4*(spin+1))) # multiple by 8 ns
    
    if spin >= 1:
        samples = arange(0,freeduration,gap*1e-9)
        samples = modify_time_point(samples, spin*8e-9)
    else:
        samples = arange(0,freeduration,gap*1e-9)
        samples = modify_time_point(samples, 4e-9)
    
    sche_func= Ramsey_sche
    sched_kwargs = dict(
        q=q,
        pi_amp={str(q):qubit_info.rxy.amp180()},
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        pi_dura=qubit_info.rxy.duration(),
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        echo_pi_num=spin,
        second_pulse_phase=second_phase
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                     f_xy='%E' %sched_kwargs['New_fxy'],
                     )
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
        )
        
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(Para_free_Du)
        meas_ctrl.setpoints(samples)
        
        
        ramsey_ds = meas_ctrl.run('Ramsey')

        # Save the raw data into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=ramsey_ds,label=exp_idx,qb=q,exp_type='T2',specific_dataFolder=data_folder)
        
        I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
        
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        
        try:
            if spin == 0:
                data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=20e-6)
                phase = round(data_fit.attrs['phase']*180/pi,1)
                T2_us[q] = data_fit.attrs['T2_fit']*1e6
                Real_detune[q] = data_fit.attrs['f']
            else:
                data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=40e-6)
                T2_us[q] = data_fit.attrs['T1_fit']*1e6
                Real_detune[q] = 0
        except:
            warning_print("T2 fitting error")
        
            data_fit=[]
            T2_us[q] = 0
            Real_detune[q] = 0

        analysis_result[q] = data_fit

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        # n_s = 2
        sweep_para= array([samples[0],samples[-1]])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    return analysis_result, T2_us, Real_detune


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


def ramsey_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,artificial_detune:float=0e6,freeDura:float=30e-6,ith:int=1,run:bool=True,specific_folder:str='',pts:int=100, avg_n:int=800, spin_echo:int=0, IF:float=250e6, second_phase:str='x'):
    if run:
        qubit_info = QD_agent.quantum_device.get_element(specific_qubits)
        ori_reset = qubit_info.reset.duration()
        qubit_info.reset.duration(qubit_info.reset.duration()+freeDura)
    
        slightly_print(f"The {ith}-th T2:")
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        
        Ramsey_results, T2_us, average_actual_detune = Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=avg_n,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=pts,run=True,exp_idx=ith,data_folder=specific_folder,spin=spin_echo,IF=IF,second_phase=second_phase)
        Fctrl[specific_qubits](0.0)
        
        cluster.reset()
        this_t2_us = T2_us[specific_qubits]
        qubit_info.reset.duration(ori_reset)
        
    else:
        Ramsey_results, _, average_actual_detune = Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=False,spin=spin_echo)
        this_t2_us = 0

    return Ramsey_results, this_t2_us, average_actual_detune




if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        "q0":{"detune":-1e6,"evoT":30e-6,"histo_counts":1},#-0.174e6
    }
    couplers = ["c0","c1"]

    """ Optional paras """
    spin_echo_pi_num:int = 0
    time_data_points = 100
    avg_n = 1000
    xy_IF = 250e6


    """ Iteration """
    for qubit in ro_elements:
        t2_us_rec = []
        for ith_histo in range(ro_elements[qubit]["histo_counts"]):
            start_time = time.time()
            """ Preparations """
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
            chip_info = cds.Chip_file(QD_agent=QD_agent)
            
            """ Running """
            Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            slightly_print(f"Ramsey with detuning = {round(ro_elements[qubit]['detune']*1e-6,2)} MHz")
            ramsey_results, this_t2_us, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=ro_elements[qubit]["detune"],freeDura=ro_elements[qubit]["evoT"],ith=ith_histo,run=execution,pts=time_data_points,spin_echo=spin_echo_pi_num,avg_n=avg_n,IF=xy_IF)
            highlight_print(f"{qubit} XYF = {round(QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()*1e-9,5)} GHz")
            if this_t2_us > 0:
                t2_us_rec.append(this_t2_us)
            else:
                warning_print("T2 fitting got minus value")

            print(this_t2_us)
            
            """ Storing """
            if ith_histo == int(ro_elements[qubit]["histo_counts"])-1:
                if execution:
                    mean_T2_us = round(mean(array(t2_us_rec)),2)
                    std_T2_us  = round(std(array(t2_us_rec)),2)
                    highlight_print(f"{qubit}: mean T2 = {mean_T2_us} 土 {std_T2_us} µs")

                    if ro_elements[qubit]["histo_counts"] == 1:
                        mean_T2_us = t2_us_rec[0]
                        sd_T2_us = 0
                        Fit_analysis_plot(ramsey_results[qubit],P_rescale=False,Dis=None,spin_echo=True if spin_echo_pi_num ==1 else False)
                    else:
                        if spin_echo_pi_num == 0:
                            Data_manager().save_histo_pic(QD_agent,{str(qubit):t2_us_rec},qubit,mode="t2*")
                        else:
                            Data_manager().save_histo_pic(QD_agent,{str(qubit):t2_us_rec},qubit,mode="t2")
                        
                        if ro_elements[qubit]["histo_counts"] >= 50:
                            QD_agent.Notewriter.save_T2_for(mean_T2_us,qubit)
                            QD_agent.QD_keeper()
                            if chip_info_restore:
                                chip_info.update_T2(qb=qubit, T2=f'{mean_T2_us} +- {std_T2_us}')

            """ Close """
            print('T2 done!')
            shut_down(cluster,Fctrl)
            end_time = time.time()
            slightly_print(f"time cost: {round(end_time-start_time,1)} secs")
            
        
        
            
            
        


    