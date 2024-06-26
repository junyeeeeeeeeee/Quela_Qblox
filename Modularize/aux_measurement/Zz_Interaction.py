import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from numpy import std, arange, array, average, mean, sign
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Zz_Interaction, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot, Fit_T2_cali_analysis_plot


def excite_Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,arti_detune:int=0,IF:int=150e6,n_avg:int=1000,points:int=101,run:bool=True,q='q1', ref_IQ:list=[0,0],Experi_info:dict={},exp_idx:int=0,data_folder:str='',excite_qubit:str=''):
    
    T2_us = {}
    analysis_result = {}
    Real_detune= {}
    
    qubit = QD_agent.quantum_device.get_element(q)
    excite = QD_agent.quantum_device.get_element(excite_qubit)
    # Manually change f01
    # f01 = qubit.clock_freqs.f01()
    # qubit.clock_freqs.f01(f01-2.47e6)
    
    New_fxy= qubit.clock_freqs.f01()+arti_detune
    
    LO= New_fxy+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    gap = (freeduration)*1e9 // points + (((freeduration)*1e9 // points) %4) # multiple by 4 ns
    samples = arange(4e-9,freeduration,gap*1e-9)
    
    sche_func= Zz_Interaction
    sched_kwargs = dict(
        q=q,
        excite_qubit=excite_qubit,
        pi_amp={str(q):qubit.rxy.amp180()},
        excite_pi_amp={str(excite_qubit):excite.rxy.amp180()},
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        pi_dura=qubit.rxy.duration(),
        excite_pi_dura=excite.rxy.duration(),
        R_amp={str(q):qubit.measure.pulse_amp()},
        R_duration={str(q):qubit.measure.pulse_duration()},
        R_integration={str(q):qubit.measure.integration_time()},
        R_inte_delay=qubit.measure.acq_delay(),
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
        
        
        ramsey_ds = meas_ctrl.run('ZzInteraction')

        # Save the raw data into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=ramsey_ds,label=exp_idx,qb=q,exp_type='T2',specific_dataFolder=data_folder)
        
        I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])
        try:
            data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=8e-6)
            T2_us[q] = data_fit.attrs['T2_fit']*1e6
            Real_detune[q] = data_fit.attrs['f']
        except:
            data_fit=[]
            T2_us[q] = 0
            Real_detune[q] = 0

        analysis_result[q] = data_fit

        show_args(exp_kwargs, title="ZzInteraction_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="ZzInteraction_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    return analysis_result, T2_us, Real_detune


def excite_ramsey_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,artificial_detune:float=0e6,freeDura:float=30e-6,histo_counts:int=1,run:bool=True,plot:bool=True,specific_folder:str='',pts:int=100, avg_n:int=800,excite_qubit:str=''):
    if run:
        T2_us_rec = []
        detune_rec = []
        
        for ith in range(histo_counts):
            start_time = time.time()
            slightly_print(f"The {ith+1}-th T2:")
            Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
            Ramsey_results, T2_us, average_actual_detune= excite_Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=avg_n,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=pts,run=True,exp_idx=ith,data_folder=specific_folder,excite_qubit=excite_qubit)
            Fctrl[specific_qubits](0.0)
            cluster.reset()
            if T2_us[specific_qubits] != 0: T2_us_rec.append(T2_us[specific_qubits]) 
            if average_actual_detune[specific_qubits] != 0: detune_rec.append(average_actual_detune[specific_qubits])
            end_time = time.time()
            slightly_print(f"time cost: {round(end_time-start_time,1)} secs")
        T2_us_rec = array(T2_us_rec)
        
        if histo_counts == 1:
            mean_T2_us = T2_us_rec[0]
            sd_T2_us = 0
            if plot:
                Fit_analysis_plot(Ramsey_results[specific_qubits],P_rescale=False,Dis=None)
        # set the histo save path
        else:
            mean_T2_us = round(mean(T2_us_rec),1)
            sd_T2_us = round(std(T2_us_rec),1)
            Data_manager().save_histo_pic(QD_agent,{str(specific_qubits):T2_us_rec},specific_qubits,mode="t2",pic_folder=specific_folder)
    else:
        Ramsey_results, T2_hist, average_actual_detune= excite_Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=False,excite_qubit=excite_qubit)
        mean_T2_us = 0

    return Ramsey_results, mean_T2_us, sd_T2_us, average_actual_detune




if __name__ == "__main__":
    
    """ Fill in """
    execution = 1
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        "q1":{"detune":0e6,"evoT":50e-6,"histo_counts":1, "excite_qubit": "q0"}
    }
    couplers = ['c0','c1', 'c2', 'c3']
    # 1 = Store
    # 0 = not store
    chip_info_restore = 1

    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    """ Running """
    ramsey_results = {}
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        freeTime = ro_elements[qubit]["evoT"]
        histo_total = ro_elements[qubit]["histo_counts"]
        detuning = ro_elements[qubit]["detune"]
        excite_qubit = ro_elements[qubit]["excite_qubit"]
        plot_result = True
            

        slightly_print(f"ZzInteraction with detuning = {round(detuning*1e-6,2)} MHz")
        ramsey_results[qubit], mean_T2_us, sd_T2_us, average_actual_detune = excite_ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=detuning,freeDura=freeTime,histo_counts=histo_total,run=execution,plot=plot_result,excite_qubit=excite_qubit)
        highlight_print(f"{qubit} XYF = {round(QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()*1e-9,5)} GHz")
            
        
       
        """ Storing """
        if execution:
            if histo_total >= 50:
                QD_agent.Notewriter.save_T2_for(mean_T2_us,qubit)
                QD_agent.QD_keeper()
                # if chip_info_restore:
                    # chip_info.update_T2(qb=qubit, T2=f'{mean_T2_us} +- {sd_T2_us}')
        
        
    """ Close """
    print('T2 done!')
    shut_down(cluster,Fctrl,Cctrl)