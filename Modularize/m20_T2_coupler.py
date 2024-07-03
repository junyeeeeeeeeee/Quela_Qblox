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
from Modularize.support.Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot, Fit_T2_cali_analysis_plot


def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,arti_detune:int=0,IF:int=150e6,n_avg:int=1000,points:int=101,run:bool=True,q='q1', ref_IQ:list=[0,0],Experi_info:dict={},exp_idx:int=0,data_folder:str=''):
    
    T2_us = {}
    analysis_result = {}
    Real_detune= {}
    
    qubit = QD_agent.quantum_device.get_element(q)

    # Manually change f01
    # qubit.clock_freqs.f01(qubit.clock_freqs.f01()+1.343e6)
    
    New_fxy= qubit.clock_freqs.f01()+arti_detune
    
    LO= New_fxy+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    gap = (freeduration)*1e9 // points + (((freeduration)*1e9 // points) %4) # multiple by 4 ns
    samples = arange(4e-9,freeduration,gap*1e-9)
    
    sche_func= Ramsey_sche
    sched_kwargs = dict(
        q=q,
        pi_amp={str(q):qubit.rxy.amp180()},
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        pi_dura=qubit.rxy.duration(),
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
        
        
        ramsey_ds = meas_ctrl.run('Ramsey')

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

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    return analysis_result, T2_us, Real_detune


def ramsey_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,artificial_detune:float=0e6,freeDura:float=30e-6,ith:int=1,run:bool=True,specific_folder:str='',pts:int=100, avg_n:int=800):
    if run:
        start_time = time.time()
        slightly_print(f"The {ith}-th T2:")
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        Ramsey_results, T2_us, average_actual_detune = Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=avg_n,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=pts,run=True,exp_idx=ith,data_folder=specific_folder)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
        this_t2_us = T2_us[specific_qubits]
        end_time = time.time()
        slightly_print(f"time cost: {round(end_time-start_time,1)} secs")
        
    else:
        Ramsey_results, _, average_actual_detune = Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=False)
        this_t2_us = 0

    return Ramsey_results, this_t2_us, average_actual_detune




if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        "q0":{"detune":-22e6,"evoT":4e-6,"histo_counts":1},
        # "q1":{"detune":0.8e6,"evoT":50e-6,"histo_counts":1},
    }
    target_coupler = 'c0'
    target_bias = 0.0800
    couplers = ['c1','c2','c3']
    pts =1000

    """ Iteration """
    for qubit in ro_elements:
        t2_us_rec = []
        for ith_histo in range(ro_elements[qubit]["histo_counts"]):
            """ Preparations """
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            chip_info = cds.Chip_file(QD_agent=QD_agent)


            """ Running """
            cp_ctrl = QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i')
            cp_ctrl[target_coupler]= target_bias
            Cctrl = coupler_zctrl(DRandIP["dr"],cluster,cp_ctrl)
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            
            slightly_print(f"Ramsey with detuning = {round(ro_elements[qubit]['detune']*1e-6,2)} MHz")
            ramsey_results, this_t2_us, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=ro_elements[qubit]["detune"],freeDura=ro_elements[qubit]["evoT"],pts=pts,ith=ith_histo,run=execution)
            highlight_print(f"{qubit} XYF = {round(QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()*1e-9,7)} GHz")
            if this_t2_us != 0:
                t2_us_rec.append(this_t2_us)
            

            """ Close """
            print('T2 done!')
            shut_down(cluster,Fctrl,Cctrl)
            
        
        """ Storing """
        if execution:
            mean_T2_us = round(mean(array(t2_us_rec)),2)
            std_T2_us  = round(std(array(t2_us_rec)),2)
            highlight_print(f"{qubit}: mean T2 = {mean_T2_us} 土 {std_T2_us} µs")
            # QD_agent.QD_keeper()# Manual keep
            if ro_elements[qubit]["histo_counts"] == 1:
                mean_T2_us = t2_us_rec[0]
                sd_T2_us = 0
                Fit_analysis_plot(ramsey_results[qubit],P_rescale=False,Dis=None)
            else:
                Data_manager().save_histo_pic(QD_agent,{str(qubit):t2_us_rec},qubit,mode="t2")
                if ro_elements[qubit]["histo_counts"] >= 50:
                    QD_agent.Notewriter.save_T2_for(mean_T2_us,qubit)
                    # QD_agent.QD_keeper()
                    if chip_info_restore:
                        chip_info.update_T2(qb=qubit, T2=f'{mean_T2_us} +- {std_T2_us}')
            
            
        


    