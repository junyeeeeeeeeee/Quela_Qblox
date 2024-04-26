import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
<<<<<<< HEAD

=======
>>>>>>> origin/RatisWu
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import std, arange, array, average, mean
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot


def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,arti_detune:int=0,IF:int=150e6,n_avg:int=1000,points:int=101,run:bool=True,q='q1', ref_IQ:list=[0,0],Experi_info:dict={},exp_idx:int=0,data_folder:str=''):
    
    T2_us = {}
    analysis_result = {}
    Real_detune= {}
    
    qubit = QD_agent.quantum_device.get_element(q)
    
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
<<<<<<< HEAD
        I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])
        data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=8e-6)
        analysis_result[q] = data_fit
        T2_us[q] = data_fit.attrs['T2_fit']*1e6
        Real_detune[q] = data_fit.attrs['f']-arti_detune
=======
        
        I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])
        try:
            data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=8e-6)
            T2_us[q] = data_fit.attrs['T2_fit']*1e6
            Real_detune[q] = data_fit.attrs['f']-arti_detune
        except:
            data_fit=[]
            T2_us[q] = 0
            Real_detune[q] = 0

        analysis_result[q] = data_fit

>>>>>>> origin/RatisWu
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


def ramsey_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,artificial_detune:float=0e6,freeDura:float=30e-6,histo_counts:int=1,run:bool=True,plot:bool=True,specific_folder:str=''):
    init_system_atte(QD_agent.quantum_device,list([specific_qubits]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'xy'))
    if run:
        T2_us_rec = []
        detune_rec = []

        for ith in range(histo_counts):
            print(f"The {ith}-th T2:")
            Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_tuneawayBiasFor(specific_qubits)))
            Ramsey_results, T2_us, average_actual_detune= Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=True,exp_idx=ith,data_folder=specific_folder)
            Fctrl[specific_qubits](0.0)
            cluster.reset()
            T2_us_rec.append(T2_us[specific_qubits])
            detune_rec.append(average_actual_detune[specific_qubits])
        T2_us = array(T2_us)
        mean_T2_us = round(mean(T2_us_rec[T2_us_rec != 0]),1)
        sd_T2_us = round(std(T2_us_rec[T2_us_rec != 0]),1)
        if histo_counts == 1:
            if plot:
                Fit_analysis_plot(Ramsey_results[specific_qubits],P_rescale=False,Dis=None)
        # set the histo save path
        else:
            Data_manager().save_histo_pic(QD_agent,{str(specific_qubits):T2_us_rec},specific_qubits,mode="t2",T1orT2=f"{mean_T2_us}+/-{sd_T2_us}_4.4G",pic_folder=specific_folder)
    else:
        Ramsey_results, T2_hist, average_actual_detune= Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=False)
        mean_T2_us = 0

    return Ramsey_results, mean_T2_us, average_actual_detune




if __name__ == "__main__":
    
    """ Fill in """
    execution = 1
    xyf_cali = 0
<<<<<<< HEAD
    QD_path = 'Modularize/QD_backup/2024_4_25/DR2#10_SumInfo.pkl'
    ro_elements = {
        "q1":{"detune":-0.2e6,"evoT":50e-6,"histo_counts":1}
=======
    QD_path = 'Modularize/QD_backup/2024_4_25/DR1#11_SumInfo.pkl'
    ro_elements = {
        "q0":{"detune":0e6,"evoT":40e-6,"histo_counts":1}
>>>>>>> origin/RatisWu
    }


    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    

    """ Running """
    Trustable = False # don't change
    ramsey_results = {}
    for qubit in ro_elements:
        freeTime = ro_elements[qubit]["evoT"]
        histo_total = ro_elements[qubit]["histo_counts"]
        if xyf_cali:
            manual_detune = [ro_elements[qubit]["detune"],-ro_elements[qubit]["detune"]]
            plot_result = False
            if histo_total != 1: raise ValueError("Histo_counts should be 1 in the XYF_calibration mode!") 
        else:
            manual_detune = [ro_elements[qubit]["detune"]]
            plot_result = True
            
        
        actual_detune = []
        for detuning in manual_detune:
            print(f"Ramsey with detuning = {round(detuning*1e-6,2)} MHz")
            ramsey_results[qubit], mean_T2_us, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=detuning,freeDura=freeTime,histo_counts=histo_total,run=execution,plot=plot_result)
            actual_detune.append(average_actual_detune[qubit])

        if xyf_cali:
            if average(array(actual_detune))<=abs(manual_detune[0]):
                Trustable = True
                original_xyf = QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()
                QD_agent.quantum_device.get_element(qubit).clock_freqs.f01(original_xyf+average(array(actual_detune)))
            else:
                print("Warning: Please set a larger detuning !")
        
        if histo_total >= 10:
            Trustable = True
            QD_agent.Notewriter.save_T2_for(mean_T2_us,qubit)

    
    """ Storing """
    if execution:
        if Trustable:
            QD_agent.QD_keeper()
        
        
    """ Close """
    print('T2 done!')
    shut_down(cluster,Fctrl)