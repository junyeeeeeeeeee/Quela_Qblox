
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import arange, array, average, mean
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot


def Ramsey(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,arti_detune:int=0,IF:int=130e6,n_avg:int=1000,points:int=101,run:bool=True,q='q1',times=1, ref_IQ:list=[0,0],Experi_info:dict={}):
    
    T2_us = {}
    analysis_result = []
    T2_us[q] = []
    Real_detune= {}
    Real_detune[q]= []
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
        
        for i in range(times):
            ramsey_ds = meas_ctrl.run('Ramsey')

            # Save the raw data into netCDF
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=ramsey_ds,histo_label=i,qb=q,exp_type='T2')

            I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
            data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])
            data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=8e-6)
            analysis_result.append(data_fit)
            T2_us[q].append(data_fit.attrs['T2_fit']*1e6)
            Real_detune[q].append(data_fit.attrs['f']-arti_detune)
        
        average_fit_detune= average(array(Real_detune[q]))
        

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
        
    return analysis_result, T2_us, average_fit_detune


def ramsey_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,artificial_detune:float=0e6,freeDura:float=30e-6,histo_counts:int=1,run:bool=True):
    init_system_atte(QD_agent.quantum_device,list([specific_qubits]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'xy'))
    linecut = 0
    if run:
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_sweetBiasFor(specific_qubits)))
        Ramsey_results, T2_hist, average_actual_detune= Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,times=histo_counts,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=True)
        Fctrl[specific_qubits](0.0)

        Fit_analysis_plot(Ramsey_results[linecut],P_rescale=False,Dis=None)
        # set the histo save path
        Data_manager().save_histo_pic(QD_agent,T2_hist,specific_qubits,mode="t2")
        mean_T2_us = mean(array(T2_hist[specific_qubits]))
    
    else:
        Ramsey_results, T2_hist, average_actual_detune= Ramsey(QD_agent,meas_ctrl,arti_detune=artificial_detune,freeduration=freeDura,n_avg=1000,q=specific_qubits,times=1,ref_IQ=QD_agent.refIQ[specific_qubits],points=100,run=False)
        mean_T2_us = 0

    return Ramsey_results, mean_T2_us, average_actual_detune




if __name__ == "__main__":
    
    """ Fill in """
    execution = True
    modify_xyf= True
    QD_path = 'Modularize/QD_backup/2024_3_28/DR2#171_SumInfo.pkl'
    ro_elements = {
        "q4":{"detune":0,"evoT":40e-6,"histo_counts":1}
    }


    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    

    """ Running """
    ramsey_results = {}
    for qubit in ro_elements:
        manual_detune = ro_elements[qubit]["detune"]
        evo_time = ro_elements[qubit]["evoT"]
        histo_total = ro_elements[qubit]["histo_counts"]

        ramsey_results[qubit], mean_T2_us, average_actual_detune = ramsey_executor(QD_agent,meas_ctrl,Fctrl,qubit,artificial_detune=manual_detune,freeDura=evo_time,histo_counts=histo_total,run=execution)
        cluster.reset()

        if modify_xyf:
            original_xyf = QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()
            QD_agent.quantum_device.get_element(qubit).clock_freqs.f01(original_xyf-average_actual_detune)

    

    """ Storing """
    if modify_xyf:
        QD_agent.QD_keeper()



    """ Close """
    print('T2 done!')
    shut_down(cluster,Fctrl)