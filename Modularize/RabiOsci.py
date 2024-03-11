"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""

from numpy import linspace, array
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Pulse_schedule_library import Rabi_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, Rabi_fit_analysis

def Rabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,XY_amp:float, XY_duration:float=20e-9, IF:int=150e6,n_avg:int=300,points:int=200,run:bool=True,XY_theta:str='X_theta',Rabi_type:str='PowerRabi',q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0]):
    analysis_result = {}
    sche_func= Rabi_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    LO= qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    if Rabi_type.lower() in ['timerabi', 'tr']:
       osci_type = "TimeRabi"
       Para_XY_amp= XY_amp
       Sweep_para=Para_XY_Du = ManualParameter(name="XY_Duration", unit="s", label="Time")
       str_Rabi= 'XY_duration'
       Sweep_para.batched = True
       samples = linspace(0, XY_duration,points)
       exp_kwargs= dict(sweep_duration=[osci_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                        Amp='%E' %XY_amp,
                        )
    elif Rabi_type.lower() in ['powerabi', 'pr']:
        osci_type = "PowerRabi"
        Sweep_para= Para_XY_amp= ManualParameter(name="XY_amp", unit="V", label="Voltage")
        str_Rabi= 'XY_amp'
        Para_XY_amp.batched = True
        Para_XY_Du = XY_duration
        samples = linspace(0,XY_amp,points) 
        exp_kwargs= dict(sweep_amp=[osci_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                         Duration='%E' %XY_duration,
                         )
    else: raise KeyError ('Typing error: Rabi_type')
    
    sched_kwargs = dict(
        q=q,
        XY_amp=Para_XY_amp,
        XY_duration=Para_XY_Du,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        XY_theta=XY_theta,
        Rabi_type=osci_type,
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
        meas_ctrl.settables(Sweep_para)
        meas_ctrl.setpoints(samples)
    
       
        rabi_ds = meas_ctrl.run(osci_type)
        # Save the raw data into netCDF
        Data_manager.save_raw_data(QD_agent=QD_agent,ds=rabi_ds,qb=q,exp_type=osci_type)
        I,Q= dataset_to_array(dataset=rabi_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        data_fit= Rabi_fit_analysis(data=data,samples=samples,Rabi_type=osci_type)
        analysis_result[q]= data_fit
        
        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sched_kwargs[str_Rabi]= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
       

        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))

    return analysis_result
    

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from Pulse_schedule_library import Fit_analysis_plot
    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_5/DR1#170_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    # Set system attenuation
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=20)
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50) 

    execute = True
    error_log = []
    for qb in Fctrl:
        print(f"{qb} are under the measurement ...")
        Rabi_results = Rabi(QD_agent,meas_ctrl,XY_amp=0.5,q=qb,ref_IQ=QD_agent.refIQ[qb])
        if Rabi_results == {}:
            error_log.append(qb)
        else:
            if execute:
                qubit = QD_agent.quantum_device.get_element(qb)
                Fit_analysis_plot(Rabi_results[qb],P_rescale=False,Dis=None)
                qubit.rxy.amp180(Rabi_results[qb].attrs['pi_2'])
    
    if error_log != []:
        print(f"Rabi Osci error qubit: {error_log}")
    if execute:
        QD_agent.refresh_log("after Rabi")
        QD_agent.QD_keeper()
    print('Rabi osci done!')
    shut_down(cluster,Fctrl)