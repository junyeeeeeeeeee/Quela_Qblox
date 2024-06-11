"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import linspace, array, arange, NaN
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Rabi_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, Rabi_fit_analysis, Fit_analysis_plot

def Rabi(QD_agent:QDmanager,meas_ctrl:MeasurementControl,XY_amp:float=0.5, XY_duration:float=20e-9, IF:int=150e6,n_avg:int=300,points:int=100,run:bool=True,XY_theta:str='X_theta',Rabi_type:str='PowerRabi',q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0]):
    analysis_result = {}
    sche_func= Rabi_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    LO= qubit_info.clock_freqs.f01()+IF

    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    if Rabi_type.lower() in ['timerabi', 'tr']:
        osci_type = "TimeRabi"
        if qubit_info.rxy.amp180() is not NaN or qubit_info.rxy.amp180() != 0:
            Para_XY_amp= qubit_info.rxy.amp180() 
        else:
            raise ValueError(f"Your rxy.amp180 is {qubit_info.rxy.amp180()} can't perform a time Rabi!")
        Sweep_para=Para_XY_Du = ManualParameter(name="XY_Duration", unit="s", label="Time")
        str_Rabi= 'XY_duration'
        Sweep_para.batched = True
        if XY_duration < 200e-9:
            XY_duration = 200e-9
            resolu = 4e-9
        else:
            resolu = 8e-9

        gap = (XY_duration)*1e9 // points + (((XY_duration)*1e9 // points) %(resolu*1e9)) # multiple by 4 ns
        samples = arange(resolu,XY_duration,gap*1e-9)
        exp_kwargs= dict(sweep_duration=[osci_type,'start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                        Amp='%E' %XY_amp,
                        )
    elif Rabi_type.lower() in ['powerrabi', 'pr']:
        osci_type = "PowerRabi"
        Sweep_para= Para_XY_amp= ManualParameter(name="XY_amp", unit="V", label="Voltage")
        str_Rabi= 'XY_amp'
        Para_XY_amp.batched = True
        Para_XY_Du = XY_duration
        samples = linspace(-XY_amp,XY_amp,points) 
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
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=rabi_ds,qb=q,exp_type=osci_type)
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
    

def rabi_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,XYamp_max:float=0.5,XYdura_max:float=20e-9,which_rabi:str='power',run:bool=True,pts:int=100):
    if which_rabi.lower() in ['p','power']:
        exp_type = 'powerRabi'
    elif which_rabi.lower() in ['t','time']:
        exp_type = 'timeRabi'
    else:
        raise ValueError (f"Can't recognized rabi type with the given = {which_rabi}")
    
    print(f"{specific_qubits} are under the measurement ...")
    trustable = False
    if run:
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        Rabi_results = Rabi(QD_agent,meas_ctrl,Rabi_type=exp_type,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=True,XY_amp=XYamp_max,XY_duration=XYdura_max,points=pts)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
        if Rabi_results == {}:
            print(f"Rabi Osci error qubit: {specific_qubits}")
        else:
            Fit_analysis_plot(Rabi_results[specific_qubits],P_rescale=False,Dis=None)
            if abs(Rabi_results[specific_qubits].attrs['pi_2']) <= 0.5:
                qubit = QD_agent.quantum_device.get_element(specific_qubits)
                qubit.rxy.amp180(Rabi_results[specific_qubits].attrs['pi_2'])
                qubit.rxy.duration(XYdura_max)
                trustable = True
        return Rabi_results[specific_qubits], trustable
    else:
        Rabi_results = Rabi(QD_agent,meas_ctrl,Rabi_type=exp_type,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=False,XY_amp=XYamp_max,XY_duration=XYdura_max)
        return 0, 0
  

if __name__ == "__main__":
    
    """ Fill in """
    execution = True
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = ['q1']
    couplers = ["c0",'c1']
    # 1 = Store
    # 0 = not store
    chip_info_restore = 1

    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    """Running """
    rabi_results = {}
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        rabi_results[qubit], trustable = rabi_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,run=execution,XYdura_max=40e-9,XYamp_max=0.6)
        cluster.reset()
    
        """ Storing """
        if trustable:   
            QD_agent.refresh_log("after Rabi")
            QD_agent.QD_keeper()
            if chip_info_restore:
                chip_info.update_RabiOsci(qb=qubit)



    """ Close """
    print('Rabi osci done!')
    shut_down(cluster,Fctrl,Cctrl)