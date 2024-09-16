"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from Modularize.support.UserFriend import eyeson_print, slightly_print
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import linspace, array, arange, NaN
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Rabi_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, Rabi_fit_analysis, Fit_analysis_plot
from Modularize.analysis.RabiChevAna import plot_chevron

def Chevron_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,detuning:float,freq_pts:int=50,XY_amp:float=0.5, XY_duration:float=20e-9, IF:int=250e6,n_avg:int=300,points:int=100,run:bool=True,XY_theta:str='X_theta',Rabi_type:str='PowerRabi',q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0],specific_data_folder:str=''):
    
    sche_func= Rabi_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    
    LO= qubit_info.clock_freqs.f01()+IF
    
    qubit_info.measure.pulse_duration(2e-6)
    qubit_info.measure.integration_time(1e-6)
    qubit_info.reset.duration(250e-6)
    print("Integration time ",qubit_info.measure.integration_time()*1e6, "µs")
    print("Reset time ", qubit_info.reset.duration()*1e6, "µs")

    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    eyeson_print(f"XYF = {round(qubit_info.clock_freqs.f01()*1e-9,3)} GHz")
    
    xyf = ManualParameter(name="xyf", unit="Hz", label="XY Frequency")
    xyf.batched = False
    f01_samples = linspace(qubit_info.clock_freqs.f01()-detuning,qubit_info.clock_freqs.f01()+detuning,freq_pts)
    
    if Rabi_type.lower() in ['timerabi', 'tr']:
        osci_type = "TimeRabi"
        Para_XY_amp = XY_amp
        
        Sweep_para=Para_XY_Du = ManualParameter(name="XY_Duration", unit="s", label="Time")
        str_Rabi= 'XY_duration'
        Sweep_para.batched = True
        if XY_duration < 200e-9:
            XY_duration = 200e-9
            resolu = 4e-9
        else:
            resolu = 4e-9

        gap = (XY_duration)*1e9 // points + (((XY_duration)*1e9 // points) %(resolu*1e9)) # multiple by 4 ns
        samples = arange(0,XY_duration,4e-9)
        eyeson_print(f"data points = {samples.shape[0]} with gap = {gap} ns")
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
        chevron=True,
        New_fxy=xyf
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
        meas_ctrl.settables([Sweep_para,xyf])
        meas_ctrl.setpoints_grid((samples,f01_samples))
    
       
        rabi_ds = meas_ctrl.run("Rabi_Chevron")
        # Save the raw data into netCDF
        
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rabi_ds,qb=q,exp_type="Chevron",specific_dataFolder=specific_data_folder,get_data_loc=True)

        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = -2
        sweep_para= array([samples[0],samples[-1]])
        sched_kwargs[str_Rabi]= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
       

        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        nc_path = ''

    return nc_path
    

def chevron_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,detuning:float,freq_pts:int=50,XYamp_max:float=0.5,XYdura_max:float=20e-9,which_rabi:str='power',run:bool=True,pts:int=100,avg_times:int=500,data_folder:str=''):
    if which_rabi.lower() in ['p','power']:
        exp_type = 'powerRabi'
    elif which_rabi.lower() in ['t','time']:
        exp_type = 'timeRabi'
    else:
        raise ValueError (f"Can't recognized rabi type with the given = {which_rabi}")
    
    print(f"{specific_qubits} are under the measurement ...")

    if run:
        
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        result_nc_path = Chevron_spec(QD_agent,meas_ctrl,detuning=detuning,freq_pts=freq_pts,Rabi_type=exp_type,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=True,XY_amp=XYamp_max,XY_duration=XYdura_max,points=pts,n_avg=avg_times,specific_data_folder=data_folder)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
    else:
        result_nc_path = Chevron_spec(QD_agent,meas_ctrl,detuning=detuning,freq_pts=freq_pts,Rabi_type=exp_type,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=False,XY_amp=XYamp_max,XY_duration=XYdura_max)
    
    return result_nc_path
  

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = ['q0']
    couplers = []


    """ Optional paras """
    pi_duration:float = 40e-9
    pi_amp_max:float = 0.9
    rabi_type:str = 'power'  # 'time' or 'power'
    data_pts = 100
    avg_n:int = 500
    xy_atte:int = 10
    span_detuning:float = 30e6
    freq_pts:int = 60
    adj_freq:float = 0e6
    
    """ Operations """
    
    for q_idx, qubit in enumerate(ro_elements):

        """ Preparations """
        start_time = time.time()
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
        QD_agent.Notewriter.save_DigiAtte_For(xy_atte,qubit,'xy')
        QD_agent.quantum_device.get_element(qubit).clock_freqs.f01(QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()+adj_freq)
        
        """Running """
        ncs = {}
        Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        ncs[qubit] = chevron_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,detuning=span_detuning,freq_pts=freq_pts,run=execution,XYdura_max=pi_duration,XYamp_max=pi_amp_max,which_rabi=rabi_type,avg_times=avg_n,pts=data_pts)
        cluster.reset()
    
        """ Storing """
        if q_idx == len(ro_elements)-1:
            for q in ncs:
                plot_chevron(QD_agent,q,ncs[q])


        """ Close """
        print('Rabi osci done!')
        shut_down(cluster,Fctrl,Cctrl)
        end_time = time.time()
        slightly_print(f"time cost: {round(end_time-start_time,1)} secs")