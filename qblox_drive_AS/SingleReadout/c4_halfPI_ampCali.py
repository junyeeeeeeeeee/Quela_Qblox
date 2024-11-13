"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from qblox_drive_AS.support.UserFriend import eyeson_print, slightly_print
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import linspace, array, arange, NaN, ndarray
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import pi_half_cali_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, Rabi_fit_analysis, Fit_analysis_plot
import matplotlib.pyplot as plt

def half_pi_amp_cali(QD_agent:QDmanager,meas_ctrl:MeasurementControl, half_pi_quadruple_num:int=3, amp_coef_span:float=0.4, IF:int=250e6,n_avg:int=300,points:int=100,run:bool=True,q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0],specific_data_folder:str=''):
    analysis_result = {}
    sche_func= pi_half_cali_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    LO= qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    Sweep_para = ManualParameter(name="halfXY_Amp_coef")
    Sweep_para.batched = True
    
    samples = linspace(1-amp_coef_span,1+amp_coef_span,points) 
    exp_kwargs= dict(sweep_amp=["halfXY_Amp_coef",'start '+'%E' %samples[0],'end '+'%E' %samples[-1]])
    
    sched_kwargs = dict(
        q=q,
        pi_half_coefs=Sweep_para,
        half_pi_quadruple_num=half_pi_quadruple_num,
        pi_amp={str(q):qubit_info.rxy.amp180()},
        XY_duration=qubit_info.rxy.duration(),
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
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
    
       
        cali_ds = meas_ctrl.run("half-Pi amp calibration")
        # Save the raw data into netCDF
        
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=cali_ds,qb=q,exp_type="xyl05cali",label=f"{half_pi_quadruple_num}HalfPi",specific_dataFolder=specific_data_folder)
        I,Q= dataset_to_array(dataset=cali_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        
        analysis_result[q]= data
        
        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        sweep_para= array([samples[0],samples[-1]])
        sched_kwargs["pi_half_coefs"]= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
       

        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))

    return analysis_result

def plot_cali_results(data:dict, samples:ndarray):
    for pi_num in data:
        plt.plot(samples,data[pi_num],label=f"{pi_num}*4 pi/2 pulses")
    plt.xlabel("pi/2 coefficient")
    plt.ylabel("Contrast (mV)")
    plt.legend()
    plt.title("Pi/2-pulse amplitude calibration")
    plt.show()
    

def half_pi_amp_calibrator(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,XYamp_coef_span:float=0.5,half_pi_quadruple_num:int=3,run:bool=True,pts:int=100,avg_times:int=500,data_folder:str='',IF:float=250e6):

    print(f"{specific_qubits} are under the measurement ...")
    if run:
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        
        Rabi_results = half_pi_amp_cali(QD_agent,meas_ctrl,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=True,half_pi_quadruple_num=half_pi_quadruple_num,amp_coef_span=XYamp_coef_span,points=pts,n_avg=avg_times,specific_data_folder=data_folder,IF=IF)
        Fctrl[specific_qubits](0.0)
        
        cluster.reset()
        
        return Rabi_results[specific_qubits]
    else:
        Rabi_results = half_pi_amp_cali(QD_agent,meas_ctrl,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=False,half_pi_quadruple_num=half_pi_quadruple_num,amp_coef_span=XYamp_coef_span,points=pts,n_avg=avg_times,specific_data_folder=data_folder)
        return 0
  

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = ['q0']
    couplers = []


    """ Optional paras """
    half_pi_quadruple_num:list = [7,9]
    pi_amp_coef_span:float = 0.1
    avg_n:int = 1000
    data_pts:int = 80
    xy_IF = 250e6
    
    """ Operations """
    
    for qubit in ro_elements:
        data = {}
        for idx, pi_num in enumerate(half_pi_quadruple_num):
            slightly_print(f"Driving with {pi_num}*4 pi/2 pulses")
            """ Preparations """
            start_time = time.time()
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
            chip_info = cds.Chip_file(QD_agent=QD_agent)
            
    
            """Running """
            rabi_results = {}
            Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            data[str(pi_num)] = half_pi_amp_calibrator(QD_agent,cluster,meas_ctrl,Fctrl,qubit,run=execution,avg_times=avg_n,pts=data_pts,XYamp_coef_span=pi_amp_coef_span,half_pi_quadruple_num=pi_num,IF=xy_IF)
            cluster.reset()
        
            """ Storing """
            if idx == len(half_pi_quadruple_num)-1:
                plot_cali_results(data,linspace(1-pi_amp_coef_span,1+pi_amp_coef_span,data_pts))
            """ Close """
            shut_down(cluster,Fctrl)
            end_time = time.time()
            slightly_print(f"time cost: {round(end_time-start_time,1)} secs")

            
        
