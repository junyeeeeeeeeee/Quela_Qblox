"""This program includes PowerRabi and TimeRabi. When it's PoweRabi, default ctrl pulse duration is 20ns."""
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from Modularize.support.UserFriend import *
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from numpy import linspace, array, arange, NaN, ndarray, sin, mean, cos
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr, meas_raw_dir
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import PI_amp_cali_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def sin_wave(x,A,k,phi,B):
    return A*sin(k*x+phi)+B

def pi_amp_cali(QD_agent:QDmanager,meas_ctrl:MeasurementControl, pi_pair_num:int=3, amp_coef_span:float=0.4, IF:int=250e6,n_avg:int=300,points:int=100,run:bool=True,q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0],specific_data_folder:str=''):
    analysis_result = {}
    sche_func= PI_amp_cali_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    # qubit_info.measure.acq_delay(0)
    # qubit_info.reset.duration(250e-6)
    LO= qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    Sweep_para = ManualParameter(name="XY_Amp_coef")
    Sweep_para.batched = True
    
    samples = linspace(1-amp_coef_span,1+amp_coef_span,points) 
    exp_kwargs= dict(sweep_amp=["XY_Amp_coef",'start '+'%E' %samples[0],'end '+'%E' %samples[-1]])
    
    sched_kwargs = dict(
        q=q,
        pi_amp_coefs=Sweep_para,
        pi_pair_num=pi_pair_num,
        XY_amp={str(q):qubit_info.rxy.amp180()},
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
    
       
        cali_ds = meas_ctrl.run("Pi amp calibration")
        # Save the raw data into netCDF
        
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=cali_ds,qb=q,exp_type="xylcali",label=f'{pi_pair_num}Pi',specific_dataFolder=specific_data_folder)
        I,Q= dataset_to_array(dataset=cali_ds,dims=1)
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        
        analysis_result[q]= data
        
        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        sweep_para= array([samples[0],samples[-1]])
        sched_kwargs["pi_amp_coefs"]= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
       

        show_args(exp_kwargs, title="Rabi_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))

    return analysis_result

def cali_fit(x:ndarray,y:ndarray):
    
    init_guess = [max(y)-min(y),1/(max(x)-min(x)),max(x),mean(y)]
    bound = [[0.8*init_guess[0],0.3*init_guess[1],min(x),min(y)],[1.2*init_guess[0],3*init_guess[1],1.5*init_guess[2],max(y)]]
    p, e = curve_fit(sin_wave,x,y,p0=init_guess,bounds=bound)
    return x, p, sin_wave(x,*p)

def plot_cali_results(data:dict, samples:ndarray):
    colors = [["orange","red"],["cyan","blue"]]
    ans = {}
    for idx, pi_num in enumerate(list(data.keys())):
        x, p, y = cali_fit(samples,data[pi_num])
        plt.plot(samples,1000*data[pi_num],c=colors[idx][0],label=f'{pi_num}*2 pi-pulses')
        # plt.plot(x,y,label=f"{pi_num}*2 pi-pulses",c=colors[idx][1])
        # plt.scatter(1,sin_wave(array([1]),*p),marker='*')
        # plt.scatter(p[2],sin_wave(array([p[2]]),*p),marker="x")
        # ans[f"pi_num={pi_num}"] = round(sin_wave(array([1]),*p)[0],5)
    plt.xlabel("Amplitude coefficient")
    plt.ylabel("Contrast (mV)")
    plt.legend()
    plt.title("Pi-pulse amplitude calibration")
    plt.show()
    # print(ans)
    # print(f"Avg coef = {round(mean(array(list(ans.values()))),4)} ")
    return ans
    

def pi_amp_calibrator(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,XYamp_coef_span:float=0.5,pi_pair_num:int=3,run:bool=True,pts:int=100,avg_times:int=500,data_folder:str='',IF:float=250e6):

    print(f"{specific_qubits} are under the measurement ...")
    if run:
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        Rabi_results = pi_amp_cali(QD_agent,meas_ctrl,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=True,pi_pair_num=pi_pair_num,amp_coef_span=XYamp_coef_span,points=pts,n_avg=avg_times,specific_data_folder=data_folder,IF=IF)
        Fctrl[specific_qubits](0.0)
    
        cluster.reset()
        
        return Rabi_results[specific_qubits]
    else:
        Rabi_results = pi_amp_cali(QD_agent,meas_ctrl,q=specific_qubits,ref_IQ=QD_agent.refIQ[specific_qubits],run=False,pi_pair_num=pi_pair_num,amp_coef_span=XYamp_coef_span,points=pts,n_avg=avg_times,specific_data_folder=data_folder)
        return 0
  

if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1 
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = ['q0']
    couplers = []


    """ Optional paras """
    pi_pair_num:list = [7,9]
    pi_amp_coef_span:float = 0.1
    avg_n:int = 800
    data_pts:int = 80
    xy_IF = 250e6
    
    """ Operations """
    
    for qubit in ro_elements:
        data = {}
        for idx, pi_num in enumerate(pi_pair_num):
            slightly_print(f"Driving with {pi_num}*2 pi pulses")
            """ Preparations """
            start_time = time.time()
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            chip_info = cds.Chip_file(QD_agent=QD_agent)
            
    
            """Running """
            rabi_results = {}
            Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            data[str(pi_num)] = pi_amp_calibrator(QD_agent,cluster,meas_ctrl,Fctrl,qubit,run=execution,avg_times=avg_n,pts=data_pts,XYamp_coef_span=pi_amp_coef_span,pi_pair_num=pi_num,IF=xy_IF)
            cluster.reset()
        
            """ Storing """
            if idx == len(pi_pair_num)-1:
                ans = plot_cali_results(data,linspace(1-pi_amp_coef_span,1+pi_amp_coef_span,data_pts))
                coef = mark_input("Input the modified coef, n to cancel: ")
                if str(coef).lower() not in ['n', 'no', '']:
                    qubit_info = QD_agent.quantum_device.get_element(qubit)
                    qubit_info.rxy.amp180(qubit_info.rxy.amp180()*float(coef))
                    QD_agent.QD_keeper()

            """ Close """
            shut_down(cluster,Fctrl,Cctrl)
            end_time = time.time()
            slightly_print(f"time cost: {round(end_time-start_time,1)} secs")

            
        
