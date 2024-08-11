import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from numpy import mean, array, arange, std, linspace
from utils.tutorial_utils import show_args
from Modularize.support.UserFriend import *
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Zgate_T1_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T1_fit_analysis, Fit_analysis_plot, Z_bias_error_bar_plot


def Zgate_T1(QD_agent:QDmanager,meas_ctrl:MeasurementControl,freeduration:float,Z_amp_min:float,Z_amp_max:float,n_avg:int=500,freeDu_points:int=101,Z_points:int=201,IF:float=150e6,run:bool=True,q='q1',Experi_info:dict={},exp_idx:int=0,no_pi_pulse:bool=False,data_folder:str=''):
    times = 1
    qubit_info = QD_agent.quantum_device.get_element(q)
    if not no_pi_pulse:
        pi_amp = qubit_info.rxy.amp180()
    else:
        pi_amp = 0
    analysis_result= {}
    T1_us = {}
    analysis_result[q]= []
    T1_us[q] = []
    LO= qubit_info.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    gap = (freeduration*1e9 // freeDu_points) + ((freeduration*1e9 // freeDu_points)%4)
    freeDu_samples = arange(4e-9,freeduration,gap*1e-9)
    
    Z_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    Z_bias.batched = False
    Z_samples = linspace(Z_amp_min,Z_amp_max,Z_points)
    sche_func= Zgate_T1_sche
    sched_kwargs = dict(
        q=q,
        pi_amp={str(q):pi_amp},
        pi_dura={str(q):qubit_info.rxy.duration()},
        freeduration=Para_free_Du,
        Z_amp=Z_bias,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %freeDu_samples[0],'end '+'%E' %freeDu_samples[-1]],
                     Z_amp=['start '+'%E' %Z_samples[0],'end '+'%E' %Z_samples[-1]],
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
        meas_ctrl.settables([Para_free_Du,Z_bias])
        meas_ctrl.setpoints_grid((freeDu_samples,Z_samples))

        
        for i in range(times):
            T1_ds = meas_ctrl.run('Zgate_T1')
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=T1_ds,label=exp_idx,qb=q,exp_type='zT1',specific_dataFolder=data_folder)
            for j in range(Z_points):
                I,Q= dataset_to_array(dataset=T1_ds,dims=2)
                data= IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[q][0],ref_Q=QD_agent.refIQ[q][1]).transpose()
                if not no_pi_pulse:
                    data_fit= T1_fit_analysis(data=data[j],freeDu=freeDu_samples,T1_guess=8e-6)
                    Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)
                    analysis_result[q].append(data_fit)
                    T1_us[q].append(data_fit.attrs['T1_fit']*1e6)
                else:
                    analysis_result[q].append(data)
                    T1_us[q].append(0)
        T1_us['plot_parameters']=[times,Z_samples]
        show_args(exp_kwargs, title="Zgate_T1_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para1= array(freeDu_samples[:n_s])
        sweep_para2= array(Z_samples[:2])
        sched_kwargs['freeduration']= sweep_para1.reshape(sweep_para1.shape or (1,))
        sched_kwargs['Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Zgate_T1_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    return analysis_result, T1_us

def zgateT1_executor(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,freeDura:float=30e-6,z_span_period_factor:int=2,run:bool=True,specific_folder:str='',tpts:int=101,zpts=201,ith:int=0,pi_pulse:bool=True):
    flux_guard_Volt = 0.4

    if run:
        offset_flux = float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits))
        qubit_info = QD_agent.quantum_device.get_element(specific_qubits)
        ori_reset = qubit_info.reset.duration()
        qubit_info.reset.duration(qubit_info.reset.duration()+freeDura)

        span_period = float(QD_agent.Fluxmanager.get_PeriodFor(specific_qubits))/z_span_period_factor
        z_from = -span_period/2
        z_end = +span_period/2

        if abs(offset_flux+z_from)>flux_guard_Volt or abs(offset_flux+z_end)>flux_guard_Volt :
            raise ValueError(f"Z span from {round(z_from,2)} to {round(z_end,2)} may not be an appropriate voltage !")

        every_start = time.time()
        slightly_print(f"The {ith}-th T1:")
        Fctrl[specific_qubits](offset_flux)
        T1_results, this_t1_us = Zgate_T1(QD_agent,meas_ctrl,q=specific_qubits,freeduration=freeDura,run=True,exp_idx=ith,data_folder=specific_folder,Z_points=zpts,freeDu_points=tpts,no_pi_pulse = not pi_pulse,Z_amp_max=z_end,Z_amp_min=z_from)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
        
        slightly_print(f"T1: {this_t1_us} µs")
        every_end = time.time()
        slightly_print(f"time cost: {round(every_end-every_start,1)} secs")
        
        qubit_info.reset.duration(ori_reset)
        
    else:
        T1_results, this_t1_us = Zgate_T1(QD_agent,meas_ctrl,q=specific_qubits,freeduration=freeDura,run=False,exp_idx=ith,data_folder=specific_folder,Z_points=zpts,freeDu_points=tpts,no_pi_pulse = not pi_pulse,Z_amp_max=0.01,Z_amp_min=-0.01)
        

    
    return T1_results, this_t1_us

if __name__ == "__main__":
    

    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 0
    prepare_excited:bool = 1
    DRandIP = {"dr":"dr3","last_ip":"13"}
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        "q4":{"evoT":40e-6,"histo_counts":2}
    }
    couplers = ['c0','c1','c2','c3']

    """ Optional paras """
    flux_span_period_factor = 2 # flux axis = [sweet-period/flux_span_period_factor, sweet+period/flux_span_period_factor]
    flux_data_points = 10
    evotime_data_points = 101
    
    

    """ Iterations """
    t1_us_rec = {}
    for qubit in ro_elements:

        t1_us_rec[qubit] = []
        for ith_histo in range(ro_elements[qubit]["histo_counts"]):
            """ Preparations """
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            chip_info = cds.Chip_file(QD_agent=QD_agent)
            
            """ Running """
            Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            T1_results, these_t1_us = zgateT1_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,freeDura=ro_elements[qubit]["evoT"],z_span_period_factor=flux_span_period_factor,tpts=evotime_data_points,zpts=flux_data_points,run=execution,ith=ith_histo,pi_pulse=prepare_excited)
            t1_us_rec[qubit]+= these_t1_us[qubit]

            """ Close """
            print('Zgate T1 done!')
            shut_down(cluster,Fctrl,Cctrl)
        
        t1_us_rec['plot_parameters'] = [ro_elements[qubit]["histo_counts"],these_t1_us['plot_parameters'][1]]
        """ Storing """
        if execution:
            if ro_elements[qubit]["histo_counts"] == 1:
                Fit_analysis_plot(T1_results[qubit][0],P_rescale=False,Dis=None)
            Z_bias_error_bar_plot(qubit,t1_us_rec,title=r"$T_{1}\  (\mu$s)")

                

            
       

    
        