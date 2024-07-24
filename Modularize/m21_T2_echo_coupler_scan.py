import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from qcodes.instrument import find_or_create_instrument
from quantify_scheduler.device_under_test.transmon_element import BasicTransmonElement
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from numpy import std, arange, array, linspace, pi, average, mean, sign, arctan
from quantify_core.measurement.control import MeasurementControl
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Ramsey_echo_coupler_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot, Fit_T2_cali_analysis_plot


def Ramsey_coupler(QD_agent:QDmanager,
           meas_ctrl:MeasurementControl,
           target_coupler:str,
           flux_center:float,
           flux_span:float,
           flux_points:int,
           freeduration:float,
           z_delay:float,
           arti_detune:int=0,
           IF:int=150e6,
           n_avg:int=1000,
           points:int=101,
           run:bool=True,
           q='q0',
           excite_q='q1',
           zz=0,
           ref_IQ:list=[0,0],
           Experi_info:dict={},
           exp_idx:int=0,
           data_folder:str=''
           ):
    
    T2_us = {}
    T2_us[q] = []
    analysis_result = {}
    analysis_result[q]= []
    Real_detune= {}
    
    qubit = QD_agent.quantum_device.get_element(q)
    excite_qubit = QD_agent.quantum_device.get_element(excite_q)
    excite_amp = {str(excite_q):excite_qubit.rxy.amp180()}

    # Manually change f01
    # qubit.clock_freqs.f01(qubit.clock_freqs.f01()+1.343e6)
    
    New_fxy= qubit.clock_freqs.f01()+arti_detune
    
    LO= New_fxy+IF
    excite_LO = excite_qubit.clock_freqs.f01()+IF
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    set_LO_frequency(QD_agent.quantum_device,q=excite_q,module_type='drive',LO_frequency=excite_LO)
    
    Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
    Para_free_Du.batched = True
    gap = (freeduration)*1e9 // points + (((freeduration)*1e9 // points) %4) # multiple by 4 ns
    samples = arange(8e-9,freeduration,gap*1e-9)

    flux_bias = ManualParameter(name="Z", unit="V", label="Z bias")
    flux_bias.batched = False
    flux_samples = linspace(flux_center-flux_span,flux_center+flux_span,flux_points)
    
    sche_func= Ramsey_echo_coupler_sche
    sched_kwargs = dict(
        q=q,
        excite_q=excite_q,
        coupler=target_coupler,
        flux_amp=flux_bias,
        z_delay=z_delay,
        pi_amp={str(q):qubit.rxy.amp180()},
        excite_pi_amp=excite_amp,
        New_fxy=New_fxy,
        freeduration=Para_free_Du,
        R_amp={str(q):qubit.measure.pulse_amp()},
        R_duration={str(q):qubit.measure.pulse_duration()},
        R_integration={str(q):qubit.measure.integration_time()},
        R_inte_delay=qubit.measure.acq_delay(),
        zz=zz,
        pi_dura=qubit.rxy.duration(),
        )
    exp_kwargs= dict(sweep_freeDu=['start '+'%E' %samples[0],'end '+'%E' %samples[-1]],
                     Flux=['start '+'%E' %flux_samples[0],'end '+'%E' %flux_samples[-1]],
                     f_xy='%E' %sched_kwargs['New_fxy'],
                     )
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=False,
            batched=True,
        )
        
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([Para_free_Du,flux_bias])
        meas_ctrl.setpoints_grid((samples, flux_samples))
        
        
        ramsey_ds = meas_ctrl.run('Ramsey')

        # Save the raw data into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=ramsey_ds,label=exp_idx,qb=q,exp_type='t2c',specific_dataFolder=data_folder)
        Basic2DAnalysis(tuid=ramsey_ds.attrs["tuid"], dataset=ramsey_ds).run()
        for j in range(flux_points):
            I,Q= dataset_to_array(dataset=ramsey_ds,dims=1)
            data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[1])
            # #phase=arctan((I-ref_IQ[0])/(Q-ref_IQ[1]))
            
            # ramsey_ds.y0.data = data
            # #ramsey_ds.y1.data = phase
            # Basic2DAnalysis(tuid=ramsey_ds.attrs["tuid"], dataset=ramsey_ds).run()

            analysis_result[q].append(data)
            T2_us[q].append(0)

        show_args(exp_kwargs, title="Ramsey_kwargs: Meas.qubit="+q)
        T2_us[q] = 0
        if Experi_info != {}:
            show_args(Experi_info(q))
    else:
        n_s = 2
        sweep_para= array(samples[:n_s])
        sweep_para2= array(flux_samples[:2])
        sched_kwargs['freeduration']= sweep_para.reshape(sweep_para.shape or (1,))
        sched_kwargs['flux_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        

        show_args(exp_kwargs, title="Ramsey_echo_coupler_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    return analysis_result, T2_us, Real_detune


def ramsey_coupler_executor(QD_agent:QDmanager,
                            cluster:Cluster,
                            meas_ctrl:MeasurementControl,
                            target_coupler:str,
                            flux_center:float,
                            flux_span:float,
                            flux_pts:int,
                            z_delay:float,
                            Fctrl:dict,
                            specific_qubits:str,
                            excite_q:str='q1',
                            zz:int=0,
                            artificial_detune:float=0e6,
                            freeDura:float=30e-6,
                            ith:int=1,
                            run:bool=True,
                            specific_folder:str='',
                            pts:int=100,
                            avg_n:int=800
                            ):
    if run:
        qubit_info = QD_agent.quantum_device.get_element(specific_qubits)
        ori_reset = qubit_info.reset.duration()
        qubit_info.reset.duration(qubit_info.reset.duration()+freeDura)
        
        start_time = time.time()
        slightly_print(f"The {ith}-th T2:")
        Fctrl[specific_qubits](float(QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)))
        Ramsey_results, T2_us, average_actual_detune = Ramsey_coupler(QD_agent,
                                                                      meas_ctrl,
                                                                      target_coupler=target_coupler,
                                                                      flux_center=flux_center,
                                                                      flux_span=flux_span,
                                                                      flux_points=flux_pts,
                                                                      z_delay=z_delay,
                                                                      arti_detune=artificial_detune,
                                                                      freeduration=freeDura,
                                                                      n_avg=avg_n,
                                                                      q=specific_qubits,
                                                                      excite_q=excite_q,
                                                                      zz=zz,
                                                                      ref_IQ=QD_agent.refIQ[specific_qubits],
                                                                      points=pts,
                                                                      run=True,
                                                                      exp_idx=ith,
                                                                      data_folder=specific_folder)
        Fctrl[specific_qubits](0.0)
        cluster.reset()
        this_t2_us = T2_us[specific_qubits]
        end_time = time.time()
        slightly_print(f"time cost: {round(end_time-start_time,1)} secs")
        
        qubit_info.reset.duration(ori_reset)

    else:
        Ramsey_results, _, average_actual_detune = Ramsey_coupler(QD_agent,
                                                                  meas_ctrl,
                                                                  target_coupler=target_coupler,
                                                                  flux_center=flux_center,
                                                                  flux_span=flux_span,
                                                                  flux_points=flux_pts,
                                                                  z_delay=z_delay,
                                                                  arti_detune=artificial_detune,
                                                                  freeduration=freeDura,
                                                                  n_avg=1000,
                                                                  q=specific_qubits,
                                                                  excite_q=excite_q,
                                                                  zz=zz,
                                                                  ref_IQ=QD_agent.refIQ[specific_qubits],
                                                                  points=100,
                                                                  run=False)
        this_t2_us = 0

    return Ramsey_results, this_t2_us, average_actual_detune




if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 1
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = {
        "q4":{"detune":0e6,"evoT":16e-6,"histo_counts":1},
        # "q1":{"detune":0.8e6,"evoT":50e-6,"histo_counts":1},
    }
    target_coupler = 'q4'
    target_bias_span = 0.2
    target_center = 0
    target_pts = 40
    couplers = ['c0','c1','c2','c3']
    pts = 100
    zz_on = 0
    z_delay = 100e-9
    excite_qubit = 'q3'


    """ Iteration """
    for qubit in ro_elements:
        t2_us_rec = []
        for ith_histo in range(ro_elements[qubit]["histo_counts"]):
            """ Preparations """
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l',dr_loc=DRandIP["dr"])
            cplr=find_or_create_instrument(BasicTransmonElement, recreate=True, name=f"q{target_coupler}")
            QD_agent.quantum_device.add_element(cplr)
            chip_info = cds.Chip_file(QD_agent=QD_agent)


            """ Running """
            cp_ctrl = QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i')
            Cctrl = coupler_zctrl(DRandIP["dr"],cluster,cp_ctrl)
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            
            slightly_print(f"Ramsey with detuning = {round(ro_elements[qubit]['detune']*1e-6,2)} MHz")
            ramsey_results, this_t2_us, average_actual_detune = ramsey_coupler_executor(QD_agent,
                                                                                        cluster,
                                                                                        meas_ctrl,
                                                                                        target_coupler=target_coupler,
                                                                                        flux_center=target_center,
                                                                                        flux_span=target_bias_span,
                                                                                        flux_pts=target_pts,
                                                                                        z_delay=z_delay,
                                                                                        Fctrl=Fctrl,
                                                                                        specific_qubits=qubit,
                                                                                        excite_q=excite_qubit,
                                                                                        zz=zz_on,
                                                                                        artificial_detune=ro_elements[qubit]["detune"],
                                                                                        freeDura=ro_elements[qubit]["evoT"],
                                                                                        ith=ith_histo,
                                                                                        pts=pts,
                                                                                        run=execution)
            highlight_print(f"{qubit} XYF = {round(QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()*1e-9,5)} GHz")
            if this_t2_us != 0:
                t2_us_rec.append(this_t2_us)
            

            """ Close """
            print('T2 done!')
            shut_down(cluster,Fctrl,Cctrl)
            
        
        """ Storing """
        if execution:
            # mean_T2_us = round(mean(array(t2_us_rec)),2)
            # std_T2_us  = round(std(array(t2_us_rec)),2)
            # highlight_print(f"{qubit}: mean T2 = {mean_T2_us} 土 {std_T2_us} µs")
            # QD_agent.QD_keeper()# Manual keep
            if ro_elements[qubit]["histo_counts"] == 1:
                mean_T2_us = 0
                std_T2_us = 0
                # Fit_analysis_plot(ramsey_results[qubit],P_rescale=False,Dis=None)
            else:
                Data_manager().save_histo_pic(QD_agent,{str(qubit):t2_us_rec},qubit,mode="t2")
                if ro_elements[qubit]["histo_counts"] >= 50:
                    QD_agent.Notewriter.save_T2_for(mean_T2_us,qubit)
                    # QD_agent.QD_keeper()
                    if chip_info_restore:
                        chip_info.update_T2(qb=qubit, T2=f'{mean_T2_us} +- {std_T2_us}')
            
            
        


    