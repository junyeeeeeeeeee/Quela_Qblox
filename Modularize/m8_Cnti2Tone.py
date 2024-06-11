
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, linspace, NaN
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support.QuFluxFit import calc_Gcoef_inFbFqFd, calc_g
from Modularize.support import init_meas, shut_down,  advise_where_fq, init_system_atte, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Two_tone_sche, set_LO_frequency, pulse_preview, IQ_data_dis, QS_fit_analysis, dataset_to_array, twotone_comp_plot

def Two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,IF:float=100e6,f01_guess:int=0,xyf_span_Hz:int=400e6,xyamp:float=0.02,n_avg:int=500,points:int=200,run:bool=True,q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0]):
    
    sche_func = Two_tone_sche   
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    # qubit_info.reset.duration(0)
    if f01_guess != 0:
        f01_high = f01_guess+IF
    else:
        f01_high = qubit_info.clock_freqs.f01()+IF
    # if xyamp == 0:
    #     xyamp = qubit_info.rxy.amp180(XYL)
    # Avoid warning
    qubit_info.clock_freqs.f01(NaN)

    f01_samples = linspace(f01_high-xyf_span_Hz,f01_high,points)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=f01_high)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True

    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        spec_amp=xyamp,
        spec_Du=48e-6,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        ref_pt="start"
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %f01_samples[0],'end '+'%E' %f01_samples[-1]],
                     spec_amp='%E' %spec_sched_kwargs['spec_amp'],
                     spec_Du='%E' %spec_sched_kwargs['spec_Du'])
    highlight_print(f"Now spec_amp = {xyamp}")
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(f01_samples)
        
        qs_ds = meas_ctrl.run("Two-tone")
        # Save the raw data into netCDF
        Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='2tone')
        I,Q= dataset_to_array(dataset=qs_ds,dims=1)
        
        data= IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1]) 
        analysis_result[q] = QS_fit_analysis(data,f=f01_samples)
        
        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        n_s = 2
        sweep_para= array(f01_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="Two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
        analysis_result[q] = {}

    return analysis_result


def update_2toneResults_for(QD_agent:QDmanager,qb:str,QS_results:dict,XYL:float):
    qubit = QD_agent.quantum_device.get_element(qb)
    Revised_f01= QS_results[qb].attrs['f01_fit']
    fb = float(QD_agent.Notewriter.get_bareFreqFor(target_q=qb))*1e-6
    fd = QD_agent.quantum_device.get_element(qb).clock_freqs.readout()*1e-6
    A = calc_Gcoef_inFbFqFd(fb,Revised_f01*1e-6,fd)
    sweet_g = calc_g(fb,Revised_f01*1e-6,A)
    qubit.clock_freqs.f01(Revised_f01)
    QD_agent.Notewriter.save_2tone_piamp_for(qb,XYL)
    QD_agent.Notewriter.save_CoefInG_for(target_q=qb,A=A)
    QD_agent.Notewriter.save_sweetG_for(target_q=qb,g_Hz=sweet_g*1e6)



def conti2tone_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,cluster:Cluster,specific_qubits:str,xyf_guess:float=0,guess_g:float=48e6,xyAmp_guess:list=[],xyf_span:float=500e6,xy_if:float=100e6,run:bool=True,V_away_from:float=0):
    
    if run:
        advised_fq = advise_where_fq(QD_agent,specific_qubits,guess_g)
        print(f"fq advice for {specific_qubits} @ {round(advised_fq*1e-9,4)} GHz")
        if xyf_guess != 0:
            guess_fq = [xyf_guess]
        else:
            guess_fq = [advised_fq-500e6, advised_fq, advised_fq+500e6]

        if len(xyAmp_guess) == 0 or (len(xyAmp_guess) == 1 and xyAmp_guess[0] == 0):
            xyAmp_guess = [0, 0.03, 0.07, 0.1]
        else:
            xyAmp_guess = xyAmp_guess
        
        for XYF in guess_fq:
            ori_data = []
            for XYL in xyAmp_guess:
                want_bias = QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)-V_away_from
                if V_away_from != 0:
                    rof = QD_agent.Fluxmanager.sin_for_cav(specific_qubits,array([want_bias]))[0]
                    QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(rof)
                    QD_agent.Fluxmanager.save_tuneawayBias_for('manual',specific_qubits,want_bias)
                Fctrl[specific_qubits](want_bias) 
                QS_results = Two_tone_spec(QD_agent,meas_ctrl,xyamp=XYL,IF=xy_if,f01_guess=XYF,q=specific_qubits,xyf_span_Hz=xyf_span,points=50,n_avg=500,run=True,ref_IQ=QD_agent.refIQ[specific_qubits]) # 
                Fctrl[specific_qubits](0.0)
                cluster.reset() # *** important
                if XYL != 0:
                    twotone_comp_plot(QS_results[specific_qubits], ori_data, True)
                else:
                    twotone_comp_plot(QS_results[specific_qubits], ori_data, False)
                    ori_data = QS_results[specific_qubits].data_vars['data']
            
        return QS_results[specific_qubits]
                
    else:
        qu = specific_qubits
        QS_results = Two_tone_spec(QD_agent,meas_ctrl,xyamp=0.1,IF=xy_if,f01_guess=4e9,q=qu,xyf_span_Hz=xyf_span,points=50,n_avg=500,run=False,ref_IQ=QD_agent.refIQ[qu])

        return 0
   
    

"""
NOTE: If you find a XYF in a tuneaway z-bias which means the `tune_bias` != 0, please go Modularize/support/meas_switch.py to save it.
"""

if __name__ == "__main__":

    """ Fill in """
    execution = True
    update = 1
    #
    DRandIP = {"dr":"dr3","last_ip":"13"}
    #
    ro_elements = {
        "q1":{"xyf_guess":[4.02e9],"xyl_guess":[0.03],"g_guess":0e6, "tune_bias":0} # g you can try a single value in  [42e6, 54e6, 62e6], higher g makes fq lower.
    }                                                                            # tune_bias is the voltage away from sweet spot. If it was given, here will calculate a ROF according to that z-bias and store it in Notebook.
    couplers = ["c0", 'c1']
    # 1 = Store
    # 0 = not store
    chip_info_restore = 1

    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    chip_info = cds.Chip_file(QD_agent=QD_agent)

    """ Running """
    tt_results = {}
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        QD_agent.Notewriter.save_DigiAtte_For(0,qubit,'xy')
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
        tune_bias = ro_elements[qubit]["tune_bias"]
        for xyf in ro_elements[qubit]["xyf_guess"]:
            g = 48e6 if ro_elements[qubit]["g_guess"] == 0 else ro_elements[qubit]["g_guess"]

            tt_results[qubit] = conti2tone_executor(QD_agent,meas_ctrl,cluster,specific_qubits=qubit,xyf_guess=xyf,xyAmp_guess=ro_elements[qubit]["xyl_guess"],run=execution,guess_g=g,xy_if=100e6,xyf_span=500e6,V_away_from=tune_bias)

            if execution and ro_elements[qubit]["xyl_guess"][0] != 0:
                print(f'update xyl={ro_elements[qubit]["xyl_guess"][0]}')
                update_2toneResults_for(QD_agent,qubit,tt_results,ro_elements[qubit]["xyl_guess"][0])
            
    """ Storing """
    if update :
        QD_agent.refresh_log("After continuous 2-tone!")
        QD_agent.QD_keeper()
    if chip_info_restore:
        chip_info.update_Cnti2Tone(tt_results)

    """ Close """
    print('2-tone done!')
    shut_down(cluster,Fctrl,Cctrl)
    