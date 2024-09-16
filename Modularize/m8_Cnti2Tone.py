
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

def Two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,IF:float=100e6,f01_guess:int=0,xyf_span_Hz:int=400e6,xyamp:float=0.02,n_avg:int=500,points:int=200,run:bool=True,q:str='q1',Experi_info:dict={},ref_IQ:list=[0,0],drive_read_overlap:bool=False):
    sche_func = Two_tone_sche   
    analysis_result = {}
    qubit_info = QD_agent.quantum_device.get_element(q)
    original_row = qubit_info.measure.pulse_duration()
    original_inteW = qubit_info.measure.integration_time()
    if not drive_read_overlap:
        drive_pulse_ref_pt = 'start'
        qubit_info.measure.pulse_duration(2e-6)
        qubit_info.measure.integration_time(1e-6)
        drive_pulse_length = 50e-6
    else:
        drive_pulse_ref_pt = 'end'
        drive_pulse_length = 100e-6
        qubit_info.measure.pulse_duration(drive_pulse_length)
        qubit_info.measure.integration_time(drive_pulse_length)

    qubit_info.reset.duration(250e-6)
    eyeson_print(f"RO amp = {qubit_info.measure.pulse_amp()}")
    
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
        spec_Du=drive_pulse_length,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        ref_pt=drive_pulse_ref_pt
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

    qubit_info.measure.pulse_duration(original_row)
    qubit_info.measure.integration_time(original_inteW)

    return analysis_result

def paras_guess_determinator(QD_path:str, specific_qubits:str, execution:bool=True, xyf_guess:list=[], guess_g:float=48e6, xyAmp_guess:list=[])->tuple[list, list]:
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()
    advised_fq = advise_where_fq(QD_agent,specific_qubits,guess_g) 
    print(f"fq advice for {specific_qubits} @ {round(advised_fq*1e-9,4)} GHz")
    if execution:
        if len(xyf_guess) == 0 or (len(xyf_guess) == 1 and xyf_guess[0] == 0):
            guess_fq = [advised_fq-500e6, advised_fq, advised_fq+500e6]
        else:
            guess_fq = xyf_guess

        if len(xyAmp_guess) == 0 or (len(xyAmp_guess) == 1 and xyAmp_guess[0] == 0):
            xyAmp_guess = [0, 0.1]
        else:
            xyAmp_guess = xyAmp_guess
    else:
        guess_fq = [4.5e9]
        xyAmp_guess = [0.05]

    return guess_fq, xyAmp_guess


def update_2toneResults_for(QD_agent:QDmanager,qb:str,QS_results:dict,XYL:float):
    qubit = QD_agent.quantum_device.get_element(qb)
    Revised_f01 = QS_results[qb].attrs['f01_fit']
    fb = float(QD_agent.Notewriter.get_bareFreqFor(target_q=qb))*1e-6
    fd = QD_agent.quantum_device.get_element(qb).clock_freqs.readout()*1e-6
    A = calc_Gcoef_inFbFqFd(fb,Revised_f01*1e-6,fd)
    sweet_g = calc_g(fb,Revised_f01*1e-6,A)
    qubit.clock_freqs.f01(Revised_f01)
    QD_agent.Notewriter.save_2tone_piamp_for(qb,XYL)
    QD_agent.Notewriter.save_CoefInG_for(target_q=qb,A=A)
    QD_agent.Notewriter.save_sweetG_for(target_q=qb,g_Hz=sweet_g*1e6)



def conti2tone_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,cluster:Cluster,specific_qubits:str,XYF:float,XYL:float,xyf_span:float=500e6,xy_if:float=100e6,run:bool=True,V_away_from:float=0,drive_read_overlap:bool=False,avg_times:int=500,fpts:int=100):
    
    if run:
        ori_data = []
        offset = QD_agent.Fluxmanager.get_proper_zbiasFor(specific_qubits)
        want_bias = offset+V_away_from
        if V_away_from != 0:
            rof = QD_agent.Fluxmanager.sin_for_cav(specific_qubits,array([want_bias]))[0]
            QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(rof)
            QD_agent.Fluxmanager.save_tuneawayBias_for('manual',specific_qubits,want_bias)
            warning_print(f"meas under flux = {round(offset,3)}+{round(V_away_from,3)} V")
        Fctrl[specific_qubits](want_bias) 
        
        QS_results = Two_tone_spec(QD_agent,meas_ctrl,xyamp=XYL,IF=xy_if,f01_guess=XYF,q=specific_qubits,xyf_span_Hz=xyf_span,points=fpts,n_avg=avg_times,run=True,ref_IQ=QD_agent.refIQ[specific_qubits],drive_read_overlap=drive_read_overlap) # 
        Fctrl[specific_qubits](0.0)
        
        cluster.reset() # *** important
        if XYL != 0:
            return QS_results[specific_qubits]
        else:
            ori_data = QS_results[specific_qubits].data_vars['data']
            return ori_data
                
    else:
        qu = specific_qubits
        QS_results = Two_tone_spec(QD_agent,meas_ctrl,xyamp=0.1,IF=xy_if,f01_guess=4e9,q=qu,xyf_span_Hz=xyf_span,points=50,n_avg=avg_times,run=False,ref_IQ=QD_agent.refIQ[qu],drive_read_overlap=drive_read_overlap)

        return 0
   
    

"""
NOTE: If you find a XYF in a tuneaway z-bias which means the `tune_bias` != 0, please go Modularize/support/meas_switch.py to save it.
"""

if __name__ == "__main__":

    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 0
    update:bool = 1
    #
    DRandIP = {"dr":"dr4","last_ip":"81"}
    #
    ro_elements = {
        "q0":{"xyf_guess":[5.11e9],"xyl_guess":[0.06],"g_guess":45e6, "tune_bias":0} # g you can try a single value about 90e6 for a 5Q4C chip.
    }                                                                            # tune_bias is the voltage away from sweet spot. If it was given, here will calculate a ROF according to that z-bias and store it in Notebook.
    couplers = []

    """ Optional paras """
    drive_read_overlap:bool = 0
    xy_IF = 100e6
    xyf_range = 500e6
    fpts:int = 100
    avg_n:int = 1000


    """ Running """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    for qubit in ro_elements:
        background = []
        g = 90e6 if ro_elements[qubit]["g_guess"] == 0 else ro_elements[qubit]["g_guess"]
        fq_guesses, xyl_guesses = paras_guess_determinator(QD_path,qubit,execution,ro_elements[qubit]["xyf_guess"],g,ro_elements[qubit]["xyl_guess"])
        for xyf in fq_guesses:
            for xyl in xyl_guesses: 
                """ Preparations """
                QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
                chip_info = cds.Chip_file(QD_agent=QD_agent)

                """ Running """
                tt_results = {}
                Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            
                QD_agent.Notewriter.save_DigiAtte_For(0,qubit,'xy')
                init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
                tune_bias = ro_elements[qubit]["tune_bias"]
                print(xyf)
                tt_results = conti2tone_executor(QD_agent,meas_ctrl,cluster,specific_qubits=qubit,XYF=xyf,XYL=xyl,run=execution,xy_if=xy_IF,xyf_span=xyf_range,V_away_from=tune_bias,drive_read_overlap=drive_read_overlap,avg_times=avg_n,fpts=fpts)

                if xyl == 0: 
                    background = tt_results
                else:
                    twotone_comp_plot(tt_results, background, True)
                
                if execution and ro_elements[qubit]["xyl_guess"][0] != 0 and update:
                    update = True
                    print(f'update xyl={ro_elements[qubit]["xyl_guess"][0]}')
                    update_2toneResults_for(QD_agent,qubit,{str(qubit):tt_results},xyl)
                else:
                    update = False
                    
                """ Storing """
                if update:
                    QD_agent.refresh_log("After continuous 2-tone!")
                    QD_agent.QD_keeper()
                    if chip_info_restore:
                        chip_info.update_Cnti2Tone({str(qubit):tt_results})

                """ Close """
                print('2-tone done!')
                shut_down(cluster,Fctrl,Cctrl)
    