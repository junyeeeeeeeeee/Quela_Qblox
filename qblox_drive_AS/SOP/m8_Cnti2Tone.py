
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import array, linspace, NaN, arange
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support.QuFluxFit import calc_Gcoef_inFbFqFd, calc_g
from qblox_drive_AS.support import init_meas, shut_down,  advise_where_fq, init_system_atte, coupler_zctrl, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import QS_fit_analysis, multi_Two_tone_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, twotone_comp_plot
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer
from qblox_drive_AS.analysis.raw_data_demolisher import Conti2tone_dataReducer

def Two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,xyf_guess:dict,xyl_guess:list,IF:float=100e6,xyf_span_Hz:int=400e6,n_avg:int=500,points:int=200,run:bool=True,drive_read_overlap:bool=False)->str:
    sche_func = multi_Two_tone_sche   
    
    original_qubit_info = {}
    for q in xyf_guess:
        original_qubit_info[q] = {}
        qubit_info = QD_agent.quantum_device.get_element(q)
        original_qubit_info[q]["ROW"] = qubit_info.measure.pulse_duration()
        original_qubit_info[q]["ITW"] = qubit_info.measure.integration_time()
        original_qubit_info[q]["XYF"] = qubit_info.clock_freqs.f01()
            
        if not drive_read_overlap:
            drive_pulse_ref_pt = 'start'
            qubit_info.measure.pulse_duration(1.5e-6)
            qubit_info.measure.integration_time(1e-6)
            drive_pulse_length = 100e-6
        else:
            drive_pulse_ref_pt = 'end'
            drive_pulse_length = 100e-6
            qubit_info.measure.pulse_duration(drive_pulse_length)
            qubit_info.measure.integration_time(drive_pulse_length)

        qubit_info.reset.duration(280e-9+drive_pulse_length+0.5e-6) 
        qubit_info.clock_freqs.f01(NaN)
    
    freqs_to_set = {}
    freq_datapoint_idx = arange(0,points)
    for q in xyf_guess:
        f01_high = xyf_guess[q][0]+IF
        freqs_to_set[q] = linspace(f01_high-xyf_span_Hz,f01_high,points)
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=f01_high)
        

    freq = ManualParameter(name="XYfreq", unit="Hz", label="Frequency")
    freq.batched = True
    xy_amp = ManualParameter(name="XYamp", unit="V", label="Amplitude")
    xy_amp.batched = False

    spec_sched_kwargs = dict(   
        frequencies=freqs_to_set,
        spec_amp=xy_amp,
        spec_Du=drive_pulse_length,
        R_amp=compose_para_for_multiplexing(QD_agent,xyf_guess,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,xyf_guess,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,xyf_guess,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,xyf_guess,'r2'),
        ref_pt=drive_pulse_ref_pt
    )

    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables((freq,xy_amp))
        meas_ctrl.setpoints_grid((freq_datapoint_idx,xyl_guess))
        
        qs_ds = meas_ctrl.run("Two-tone")
        qs_ds.attrs["RO_qs"] = ""
        qs_ds.attrs["execution_time"] = Data_manager().get_time_now()
        for idx, q in enumerate(freqs_to_set):
            attr_0 = qs_ds['x0'].attrs
            qs_ds[f'x{2*idx}'] = array(freqs_to_set[q])
            qs_ds[f'x{2*idx}'].attrs = attr_0
            
            attr_1 = qs_ds['x1'].attrs
            qs_ds[f'x{2*idx+1}'] = qs_ds['x1']
            qs_ds[f'x{2*idx+1}'].attrs = attr_1

            qs_ds.attrs["RO_qs"] += f" {q}"

        # Save the raw data into netCDF
        nc_path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb="multiQ",exp_type='2tone',get_data_loc=True)
        
        
    else:
        n_s = 2
        preview_para = {}
        for q in xyf_guess:
            preview_para[q] = xyf_guess[q][:n_s]
        sweep_para2 = array(xyl_guess[:2])
        spec_sched_kwargs['frequencies']= preview_para
        spec_sched_kwargs['spec_amp']= {q:sweep_para2.reshape(sweep_para2.shape or (1,))[0]}
        pulse_preview(QD_agent.quantum_device,sche_func,spec_sched_kwargs)
        nc_path = ""
    
    # restore
    for q in original_qubit_info:
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.measure.pulse_duration(original_qubit_info[q]["ROW"])
        qubit_info.measure.integration_time(original_qubit_info[q]["ITW"])
        qubit_info.clock_freqs.f01(original_qubit_info[q]["XYF"])

    return nc_path

def paras_guess_determinator(QD_path:str, ro_elements:dict, execution:bool=True)->dict:
    """ 
    ro_elements should be a dict with q as the key, each value also should be a dict with the keys : 'xyf_guess', 'xyl_guess', 'g_guess', 'tune_bias'\n
    for example: ro_elements = {'q0':{'xyf_guess': [ ], 'g_guess': 0, 'tune_bias': 0}, ...}
    ### Return the dict contains the xyf_guess for the qubits in ro_elements
    """
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()
    
    xy_freq_settings = {}

    for specific_qubit in ro_elements:
        guess_g = ro_elements[specific_qubit]["g_guess"] if ro_elements[specific_qubit]["g_guess"] != 0 else 90e6
        xyf_guess = list(ro_elements[specific_qubit]["xyf_guess"])

        advised_fq = advise_where_fq(QD_agent,specific_qubit,guess_g) 
        eyeson_print(f"fq advice for {specific_qubit} @ {round(advised_fq*1e-9,4)} GHz")
        if execution:
            if len(xyf_guess) == 0 or (len(xyf_guess) == 1 and xyf_guess[0] == 0):  # xy-freq guess is empty or 0 only
                guess_fq = [advised_fq]
            else:
                guess_fq:list = xyf_guess
        else:
            guess_fq = [4.5e9]
        
        xy_freq_settings[specific_qubit] = guess_fq

    return xy_freq_settings


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

def tune_away_setup(QD_agent:QDmanager,Fctrl:dict,bias_setup:dict,ro_qs:list,zero_mode:bool=False):
    for q in ro_qs:
        if zero_mode:
            Fctrl[q](0.0)
        else:
            offset = QD_agent.Fluxmanager.get_proper_zbiasFor(q)
            if q not in list(bias_setup.keys()):
                bias_setup[q] = 0
            want_bias = offset+bias_setup[q]
            if bias_setup[q] != 0:
                rof = QD_agent.Fluxmanager.sin_for_cav(q,array([want_bias]))[0]
                QD_agent.quantum_device.get_element(q).clock_freqs.readout(rof)
                QD_agent.Fluxmanager.save_tuneawayBias_for('manual',q,want_bias)
                warning_print(f"meas under flux = {round(offset,3)}+{round(bias_setup[q],3)} V")
            Fctrl[q](want_bias) 




def conti2tone_executor(QD_agent:QDmanager,Fctrl:dict,meas_ctrl:MeasurementControl,cluster:Cluster,XYFs:dict,XYLs:list,xyf_span:float=500e6,xy_if:float=100e6,run:bool=True,V_away_from:dict={},drive_read_overlap:bool=False,avg_times:int=500,fpts:int=100)->str:
    
    if run:
        tune_away_setup(QD_agent,Fctrl,V_away_from,list(XYFs.keys()))
        nc_path = Two_tone_spec(QD_agent,meas_ctrl,xyl_guess=XYLs,IF=xy_if,xyf_guess=XYFs,xyf_span_Hz=xyf_span,points=fpts,n_avg=avg_times,run=True,drive_read_overlap=drive_read_overlap) # 
        tune_away_setup(QD_agent,Fctrl,V_away_from,list(XYFs.keys()),True)
        cluster.reset() # *** important       
    else:
        nc_path = Two_tone_spec(QD_agent,meas_ctrl,xyl_guess=XYLs,IF=xy_if,xyf_guess=XYFs,xyf_span_Hz=xyf_span,points=50,n_avg=avg_times,run=False,drive_read_overlap=drive_read_overlap)

    return nc_path
   
    

"""
NOTE: If you find a XYF in a tuneaway z-bias which means the `tune_bias` != 0, please go Modularize/support/meas_switch.py to save it.
"""

if __name__ == "__main__":

    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 0
    update:bool = 1
    #
    DRandIP = {"dr":"dr2","last_ip":"10"}
    #
    xyl_guess:list = linspace(0,0.1,10)
    xyf_setup = {
        "q0":{"xyf_guess":[4.74e9],"g_guess":95e6},   # g you can try a single value about 90e6 for a 5Q4C chip.
        # "q1":{"xyf_guess":[4.21e9],"g_guess":95e6}
    }                                                # tune_bias is the voltage away from sweet spot. If it was given, here will calculate a ROF according to that z-bias and store it in Notebook.
    tunebias_setup = {}
    couplers = []

    """ Optional paras """
    drive_read_overlap:bool = 0
    xy_IF = 100e6
    xyf_range = 500e6
    fpts:int = 100
    avg_n:int = 100


    """ Running """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    ro_elements = paras_guess_determinator(QD_path,xyf_setup,execution)
    
        
    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in xyf_setup:
        QD_agent.Notewriter.save_DigiAtte_For(0,qubit,'xy')
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))


    """ Running """
    tt_results = {}
    nc_path = conti2tone_executor(QD_agent,Fctrl,meas_ctrl,cluster,XYFs=ro_elements,XYLs=xyl_guess,run=execution,xy_if=xy_IF,xyf_span=xyf_range,V_away_from=tunebias_setup,drive_read_overlap=drive_read_overlap,avg_times=avg_n,fpts=fpts)
    slightly_print(f"Raw data loc:\n{nc_path}")


    """ Analysis """
    dss = Conti2tone_dataReducer(nc_path)
    ANA = Multiplex_analyzer("m8")      
    for q in dss:
        ANA._import_data(dss[q],2,QD_agent.refIQ[q],QS_fit_analysis)
        ANA._start_analysis()
        ANA._export_result(Data_manager().get_today_picFolder())

        if ANA.fit_packs != {}:
            analysis_result = QS_fit_analysis(ANA.fit_packs[q]["contrast"],f=ANA.fit_packs[q]["xyf_data"])
            twotone_comp_plot(analysis_result,[],True)
            update_2toneResults_for(QD_agent,q,{str(q):analysis_result},xyl_guess[0])
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
    shut_down(cluster,Fctrl)
