import datetime, os, json
from numpy import array, ndarray, mean, std
from Modularize.RefIQ import Single_shot_ref_spec
from Modularize.FluxQubit import Zgate_two_tone_spec
from Modularize.support import QDmanager, Data_manager, init_system_atte, reset_offset
from quantify_core.measurement.control import MeasurementControl
from Modularize.QuFluxFit import calc_fq_g_excluded, convert_netCDF_2_arrays, data2plot, fq_fit

def sweepZ_arranger(QD_agent:QDmanager,qb:str,Z_guard:float=0.4):
    """
    Split the half flux period in to 6 slices, and composing 4 flux axis.\n
    These 4 flux axis will be return in [[start flux, center flux, end flux], [...], ...]
    """
    start = QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qb)
    flux_step = (QD_agent.Fluxmanager.get_PeriodFor(target_q=qb)/2)/6
    windows = []
    for step_idx in range(0,5,2): # 0,2,4  reject bottom at flux 
        for sign in [1,-1] if step_idx != 0 else [1]:
            center = start + sign*step_idx*flux_step
            window_from = center - flux_step  
            window_end = center + flux_step
            sweet_mark = 0 if step_idx != 0 else 1
            if abs(window_from) <= Z_guard and abs(window_end) <= Z_guard:
                windows.append([window_from,sweet_mark,window_end])
            
    return windows

def Fq_z_coordinator(QD_agent:QDmanager,qb:str,IF:float,span:float,Z_guard:float=0.4):
    windows = sweepZ_arranger(QD_agent,qb,Z_guard)
    meas_var = []
    for window in windows:
        start_rof = QD_agent.Fluxmanager.sin_for_cav(qb,array([window[0]]))[0]
        end_rof = QD_agent.Fluxmanager.sin_for_cav(qb,array([window[-1]]))[0]
        sweet_rof = QD_agent.Fluxmanager.sin_for_cav(qb,array([QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qb)]))[0]
        A = QD_agent.Notewriter.get_CoefInGFor(target_q=qb)
        fb = QD_agent.Notewriter.get_bareFreqFor(target_q=qb)*1e-6
        if window[1]: 
            fq_start = calc_fq_g_excluded(A,sweet_rof*1e-6,fb)*1e6 + IF
        else:
            fq_start = calc_fq_g_excluded(A,start_rof*1e-6,fb)*1e6
        fq_end = calc_fq_g_excluded(A,end_rof*1e-6,fb)*1e6
        fq_range = abs(fq_start-fq_end)
        center_fq = (fq_start+fq_end)/2
        # check fq range > span or not
        if fq_range > span:
            pices = round(fq_range / span)
            if fq_start > fq_end:
                for step in range(pices):
                    fq_pred = fq_start - step*span - IF
                    if center_fq > 2.5e9 and center_fq < 5.5e9:
                        meas_var.append({"z_start":window[0],"fq_predict":fq_pred,"z_end":window[-1]})
                    else:
                        print("Center of the fq range is not in the normal interval [2.5 , 5.5] GHz.")
            else:
                for step in range(pices):
                    fq_pred = fq_end - step*span - IF
                    if center_fq > 2.5e9 and center_fq < 5.5e9:
                        meas_var.append({"z_start":window[0],"fq_predict":fq_pred,"z_end":window[-1]})
                    else:
                        print("Center of the fq range is not in the normal interval [2.5 , 5.5] GHz.")

        else:
            # check center of this span higher than 2.5 GHz or not
            if center_fq > 2.5e9 and center_fq < 5.5e9:
                fq_pred = fq_start-IF if fq_start > fq_end else fq_end-IF
                meas_var.append({"z_start":window[0],"fq_predict":fq_pred,"z_end":window[-1]})
            else:
                print("Center of the fq range is not in the normal interval [2.5 , 5.5] GHz.")

    if len(meas_var) <= 3:
        z_pts = 20
    else:
        z_pts = 10
    print(f"Total {len(meas_var)} pics for Flux vs. Fq exp.")
    return meas_var, z_pts

    
def mag_static_filter(x_array:ndarray,y_array:ndarray,mag_array:ndarray,threshold:float=2.0):
    x_choosen = []
    y_choosen = []
    while True:
        bottom_line = mean(mag_array) - threshold*std(mag_array)
        for idx in range(x_array.shape[0]):
            if mag_array[idx] > bottom_line:
                x_choosen.append(x_array[idx])
                y_choosen.append(y_array[idx])
        if len(x_choosen) > 0 :
            break
        else:
            if threshold != 0.5:
                threshold -= 0.5
            else:
                print(f"No peaks are found in this Z interval: {x_array[0]}~{x_array[-1]} V")
    
    return x_choosen, y_choosen

# main execute program
def FluxFqFit_execution(QD_agent:QDmanager, meas_ctrl:MeasurementControl, Fctrl:dict, target_q:str, execute:bool=True, f_pts:int=30, re_n:int=500):
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(target_q,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(target_q,'xy'))
    original_rof = QD_agent.quantum_device.get_element(target_q).clock_freqs.readout()
    # ROF is at zero bias
    rof = QD_agent.Fluxmanager.sin_for_cav(target_q,array([0]))[0]
    QD_agent.quantum_device.get_element(target_q).clock_freqs.readout(rof)
    # get the ref IQ according to the ROF (zero bias)
    analysis_result = Single_shot_ref_spec(QD_agent,q=target_q,want_state='g',shots=10000)
    I_ref, Q_ref= analysis_result[target_q]['fit_pack'][0],analysis_result[target_q]['fit_pack'][1]
    ref = [I_ref,Q_ref]
    f_span = 500e6
    xy_if = 100e6
    meas_vars, z_pts = Fq_z_coordinator(QD_agent,target_q,xy_if,f_span)
    x2static = []
    y2static = []
    mag2static=[]
    for window in meas_vars:
        # main execution
        results, raw_path, _ = Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=window["z_start"],Z_amp_end=window["z_end"],q=target_q,run=execute,get_data_path=True,Z_points=z_pts,f_points=f_pts,n_avg=re_n,xyf=window["fq_predict"],xyf_span_Hz=f_span,IF=xy_if)
        reset_offset(Fctrl)
        xyf,z,ii,qq = convert_netCDF_2_arrays(raw_path)
        # below can be switch into QM system with the arg `qblox=False`
        x, y, mag = data2plot(xyf,z,ii,qq,specified_refIQ=ref,qblox=True,q=target_q)
        x2static += x.tolist()
        y2static += y.tolist()
        mag2static+=mag.tolist()
    
    x2fit, y2fit = mag_static_filter(array(x2static),array(y2static),array(mag2static))
    data2fit = {"x":x2fit,"y":y2fit}
    json_path = Data_manager().save_dict2json(QD_agent,data2fit,target_q,True)

    pic_parentpath = os.path.join(Data_manager().get_today_picFolder())
    fq_fit(QD_agent,json_path,target_q,savefig_path=pic_parentpath,saveParas=True) 
    QD_agent.quantum_device.get_element(target_q).clock_freqs.readout(original_rof)
    QD_agent.QD_keeper()

    

if __name__ == "__main__":
    from Modularize.support import init_meas, shut_down, reset_offset
    from numpy import pi, absolute
    import matplotlib.pyplot as plt
   
    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_21/DR2#171_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    
    # Set system attenuation
    # init_system_atte(QDmanager.quantum_device,list(Fctrl.keys()))
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)
    
    execution = True
    for qb in ["q0"]:
        FluxFqFit_execution(QD_agent, meas_ctrl, Fctrl, target_q=qb, execute=execution)

    print('done!')
    shut_down(cluster,Fctrl)
    
        
        

