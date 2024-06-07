import datetime, json
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')) 
from numpy import array, ndarray, mean, std, where, sort
from Modularize.m7_RefIQ import Single_shot_ref_spec
from Modularize.support.UserFriend import *
from Modularize.m9_FluxQubit import Zgate_two_tone_spec
from qblox_instruments import Cluster
from Modularize.support import QDmanager, Data_manager, init_system_atte, reset_offset, shut_down, init_meas, coupler_zctrl
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support.QuFluxFit import calc_fq_g_excluded, convert_netCDF_2_arrays, data2plot, fq_fit

def sweepZ_arranger(QD_agent:QDmanager,qb:str,Z_guard:float=0.4):
    """
    Split the half flux period in to 6 slices, and composing 4 flux axis.\n
    These 4 flux axis will be return in [[start flux, center flux, end flux], [...], ...]
    """
    start = QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qb)
    flux_step = (QD_agent.Fluxmanager.get_PeriodFor(target_q=qb)/2)/6
    windows = []
    for step_idx in range(0,3,2): # 0,2,4  reject bottom at flux 

        for sign in [1,-1] if step_idx != 0 else [1]:
            center = sign*step_idx*flux_step
            window_from = center - flux_step  
            window_end = center + flux_step
            sweet_mark = 0 if step_idx != 0 else 1
            if abs(window_from) <= Z_guard and abs(window_end) <= Z_guard:
                windows.append([window_from,sweet_mark,window_end])
            
    return windows

def window_picker(meas_var:list,maxi_window_num:int=5)->list:
    """
    pick the `maxi_window_num` highest fq windows
    """
    tran_freq = []
    for item in meas_var:
        tran_freq.append(item["fq_predict"])
    
    max_idx = where(array(tran_freq)>sort(array(tran_freq))[-(maxi_window_num+1)])[0]

    picked_item = []
    
    for i in max_idx:
        picked_item.append(meas_var[i])

    return picked_item



def Fq_z_coordinator(QD_agent:QDmanager,qb:str,IF:float,span:float,Z_guard:float=0.4):
    windows = sweepZ_arranger(QD_agent,qb,Z_guard)
    sweet_spot_v = QD_agent.Fluxmanager.get_sweetBiasFor(qb)
    meas_var = []
    for window in windows:
        start_rof = QD_agent.Fluxmanager.sin_for_cav(qb,array([window[0]+sweet_spot_v]))[0]
        end_rof = QD_agent.Fluxmanager.sin_for_cav(qb,array([window[-1]+sweet_spot_v]))[0]
        sweet_rof = QD_agent.Fluxmanager.sin_for_cav(qb,array([sweet_spot_v]))[0]
        A = QD_agent.Notewriter.get_CoefInGFor(target_q=qb)
        fb = QD_agent.Notewriter.get_bareFreqFor(target_q=qb)*1e-6
        if window[1]: 
            fq_start = calc_fq_g_excluded(A,sweet_rof*1e-6,fb)*1e6 + IF
        else:
            fq_start = calc_fq_g_excluded(A,start_rof*1e-6,fb)*1e6
        fq_end = calc_fq_g_excluded(A,end_rof*1e-6,fb)*1e6
        center_fq = (fq_start+fq_end)/2
        
        # check fq range > span or not
        # fq_range = abs(fq_start-fq_end)
        # if fq_range > span:
        #     pices = round(fq_range / span)
        #     if fq_start > fq_end:
        #         for step in range(pices):
        #             fq_pred = fq_start - step*span - IF
        #             if center_fq > 2.5e9 and center_fq < 5.5e9:
        #                 meas_var.append({"z_start":window[0],"fq_predict":fq_pred,"z_end":window[-1]})
        #             else:
        #                 print("Center of the fq range is not in the normal interval [2.5 , 5.5] GHz.")
        #     else:
        #         for step in range(pices):
        #             fq_pred = fq_end - step*span - IF
        #             if center_fq > 2.5e9 and center_fq < 5.5e9:
        #                 meas_var.append({"z_start":window[0],"fq_predict":fq_pred,"z_end":window[-1]})
        #             else:
        #                 print("Center of the fq range is not in the normal interval [2.5 , 5.5] GHz.")

        # else:
        #     # check center of this span higher than 2.5 GHz or not
        #     if center_fq > 2.5e9 and center_fq < 5.5e9:
        #         fq_pred = fq_start-IF if fq_start > fq_end else fq_end-IF
        #         meas_var.append({"z_start":window[0],"fq_predict":fq_pred,"z_end":window[-1]})
        #     else:
        #         print("Center of the fq range is not in the normal interval [2.5 , 5.5] GHz.")

        if center_fq > 2.5e9 and center_fq < 5.5e9:
            meas_var.append({"z_start":window[0],"fq_predict":center_fq,"z_end":window[-1]})
    
    # if len(meas_var) <= 3:
    #     z_pts = 20
    # else:
    #     z_pts = 10
    z_pts = 10
        

    eyeson_print(f"Total {len(meas_var)} pics for Flux vs. Fq exp.")
    return meas_var, z_pts

    
def mag_static_filter(x_array:ndarray,y_array:ndarray,mag_array:ndarray,threshold:float=3.0):
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
def FluxFqFit_execution(QD_agent:QDmanager, meas_ctrl:MeasurementControl, Fctrl:dict, cluster:Cluster,target_q:str, run:bool=True, f_pts:int=30, re_n:int=500, peak_threshold:float=3):
    f_span = 500e6
    xy_if = 100e6
    failed = False
    original_rof = QD_agent.quantum_device.get_element(target_q).clock_freqs.readout()
    if run:
        meas_vars, z_pts = Fq_z_coordinator(QD_agent,target_q,xy_if,f_span)
        x2static = []
        y2static = []
        mag2static=[]
        n = 0
        for window in meas_vars:
            # main execution
            ref_z = QD_agent.Fluxmanager.get_sweetBiasFor(target_q)
            cluster.reset()
            Fctrl[target_q](ref_z)
            _, raw_path, _ = Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=window["z_start"],Z_amp_end=window["z_end"],q=target_q,run=True,get_data_path=True,Z_points=z_pts,f_points=f_pts,n_avg=re_n,xyf=window["fq_predict"],xyf_span_Hz=f_span,IF=xy_if,analysis=False)
            reset_offset(Fctrl)
            
            xyf,z,ii,qq = convert_netCDF_2_arrays(raw_path)
            # below can be switch into QM system with the arg `qblox=False`
            x, y, mag = data2plot(xyf,z+ref_z,ii,qq,specified_refIQ=QD_agent.refIQ[target_q],qblox=True,q=target_q,filter2D_threshold=peak_threshold,plot_scatter=1)
            x2static += x.tolist()
            y2static += y.tolist()
            mag2static+=mag.tolist()
            
        
        x2fit, y2fit = mag_static_filter(array(x2static),array(y2static),array(mag2static))
        data2fit = {"x":x2fit,"y":y2fit}
        json_path = Data_manager().save_dict2json(QD_agent,data2fit,target_q,True)

        pic_parentpath = os.path.join(Data_manager().get_today_picFolder())
        try:
            fq_fit(QD_agent,json_path,target_q,savefig_path=pic_parentpath,saveParas=True,plot=False,FitFilter_threshold=2.5) 
        except:
            failed = True
        QD_agent.quantum_device.get_element(target_q).clock_freqs.readout(original_rof)
    else:
        meas_vars, z_pts = Fq_z_coordinator(QD_agent,target_q,xy_if,f_span)
        window = meas_vars[0]
        # main execution
        _, raw_path, _ = Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=window["z_start"],Z_amp_end=window["z_end"],q=target_q,run=False,get_data_path=True,Z_points=z_pts,f_points=f_pts,n_avg=re_n,xyf=window["fq_predict"],xyf_span_Hz=f_span,IF=xy_if)
        reset_offset(Fctrl)
        
    return failed

    

if __name__ == "__main__":

   
    """ Fill in """
    execution = True
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = ['q0']
    couplers = ["c0"]


    """ Preparations"""
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    if ro_elements == 'all':
        ro_elements = list(Fctrl.keys())
   
    
    """ Running """
    fit_error = []
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        if QD_agent.Fluxmanager.get_offsweetspot_button(qubit): raise ValueError("m10 should be performed at sweet spot, now is deteced in off-sweetspot mode!")
        init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
    
        error = FluxFqFit_execution(QD_agent, meas_ctrl, Fctrl, cluster, target_q=qubit, run=execution,peak_threshold=2)
        if error:
            fit_error.append(qubit)


    
    """ Storing (Future) """
    if execution:
        pass # QD_agent.QD_keeper()


    """ Close """
    print('done!')
    print(f"fitting error occured at {fit_error}")
    shut_down(cluster,Fctrl,Cctrl)
    
        
        

