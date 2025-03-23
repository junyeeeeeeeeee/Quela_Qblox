import os, sys, json 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager
import matplotlib.pyplot as plt
from numpy import ndarray, array, median, std, mean
from datetime import datetime 
from qblox_drive_AS.support.UserFriend import *
from matplotlib.gridspec import GridSpec as GS
from qblox_drive_AS.support.ExpFrames import nSingleShot, SpinEcho, CPMG, Ramsey, EnergyRelaxation


def time_label_sort(nc_file_name:str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0],"H%HM%MS%S")

def plot_coherence_timetrace(raw_data:ndarray, time_samples:ndarray, ans:ndarray, q:str, raw_data_folder:str, exp:str, detunings:ndarray=[]):
    """
    Plot color map, histogram, and detuning if ramsey in one figure
    """
    time_json_path = [os.path.join(raw_data_folder,name) for name in os.listdir(raw_data_folder) if (os.path.isfile(os.path.join(raw_data_folder,name)) and name.split(".")[-1] == "json")][0]
    with open(time_json_path) as time_record_file:
        time_past_dict:dict = json.load(time_record_file)
    time_array = array(list(time_past_dict.values())[0])

    median_ans, std_ans = median(ans), std(ans)
    fig = plt.figure(dpi=100, figsize=(12,9))
    
    gs = GS(2,2, width_ratios=[2,1],height_ratios=[1,5])
    if exp.lower() == 't2':
        ax0 = fig.add_subplot(gs[0,0])
        d = ax0.plot(time_array,detunings,c='magenta')
        ax0.set_xlim(min(time_array),max(time_array))
        ax0.set_title("Transition detuning (MHz)",fontsize=20)
        ax0.grid()
        ax0.xaxis.set_tick_params(labelsize=16)
        ax0.yaxis.set_tick_params(labelsize=16)

    ax1 = fig.add_subplot(gs[:,1])
    ax1.hist(ans, bins='auto', density=False)
    ax1.axvline(median_ans,c='k',ls='--',lw='1')
    ax1.set_ylabel("Counts",fontsize=20)
    ax1.set_xlabel(f"{exp.upper()} (us)", fontsize=20)
    ax1.set_title(f"{exp.upper()} = {round(median_ans,1)} $\pm$ {round(std_ans,1)} us",fontsize=20)

    
    if exp.lower() in ['t1', 'os']:
        ax2 = fig.add_subplot(gs[:,0])
    elif exp.lower() == 't2':
        ax2 = fig.add_subplot(gs[1,0])
    else:
        pass
    c = ax2.pcolormesh(time_array, time_samples*1e6, raw_data.transpose()*1000, cmap='RdBu')
    fig.colorbar(c, ax=ax2, label='Contrast (mV)',location='bottom',)
    ax2.plot(time_array,ans, c='green')
    ax2.axhline(median_ans+std_ans,linestyle='--',c='orange')
    ax2.axhline(median_ans,linestyle='-',c='orange')
    ax2.axhline(median_ans-std_ans,linestyle='--',c='orange')
    ax2.set_xlabel('Time past (min)',fontsize=20)
    ax2.set_ylabel("Free evolution time (us)",fontsize=20)
    

    for Ax in [ax1, ax2]:
        Ax:plt.Axes
        Ax.xaxis.set_tick_params(labelsize=16)
        Ax.yaxis.set_tick_params(labelsize=16)

    plt.title(f"Time dependent {exp.upper()}",fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(raw_data_folder,f"{exp.upper()}_{q}_timeDep.png"))
    plt.close()

def colormap(x:ndarray, y:ndarray, z:ndarray, fit_values:ndarray, ax:plt.Axes=None, fig_path:str=None):
    if ax is None:
        fig, ax = plt.subplots()
    c = ax.pcolormesh(x, y, z.transpose(), cmap='RdBu')
    ax.plot(x,fit_values, c='green')
    ax.axhline(median(fit_values)+std(fit_values),linestyle='--',c='orange')
    ax.axhline(median(fit_values),linestyle='-',c='orange')
    ax.axhline(median(fit_values)-std(fit_values),linestyle='--',c='orange')
    ax.set_xlabel('Time past (min)',fontsize=20)
    ax.set_ylabel("Free evolution time (us)",fontsize=20)
    plt.colorbar(c, ax=ax, label='I channel (mV)')
    plt.title(os.path.split(fig_path)[-1].split(".")[0])
    plt.tight_layout()
    if fig_path is not None:
        plt.savefig(fig_path)
    plt.close()

def plot_timeDepCohe(time_values:ndarray, y_values:ndarray, exp:str, fig_path:str=None, units:dict={"x":"min","y":"µs"}):
    fig, axs = plt.subplots(1,3,figsize=(16,8))
    ax:plt.Axes = axs[0]
    ax.grid()
    ax.plot(time_values, y_values)
    ax.axhline(median(y_values)+std(y_values),linestyle='--',c='orange',label='1σ')
    ax.axhline(mean(y_values),linestyle='-',c='orange',label='mean')
    ax.axhline(median(y_values)-std(y_values),linestyle='--',c='orange')
    ax.axhline(median(y_values),linestyle='-.',c='red',label='median')
    ax.set_xlabel(f"time past ({units['x']})",fontsize=16)
    ax.set_ylabel(f"{exp.upper()} ({units['y']})",fontsize=16)
    ax.set_title(f"Time dependent {exp.upper()}",fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.legend()

    ax:plt.Axes=axs[1]
    ax.grid()
    ax.hist(y_values, bins='auto', density=False,orientation="horizontal")
    ax.axhline(median(y_values),c='k',ls='--',lw='1')
    ax.set_xlabel("Counts",fontsize=16)
    ax.set_ylabel(f"{exp.upper()} ({units['y']})", fontsize=16)
    ax.set_title(f"{exp.upper()} = {round(median(y_values),2 if exp == 'δf' else 1)} $\pm$ {round(std(y_values),2 if exp == 'δf' else 1)} {units['y']}",fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)

    ax:plt.Axes=axs[2]
    ax.hist(y_values, bins=5000, orientation="horizontal", cumulative=True, density=True,histtype="step")
    ax.axhline(median(y_values),c='k',ls='--',lw='1')
    ax.axvline(0.5,c='red',ls='--',lw='1',label='density=0.5')
    ax.set_xlabel("Density",fontsize=16)
    ax.set_ylabel(f"{exp.upper()} ({units['y']})", fontsize=16)
    ax.set_title("CDF",fontsize=20)
    ax.legend()
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)

    plt.tight_layout(pad=3.0)
    if fig_path is not None:
        plt.savefig(fig_path)
    plt.close()

def time_monitor_data_ana(QD_agent:QDmanager,folder_path:str,save_every_fit_pic:bool=False):
    files = [name for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[-1] == "nc")]

    # T1, T2, SS, T2s, T2cpmg = 0, 0 , 0, 0, 0 # counter 
    counters = {}
    cpmg_pi_nums = []
    T1_rec, detu_rec, T2_rec, SS_rec = {}, {}, {}, {}
    T1_raw, T2_raw = {}, {}
    T1_evo_time, T2_evo_time = {}, {}
    for idx, file in enumerate(files) :
        slightly_print(f"Analysis for the {idx}-th files ...")
        exp_type:str = file.split("_")[0]
        path = os.path.join(folder_path,file)
        ds = open_dataset(path)
        match exp_type.lower():
            case "t1":
                counters["T1"] = 1 if "T1" not in list(counters.keys()) else counters["T1"]+1
                for var in [ var for var in ds.data_vars if var.split("_")[-1] != 'x'] :
                    if counters["T1"] == 1:
                        T1_picsave_folder = os.path.join(folder_path,"T1_pics") if save_every_fit_pic else None
                        T1_rec[var], T1_raw[var] = {}, {}
                        if save_every_fit_pic:
                            if not os.path.exists(T1_picsave_folder):
                                os.mkdir(T1_picsave_folder)
        
                    T1_evo_time[var] = array(ds[f"{var}_x"])[0][0]
                this_exp_time = ds.attrs["end_time"]
                ds.close()

                Anaer = EnergyRelaxation(QD_path="")
                Anaer.keep_QD = False
                Anaer.save_pics = False
                Anaer.execution = True
                Anaer.histos = 1
                Anaer.RunAnalysis(new_file_path=path,new_QDagent=QD_agent,new_pic_save_place=T1_picsave_folder)

                for var in Anaer.sum_dict:
                    T1_raw[var][this_exp_time] = Anaer.sum_dict[var]["plot_item"]["data"]
                    T1_rec[var][this_exp_time] = Anaer.sum_dict[var]["median_T1"]
            
            case "singleshot":
                counters["SS"] = 1 if "SS" not in list(counters.keys()) else counters["SS"]+1
                for var in ds.data_vars:
                    if counters["SS"] == 1:
                        SS_picsave_folder = os.path.join(folder_path,f"SingleShot_pics") if save_every_fit_pic else None
                        SS_rec[var] = {}
                        if save_every_fit_pic:
                            if not os.path.exists(SS_picsave_folder):
                                os.mkdir(SS_picsave_folder)
                this_exp_time = ds.attrs["end_time"]
                ds.close()
                Anaer = nSingleShot(QD_path="")
                Anaer.keep_QD = False
                Anaer.save_pics = False
                Anaer.execution = True
                Anaer.save_os_model = False
                Anaer.histos = 1
                Anaer.RunAnalysis(new_file_path=path,new_QDagent=QD_agent,new_pic_save_place=SS_picsave_folder)
                
                for var in Anaer.sum_dict:
                    SS_rec[var][this_exp_time] = Anaer.sum_dict[var]["effT_mK"][0]
            
            case "ramsey":
                counters["T2s"] = 1 if "T2s" not in list(counters.keys()) else counters["T2s"]+1

                for var in [ var for var in ds.data_vars if var.split("_")[-1] != 'x']:
                    # create raw data fitting folder
                    if counters["T2s"] == 1:
                        T2_picsave_folder = os.path.join(folder_path,f"ramsey_pics") if save_every_fit_pic else None
                        if "ramsey" not in list(T2_rec.keys()):
                            T2_rec["ramsey"], T2_raw["ramsey"] = {}, {}
                        T2_rec["ramsey"][var], detu_rec[var], T2_raw["ramsey"][var] = {}, {}, {}
                        if save_every_fit_pic:
                            if not os.path.exists(T2_picsave_folder):
                                os.mkdir(T2_picsave_folder)
                    # start analysis 
                    T2_evo_time[var] = array(ds[f"{var}_x"])[0][0]
                this_exp_time = ds.attrs["end_time"]
                ds.close()
                
                Anaer = Ramsey(QD_path="")
                Anaer.sec_phase = 'y'
                Anaer.keep_QD = False
                Anaer.save_pics = False
                Anaer.execution = True
                Anaer.histos = 1
                Anaer.RunAnalysis(new_file_path=path,new_QDagent=QD_agent,new_pic_save_place=T2_picsave_folder)
                
                # keep values
                for var in Anaer.sum_dict:
                    T2_raw["ramsey"][var][this_exp_time] = Anaer.sum_dict[var]["plot_item"]["data"]
                    detu_rec[var][this_exp_time] = Anaer.corrected_detune[var]*1e-6
                    T2_rec["ramsey"][var][this_exp_time] = Anaer.sum_dict[var]["median_T2"]

            case "spinecho":
                counters["T2"] = 1 if "T2" not in list(counters.keys()) else counters["T2"]+1
                for var in [ var for var in ds.data_vars if var.split("_")[-1] != 'x']:
                    # create raw data fitting folder
                    
                    if counters["T2"] == 1:
                        T2_picsave_folder = os.path.join(folder_path,f"SpinEcho_pics") if save_every_fit_pic else None
                        if "spinecho" not in list(T2_rec.keys()):
                            T2_rec["spinecho"], T2_raw["spinecho"] = {}, {}
                        T2_rec["spinecho"][var], T2_raw["spinecho"][var] = {}, {}
                        if save_every_fit_pic:
                            if not os.path.exists(T2_picsave_folder):
                                os.mkdir(T2_picsave_folder)
                    # start analysis 
                    T2_evo_time[var] = array(ds[f"{var}_x"])[0][0]
                this_exp_time = ds.attrs["end_time"]
                ds.close()

                Anaer = SpinEcho(QD_path="")
                Anaer.keep_QD = False
                Anaer.save_pics = False
                Anaer.execution = True
                Anaer.histos = 1
                Anaer.RunAnalysis(new_file_path=path,new_QDagent=QD_agent,new_pic_save_place=T2_picsave_folder)

                # keep values
                for var in Anaer.sum_dict:
                    T2_raw["spinecho"][var][this_exp_time] = Anaer.sum_dict[var]["plot_item"]["data"]
                    T2_rec["spinecho"][var][this_exp_time] = Anaer.sum_dict[var]["median_T2"]
            
            case "cpmg": # There may do different pi_num CPMG
                a_q = [q for q in list(ds.data_vars) if q[0]=='q' and "_" not in q][0]
                pi_num =  ds[a_q].attrs["spin_num"]
                if pi_num not in cpmg_pi_nums: cpmg_pi_nums.append(pi_num) 
                counters[f"T2cpmg_{pi_num}"] = 1 if f"T2cpmg_{pi_num}" not in list(counters.keys()) else counters[f"T2cpmg_{pi_num}"]+1
                
                for var in [ var for var in ds.data_vars if var.split("_")[-1] != 'x']:
                    # create raw data fitting folder
                    if counters[f"T2cpmg_{pi_num}"] == 1:
                        T2_picsave_folder = os.path.join(folder_path,f"CPMG{pi_num}_pics") if save_every_fit_pic else None
                        if "T2cpmg" not in list(T2_rec.keys()):
                            T2_rec["T2cpmg"], T2_raw["T2cpmg"] = {}, {}
                        if pi_num not in list(T2_rec["T2cpmg"].keys()):
                            T2_rec["T2cpmg"][pi_num], T2_raw["T2cpmg"][pi_num] = {}, {}
                        T2_rec["T2cpmg"][pi_num][var], T2_raw["T2cpmg"][pi_num][var] = {}, {}
                        if save_every_fit_pic:
                            if not os.path.exists(T2_picsave_folder):
                                os.mkdir(T2_picsave_folder)
                        
                    # start analysis 
                    T2_evo_time[var] = array(ds[f"{var}_x"])[0][0]
                this_exp_time = ds.attrs["end_time"]
                ds.close()

                Anaer = CPMG(QD_path="")
                Anaer.keep_QD = False
                Anaer.save_pics = False
                Anaer.execution = True
                Anaer.histos = 1
                Anaer.RunAnalysis(new_file_path=path,new_QDagent=QD_agent,new_pic_save_place=T2_picsave_folder)
                    
                # keep values
                for var in Anaer.sum_dict:
                    T2_raw["T2cpmg"][pi_num][var][this_exp_time] = Anaer.sum_dict[var]["plot_item"]["data"]
                    T2_rec["T2cpmg"][pi_num][var][this_exp_time] = Anaer.sum_dict[var]["median_T2"]
            
    
    items_2_plot = list(counters.keys())
    pic_folder = os.path.join(folder_path, "Pics")
    if not os.path.exists(pic_folder):
        os.mkdir(pic_folder)

    slightly_print(f"\nPlotting... ")
    if "T1" in items_2_plot:
        for q in T1_rec:
            
            sorted_item_ans = sorted(T1_rec[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
            sorted_item_raw = sorted(T1_raw[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
            earliest_time = datetime.strptime(sorted_item_ans[0][0], "%Y-%m-%d %H:%M:%S")
            time_diffs = []
            sorted_values_ans, sorted_values_raw = [], []
            for idx, item in enumerate(sorted_item_ans):
                key = item[0]
                value = item[1]
                current_time = datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
                time_diff = round((current_time - earliest_time).total_seconds()/60,1)
                time_diffs.append(time_diff)
                sorted_values_ans.append(value)
                sorted_values_raw.append(sorted_item_raw[idx][1])

            
            colormap(array(time_diffs),array(T1_evo_time[q])*1e6,array(sorted_values_raw),array(sorted_values_ans),fig_path=os.path.join(pic_folder,f"{q}_T1_timeDep_colormap.png"))
            plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "t1", units={"x":"min","y":"µs"}, fig_path=os.path.join(pic_folder,f"{q}_T1_timeDep.png"))
        items_2_plot.remove("T1")
    if "SS" in list(counters.keys()):
        for q in SS_rec:
            sorted_item_ans = sorted(SS_rec[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
            earliest_time = datetime.strptime(sorted_item_ans[0][0], "%Y-%m-%d %H:%M:%S")
            
            time_diffs = []
            sorted_values_ans = []
            for key, value in sorted_item_ans:
                current_time = datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
                time_diff = round((current_time - earliest_time).total_seconds()/60,1)
                time_diffs.append(time_diff)
                sorted_values_ans.append(value)


            plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "eff_Temp.", units={"x":"min","y":"mK"}, fig_path=os.path.join(pic_folder,f"{q}_effT_timeDep.png"))
        items_2_plot.remove("SS")
    
    if  len(items_2_plot) > 0:
        for exp in T2_rec:
            match exp:
                case "ramsey":
                    for q in T2_rec["ramsey"]:
                        sorted_item_ans = sorted(T2_rec["ramsey"][q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                        sorted_item_detu = sorted(detu_rec[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                        sorted_item_raw = sorted(T2_raw["ramsey"][q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                        earliest_time = datetime.strptime(sorted_item_ans[0][0], "%Y-%m-%d %H:%M:%S")
                        
                        time_diffs = []
                        sorted_values_ans, sorted_values_detu, sorted_values_raw = [], [], []
                        for idx, item in enumerate(sorted_item_ans):
                            key = item[0]
                            value = item[1]
                            current_time = datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
                            time_diff = round((current_time - earliest_time).total_seconds()/60,1)
                            time_diffs.append(time_diff)
                            sorted_values_ans.append(value)
                            sorted_values_detu.append(sorted_item_detu[idx][1])
                            sorted_values_raw.append(sorted_item_raw[idx][1])

                        
                        colormap(array(time_diffs),array(T2_evo_time[q])*1e6,array(sorted_values_raw),array(sorted_values_ans),fig_path=os.path.join(pic_folder,f"{q}_Ramsey_timeDep_colormap.png"))
                        plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "t2*", units={"x":"min","y":"µs"}, fig_path=os.path.join(pic_folder,f"{q}_Ramsey_timeDep.png"))
                        plot_timeDepCohe(array(time_diffs), array(sorted_values_detu)-array(sorted_values_detu)[0], "δf", units={"x":"min","y":"MHz"}, fig_path=os.path.join(pic_folder,f"{q}_Detune_timeDep.png"))
                case "spinecho":
                   for q in T2_rec["spinecho"]:
                        sorted_item_ans = sorted(T2_rec["spinecho"][q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                        sorted_item_raw = sorted(T2_raw["spinecho"][q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                        earliest_time = datetime.strptime(sorted_item_ans[0][0], "%Y-%m-%d %H:%M:%S")
                        
                        time_diffs = []
                        sorted_values_ans, sorted_values_raw = [], []
                        for idx, item in enumerate(sorted_item_ans):
                            key = item[0]
                            value = item[1]
                            current_time = datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
                            time_diff = round((current_time - earliest_time).total_seconds()/60,1)
                            time_diffs.append(time_diff)
                            sorted_values_ans.append(value)
                            sorted_values_raw.append(sorted_item_raw[idx][1])

                        colormap(array(time_diffs),array(T2_evo_time[q])*1e6,array(sorted_values_raw),array(sorted_values_ans),fig_path=os.path.join(pic_folder,f"{q}_SpinEcho_timeDep_colormap.png"))
                        plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "t2", units={"x":"min","y":"µs"}, fig_path=os.path.join(pic_folder,f"{q}_SpinEcho_timeDep.png"))
                
                case "T2cpmg":
                    for pi_num in T2_rec["T2cpmg"]:
                        for q in T2_rec["T2cpmg"][pi_num]:
                            sorted_item_ans = sorted(T2_rec["T2cpmg"][pi_num][q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                            sorted_item_raw = sorted(T2_raw["T2cpmg"][pi_num][q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
                            earliest_time = datetime.strptime(sorted_item_ans[0][0], "%Y-%m-%d %H:%M:%S")
                            
                            time_diffs = []
                            sorted_values_ans, sorted_values_raw = [], []
                            for idx, item in enumerate(sorted_item_ans):
                                key = item[0]
                                value = item[1]
                                current_time = datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
                                time_diff = round((current_time - earliest_time).total_seconds()/60,1)
                                time_diffs.append(time_diff)
                                sorted_values_ans.append(value)
                                sorted_values_raw.append(sorted_item_raw[idx][1])

                            colormap(array(time_diffs),array(T2_evo_time[q])*1e6,array(sorted_values_raw),array(sorted_values_ans),fig_path=os.path.join(pic_folder,f"{q}_CPMG_{pi_num}pi_timeDep_colormap.png"))
                            plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "cpmg", units={"x":"min","y":"µs"}, fig_path=os.path.join(pic_folder,f"{q}_CPMG_{pi_num}pi_timeDep.png"))

    eyeson_print(f"\nProcedures done ! ")


if __name__ == "__main__":
    QD_path = "qblox_drive_AS/QD_backup/20250313/DR1#11_SumInfo.pkl"
    folder_path = "qblox_drive_AS/Meas_raw/20250313/H17M03S52"
    save_every_fit_fig:bool = False

    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    time_monitor_data_ana(QD_agent, folder_path, save_every_fit_fig)