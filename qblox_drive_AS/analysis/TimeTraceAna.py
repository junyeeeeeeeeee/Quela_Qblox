import os, sys, json 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import open_dataset
from qblox_drive_AS.support.QDmanager import QDmanager,Data_manager
import matplotlib.pyplot as plt
from numpy import ndarray, array, median, std, mean
from datetime import datetime 
from qblox_drive_AS.support.UserFriend import *
from matplotlib.gridspec import GridSpec as GS
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer

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
    plt.colorbar(c, ax=ax, label='Contrast (V)')
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

    plt.tight_layout(pad=5.0)
    if fig_path is not None:
        plt.savefig(fig_path)
    plt.close()



if __name__ == "__main__":
    folder_paths = "Modularize/Meas_raw/20241107/TimeMonitor_MultiQ_H16M32S52"
    QD_file_path = 'Modularize/QD_backup/20241107/DR2#10_SumInfo.pkl'
    save_every_fit_pic:bool=False
    QD_agent = QDmanager(QD_file_path)
    QD_agent.QD_loader()

    files = [name for name in os.listdir(folder_paths) if (os.path.isfile(os.path.join(folder_paths,name)) and name.split(".")[-1] == "nc")]

    T1, T2, SS = 0, 0 , 0
    T1_rec, detu_rec, T2_rec, SS_rec = {}, {}, {}, {}
    T1_raw, T2_raw = {}, {}
    T1_evo_time, T2_evo_time = {}, {}
    for idx, file in enumerate(files) :
        slightly_print(f"Analysis for the {idx}-th files ...")
        exp_type:str = file.split("_")[1].split("(")[0]
        path = os.path.join(folder_paths,file)
        ds = open_dataset(path)
        match exp_type.lower():
            case "t1":
                T1 += 1
                for var in [ var for var in ds.data_vars if var.split("_")[-1] != 'x']:
                    T1_picsave_folder = os.path.join(folder_paths,f"{var}_T1_pics")
                    if T1 == 1:
                        T1_rec[var], T1_raw[var] = {}, {}
                        if save_every_fit_pic:
                            if not os.path.exists(T1_picsave_folder):
                                os.mkdir(T1_picsave_folder)
                    time_data = array(ds[f"{var}_x"])[0][0]
                    T1_evo_time[var] = time_data
                    ds[var].attrs['end_time'] = ds.attrs["end_time"]
                    ANA = Multiplex_analyzer("m13")
                    ANA._import_data(ds[var],var_dimension=2,refIQ=QD_agent.refIQ[var])
                    ANA._import_2nddata(time_data)
                    ANA._start_analysis()
                    if save_every_fit_pic:
                        ANA._export_result(T1_picsave_folder)
                    T1_raw[var][ds.attrs["end_time"]] = ANA.plot_item["data"]
                    T1_rec[var][ds.attrs["end_time"]] = ANA.fit_packs["median_T1"]
            case "t2":
                T2 += 1
                for var in [ var for var in ds.data_vars if var.split("_")[-1] != 'x']:
                    T2_picsave_folder = os.path.join(folder_paths,f"{var}_T2_pics")
                    if T2 == 1:
                        T2_rec[var], detu_rec[var], T2_raw[var] = {}, {}, {}
                        if save_every_fit_pic:
                            if not os.path.exists(T2_picsave_folder):
                                os.mkdir(T2_picsave_folder)
                    time_data = array(ds[f"{var}_x"])[0][0]
                    T2_evo_time[var] = time_data
                    ds[var].attrs['end_time'] = ds.attrs["end_time"]
                    ANA = Multiplex_analyzer("m12")
                    ANA._import_data(ds[var],var_dimension=2,refIQ=QD_agent.refIQ[var])
                    ANA._import_2nddata(time_data)
                    ANA._start_analysis()
                    if save_every_fit_pic:
                        ANA._export_result(T2_picsave_folder)
                    T2_raw[var][ds.attrs["end_time"]] = ANA.data_n
                    detu_rec[var][ds.attrs["end_time"]] = ANA.fit_packs["freq"]*1e-6
                    T2_rec[var][ds.attrs["end_time"]] = ANA.fit_packs["median_T2"]
            case "singleshot":
                SS += 1
                for var in ds.data_vars:
                    SS_picsave_folder = os.path.join(folder_paths,f"{var}_SingleShot_pics")
                    if SS == 1:
                        SS_rec[var] = {}
                        if save_every_fit_pic:
                            if not os.path.exists(SS_picsave_folder):
                                os.mkdir(SS_picsave_folder)
                    ANA = Multiplex_analyzer("m14")
                    ANA._import_data(ds[var]*1000,var_dimension=0,fq_Hz=QD_agent.quantum_device.get_element(var).clock_freqs.f01())
                    ANA._start_analysis()
                    if save_every_fit_pic:
                        pic_path = os.path.join(SS_picsave_folder,f"{var}_SingleShot_{ds.attrs['end_time'].replace(' ', '_')}")
                        ANA._export_result(pic_path)
                    SS_rec[var][ds.attrs["end_time"]] = ANA.fit_packs["effT_mK"]

    slightly_print(f"\nPlotting... ")
    if T1 != 0:
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

            
            colormap(array(time_diffs),array(T1_evo_time[q])*1e6,array(sorted_values_raw),array(sorted_values_ans),fig_path=os.path.join(folder_paths,f"{q}_T1_timeDep_colormap.png"))
            plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "t1", units={"x":"min","y":"µs"}, fig_path=os.path.join(folder_paths,f"{q}_T1_timeDep.png"))
    if T2 != 0:
        for q in T2_rec:
            sorted_item_ans = sorted(T2_rec[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
            sorted_item_detu = sorted(detu_rec[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
            sorted_item_raw = sorted(T2_raw[q].items(), key=lambda item: datetime.strptime(item[0], "%Y-%m-%d %H:%M:%S"))
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

            
            colormap(array(time_diffs),array(T2_evo_time[q])*1e6,array(sorted_values_raw),array(sorted_values_ans),fig_path=os.path.join(folder_paths,f"{q}_T2_timeDep_colormap.png"))
            plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "t2", units={"x":"min","y":"µs"}, fig_path=os.path.join(folder_paths,f"{q}_T2_timeDep.png"))
            plot_timeDepCohe(array(time_diffs), array(sorted_values_detu)-array(sorted_values_detu)[0], "δf", units={"x":"min","y":"MHz"}, fig_path=os.path.join(folder_paths,f"{q}_Detune_timeDep.png"))
    if SS != 0:
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


            plot_timeDepCohe(array(time_diffs), array(sorted_values_ans), "eff_Temp.", units={"x":"min","y":"mK"}, fig_path=os.path.join(folder_paths,f"{q}_effT_timeDep.png"))
    
    eyeson_print(f"\nProcedures done ! ")