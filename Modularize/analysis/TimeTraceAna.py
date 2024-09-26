import os, sys, json 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import Dataset, open_dataset
from Modularize.support import QDmanager,Data_manager
import matplotlib.pyplot as plt
from numpy import ndarray, array, median, std
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T2_fit_analysis, T1_fit_analysis
from datetime import datetime 
from matplotlib.gridspec import GridSpec as GS
from Modularize.analysis.Radiator.RadiatorSetAna import sort_set

def time_label_sort(nc_file_name:str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0],"H%HM%MS%S")

def plot_coherence_timetrace(raw_data:ndarray, time_samples:ndarray, ans:ndarray, q:str, raw_data_folder:str, exp:str, detunings:ndarray=[]):
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
        print("**********************************Yes")
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


if __name__ == "__main__":
    folder_paths = {"T1_folder_path":"Modularize/Meas_raw/T1_timeDep",
                    "T2_folder_path":"Modularize/Meas_raw/T2_timeDep",
                    "OS_folder_path":""
                    }
    QD_file_path = 'Modularize/QD_backup/2024_9_25/DR4#81_SumInfo.pkl'
    qs = ['q4']
    sort_mode = 'idx' # 'idx' or 'time'

    QD_agent = QDmanager(QD_file_path)
    QD_agent.QD_loader()

    for folder_name in folder_paths: 
        folder = folder_paths[folder_name]
        if folder_paths[folder_name] != '':
            if sort_mode == 'time':
                files = sorted([name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split(".")[-1] == "nc")],key=lambda name:time_label_sort(name))
            elif sort_mode == 'idx':
                files = [name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split(".")[-1] == "nc")]
                sort_set(files,3)
            else:
                raise KeyError(f"Unsupported sort mode was given = '{sort_mode}'")
            raw_data = []
            ans = []
            detu = []
            for idx, file in enumerate(files) :
                path = os.path.join(folder,file)
                nc = open_dataset(path)
                if folder_name.split("_")[0] in ["T1", "T2"]:
                    samples = array(nc.x0)
                    I,Q= dataset_to_array(dataset=nc,dims=1)
                    data = IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[qs[0]][0],ref_Q=QD_agent.refIQ[qs[0]][-1])
                    raw_data.append(data)
                    if folder_name.split("_")[0] == "T1":   
                        data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=25e-6)
                        ans.append(data_fit.attrs['T1_fit']*1e6)
                    else:
                        data_fit= T2_fit_analysis(data=data,freeDu=samples,T2_guess=10e-6)
                        ans.append(data_fit.attrs['T2_fit']*1e6)
                        detu.append(data_fit.attrs['f']*1e-6)
                else:
                    #　OS to build
                    pass
                    

            Data_manager().save_histo_pic(QD_agent,{str(qs[0]):ans},qs[0],mode=folder_name.split("_")[0])
            
            plot_coherence_timetrace(array(raw_data),samples,array(ans),qs[0],folder,folder_name.split("_")[0],array(detu))
