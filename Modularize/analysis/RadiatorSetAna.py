import os, sys, time, json
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis, T2_fit_analysis, plot_textbox
from xarray import Dataset, open_dataset # Dataset.to_dict(SS_ds)
from numpy import array, ndarray, mean, std, round, arange
# from Modularize.support.Pulse_schedule_library import hist_plot
import matplotlib.pyplot as plt
from Modularize.support.Path_Book import meas_raw_dir
import re

def hist_plot(q:str,data:dict,title:str, save_path:str='', show:bool=True):
    fig, ax = plt.subplots(nrows =1,figsize =(2.5,2),dpi =250) 
    m, bins, patches = ax.hist(array(data[q]))
    ax.axvline(mean(array(data[q])), color = "k", ls = "--",lw=1)
    ax.set_xlabel(title)
    ax.set_ylabel('Counts')
    fig.tight_layout()
    if save_path != '':
        fig.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close()

def create_results_dict():
    exp_type = ["T1", "T2", "eff_T"]
    recor_cata =["avg", "std"]

    info_dict = {}
    for exp in exp_type:
        info_dict[exp] = {}
        for reco_item in recor_cata:
            info_dict[exp][reco_item]=0 # initialize
    
    return info_dict

def creat_weired_pic_folder(temperature_folder):
    weired_pic_path = os.path.join(temperature_folder,"weired_data")
    if not os.path.isdir(weired_pic_path):
        os.mkdir(weired_pic_path)
    return weired_pic_path

def create_result_folder(temperature_folder):
    results_path = os.path.join(temperature_folder,"results")
    if not os.path.isdir(results_path):
        os.mkdir(results_path)
    return results_path

def create_json_folder(results_folder):
    json_path = os.path.join(results_folder,"jsons")
    if not os.path.isdir(json_path):
        os.mkdir(json_path)
    return json_path

def save_weired_data_pic(x_s:ndarray,y:ndarray,exp_type:str,exp_idx:str,set_idx:str,T_folder_path:str, fit_res:dict={}):
    save_folder = creat_weired_pic_folder(T_folder_path)
    if fit_res == {}:
        plt.plot(round(x_s*1e6,1),y)
        plt.ylabel("Contrast (V)")
        plt.xlabel("Free Time (µs)")
        name = f"{exp_type} S{set_idx}-{exp_idx}"
        plt.title(f"{name} weired can't fit ")
        save_place = os.path.join(save_folder,f"{name}.png")
        plt.savefig(save_place)
        plt.close()
    else:
        fig, ax = plt.subplots(nrows =1,figsize =(6,4),dpi =250)
        text_msg = "Fit results\n"
        if fit_res.attrs['exper'] == 'T1':  
            title= 'T1 relaxation'
            x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
            x= fit_res.coords['freeDu']*1e6
            x_fit= fit_res.coords['para_fit']*1e6
            text_msg += r"$T_{1}= %.3f $"%(fit_res.attrs['T1_fit']*1e6) +r"$\ [\mu$s]"+'\n'
        elif fit_res.attrs['exper'] == 'T2':  
            title= 'Ramsey'
            x_label= r"$t_{f}$"+r"$\ [\mu$s]" 
            x= fit_res.coords['freeDu']*1e6
            x_fit= fit_res.coords['para_fit']*1e6
            text_msg += r"$T_{2}= %.3f $"%(fit_res.attrs['T2_fit']*1e6) +r"$\ [\mu$s]"+'\n'        
            text_msg += r"$detuning= %.3f $"%(fit_res.attrs['f']*1e-6) +' MHz\n'
        else:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            pass

        ax.plot(x,fit_res.data_vars['data']/(1/1000),'-', color="blue",label=r"$data$", alpha=0.5, ms=4)
        ax.plot(x_fit,fit_res.data_vars['fitting']/(1/1000),'-', color="red",label=r"$fit$", alpha=0.8, lw=1)       
        ax.set_xlabel(x_label)
        ax.set_title(title)
        ax.set_ylabel('Contrast [mV]')
        plot_textbox(ax,text_msg)
        fig.tight_layout()
        name = f"{exp_type} S{set_idx}-{exp_idx}-fit-O"
        plt.title(f"{name} weired")
        save_place = os.path.join(save_folder,f"{name}.png")
        plt.savefig(save_place)
        plt.close()

def sort_set(name_list:list, by_which_num_idx:str):
    name_list.sort(key=lambda l: int(re.findall('\d+',l)[by_which_num_idx]))
    

def sort_files(file_name_list:list):
    T1_file, T2_file, SS_file = [], [], []
    for file_name in file_name_list:
        if file_name.split("_")[1].split("(")[0] == 'T1':
            T1_file.append(file_name)
        elif file_name.split("_")[1].split("(")[0] == 'T2':
            T2_file.append(file_name)
        elif file_name.split("_")[1].split("(")[0] == 'SingleShot':
            SS_file.append(file_name)
    
    sort_set(T1_file,3)
    sort_set(T2_file,3)
    sort_set(SS_file,2)
    
    return T1_file+T2_file+SS_file


def plot_time_behavior(json_files:list, temperature_folder:str, time_axis:ndarray=array([])):
    radiator_temp = temperature_folder.split("/")[-1] if temperature_folder.split("/")[-1][-1].lower() == "k" else temperature_folder.split("/")[-2]
    a_set_time = 7 # min
    avg_t1, std_t1 = [], []
    avg_t2, std_t2 = [], []
    avg_eff_T, std_eff_T = [], []
    set_n = 1
    for a_json_file in json_files:
        with open(a_json_file) as J:
            info_dict = json.load(J)
            avg_t1.append(float(info_dict["T1"]["avg"]))
            std_t1.append(float(info_dict["T1"]["std"]))
            avg_t2.append(float(info_dict["T2"]["avg"]))
            std_t2.append(float(info_dict["T2"]["std"]))
            avg_eff_T.append(float(info_dict["eff_T"]["avg"]))
            std_eff_T.append(float(info_dict["eff_T"]["std"]))
        set_n += 1
    print(avg_t1)
    time_axis_min = arange(1,set_n)*a_set_time if time_axis.shape[0]==0 else round(time_axis/60,1)
    fig, ax = plt.subplots()
    ax.errorbar(time_axis_min,avg_t1,yerr=std_t1,fmt="o-",color='red',label='T1')
    ax.errorbar(time_axis_min,avg_t2,yerr=std_t2,fmt="o-",color='blue',label='T2')
    ax.set_ylabel("Time (µs)")
    ax.set_xlabel("Time after radiator ON (min)")
    ax.set_ylim(0,100)
    ax.legend(loc='upper left')
    axT = ax.twinx()
    axT.errorbar(time_axis_min,avg_eff_T,yerr=std_eff_T,fmt="o-",color='orange',label='eff_T')
    axT.spines['right'].set_color("orange")
    axT.yaxis.label.set_color('orange')
    axT.set_ylabel("Effective Temp. (mK)")
    axT.set_ylim(0,10)
    axT.legend(loc='upper right')
    plt.title(f"Behavior after radiator ON (T={radiator_temp}) by Qblox")
    behavior_path = os.path.join(temperature_folder,f"T_{radiator_temp}_behavior.png")
    plt.savefig(behavior_path)
    plt.close()
    

# SS_ds = open_dataset("Modularize/Meas_raw/10K/Radiator(0)/DR1q0_SingleShot(1)_H16M56S29.nc")
# ss_dict = Dataset.to_dict(SS_ds)
# print(ss_dict['data_vars']['e']['data'][0])

other_info_dict = {}
target_q = 'q0'
temperature = '10K'
parent_path = os.path.join(meas_raw_dir,temperature)

sub_folders = [name for name in os.listdir(parent_path) if (os.path.isdir(os.path.join(parent_path,name)) and name[:8]=='Radiator')] # Radiator(idx)
other_info_ = [name for name in os.listdir(parent_path) if (os.path.isfile(os.path.join(parent_path,name)) and name.split(".")[-1]=='json')]
print(other_info_)
with open(os.path.join(parent_path,other_info_[0])) as JJ:
    other_info_dict = json.load(JJ)

sort_set(sub_folders,0) # sort sets start from 0 
ref_iq = other_info_dict[target_q]["refIQ"]
transi_freq = other_info_dict[target_q]["f01"] # Hz
time_past_sec_array = array(other_info_dict[target_q]["time_past"]) # shold be the same length with sets, units in second

start = time.time()
for folder_name in sub_folders:
    info_dict = create_results_dict()
    set_idx = folder_name.split("(")[-1].split(")")[0]
    folder_path = os.path.join(parent_path,folder_name)
    print(f"==================================================== Set-{set_idx} start")
    files = [name for name in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path,name))] # DR1q0_{T1/T2/SingleShot}(exp_idx)_H17M23S19.nc
    files = sort_files(files) # sort files start from 0 
    T1_us = []
    T2_us = []
    effT_mK = []

    pgI_collection = []
    pgQ_collection = []
    peI_collection = []
    peQ_collection = []
    for file_name in files: # in a single set
        exp_idx = file_name.split("(")[-1].split(")")[0]  # histo_counts
        exp_type = file_name.split("(")[0].split("_")[-1] # T1/ T2/ SingleShot
        file_path = os.path.join(folder_path,file_name)
        print(f"{exp_type}-{exp_idx}")
        if exp_type == "T1":
            T1_ds = open_dataset(file_path)
            
            times = array(Dataset.to_dict(T1_ds)["coords"]['x0']['data']) # s
            I,Q= dataset_to_array(dataset=T1_ds,dims=1)
            data= array(IQ_data_dis(I,Q,ref_I=ref_iq[0],ref_Q=ref_iq[-1]))
            try:
                data_fit= T1_fit_analysis(data=data,freeDu=times,T1_guess=8e-6)
                T1_us.append(data_fit.attrs['T1_fit']*1e6)
                if data_fit.attrs['T1_fit']*1e6 > 50: 
                    save_weired_data_pic(times, data, "T1", exp_idx, set_idx, parent_path,data_fit)   
            except:
                save_weired_data_pic(times, data, "T1", exp_idx, set_idx, parent_path)
                T1_us.append(0)
        elif exp_type == "T2":
            T2_ds = open_dataset(file_path)
            times = array(Dataset.to_dict(T2_ds)["coords"]['x0']['data']) # s
            I,Q= dataset_to_array(dataset=T2_ds,dims=1)
            data= (IQ_data_dis(I,Q,ref_I=ref_iq[0],ref_Q=ref_iq[1]))
            try:
                data_fit= T2_fit_analysis(data=data,freeDu=times,T2_guess=8e-6)
                T2_us.append(data_fit.attrs['T2_fit']*1e6)
                if data_fit.attrs['T2_fit']*1e6 > 30:
                    save_weired_data_pic(times, data, "T2", exp_idx, set_idx, parent_path, data_fit)
            except:
                save_weired_data_pic(times, data, "T2", exp_idx, set_idx, parent_path)
                T2_us.append(0)
            
            
        elif exp_type == "SingleShot":
            # collect data to choose training and predict
            
            SS_ds = open_dataset(file_path)
            ss_dict = Dataset.to_dict(SS_ds)
            # print(ss_dict)
            pe_I, pe_Q = ss_dict['data_vars']['e']['data']
            pg_I, pg_Q = ss_dict['data_vars']['g']['data']
            pgI_collection.append(pg_I)
            pgQ_collection.append(pg_Q)
            peI_collection.append(pe_I)
            peQ_collection.append(pe_Q)
        else:
            pass

        # TODO: train GMM
        # TODO: predict all collection to calculate eff_T for every exp_idx

    print(f"First stage analysis complete for set-{set_idx}!")

    result_folder = create_result_folder(parent_path)
    json_folder = create_json_folder(result_folder)

    # calc T1
    T1_us = array(T1_us)
    mean_T1_us = round(mean(T1_us[T1_us != 0]),1)
    sd_T1_us = round(std(T1_us[T1_us != 0]),1)
    info_dict["T1"]["avg"], info_dict["T1"]["std"] = mean_T1_us, sd_T1_us
    histo_path = os.path.join(result_folder,f"T1-histo-S{set_idx}.png")
    hist_plot("nobu",{"nobu":T1_us},f"S{set_idx}, T1={mean_T1_us}+/-{sd_T1_us} us",histo_path, False)
    # calc T2
    T2_us = array(T2_us)
    mean_T2_us = round(mean(T2_us[T2_us != 0]),1)
    sd_T2_us = round(std(T2_us[T2_us != 0]),1)
    info_dict["T2"]["avg"], info_dict["T2"]["std"] = mean_T2_us, sd_T2_us
    histo_path = os.path.join(result_folder,f"T2-histo-S{set_idx}.png")
    hist_plot("nobu",{"nobu":T2_us[T2_us != 0]},f"S{set_idx}, T2={mean_T2_us}+/-{sd_T2_us} us",histo_path, False)
    
    # calc OnsShot

    # save the info to plt scatter
    with open(f"{json_folder}/setInfo({set_idx}).json", "w") as record_file:
        json.dump(info_dict,record_file)

end = time.time()
print(f"Analysis time cost: {round((end-start)/60,1)} mins")



""" plot behavior """
temp = '10K'
parent_folder = os.path.join(meas_raw_dir,temp)
jsons_folder = os.path.join(parent_folder,"results/jsons")
static_info = [name for name in os.listdir(jsons_folder) if (os.path.isfile(os.path.join(jsons_folder,name)) and name.split(".")[-1]=='json')]
sort_set(static_info,0)

j_paths = []
print(static_info)
for a_json in static_info:
    j_paths.append(os.path.join(jsons_folder,a_json))


plot_time_behavior(j_paths,parent_folder,time_past_sec_array)
