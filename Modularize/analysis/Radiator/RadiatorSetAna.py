"""
This program is only for analyzing a series of radiation tests like 0K, 20K 40K 60K and re0K with/without shielding. This didn't support comparison between with and without shielding
"""
import os, sys, time, json, pickle
from pandas import DataFrame as DF
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis, T2_fit_analysis, plot_textbox, Fit_analysis_plot
from xarray import Dataset, open_dataset # Dataset.to_dict(SS_ds)
from numpy import array, ndarray, mean, std, round, arange, moveaxis, any, zeros, delete, average, sqrt
# from Modularize.support.Pulse_schedule_library import hist_plot
import matplotlib.pyplot as plt
from xarray import DataArray
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.analysis.DRtemp import Kelvin_collector
from qcat.analysis.state_discrimination.discriminator import train_GMModel  # type: ignore
from qcat.visualization.readout_fidelity import plot_readout_fidelity

exp_items = {"1":"T1","2":"T2","3":"effT","4":"gamma1","5":"gamma2","6":"thermalPop","7":"gammaPhi"}
units = {"T1":"us","T2":"us","effT":"mK","gamma1":"MHz","gamma2":"MHz","thermalPop":"%","gammaPhi":"MHz"}
# ============================= figure settings ================================================

fig_size = (15,10)
label_font_size = 26
tick_num_size = 26
legend_font_size = 30

def ret_error_bound(errors:list,threshold:float=1)->tuple[float, float]:
    all_err = array(errors).reshape(-1)
    err_mean = mean(all_err)
    up_lim = err_mean + threshold*std(all_err)
    bt_lim = err_mean - threshold*std(all_err)
    return up_lim, bt_lim

def ax_set_y_label(ax:plt.Axes,exp:str, font_size:int=26)->plt.Axes:
    if exp.lower() == 't1':
        ax.set_ylabel("$T_{1}$ "+f"({units[exp]})",fontsize=font_size)
    #---------------------------------------------------------------------------------------------------------------
    elif exp.lower() == 'efft':
        ax.set_ylabel("Effective Temp. "+f"({units[exp]})",fontsize=font_size)
    #---------------------------------------------------------------------------------------------------------------
    elif exp.lower() == 't2':
        ax.set_ylabel("$T_{2}$ "+f"({units[exp]})",fontsize=font_size)
    #---------------------------------------------------------------------------------------------------------------
    elif exp.lower() == "gamma1":
        ax.set_ylabel("$\Gamma_{1}$ "+f"({units[exp]})",fontsize=font_size)
    #---------------------------------------------------------------------------------------------------------------
    elif exp.lower() == "gamma2":
        ax.set_ylabel("$\Gamma_{2}$ "+f"({units[exp]})",fontsize=font_size)
    #---------------------------------------------------------------------------------------------------------------
    elif exp.lower() == "gammaphi":
        ax.set_ylabel("$\Gamma_{\phi}$ "+f"({units[exp]})",fontsize=font_size)
    #---------------------------------------------------------------------------------------------------------------
    else:
        ax.set_ylabel("Thermal populations "+f"({units[exp]})",fontsize=font_size)
    return ax

# ================================ Functional =========================================================================================
def exp_item_translator(options:list)->list:
    x = []
    for option_label in options:
        x.append(exp_items[str(option_label)])
    return x

def get_ref_from_json(target_q:str,sample_folder_name:str, condition_name:str)->dict:
    """
    mode is "before" | "recover", \n
    ### Return 
    {"ref_before":{"T1":35.9,"T1_sd":2,...}, "ref_recove":{...}}
    """
    ref_path = os.path.join(meas_raw_dir,sample_folder_name,f"{target_q}_references.json")
    old_ref:dict = {}
    if os.path.isfile(ref_path):
        with open(ref_path) as record_file: 
            old_ref = json.load(record_file)
    else:
        old_ref[condition_name] = {"ref_before":{},"ref_recove":{}}
        print("Here is no reference rec !")
    
    return old_ref[condition_name]

def find_nearest(ary:ndarray, value:float)->tuple[int,float]:
    """ find the element  which is closest to the given target_value in the given array"""
    idx = (abs(ary - value)).argmin()
    return idx, float(ary[idx])

def delet_incomplete_set(temp_folder_path)->list:
    file_number = []
    set_folders = [os.path.join(temp_folder_path,name) for name in os.listdir(temp_folder_path) if (os.path.isdir(os.path.join(temp_folder_path,name)) and name[:8]=='Radiator')] # Radiator(idx)
    for folder in set_folders:
        file_number.append(len([name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split(".")[-1] == "nc")]))
    complete_file_number = max(set(file_number), key=file_number.count)
    print("***",file_number,":",complete_file_number)
    idx_to_delete = []
    for idx in range(len(file_number)):
        if file_number[idx] < complete_file_number:
            idx_to_delete.append(idx)
    
    complete_set_folder = delete(set_folders, idx_to_delete).tolist()
    return [os.path.split(path)[-1] for path in complete_set_folder]   # return the set folder names list

def save_ref(ref_dict:dict, target_q:str, sample_name:str="",conditional_name:str="", ref_type:str="before"):
        """ ref_type in ['before', 'recover']"""
    
        if conditional_name == "":
            from datetime import datetime as d
            ref_name = str(d.now().date())
        else:
            ref_name = conditional_name

        ref_path = os.path.join(meas_raw_dir,sample_name,f"{target_q}_references.json")
        if os.path.exists(ref_path):
            old_ref:dict = {}
            with open(ref_path) as record_file: 
                old_ref = json.load(record_file)
            if ref_name in list(old_ref.keys()):
                if ref_type.lower() == 'before':
                    if "ref_before" in list(old_ref[ref_name].keys()):
                        permission = input(f"There is an before reference inside: {old_ref[ref_name]['ref_before']}, replace it? [y/n]...")
                        if permission.lower() == 'y':
                            old_ref[ref_name]["ref_before"] = ref_dict
                            print("REF json refreshed!")
                        else:
                            print("Cancelled!")
                    else:
                        old_ref[ref_name]["ref_before"] = ref_dict
                elif  ref_type.lower() == 'recover':
                    if "ref_recove" in list(old_ref[ref_name].keys()):
                        permission = input(f"There is an recover reference inside: {old_ref[ref_name]['ref_recove']}, replace it? [y/n]...")
                        if permission.lower() == 'y':
                            old_ref[ref_name]["ref_recove"] = ref_dict
                            print("REF json refreshed!")
                        else:
                            print("Cancelled!")
                    else:
                        old_ref[ref_name]["ref_recove"] = ref_dict
                else:
                    raise KeyError("Unsupported ref_type was given!")
            else:
                dict_name = 'ref_before' if ref_type.lower() == 'before' else 'ref_recove'
                if ref_type.lower() not in ['before', 'recover']: raise KeyError("Unsupported ref_type was given!")
                old_ref[ref_name] = {}
                old_ref[ref_name][dict_name] = ref_dict
                    
            with open(ref_path,'w') as record_file: 
                json.dump(old_ref,record_file)
        else:
            dict_name = 'ref_before' if ref_type.lower() == 'before' else 'ref_recove'
            if ref_type.lower() not in ['before', 'recover']: raise KeyError("Unsupported ref_type was given!")
            new_ref = {str(ref_name):{dict_name:ref_dict}}
            with open(ref_path,'w') as record_file: 
                json.dump(new_ref,record_file)
    
def timelabelfolder_creator(folder_path:str,additional_folder_name:str='')->str:
    import datetime
    current_time = datetime.datetime.now()
    temp_dep_folder_name =  f"{additional_folder_name}_H{current_time.hour}M{current_time.minute}S{current_time.second}"
    temp_dep_folder_path = os.path.join(folder_path,temp_dep_folder_name)
    if not os.path.exists(temp_dep_folder_path):
        os.mkdir(temp_dep_folder_path)
    return temp_dep_folder_path

def pic_values_saver(folder_path:str,mode:str,*args)->str:
    """ 
        mode indicates the axis in this pic is time with the mode='time', or temperature with the mode='temp'.
    """
    exp_type = [data_dict["exp"] for data_dict in args]
    if "T1" in exp_type:
        title=["T1",units["T1"]]
    elif "T2" in exp_type:
        title=["T2",units["T2"]]
    elif "effT" in exp_type:
        title=["effT",units["effT"]]
    elif "gamma1" in exp_type:
        title=["gamma1",units["gamma1"]]
    elif "gamma2" in exp_type:
        title=["gamma2",units["gamma2"]]    
    elif "gammaPhi" in exp_type:
        title=["gammaPhi",units["gammaPhi"]]   
    elif "thermalPop" in exp_type:
        title=["thermalPop",units["thermalPop"]]   
    else:
        title=["EXP","a.u."]

    if mode.lower() in ["time"]:
        x_axis_unit = "min"
    elif mode.lower() in ["temp","temperature","tempera"]:
        x_axis_unit = "K"
    else:
        raise KeyError("Unsupported mode given!")

    to_csv_dict = {f"exp_x_({x_axis_unit})":[],f"exp_y_({title[-1]})":[],f"exp_yerr_({title[-1]})":[],
                   f"exp_dr_temp_x_({x_axis_unit})":[],"exp_dr_temp_(K)":[],"exp_dr_temp_err_(K)":[],
                   f"ref_x_(K)":[],f"ref_y_({title[-1]})":[],f"ref_yerr_({title[-1]})":[]
                   }
    for data_dict in args:
        if data_dict["exp"] in list(exp_items.values()):
            to_csv_dict[f"exp_x_({x_axis_unit})"] = data_dict["x"]
            to_csv_dict[f"exp_y_({title[-1]})"] += data_dict["y"]
            to_csv_dict[f"exp_yerr_({title[-1]})"] += data_dict["yerr"]

        elif data_dict["exp"].lower() in ["ref-t1", "ref-t2", "ref-efft", "ref-gamma1", "ref-gamma2", "ref-gammaphi", "ref-thermalpop"]:
            to_csv_dict[f"ref_x_(K)"] += data_dict["x"] # reference measured at 4K
            to_csv_dict[f"ref_y_({title[-1]})"] += data_dict["y"]
            to_csv_dict[f"ref_yerr_({title[-1]})"] += data_dict["yerr"]
        else:
            if data_dict["exp"].split("-")[-1] == '6':
                n = 1/1000
            else:
                n = 1
            to_csv_dict[f"exp_dr_temp_x_({x_axis_unit})"] = data_dict["x"]
            to_csv_dict["exp_dr_temp_(K)"] = array(data_dict["y"])*n
            to_csv_dict["exp_dr_temp_err_(K)"] = array(data_dict["yerr"])*n

    df = DF.from_dict(to_csv_dict, orient='index').transpose()
    DF.to_csv(df, os.path.join(folder_path,f"{title[0]}_dataValues.csv"))

    return folder_path
    
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
    exp_type = list(exp_items.values())
    recor_cata =["avg", "std"]

    info_dict = {}
    for exp in exp_type:
        info_dict[exp] = {}
        for reco_item in recor_cata:
            info_dict[exp][reco_item]=0 # initialize
    
    return info_dict

def create_T1T2_folder(results_folder:str, mode:str):
    if mode.lower() == 't1':
        folder = os.path.join(results_folder,"T1")
    else:
        folder = os.path.join(results_folder,"T2")
    if not os.path.isdir(folder):
        os.mkdir(folder)
    return folder

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

def creat_oneshot_folder(result_folder):
    results_path = os.path.join(result_folder,"SingleShot")
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
    import re
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

def OSdata_arranger(total_array:ndarray, want_IQset_num:int=1)->tuple[list, list]:
    """
    total_array shape: (2,2,M,N) which is (IQ, state, histos, shots)\n
    return traning and predict list
    """
    from numpy.random import randint
    all_datasets = []
    train_dataset  = []
    
    for histo in range(total_array.shape[2]):
        IQ_folder = []
        for iq in range(total_array.shape[0]):
            state_folder = []
            for state in range(total_array.shape[1]):
                state_folder.append(total_array[iq][state][histo])
            IQ_folder.append(state_folder)
        ds = DataArray(array(IQ_folder), coords= [("mixer",["I","Q"]), ("prepared_state",[0,1]), ("index",arange(array(IQ_folder).shape[2]))] )
        all_datasets.append(ds)
    # total_sets.shape shape = (histo, IQ, State, shots)

    for pick_number in range(want_IQset_num):
        rand_pick_set_idx = randint(0 ,len(all_datasets))
        train_dataset.append(all_datasets[rand_pick_set_idx])
    
    return train_dataset, all_datasets

def collect_allSets_inTempera(temperature_folder_path:str, refresh:bool=False)->dict:
    """ 
    ## Return \n
    info_dict = {"T1":{"avg","std"},"gamma1":{...},"T2":{...},"gamma2":{...},"gammaPhi":{...},"effT":{...},"thermalPop":{...}}
    """
    json_path = os.path.join(f"{temperature_folder_path}","results","jsons")
    info_dict = {}
    if not refresh:
        if os.path.exists(os.path.join(json_path,"temperatureInfo.json")):
            print("read old json")
            with open(os.path.join(json_path,"temperatureInfo.json")) as J:
                info_dict = json.load(J)
    if refresh or not os.path.exists(os.path.join(json_path,"temperatureInfo.json")):
        avg_t1, std_t1 = [], []
        avg_gamma1, std_gamma1 = [], []
        avg_t2, std_t2 = [], []
        avg_gamma2, std_gamma2 = [], []
        avg_eff_T, std_eff_T = [], []
        avg_thermalPop, std_thermalPop = [], []
        avg_gammaPhi, std_gammaPhi = [], []
        json_files = [name for name in os.listdir(json_path) if (os.path.isfile(os.path.join(json_path,name)) and name.split(".")[-1]=='json' and name.split("(")[0]=='setInfo')]
        sort_set(json_files,0)
        j_paths = []
        for a_json in json_files:
            j_paths.append(os.path.join(json_path,a_json))
        for a_json_file in j_paths:
            with open(a_json_file) as J:
                set_dict = json.load(J)
                avg_t1.append(float(set_dict["T1"]["avg"]))
                std_t1.append(float(set_dict["T1"]["std"]))
                avg_t2.append(float(set_dict["T2"]["avg"]))
                std_t2.append(float(set_dict["T2"]["std"]))
                avg_eff_T.append(float(set_dict["effT"]["avg"]))
                std_eff_T.append(float(set_dict["effT"]["std"]))
                try: # added after the first run exp at 2024/05/28
                    avg_gamma1.append(float(set_dict["gamma1"]["avg"]))
                    std_gamma1.append(float(set_dict["gamma1"]["std"]))
                    avg_gamma2.append(float(set_dict["gamma2"]["avg"]))
                    std_gamma2.append(float(set_dict["gamma2"]["std"]))
                    avg_gammaPhi.append(float(set_dict["gammaPhi"]["avg"]))
                    std_gammaPhi.append(float(set_dict["gammaPhi"]["std"]))
                    avg_thermalPop.append(float(set_dict["thermalPop"]["avg"]))
                    std_thermalPop.append(float(set_dict["thermalPop"]["std"]))
                except:
                    pass
        info_dict = {"T1":{"avg":avg_t1,"std":std_t1},"gamma1":{"avg":avg_gamma1,"std":std_gamma1},"T2":{"avg":avg_t2,"std":std_t2},"gamma2":{"avg":avg_gamma2,"std":std_gamma2},"effT":{"avg":avg_eff_T,"std":std_eff_T},"thermalPop":{"avg":avg_thermalPop,"std":std_thermalPop},"gammaPhi":{"avg":avg_gammaPhi,"std":std_gammaPhi}}  # contains the AVG and SG for every set
        with open(os.path.join(json_path,"temperatureInfo.json"), "w") as record_file:
            json.dump(info_dict,record_file)
    
    return info_dict

def get_time_axis(target_q:str, folder_path:str)->ndarray:
    other_info_dict = {}
    other_info_ = [name for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[0]=='otherInfo')]
    with open(os.path.join(folder_path,other_info_[0])) as JJ:
        other_info_dict = json.load(JJ)
    # extract every set time
    return array(other_info_dict[target_q]["time_past"]) # shold be the same length with sets, units in min

def filter_zero(exp_values:ndarray, time_axis:ndarray, exp_stds:ndarray=[])->tuple[ndarray, ndarray, ndarray]:
    del_idx = []
    for idx, value in enumerate(exp_values):
        if value == 0:
            del_idx.append(idx)
    
    new_time_axis = delete(array(time_axis), del_idx)
    new_exp_values = delete(array(exp_values), del_idx)
    if array(exp_stds).shape[0] != 0:
        new_std_values = delete(array(exp_stds), del_idx)
    else:
        new_std_values = array([])
    return new_time_axis, new_exp_values, new_std_values

# ============================================ Analysis ================================================================================
def get_references_from_ResultJson(temperature_folder_path:str,minute_from_the_end:int=0,target_q:str=""):
    json_file = os.path.join(temperature_folder_path,"results","jsons","every_values.json")
    if minute_from_the_end != 0 and target_q != "":
        time_past_min_array = get_time_axis(target_q,temperature_folder_path)/60
        slice_from_this_min = int(time_past_min_array[-1]-minute_from_the_end)
        set_number = time_past_min_array.shape[0]
    ref_dict = {}
    every_values_dict = {}
    with open(json_file) as JJ:
        every_values_dict = json.load(JJ)
    
    for keyname in every_values_dict: # "T1s", "T2s", ...
        values = array(every_values_dict[keyname])
        if minute_from_the_end != 0 and target_q != "":
            histo_count_in_set = int(values.shape[0]/set_number) # data number of a exp in every_value_dict = set_number * histo_count_in_set
            this_min_idx_in_set, _ = find_nearest(time_past_min_array, slice_from_this_min) # this index is equal to the set index
            histo_idx_this_min = (this_min_idx_in_set+1)*histo_count_in_set 
            values = values[histo_idx_this_min:]
        if keyname[:-1] in ["T1", "T2", "effT"]:
            if  keyname[:-1] == "effT":
                ref_dict[f"{keyname[:-1]}"] = round(mean(values[values != 0]),1)
            else:
                errors = array(every_values_dict[keyname[:-1]+"ers"])[histo_idx_this_min:]
                ref_dict[f"{keyname[:-1]}"] = round(average(values[values != 0], weights=1/(errors[errors!=0])),1)
            ref_dict[f"{keyname[:-1]}_sd"] = round(std(values[values != 0]),1)
        elif keyname[:-1] in ["gamma1", "gamma2", "gammaPhi"]:
            errors = array(every_values_dict[keyname[:-1]+"ers"])[histo_idx_this_min:]
            ref_dict[f"{keyname[:-1]}"] = round(average(values[values != 0], weights=1/(errors[errors!=0])),3)
            ref_dict[f"{keyname[:-1]}_sd"] = round(std(values[values != 0]),3)
        elif keyname[:-1] in ["thermalPop"]:
            ref_dict[f"{keyname[:-1]}"] = round(mean(values[values != 0]),2)
            ref_dict[f"{keyname[:-1]}_sd"] = round(std(values[values != 0]),2)
        else:
            pass
    
    return ref_dict

def a_set_analysis(set_folder_path:str, old_monitor_dict:dict, ref_iq:list, transi_freq:float)->dict:
    """ old_monitor_dict = {"x_minutes":[],"T1":[],"T1_sd":[],"T2":[],"T2_sd":[],"effT":[],"effT_sd":[],"thermalPop":[],"thermalPop_sd":[]} """
    if old_monitor_dict == {}:
        old_monitor_dict = {"x_minutes":[],"T1":[],"T1_sd":[],"T2":[],"T2_sd":[],"effT":[],"effT_sd":[],"thermalPop":[],"thermalPop_sd":[]}
    folder_name = os.path.split(set_folder_path)[-1]
    set_idx = folder_name.split("(")[-1].split(")")[0]
    folder_path = set_folder_path
    print(f"==================================================== Set-{set_idx} start")
    files = [name for name in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path,name))] # DR1q0_{T1/T2/SingleShot}(exp_idx)_H17M23S19.nc
    files = sort_files(files) # sort files start from 0 
    T1_us = []
    T2_us = []
    effT_mK = []
    thermal_pop = []

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
            except:
                T1_us.append(0)
        elif exp_type == "T2":
            T2_ds = open_dataset(file_path)
            times = array(Dataset.to_dict(T2_ds)["coords"]['x0']['data']) # s
            I,Q= dataset_to_array(dataset=T2_ds,dims=1)
            data= (IQ_data_dis(I,Q,ref_I=ref_iq[0],ref_Q=ref_iq[1]))
            try:
                data_fit= T2_fit_analysis(data=data,freeDu=times,T2_guess=8e-6)
                T2_us.append(data_fit.attrs['T2_fit']*1e6)
            except:
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
        
    
    # reshape data to (I,Q)*(g,e)*shots       
    OS_data = 1000*array([[pgI_collection,peI_collection],[pgQ_collection,peQ_collection]]) # can train or predict 2*2*histo_counts*shot
    tarin_data, fit_arrays = OSdata_arranger(OS_data)
    # train GMM
    dist_model = train_GMModel (tarin_data[0])
    dist_model.relabel_model(array([ref_iq]).transpose())
    # predict all collection to calculate eff_T for every exp_idx
    for histo_i in range(fit_arrays.shape[0]):
        analysis_data = fit_arrays[histo_i] #your (2,2,N) data to analysis

        new_data = moveaxis( analysis_data ,1,0)
        p0_pop = dist_model.get_state_population(new_data[0].transpose()) # [p0in0, p0in1]
        p1_pop = dist_model.get_state_population(new_data[1].transpose()) # [p1in0, p1in1]
        if transi_freq != 0 and transi_freq is not None:
            fig , eff_t, snr = plot_readout_fidelity(analysis_data, transi_freq, None, False)
            effT_mK.append(eff_t)
        else:
            thermal_pop.append(p0_pop[1])
        
        plt.close()
    
    # calc T1
    T1_us = array(T1_us)
    old_monitor_dict["T1"].append(round(mean(T1_us[T1_us != 0]),1))
    old_monitor_dict["T1_sd"].append(round(std(T1_us[T1_us != 0]),1))

    # calc T2
    T2_us = array(T2_us)
    old_monitor_dict["T2"].append(round(mean(T2_us[T2_us != 0]),1))
    old_monitor_dict["T2_sd"].append(round(std(T2_us[T2_us != 0]),1))

    # calc OnsShot
    if len(effT_mK) != 0:
        effT_mK = array(effT_mK)
        old_monitor_dict["effT"].append(round(mean(effT_mK),1))
        old_monitor_dict["effT_sd"].append(round(std(effT_mK),1))
    else:
        thermal_pop = array(thermal_pop)
        old_monitor_dict["thermalPop"].append(round(mean(thermal_pop),3))
        old_monitor_dict["thermalPop"].append(round(std(thermal_pop),3))

    return old_monitor_dict
    
def main_analysis(target_q:str, temperature_folder_path:str, mode:str='all' ):
    """
    target_q: 'q0'\n
    temperature: '10K'\n
    mode: 'all' or 'jump', 'all' for analyze all set again. 'jump' analyze those haven't been analyzed set.\n
    """
    other_info_dict = {}
    temperature = os.path.split(temperature_folder_path)[-1]
    parent_path = temperature_folder_path

    sub_folders = delet_incomplete_set(temperature_folder_path)
    print(sub_folders)
    other_info_ = [name for name in os.listdir(temperature_folder_path) if (os.path.isfile(os.path.join(temperature_folder_path,name)) and name.split(".")[0]=='otherInfo')]
    with open(os.path.join(temperature_folder_path,other_info_[0])) as JJ:
        other_info_dict = json.load(JJ)

    sort_set(sub_folders,0) # sort sets start from 0 
    ref_iq = other_info_dict[target_q]["refIQ"]
    transi_freq = other_info_dict[target_q]["f01"] # Hz
    result_folder = create_result_folder(temperature_folder_path)
    json_folder = create_json_folder(result_folder)
    SS_folder = creat_oneshot_folder(result_folder)
    T1_pic_folder = create_T1T2_folder(result_folder,"t1")
    T2_pic_folder = create_T1T2_folder(result_folder,"t2")
    start = time.time()
    T1s = {}
    T1ers = {}
    gamma1s = {}
    gamma1ers = {}
    T2s = {}
    T2ers = {}
    gamma2ers = {}
    gamma2s = {}
    effTs = {}
    thermalpops = {}
    gammaphis = {}
    gammaphiers = {}

    if mode.lower() == 'jump':
        set_idx_done = [int(name.split("(")[-1].split(")")[0]) for name in os.listdir(json_folder) if (os.path.isfile(os.path.join(json_folder,name)) and name.split(".")[0].split("(")[0]=='setInfo')]
        final_complete_setidx = max(set_idx_done)
        print(f"start set-idx = {final_complete_setidx+1}")
        sub_folders = sub_folders[final_complete_setidx+1:]

    for folder_name in sub_folders:
        set_idx = folder_name.split("(")[-1].split(")")[0]
        folder_path = os.path.join(temperature_folder_path,folder_name)
        print(f"==================================================== Set-{set_idx} start")
        files = [name for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[-1] == "nc")] # DR1q0_{T1/T2/SingleShot}(exp_idx)_H17M23S19.nc
        files = sort_files(files) # sort files start from 0 
        T1_us = []
        T1_err = []
        gamma1_MHz = []
        gamm1_err = []
        T2_us = []
        T2_err = []
        gamma2_MHz = []
        gamma2_err = []
        effT_mK = []
        therm_pop = []

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
                    data_fit, fit_error = T1_fit_analysis(data=data,freeDu=times,T1_guess=8e-6, return_error=True)
                    t1 = data_fit.attrs['T1_fit']*1e6
                    print(f"T1={t1}, error={fit_error}, ~ {round(fit_error*100/t1,1)} %")
                    if t1 < 50 and t1 > 0: 
                        T1_us.append(t1)
                        T1_err.append(fit_error)
                        gamma1_MHz.append(1/t1)
                        if fit_error*100/t1 > 50:
                            save_weired_data_pic(times, data, "T1", exp_idx, set_idx, temperature_folder_path,data_fit) 
                    else: 
                        T1_us.append(0)
                        gamma1_MHz.append(0)
                        T1_err.append(0)
                        save_weired_data_pic(times, data, "T1", exp_idx, set_idx, temperature_folder_path,data_fit)   
                except:
                    save_weired_data_pic(times, data, "T1", exp_idx, set_idx, temperature_folder_path)
                    T1_us.append(0)
                    gamma1_MHz.append(0)
                    T1_err.append(0)
            elif exp_type == "T2":
                T2_ds = open_dataset(file_path)
                times = array(Dataset.to_dict(T2_ds)["coords"]['x0']['data']) # s
                I,Q= dataset_to_array(dataset=T2_ds,dims=1)
                data= (IQ_data_dis(I,Q,ref_I=ref_iq[0],ref_Q=ref_iq[1]))
                try:
                    data_fit, fit_error = T2_fit_analysis(data=data,freeDu=times,T2_guess=8e-6, return_error=True)
                    t2 = data_fit.attrs['T2_fit']*1e6
                    print(f"T2={t2}, error={fit_error}, ~ {round(fit_error*100/t2,1)} %")
                    if t2 < 50 and t2 > 0: 
                        T2_us.append(t2)
                        gamma2_MHz.append(1/t2)
                        T2_err.append(fit_error)
                        if fit_error*100/t2 > 50 :
                            save_weired_data_pic(times, data, "T2", exp_idx, set_idx, temperature_folder_path, data_fit)
                    else:
                        T2_us.append(0)
                        gamma2_MHz.append(0)
                        T2_err.append(0)
                        save_weired_data_pic(times, data, "T2", exp_idx, set_idx, temperature_folder_path, data_fit)
                except:
                    save_weired_data_pic(times, data, "T2", exp_idx, set_idx, temperature_folder_path)
                    T2_us.append(0)
                    gamma2_MHz.append(0)
                    T2_err.append(0)
                
                
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
            
            
        print(f"First stage analysis complete for set-{set_idx}!")
            
        
        # reshape data to (I,Q)*(g,e)*shots       
        OS_data = 1000*array([[pgI_collection,peI_collection],[pgQ_collection,peQ_collection]]) # can train or predict 2*2*histo_counts*shot
        tarin_data, fit_arrays = OSdata_arranger(OS_data)
        # train GMM
        dist_model = train_GMModel (tarin_data[0])
        dist_model.relabel_model(array([ref_iq]).transpose())
        # predict all collection to calculate eff_T for every exp_idx
        for histo_i in range(fit_arrays.shape[0]):
            analysis_data = fit_arrays[histo_i] #your (2,2,N) data to analysis

            new_data = moveaxis( analysis_data ,1,0)
            p0_pop = dist_model.get_state_population(new_data[0].transpose()) # [p0in0, p0in1]
            p1_pop = dist_model.get_state_population(new_data[1].transpose()) # [p1in0, p1in1]
            OneShot_pic_path = os.path.join(SS_folder,f"SingleShot-S{set_idx}-{histo_i}")
            fig , eff_t, snr = plot_readout_fidelity(analysis_data, transi_freq, None, False)
            therm_pop.append(100-(p0_pop[0]*100/(p0_pop[0]+p0_pop[1])))
            effT_mK.append(eff_t)
            plt.close()
        
        # T1
        T1s[set_idx] = T1_us
        T1ers[set_idx] = T1_err 
        gamma1s[set_idx] = gamma1_MHz
    
        # T2
        T2s[set_idx] = T2_us
        T2ers[set_idx] = T2_err 
        gamma2s[set_idx] = gamma2_MHz

        # calc OnsShot
        effTs[set_idx] = effT_mK
        thermalpops[set_idx] = therm_pop

    T1_err_up_lim, T1_err_bt_lim = ret_error_bound([ T1ers[set_id] for set_id in T1ers],2)
    T2_err_up_lim, T2_err_bt_lim = ret_error_bound([ T2ers[set_id] for set_id in T2ers],2)

    # # Filter outlier
    for set_idx in T1s:
        info_dict = create_results_dict()
        # # cala T1
        gamma1ers[set_idx] = []
        for histo_idx in range(len(T1ers[set_idx])):
            error_value = T1ers[set_idx][histo_idx]
            if T1s[set_idx][histo_idx] != 0:
                if error_value > T1_err_up_lim or error_value < T1_err_bt_lim:
                    T1s[set_idx][histo_idx] = 0
                    T1ers[set_idx][histo_idx] = 0
                    gamma1s[set_idx][histo_idx] = 0
                    gamma1ers[set_idx].append(0)
                else:
                    err_ratio = T1ers[set_idx][histo_idx]/T1s[set_idx][histo_idx]
                    gamma1ers[set_idx].append(gamma1s[set_idx][histo_idx]*err_ratio)
            else:
                gamma1ers[set_idx].append(0)
        
        t1_aray, t1_err_aray, g1_aray, g1_err_aray = array(T1s[set_idx]), array(T1ers[set_idx]), array(gamma1s[set_idx]), array(gamma1ers[set_idx])
        if sum(t1_err_aray) == 0:
            mean_T1_us = 0
            sd_T1_us = 0
        else:
            mean_T1_us = round(average(t1_aray[t1_aray != 0],weights=1/(t1_err_aray[t1_err_aray != 0])),1)
            sd_T1_us = round(std(t1_aray[t1_aray != 0]),1)
        if sum(g1_err_aray[g1_err_aray != 0]) != 0:
            mean_gamma1_MHz = round(average(g1_aray[g1_aray != 0],weights=1/(g1_err_aray[g1_err_aray != 0])),2)
        else:
            mean_gamma1_MHz = round(mean(g1_aray[g1_aray != 0]),2)
        sd_gamma1_MHz = round(std(g1_aray[g1_aray != 0]),2)
        info_dict["T1"]["avg"], info_dict["T1"]["std"] = mean_T1_us, sd_T1_us
        info_dict["gamma1"]["avg"], info_dict["gamma1"]["std"] = mean_gamma1_MHz, sd_gamma1_MHz
        # histo_path = os.path.join(result_folder,f"T1-histo-S{set_idx}.png")
        # hist_plot("nobu",{"nobu":t1_aray},f"S{set_idx}, T1={mean_T1_us}+/-{sd_T1_us} us",histo_path, False)

        # # cala T2
        gamma2ers[set_idx] = []
        for histo_idx in range(len(T2ers[set_idx])):
            error_value = T2ers[set_idx][histo_idx]
            if T2s[set_idx][histo_idx] != 0:
                if error_value > T2_err_up_lim or error_value < T2_err_bt_lim:
                    T2s[set_idx][histo_idx] = 0
                    T2ers[set_idx][histo_idx] = 0
                    gamma2s[set_idx][histo_idx] = 0
                    gamma2ers[set_idx].append(0)
                else:
                    err_ratio = T2ers[set_idx][histo_idx]/T2s[set_idx][histo_idx]
                    gamma2ers[set_idx].append(gamma2s[set_idx][histo_idx]*err_ratio)
            else:
                gamma2ers[set_idx].append(0)
        
        t2_aray, t2_err_aray, g2_aray, g2_err_aray = array(T2s[set_idx]), array(T2ers[set_idx]), array(gamma2s[set_idx]), array(gamma2ers[set_idx])
        if sum(t2_err_aray) == 0:
            mean_T2_us = 0
            sd_T2_us = 0
        else:
            mean_T2_us = round(average(t2_aray[t2_aray != 0],weights=1/(t2_err_aray[t2_err_aray != 0])),1)
            sd_T2_us = round(std(t2_aray[t2_aray != 0]),1)
        if sum(g2_err_aray) != 0:
            mean_gamma2_MHz = round(average(g2_aray[g2_aray != 0],weights=1/(g2_err_aray[g2_err_aray != 0])),2)
        else:
            mean_gamma2_MHz = round(mean(g2_aray[g2_aray != 0]),2)
        sd_gamma2_MHz = round(std(g2_aray[g2_aray != 0]),2)
        info_dict["T2"]["avg"], info_dict["T2"]["std"] = mean_T2_us, sd_T2_us
        info_dict["gamma2"]["avg"], info_dict["gamma2"]["std"] = mean_gamma2_MHz, sd_gamma2_MHz
        # histo_path = os.path.join(result_folder,f"T2-histo-S{set_idx}.png")
        # hist_plot("nobu",{"nobu":t2_aray},f"S{set_idx}, T2={mean_T2_us}+/-{sd_T2_us} us",histo_path, False)

        # # calc T-phi
        gammaphis[set_idx] = []
        gammaphiers[set_idx] = []
        for histo_idx in range(len(T1s[set_idx])):
            if gamma1s[set_idx][histo_idx] != 0 and gamma2s[set_idx][histo_idx] != 0:
                gammaphis[set_idx].append(gamma2s[set_idx][histo_idx] - 0.5*(gamma1s[set_idx][histo_idx]))
                error_gamma_phi = sqrt((gamma2ers[set_idx][histo_idx])**2 + (0.5*(gamma1ers[set_idx][histo_idx]))**2)
                gammaphiers[set_idx].append(error_gamma_phi)
            else:
                gammaphis[set_idx].append(0)
                gammaphiers[set_idx].append(0)

        gphi_aray, aphi_err_aray = array(gammaphis[set_idx]), array(gammaphiers[set_idx])
        if sum(aphi_err_aray) == 0 :
            mean_gammaphi_MHz = 0
            sd_gammaphi_MHz = 0
        else:
            mean_gammaphi_MHz = round(average(gphi_aray[gphi_aray != 0],weights=1/(aphi_err_aray[aphi_err_aray != 0])),2)
            sd_gammaphi_MHz = round(std(gphi_aray[gphi_aray != 0]),2)
        info_dict["gammaPhi"]["avg"], info_dict["gammaPhi"]["std"] = mean_gammaphi_MHz, sd_gammaphi_MHz

        # # calc effT (no filters now)
        effT_mK = array(effTs[set_idx])
        therm_pop = array(thermalpops[set_idx])
        mean_effT_mK = round(mean(effT_mK),1)
        sd_effT_mK = round(std(effT_mK),1)
        info_dict["effT"]["avg"], info_dict["effT"]["std"] = mean_effT_mK, sd_effT_mK
        info_dict["thermalPop"]["avg"], info_dict["thermalPop"]["std"] = round(mean(therm_pop),2), round(std(therm_pop),2)

        # save the info to plt scatter
        with open(f"{json_folder}/setInfo({set_idx}).json", "w") as record_file: 
            json.dump(info_dict,record_file)
    
    every_value = {}
    parts = {"T1s":T1s,"T1ers":T1ers,"gamma1s":gamma1s,"gamma1ers":gamma1ers,"T2s":T2s,"T2ers":T2ers,"gamma2ers":gamma2ers,"gamma2s":gamma2s,"effTs":effTs,"thermalPops":thermalpops,"gammaPhis":gammaphis,"gammaPhiers":gammaphiers}
    for part in parts:
        infos = parts[part]
        myKeys = list(infos.keys())
        myKeys.sort()
        sorted_infos = {i: infos[i] for i in myKeys}
        every_value[part] = array([sorted_infos[set_idx] for set_idx in sorted_infos]).reshape(-1).tolist()
    

    with open(f"{json_folder}/every_values.json", "w") as record_file:   
        json.dump(every_value,record_file)

    end = time.time()
    print(f"Analysis time cost: {round((end-start)/60,1)} mins")

# ============================================ Plot ================================================================================
def live_time_monitoring_plot(monitor_dict:dict,pic_save_folder_path:str):
    plot_catas = ["T1", "T2"]
    if len(monitor_dict["effT"]) != 0:
        plot_catas.append("effT")
    else:
        plot_catas.append("thermalPop")
    
    for exp in plot_catas:
        fig, ax = plt.subplots(1,1,figsize=fig_size)
        ax:plt.Axes
        ax.scatter( monitor_dict["x_minutes"],monitor_dict[exp], s=150)
        ax.errorbar(monitor_dict["x_minutes"],monitor_dict[exp],yerr=monitor_dict[f"{exp}_sd"])
        ax.xaxis.set_tick_params(labelsize=tick_num_size)
        ax.yaxis.set_tick_params(labelsize=tick_num_size)
        ylabel = f"{exp} (us)" if exp in ["T1","T2"] else ("Effective temp. (mK)" if exp=="effT" else "Thermal Populations")
        ax.set_ylabel(ylabel, fontsize=label_font_size)
        ax.set_xlabel(f"time past (min)", fontsize=label_font_size)
        ax.set_ylim(min(array(monitor_dict[exp]))-std(array(monitor_dict[exp])), max(array(monitor_dict[exp]))+std(array(monitor_dict[exp])))
        plt.title(f"{exp} live monitor",fontsize=30)
        plt.tight_layout()
        plt.savefig(os.path.join(pic_save_folder_path,f"{exp}_TimeMonitor.png"))
        plt.close()

def plot_temp_compa_timeMonitor(temp_name_list:list, target_q:str, sample_folder_name:str="", conditional_folder_name:str="", exp_type="T1")->plt.Axes:
    """
    plot T1/effT/T2/... with different temperature in a fiqure along time axis.\n
    exp_type: "T1",... refewr to exp_items in RadiatorSetana.py .
    """
    if exp_type not in list(exp_items.values()): raise KeyError(f"Un-supported exp_type={exp_type} was given!")
    
    fig,ax = plt.subplots(1,1,figsize=fig_size)
    plt.grid()
    ax:plt.Axes
    
     # Ex. 20K-2: {"T1":{"avg","std"},"T2":{"avg","std"},"eff_T":{"avg","std"}}
    for temperature in temp_name_list:
        info_recorder = {}
        temperature_folder = os.path.join(meas_raw_dir,sample_folder_name,conditional_folder_name,temperature)
        json_folder = os.path.join(temperature_folder,"results","jsons")
        with open(os.path.join(json_folder,"temperatureInfo.json")) as J:
            info_recorder = json.load(J)
        
        time_past_sec_array = round(get_time_axis(target_q,temperature_folder)/60,1)
        ax.errorbar(time_past_sec_array,info_recorder[exp_type]["avg"],yerr=info_recorder[exp_type]["std"],fmt='o-',label=f"{temperature}")
        
    ax.set_title(f"Temperature comparisom on {exp_type} by Qblox", fontsize=label_font_size)
    ax.set_xlabel(f"Time past (min)",fontsize=label_font_size) 
    ax.xaxis.set_tick_params(labelsize=tick_num_size)
    ax.yaxis.set_tick_params(labelsize=tick_num_size)   
    ax = ax_set_y_label(ax,exp_type,label_font_size)
    
    return ax

def plot_time_behavior_sep(temperature_folder:str, mode:str, time_axis:ndarray=array([]), radiator_act:str="ON", refresh_tempera_info:bool=False)->plt.Axes:
    """
    Plot T1/T2/eff_T trend along time axis.\n
    arg `mode` assign which value to plot : 'T1', 'T2', 'effT',  'gamma1', 'gamma2', 'gammaPhi', 'thermalPop'
    """
    tempera_into = collect_allSets_inTempera(temperature_folder,refresh_tempera_info)
    a_set_time = 7 # mins
    time_axis_min = arange(1,len(tempera_into["T1"]["avg"]))*a_set_time if time_axis.shape[0]==0 else round(time_axis/60,1)  
    
    ax_info = {"exp":mode, "x":time_axis_min.tolist()}

    fig, ax = plt.subplots(1,1,figsize=fig_size)
    ax:plt.Axes
    if mode in list(exp_items.values()):
        times, ys, yerrs = filter_zero(array(tempera_into[mode]["avg"]), array(time_axis_min), array(tempera_into[mode]["std"]))
        ax.errorbar(times,ys,yerr=yerrs,fmt="o-",color='red')
        ax_info["y"], ax_info["yerr"] = tempera_into[mode]["avg"], tempera_into[mode]["std"]
        ax.yaxis.label.set_color('red')

    ax = ax_set_y_label(ax, mode, label_font_size)
        
    ax.xaxis.set_tick_params(labelsize=tick_num_size)
    ax.yaxis.set_tick_params(labelsize=tick_num_size)
    ax.set_xlabel(f"time after radiator {radiator_act} (min)", fontsize=label_font_size)
    
    return ax, ax_info

def plot_ref_onAxes(ax:plt.Axes,ref_dict_b4:dict={},ref_dict_recover:dict={}, which:str="T1")->tuple[plt.Axes,dict]:
    """ #### ref_dict keys should follow the form which return from get_ref_from_json in a given condition\n
        which: see exp_items in RadiatorSetAna.py\n
        If the dict is empty, we don't plot it on Axes. 
    """
    ax_info = {"exp":f"ref-{which}","x":[],"y":[],"yerr":[]}
    if ref_dict_b4 != {}:
        key_names_b4 = list(ref_dict_b4.keys())
        if which in list(exp_items.values()):
            if (which in key_names_b4) and (f"{which}_sd" in key_names_b4):
                ax.axhline(y=ref_dict_b4[which],color='#000080',label=f"{which}-Before turn ON radiator",lw=3.5)
                ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_b4[which]-ref_dict_b4[f"{which}_sd"],y2=ref_dict_b4[which]+ref_dict_b4[f"{which}_sd"],color='#1E90FF',alpha=0.3)
                ax_info["x"].append(4), ax_info["y"].append(ref_dict_b4[which]), ax_info["yerr"].append(ref_dict_b4[f"{which}_sd"])

    if ref_dict_recover != {}:
        key_names_rc = list(ref_dict_recover.keys())
        if which in list(exp_items.values()):
            if (which in key_names_rc) and (f"{which}_sd" in key_names_rc):
                ax.axhline(y=ref_dict_recover[which],color='#FA8072',label=f"{which}-LONG enough After turn OFF radiator",lw=3.5)
                ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_recover[which]-ref_dict_recover[f"{which}_sd"],y2=ref_dict_recover[which]+ref_dict_recover[f"{which}_sd"],color='#FFA07A',alpha=0.3)
                ax_info["x"].append(4), ax_info["y"].append(ref_dict_recover[which]), ax_info["yerr"].append(ref_dict_recover[f"{which}_sd"])  
    
    return ax, ax_info
    
def plot_DR_tempera(start_date:str="",start_time:str="",temp_chennel:int=6,DR_log_folder_path:str="",ax:plt.Axes=None,time_axis:ndarray=array([]))-> plt.Axes:
    """
    ## Plot only when start_date, start_time, and DR_log_folder_path are ALL given.
    """
    dr_info = {"exp":f"Chennel-{temp_chennel}","x":[],"y":[],"yerr":[]} # no average here
    if ax is None:
        fig, axDR = plt.subplots(1,1,figsize=fig_size)
        axDR.set_xlabel("Time past (min)", fontsize=label_font_size)
    else:
        if start_date == "" or start_time == "" or DR_log_folder_path == "":
            axDR = ax
        else:
            axDR:plt.Axes = ax.twinx()
            axDR.spines['right'].set_color("green")
            axDR.yaxis.label.set_color("green")
    # call DR temperature coresponding to temp. chennel
    if start_date == "" or start_time == "" or DR_log_folder_path == "":
        pass
    else:
        dr_time_array, dr_temp_array = Kelvin_collector(DR_log_folder_path,start_date,start_time,int((time_axis[-1]-time_axis[0])/60),temp_chennel)
        
        if temp_chennel in [6, "6"]:
            y_name = 'MXC temp. (mK)'
            dr_temp_array *= 1000
        elif temp_chennel in [2, "2"]:
            y_name = '4K-plate temp. (K)'
        else:
            y_name = f"Chennel-{temp_chennel} temp. (K)"
        
        axDR.plot(dr_time_array, dr_temp_array, c="cyan",lw=3.5)
        axDR.plot(dr_time_array, dr_temp_array, c="green",lw=2.8)
        axDR.set_ylabel(y_name,fontsize=label_font_size)
        axDR.set_ylim(min(dr_temp_array)-std(dr_temp_array), max(dr_temp_array)+std(dr_temp_array))
        axDR.xaxis.set_tick_params(labelsize=tick_num_size)
        axDR.yaxis.set_tick_params(labelsize=tick_num_size)
        dr_info["x"] = dr_time_array.tolist()
        dr_info["y"] = dr_temp_array.tolist()
        dr_info["yerr"] = zeros(dr_temp_array.shape[0]).tolist()
    return axDR, dr_info

def scat_DR_avg_temp(need_log_info:dict,sample_folder_name:str="",conditional_folder_name:str="",DR_log_folder_path:str="",ax:plt.Axes=None,temp_chennel:int=6)->tuple[plt.Axes, plt.Axes]:
    """ 
    `need_log_info` follows the form: {temp:{"keep_time_min", "start_date", "start_time", "avg_min_from_the_end"}}
    """
    button:bool = 1
    dr_info = {"exp":f"Chennel-{temp_chennel}"}
    if ax is None :
        plt.close()
        fig, axDR = plt.subplots(1,1,figsize=fig_size)
        axDR.set_xlabel("Radiator temp. (K)", fontsize=label_font_size)
        x_axis = []
        x_axis += [int(temp.split("K")[0]) for temp in need_log_info]
   

    dr_temp_avg = []
    dr_temp_std = []
    x_axis = []
    for temperature in need_log_info:
        exp_keep_time_min:int = need_log_info[temperature]["keep_time_min"]
        avg_min_from_the_end:int = need_log_info[temperature]["avg_min_from_the_end"] if "avg_min_from_the_end" in  need_log_info[temperature] else 60
        try:
            other_info = {}
            with open(os.path.join(meas_raw_dir,sample_folder_name,conditional_folder_name,temperature,"otherInfo.json")) as JJ:
                other_info = json.load(JJ)
            start_date:str = other_info[target_q]["start_time"].split(" ")[0]
            start_time:str = other_info[target_q]["start_time"].split(" ")[-1]
        except:
            try:
                start_date:str = need_log_info[temperature]["start_date"]
                start_time:str = need_log_info[temperature]["start_time"]
            except:
                start_date = ""
                start_time = ""

        if DR_log_folder_path == "" or start_date == "" or start_time == "":
            pass
        else:
            if button:
                axDR:plt.Axes = ax.twinx()
                axDR.spines['right'].set_color("green")
                axDR.yaxis.label.set_color("green")
                button = 0

            dr_time_array, dr_temp_array = Kelvin_collector(DR_log_folder_path,start_date,start_time,exp_keep_time_min,temp_chennel)
            avg_start_min = int(dr_time_array[-1] - avg_min_from_the_end)
            idx, _ = find_nearest(dr_time_array, avg_start_min)
            dr_temp_avg.append(mean(dr_temp_array[idx:]))
            dr_temp_std.append(std(dr_temp_array[idx:]))
            x_axis.append(int(temperature.split("K")[0]))
    dr_temp_avg = array(dr_temp_avg)
    dr_temp_std = array(dr_temp_std)

    if int(temp_chennel) == 6 :
        y_name = 'MXC temp. (mK)'
        dr_temp_avg *= 1000
        dr_temp_std *= 1000
    elif int(temp_chennel) == 2 :
        y_name = '4K-plate temp. (K)'
    else:
        y_name = f"Chennel-{temp_chennel} temp. (K)"

    if button:
        axDR = ax
        dr_info["x"], dr_info["y"], dr_info["yerr"] = [], [], []
    else:
        dr_info["x"], dr_info["y"], dr_info["yerr"] = x_axis, dr_temp_avg, dr_temp_std
        axDR.errorbar(x_axis, dr_temp_avg, yerr=dr_temp_std, fmt='o-', c="green")
        axDR.set_ylabel(y_name,fontsize=label_font_size)
        axDR.set_ylim(min(dr_temp_avg)-std(dr_temp_avg), max(dr_temp_avg)+std(dr_temp_avg))
        axDR.xaxis.set_tick_params(labelsize=tick_num_size) 
        axDR.yaxis.set_tick_params(labelsize=tick_num_size)    
    return axDR, dr_info

def plot_stable_temp_dep(temp_folder_names:list, exp_type:str, slice_min_from_the_end:list=[0], sameple_folder:str="", conditional_folder:str="")->tuple[plt.Axes, dict, str, dict]:
    """exp_type should be one of the exp_items"""
    keep_time = {}
    x_axis = []
    exp_value = {"avgs":[],"stds":[]}
    if exp_type in list(exp_items.values()):
        for idx, temperature in enumerate(temp_folder_names):
            every_value_dict = {}
            tempera_folder = os.path.join(meas_raw_dir,sameple_folder,conditional_folder,temperature)
            time_past_min_array = get_time_axis(target_q,tempera_folder)/60
            keep_time[temperature] = int(time_past_min_array[-1]-time_past_min_array[0])
            set_number = time_past_min_array.shape[0]
            json_file = os.path.join(tempera_folder,"results","jsons","every_values.json")
            with open(json_file) as JJ:
                every_value_dict = json.load(JJ)
                x_axis.append(int(temperature.split("K")[0]))

            slice_from_this_min = int(time_past_min_array[-1]-slice_min_from_the_end[idx])
            histo_count_in_set = int(len(every_value_dict[f"{exp_type}s"])/set_number) # data number of a exp in every_value_dict = set_number * histo_count_in_set
            this_min_idx_in_set, _ = find_nearest(time_past_min_array, slice_from_this_min) # this index is equal to the set index
            histo_idx_this_min = (this_min_idx_in_set+1)*histo_count_in_set 
            target_data = array(every_value_dict[f"{exp_type}s"][histo_idx_this_min:]) 
            this_std = std(target_data[target_data != 0])
            if exp_type in ["effT","thermalPop"]:
                exp_value["avgs"].append(average(target_data[target_data != 0]))
            else:
                target_erro = array(every_value_dict[f"{exp_type}ers"][histo_idx_this_min:]) 
                exp_value["avgs"].append(average(target_data[target_data != 0], weights=1/(target_erro[target_erro != 0])))
            exp_value["stds"].append(this_std)
    else:
        raise ValueError (f" The given exp_type={exp_type} can't be recognized!")

    fig, ax = plt.subplots(1,1,figsize=fig_size)   
    plt.grid()  
    ax: plt.Axes
    ax.yaxis.label.set_color("red")
    ax.errorbar(x_axis,exp_value["avgs"],yerr=exp_value["stds"],fmt="o-",c='red')
    ax.set_xlabel("Radiator temp. (K)",fontsize=label_font_size)
    
    ax = ax_set_y_label(ax, exp_type, label_font_size)
    
    ax.xaxis.set_tick_params(labelsize=tick_num_size)
    ax.yaxis.set_tick_params(labelsize=tick_num_size)


    return ax, keep_time, os.path.split(tempera_folder)[0], {"exp":exp_type,"x":x_axis,"y":exp_value["avgs"],"yerr":exp_value["stds"]}

def scat_ref_temp_dep(ax:plt.Axes,ref_dict_b4:dict={},ref_dict_recover:dict={}, which:str="T1")->tuple[plt.Axes, list, list, dict]:
    """ #### ref_dict keys should follow the form: `{"T1","T1_sd","T2","T2_sd","effT","effT_sd"}`\n
        the element in dict must be given in value-sd pair, which kinds of pair is optional.\n
        If the dict is empty, we don't plot it on Axes. 
    """
    x = []
    y = []
    yerr = []
    if ref_dict_b4 != {}:
        key_names_b4 = list(ref_dict_b4.keys())
        if which in list(exp_items.values()):
            if (which in key_names_b4) and (f"{which}_sd" in key_names_b4):
                ax.scatter([4],ref_dict_b4[which], s=150, marker='X')
                ax.errorbar([4],ref_dict_b4[which],yerr=ref_dict_b4[f"{which}_sd"],fmt="X",label=f"{which}-Before turn ON radiator")
                x.append(4)
                y.append(ref_dict_b4[which])
                yerr.append(ref_dict_b4[f"{which}_sd"])
    if ref_dict_recover != {}:
        key_names_rc = list(ref_dict_recover.keys())
        if which in list(exp_items.values()):
            if (which in key_names_rc) and (f"{which}_sd" in key_names_rc):
                ax.scatter([4],ref_dict_recover[which], s=150, marker='D')
                ax.errorbar([4],ref_dict_recover[which],yerr=ref_dict_recover[f"{which}_sd"],fmt="D",label=f"{which}-LONG enough After turn OFF radiator")
                x.append(4)
                y.append(ref_dict_recover[which])
                yerr.append(ref_dict_recover[f"{which}_sd"])

    handles, labels = ax.get_legend_handles_labels()
    return ax, handles, labels, {"exp":f"ref-{which}","x":x,"y":y,"yerr":yerr}
    
#                         ==============
#                         =            =
# ========================= user layer ================================
#                         =            =
#                         ==============
def time_trend_artist(tempera_folder:str, target_q:str, exp_catas:list, time_past_sec_array:ndarray, ref_before:dict, ref_recove:dict, DR_time_info:dict, log_folder:str, DRtemp_che:int=6, show:bool=0, refresh_tempera_info:bool=False):
    try:
        other_info = {}
        with open(os.path.join(tempera_folder,"otherInfo.json")) as JJ:
            other_info = json.load(JJ)
        start_date = other_info[target_q]["start_time"].split(" ")[0]
        start_time = other_info[target_q]["start_time"].split(" ")[-1]
    except:
        print("OtherInfo.json didn't work!")
        try:
            start_date = DR_time_info["start_date"]
            start_time = DR_time_info["start_time"]
        except:
            start_date = ""
            start_time = ""
    temperature = os.path.split(tempera_folder)[-1]
    if temperature[:2] == "re":
        action = "OFF"
    else:
        action = "ON"
    for exp in exp_item_translator(exp_catas):
        ax, trend_info = plot_time_behavior_sep(tempera_folder,exp,time_axis=time_past_sec_array,radiator_act=action,refresh_tempera_info=refresh_tempera_info)
        ax, ref_info = plot_ref_onAxes(ax,ref_before,ref_recove,exp)
        handles, labels = ax.get_legend_handles_labels()
        plt.grid(axis='both')
        ax, dr_info = plot_DR_tempera(start_date,start_time,DRtemp_che,log_folder,ax,time_axis=time_past_sec_array)
        new_folder_path = timelabelfolder_creator(tempera_folder,f"{target_q}_{temperature}TIME_{exp}")
        pic_values_saver(new_folder_path,'time',trend_info, ref_info, dr_info)
        plt.legend(handles, labels, fontsize=legend_font_size,bbox_to_anchor=(1,1.2))
        plt.tight_layout()
        if not show:
            plt.savefig(os.path.join(new_folder_path,f"{target_q}_{exp}_timeMonitor.png"))
            plt.close()
        else:
            plt.show()

def temp_depen_artist(temperature_list:list, target_q:str, sample_folder:str, conditional_folder:str, log_info_dict:dict, exp_catas:list, ref_before:dict={}, ref_recove:dict={}, log_folder="", tempera_che:int=6):
    avg_time_list = []
    for temp in temperature_list: 
        try:
            avg_time_list.append(log_info_dict[temp]["avg_min_from_the_end"])
        except:
            avg_time_list.append(60)
    for exp in exp_item_translator(exp_catas):
        ax, keep_time_info, conditional_folder_path, axes_info = plot_stable_temp_dep(temperature_list, exp, avg_time_list, sample_folder, conditional_folder)
        for tempera in log_info_dict:
            log_info_dict[tempera]["keep_time_min"] = keep_time_info[tempera]
        ax, handles, labels, ref_info = scat_ref_temp_dep(ax,ref_before,ref_recove,exp)
        
        ax, dr_info = scat_DR_avg_temp(log_info_dict,sample_folder,conditional_folder,DR_log_folder_path=log_folder,ax=ax,temp_chennel=tempera_che)
        new_folder_path = timelabelfolder_creator(conditional_folder_path,f"{target_q}_TEMP_{exp}")
        pic_values_saver(new_folder_path,'temp',axes_info, ref_info, dr_info)
        plt.legend(handles, labels, fontsize=legend_font_size,bbox_to_anchor=(1,1.15))
        plt.tight_layout()
        plt.savefig(os.path.join(new_folder_path,f"{target_q}_{exp}_tempDependence.png"))
        plt.close()
        # plt.show()
    
def TimeMonitor_tempCompa(temp_name_list:list, target_q, exp_catas, sample_folder_name:str, conditional_folder_name:str):
    temp_compa_path = os.path.join(meas_raw_dir,"Temperature Compa")
    if not os.path.isdir(temp_compa_path):
        os.mkdir(temp_compa_path)
    for exp in exp_item_translator(exp_catas):
        ax = plot_temp_compa_timeMonitor(temp_name_list, target_q, sample_folder_name, conditional_folder_name, exp)
        ax.legend(fontsize=legend_font_size)
        plt.tight_layout()
        plt.savefig(os.path.join(temp_compa_path,f"{exp}_temp_compa.png"))
        plt.close()
        # plt.show()

if __name__ == '__main__':
    # // If you wanna plot DR temperature, "start_date" and "start_time" in log_info_dict and also log_folder ALL should be given !
    # *********** Manully settings ***********
    analysis:bool = 1   # analysis or not. Once u had analyzed, u won't need it again
    plot_time_trend:bool = 1   # If there is only one element in `log_info_dict`, polt the time monitoring, which is time as the x-axis.
    always_plot_timeTrend:bool = 1  # if there are more than one element in `log_info_dict`, turn on this will plot time monitoring for all elements. 
    plot_temp_dependence: bool = 0 # if there are more than one element in `log_info_dict`, plot the radiator-temp dependence, radiator_temp as x-asix. 
    plot_tempComp_alongTime:bool = 0   # plot data of different radiator temp along time as x-axis.
    
    target_q = 'q0'
    exp_catas = [1,2,6]       # {"1":"T1","2":"T2","3":"effT","4":"gamma1","5":"gamma2","6":"thermal_Pop","7":"gammaPhi"}
    conditional_folder = "ScalinQ_q4"          # the previous folder from temperature folder     
    sample_folder = "MXC_heater"     # the previous folder from conditional_folder
    
    log_info_dict = {"re0K":{}} # if keyname 'avg_min_from_the_end' is not inside, the default is 60 minutes
    
    # C:\Users\ASQUM\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\MXC_heater\ScalinQ_q4\20K
    # ? For references.
    ref_before = get_ref_from_json(target_q,sample_folder,conditional_folder)["ref_before"] if "ref_before" in get_ref_from_json(target_q,sample_folder,conditional_folder).keys() else {}
    ref_recove = get_ref_from_json(target_q,sample_folder,conditional_folder)["ref_recove"] if "ref_recove" in get_ref_from_json(target_q,sample_folder,conditional_folder).keys() else {}


    # ? If this folder is "", it won't plot the MXC temperature. 
    log_folder = ""

    # # ! Don't touch below !
    temperature_list = list(log_info_dict.keys())
    if always_plot_timeTrend :
        for tempera in temperature_list:
            tempera_folder = os.path.join(meas_raw_dir,sample_folder,conditional_folder,tempera)
            if analysis:
                print("y")
                main_analysis(target_q, tempera_folder)
            time_past_sec_array = get_time_axis(target_q,tempera_folder)
            # plot time trend with a given temperature (time as x-axis, in min)
            if plot_time_trend:
                time_trend_artist(tempera_folder, target_q, exp_catas, time_past_sec_array, ref_before, ref_recove, log_info_dict[tempera], log_folder, DRtemp_che=6,refresh_tempera_info=True)
    else:
        if plot_temp_dependence:
            temp_depen_artist(temperature_list, target_q, sample_folder, conditional_folder, log_info_dict, exp_catas, ref_before, ref_recove, log_folder, tempera_che=6)

    if plot_tempComp_alongTime:
        TimeMonitor_tempCompa(temperature_list, target_q, exp_catas, sample_folder, conditional_folder)


                


            
        
        
        
   