"""
This program is only for analyzing a series of radiation tests like 0K, 20K 40K 60K and re0K with/without shielding. This didn't support comparison between with and without shielding
"""
import os, sys, time, json, pickle
from pandas import DataFrame as DF
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis, T2_fit_analysis, plot_textbox, Fit_analysis_plot
from xarray import Dataset, open_dataset # Dataset.to_dict(SS_ds)
from numpy import array, ndarray, mean, std, round, arange, moveaxis, any, zeros
# from Modularize.support.Pulse_schedule_library import hist_plot
import matplotlib.pyplot as plt
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.analysis.DRtemp import Kelvin_collector
from qcat.state_discrimination.discriminator import train_GMModel  # type: ignore
from qcat.visualization.readout_fidelity import plot_readout_fidelity

# ================================ Functional =========================================================================================
def find_nearest(ary:ndarray, value:float)->tuple[int,float]:
    """ find the element  which is closest to the given target_value in the given array"""
    idx = (abs(ary - value)).argmin()
    return idx, float(ary[idx])

def timelabelfolder_creator(folder_path:str,additional_folder_name:str='')->str:
    import datetime
    current_time = datetime.datetime.now()
    temp_dep_folder_name =  f"{additional_folder_name}_H{current_time.hour}M{current_time.minute}S{current_time.second}"
    temp_dep_folder_path = os.path.join(folder_path,temp_dep_folder_name)
    if not os.path.exists(temp_dep_folder_path):
        os.mkdir(temp_dep_folder_path)
    return temp_dep_folder_path

def pic_values_saver(folder_path:str,mode:str,*args):
    """ 
        mode indicates the axis in this pic is time with the mode='time', or temperature with the mode='temp'.
    """
    exp_type = [data_dict["exp"].upper() for data_dict in args]
    if "T1" in exp_type:
        title=["T1","us"]
    elif "T2" in exp_type:
        title=["T2","us"]
    elif  any(array(["EFFT" in exp_type, "EFF_T" in exp_type])):
        title=["effT","mK"]
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
        if data_dict["exp"].lower() in ["t1", "t2", "efft", "eff_t"]:
            to_csv_dict[f"exp_x_({x_axis_unit})"] = data_dict["x"]
            to_csv_dict[f"exp_y_({title[-1]})"] += data_dict["y"]
            to_csv_dict[f"exp_yerr_({title[-1]})"] += data_dict["yerr"]

        elif data_dict["exp"].lower() in ["ref-t1", "ref-t2", "ref-efft", "ref-eff_t"]:
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
    exp_type = ["T1", "T2", "eff_T"]
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


def OSdata_arranger(total_array:ndarray, want_IQset_num:int=1):
    """
    total_array shape: (2,2,M,N) which is (IQ, state, histos, shots)\n
    return traning, predict, label arrays:\n 1) (want_IQset_num, IQ, state, shots)\n 2) (histos, IQ, state, shots)\n
    3) (g, IQ, shot)
    """
    from numpy.random import randint
    total_sets = []
    train_set  = []
    for_label = []
    for histo in range(total_array.shape[2]):
        IQ_folder = []
        for iq in range(total_array.shape[0]):
            state_folder = []
            for state in range(total_array.shape[1]):
                state_folder.append(total_array[iq][state][histo])
            IQ_folder.append(state_folder)
        total_sets.append(IQ_folder)
    print(array(total_sets).shape) # shape = (histo, IQ, State, shots)

    for pick_number in range(want_IQset_num):
        rand_pick_set_idx = randint(0 ,len(total_sets))
        train_set.append(total_sets[rand_pick_set_idx])
    
    return array(train_set), array(total_sets)

def collect_allSets_inTempera(temperature_folder_path:str)->dict:
    json_path = os.path.join(f"{temperature_folder_path}","results","jsons")
    info_dict = {}
    if os.path.exists(os.path.join(json_path,"temperatureInfo.json")):
        print("read old json")
        with open(os.path.join(json_path,"temperatureInfo.json")) as J:
            info_dict = json.load(J)
    else:
        avg_t1, std_t1 = [], []
        avg_t2, std_t2 = [], []
        avg_eff_T, std_eff_T = [], []
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
                avg_eff_T.append(float(set_dict["eff_T"]["avg"]))
                std_eff_T.append(float(set_dict["eff_T"]["std"]))
        info_dict = {"T1":{"avg":avg_t1,"std":std_t1},"T2":{"avg":avg_t2,"std":std_t2},"eff_T":{"avg":avg_eff_T,"std":std_eff_T}}  # contains the AVG and SG for every set
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

# ============================================ Analysis ================================================================================

def main_analysis(target_q:str, temperature_folder_path:str, mode:str='quick'):
    """
    target_q: 'q0'\n
    temperature: '10K'\n
    mode: 'quick' or 'detail', represented save all the results pic or just histogram.\n
    """
    other_info_dict = {}
    temperature = os.path.split(temperature_folder_path)[-1]
    parent_path = temperature_folder_path

    sub_folders = [name for name in os.listdir(temperature_folder_path) if (os.path.isdir(os.path.join(temperature_folder_path,name)) and name[:8]=='Radiator')] # Radiator(idx)
    other_info_ = [name for name in os.listdir(temperature_folder_path) if (os.path.isfile(os.path.join(temperature_folder_path,name)) and name.split(".")[-1]=='json')]
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
    T1s = []
    T2s = []
    effTs = []


    for folder_name in sub_folders:
        info_dict = create_results_dict()
        set_idx = folder_name.split("(")[-1].split(")")[0]
        folder_path = os.path.join(temperature_folder_path,folder_name)
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
                    if mode == 'detail':
                        Fit_analysis_plot(data_fit,P_rescale=False,Dis=None,save_folder=os.path.join(T1_pic_folder,f"T1_S{set_idx}-{exp_idx}.png"))
                    T1_us.append(data_fit.attrs['T1_fit']*1e6)
                    if data_fit.attrs['T1_fit']*1e6 > 50: 
                        save_weired_data_pic(times, data, "T1", exp_idx, set_idx, temperature_folder_path,data_fit)   
                except:
                    save_weired_data_pic(times, data, "T1", exp_idx, set_idx, temperature_folder_path)
                    T1_us.append(0)
            elif exp_type == "T2":
                T2_ds = open_dataset(file_path)
                times = array(Dataset.to_dict(T2_ds)["coords"]['x0']['data']) # s
                I,Q= dataset_to_array(dataset=T2_ds,dims=1)
                data= (IQ_data_dis(I,Q,ref_I=ref_iq[0],ref_Q=ref_iq[1]))
                try:
                    data_fit= T2_fit_analysis(data=data,freeDu=times,T2_guess=8e-6)
                    if mode == 'detail':
                        Fit_analysis_plot(data_fit,P_rescale=False,Dis=None,save_folder=os.path.join(T2_pic_folder,f"T2_S{set_idx}-{exp_idx}.png"))
                    T2_us.append(data_fit.attrs['T2_fit']*1e6)
                    if data_fit.attrs['T2_fit']*1e6 > 30:
                        save_weired_data_pic(times, data, "T2", exp_idx, set_idx, temperature_folder_path, data_fit)
                except:
                    save_weired_data_pic(times, data, "T2", exp_idx, set_idx, temperature_folder_path)
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
            p0_pop = dist_model.get_state_population(new_data[0].transpose())
            p1_pop = dist_model.get_state_population(new_data[1].transpose())
            OneShot_pic_path = os.path.join(SS_folder,f"SingleShot-S{set_idx}-{histo_i}")
            fig , eff_t, snr = plot_readout_fidelity(analysis_data, transi_freq, OneShot_pic_path, False)
            effT_mK.append(eff_t)
            plt.close()
        
        # calc T1
        T1s.append(T1_us)
        T1_us = array(T1_us)
        mean_T1_us = round(mean(T1_us[T1_us != 0]),1)
        sd_T1_us = round(std(T1_us[T1_us != 0]),1)
        info_dict["T1"]["avg"], info_dict["T1"]["std"] = mean_T1_us, sd_T1_us
        histo_path = os.path.join(result_folder,f"T1-histo-S{set_idx}.png")
        hist_plot("nobu",{"nobu":T1_us},f"S{set_idx}, T1={mean_T1_us}+/-{sd_T1_us} us",histo_path, False)
        # calc T2
        T2s.append(T2_us)
        T2_us = array(T2_us)
        mean_T2_us = round(mean(T2_us[T2_us != 0]),1)
        sd_T2_us = round(std(T2_us[T2_us != 0]),1)
        info_dict["T2"]["avg"], info_dict["T2"]["std"] = mean_T2_us, sd_T2_us
        histo_path = os.path.join(result_folder,f"T2-histo-S{set_idx}.png")
        hist_plot("nobu",{"nobu":T2_us[T2_us != 0]},f"S{set_idx}, T2={mean_T2_us}+/-{sd_T2_us} us",histo_path, False)
        # calc OnsShot
        effTs.append(effT_mK)
        effT_mK = array(effT_mK)
        mean_effT_mK = round(mean(effT_mK),1)
        sd_effT_mK = round(std(effT_mK),1)
        info_dict["eff_T"]["avg"], info_dict["eff_T"]["std"] = mean_effT_mK, sd_effT_mK
        # save the info to plt scatter
        with open(f"{json_folder}/setInfo({set_idx}).json", "w") as record_file: 
            json.dump(info_dict,record_file)

    every_value = {"T1s":list(array(T1s).reshape(-1)),"T2s":list(array(T2s).reshape(-1)),"effTs":list(array(effTs).reshape(-1))} # This contains all the data point to calc a AVG and SD in a temp
    with open(f"{json_folder}/every_values.json", "w") as record_file:   
            json.dump(every_value,record_file)

    end = time.time()
    print(f"Analysis time cost: {round((end-start)/60,1)} mins")



# ============================================ Plot ================================================================================

def plot_temp_compa(mode:str="all"):
    """
    plot T1/effT/T2 with different temperature in a fiqure along time axis.\n
    mode: 'all', 'part1' or 'part2'. Default is 'all'.
    ## this will be modified in the future
    """
    temp_compa_path = os.path.join(meas_raw_dir,"Temperature Compa")
    if not os.path.isdir(temp_compa_path):
        os.mkdir(temp_compa_path)

    temp_folders = [name for name in os.listdir(meas_raw_dir) if (os.path.isdir(os.path.join(meas_raw_dir,name)) and (len(name.split("K"))==2 and name.split("K")[-1][0]=='-'))]
    
    include_mode = []
    for temp in temp_folders:
        if mode.lower() == 'part1':
            if temp[-1] == '1':
                include_mode.append(temp)
                sub_title = "after radiator ON"
        elif mode.lower() == 'part2':
            if temp[-1] == '2':
                include_mode.append(temp)
                sub_title = "during stable env"
        else:
            include_mode.append(temp)
            sub_title = ""

    info_recorder = {} # Ex. 20K-2: {"T1":{"avg","std"},"T2":{"avg","std"},"eff_T":{"avg","std"}}
    for temperature in include_mode:
        temperature_folder = os.path.join(meas_raw_dir,temperature)
        json_folder = os.path.join(temperature_folder,"results","jsons")
        with open(os.path.join(json_folder,"temperatureInfo.json")) as J:
            info_recorder[str(temperature)] = json.load(J)

    plot_item = ["T1","eff_T"] #list(info_recorder[list(info_recorder.keys())[0]].keys())
    fig,ax = plt.subplots(len(plot_item),1,figsize=(15,12))
    for exp_idx in range(len(plot_item)):
        exp_type = plot_item[exp_idx]
        if exp_type != "eff_T":
            y_label = "Time (µs)"
        else:
            y_label = "Effective temperature (mK)"
        if exp_type != "Time":
            for temperature in info_recorder:
                ax[exp_idx].errorbar(info_recorder[temperature]['Time'],info_recorder[temperature][exp_type]["avg"],yerr=info_recorder[temperature][exp_type]["std"],fmt='o-',label=f"{temperature}")
                ax[exp_idx].set_xlabel(f"Time {sub_title} (min)")
                ax[exp_idx].set_ylabel(y_label)
                ax[exp_idx].set_title(f"Temperature comparisom on {exp_type} by Qblox")
            ax[exp_idx].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(temp_compa_path,f"{mode}_temp_compa.png"))
    plt.close()

def plot_time_behavior_sep(temperature_folder:str, mode:str, time_axis:ndarray=array([]), radiator_act:str="ON")->plt.Axes:
    """
    Plot T1/T2/eff_T trend along time axis.\n
    arg `mode` assign which value to plot : 'T1', 'T2' or 'effT' 
    """
    
    tempera_into = collect_allSets_inTempera(temperature_folder)
    a_set_time = 7 # mins
    time_axis_min = arange(1,len(tempera_into["T1"]["avg"]))*a_set_time if time_axis.shape[0]==0 else round(time_axis/60,1)  
    
    ax_info = {"exp":mode, "x":time_axis_min.tolist()}

    fig, ax = plt.subplots(1,1,figsize=(15,10))
    ax:plt.Axes
    if mode.lower() == 't1':
        ax.errorbar(time_axis_min,tempera_into["T1"]["avg"],yerr=tempera_into["T1"]["std"],fmt="o-",color='red')
        ax_info["y"], ax_info["yerr"] = tempera_into["T1"]["avg"], tempera_into["T1"]["std"]
        ax.set_ylabel("$T_{1}$ (µs)",fontsize=26)
        ax.yaxis.label.set_color("red")
    
    elif mode.lower() == 'efft':
        ax.errorbar(time_axis_min,tempera_into["eff_T"]["avg"],yerr=tempera_into["eff_T"]["std"],fmt="o-",color='orange')
        ax_info["y"], ax_info["yerr"] = tempera_into["eff_T"]["avg"], tempera_into["eff_T"]["std"]
        ax.set_ylabel("Effective Temp. (mK)",fontsize=26)
        ax.yaxis.label.set_color('orange')
        
    
    else:
        ax.errorbar(time_axis_min,tempera_into["T2"]["avg"],yerr=tempera_into["T2"]["std"],fmt="o-",color='blue')
        ax_info["y"], ax_info["yerr"] = tempera_into["T2"]["avg"], tempera_into["T2"]["std"]
        ax.set_ylabel("$T_{2}$ (µs)",fontsize=26)
        ax.yaxis.label.set_color('blue')
        
    ax.xaxis.set_tick_params(labelsize=26)
    ax.yaxis.set_tick_params(labelsize=26)
    ax.set_xlabel(f"time after radiator {radiator_act} (min)", fontsize=26)
    
    return ax, ax_info

    

def plot_ref_onAxes(ax:plt.Axes,ref_dict_b4:dict={},ref_dict_recover:dict={}, which:str="T1")->plt.Axes:
    """ #### ref_dict keys should follow the form: `{"T1","T1_sd","T2","T2_sd","effT","effT_sd"}`\n
        the element in dict must be given in value-sd pair, which kinds of pair is optional.\n
        If the dict is empty, we don't plot it on Axes. 
    """
    ax_info = {"exp":f"ref-{which}","x":[],"y":[],"yerr":[]}
    if ref_dict_b4 != {}:
        key_names_b4 = list(ref_dict_b4.keys())
        if ("T1" in key_names_b4) and ("T1_sd" in key_names_b4) and (which.lower() == 't1'):
            ax.axhline(y=ref_dict_b4["T1"],color='#000080',label=f"{which.upper()}-Before turn ON radiator",lw=3.5)
            ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_b4["T1"]-ref_dict_b4["T1_sd"],y2=ref_dict_b4["T1"]+ref_dict_b4["T1_sd"],color='#1E90FF',alpha=0.3)
            ax_info["x"].append(4), ax_info["y"].append(ref_dict_b4["T1"]), ax_info["yerr"].append(ref_dict_b4["T1_sd"])
        if  ("T2" in key_names_b4) and ("T2_sd" in key_names_b4) and (which.lower() == 't2'):
            ax.axhline(y=ref_dict_b4["T2"],color='#000080',label=f"{which.upper()}-Before turn ON radiator",lw=3.5)
            ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_b4["T2"]-ref_dict_b4["T2_sd"],y2=ref_dict_b4["T2"]+ref_dict_b4["T2_sd"],color='#1E90FF',alpha=0.3)
            ax_info["x"].append(4), ax_info["y"].append(ref_dict_b4["T2"]), ax_info["yerr"].append(ref_dict_b4["T2_sd"])
        if ("effT" in key_names_b4) and ("effT_sd" in key_names_b4) and (which.lower() == 'efft'):
            ax.axhline(y=ref_dict_b4["effT"],color='#000080',label=f"{which.upper()}-Before turn ON radiator",lw=3.5)
            ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_b4["effT"]-ref_dict_b4["effT_sd"],y2=ref_dict_b4["effT"]+ref_dict_b4["effT_sd"],color='#1E90FF',alpha=0.3)
            ax_info["x"].append(4), ax_info["y"].append(ref_dict_b4["effT"]), ax_info["yerr"].append(ref_dict_b4["effT_sd"])
    if ref_dict_recover != {}:
        key_names_rc = list(ref_dict_recover.keys())
        if ("T1" in key_names_rc) and ("T1_sd" in key_names_rc) and (which.lower() == 't1'):
            ax.axhline(y=ref_dict_recover["T1"],color='#FA8072',label=f"{which.upper()}-LONG enough After turn OFF radiator",lw=3.5)
            ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_recover["T1"]-ref_dict_recover["T1_sd"],y2=ref_dict_recover["T1"]+ref_dict_recover["T1_sd"],color='#FFA07A',alpha=0.3)
            ax_info["x"].append(4), ax_info["y"].append(ref_dict_recover["T1"]), ax_info["yerr"].append(ref_dict_recover["T1_sd"])
        if  ("T2" in key_names_rc) and ("T2_sd" in key_names_rc) and (which.lower() == 't2'):
            ax.axhline(y=ref_dict_recover["T2"],color='#FA8072',label=f"{which.upper()}-LONG enough After turn OFF radiator",lw=3.5)
            ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_recover["T2"]-ref_dict_recover["T2_sd"],y2=ref_dict_recover["T2"]+ref_dict_recover["T2_sd"],color='#FFA07A',alpha=0.3)
            ax_info["x"].append(4), ax_info["y"].append(ref_dict_recover["T2"]), ax_info["yerr"].append(ref_dict_recover["T2_sd"])
        if ("effT" in key_names_rc) and ("effT_sd" in key_names_rc) and (which.lower() == 'efft'):
            ax.axhline(y=ref_dict_recover["effT"],color='#FA8072',label=f"{which.upper()}-LONG enough After turn OFF radiator",lw=3.5)
            ax.fill_between(ax.lines[0].get_xdata(),y1=ref_dict_recover["effT"]-ref_dict_recover["effT_sd"],y2=ref_dict_recover["effT"]+ref_dict_recover["effT_sd"],color='#FFA07A',alpha=0.3)
            ax_info["x"].append(4), ax_info["y"].append(ref_dict_recover["effT"]), ax_info["yerr"].append(ref_dict_recover["effT_sd"])
    
    return ax, ax_info
    

def plot_DR_tempera(start_date:str="",start_time:str="",temp_chennel:int=6,DR_log_folder_path:str="",ax:plt.Axes=None,time_axis:ndarray=array([]))-> plt.Axes:
    """
    ## Plot only when start_date, start_time, and DR_log_folder_path are ALL given.
    """
    dr_info = {"exp":f"Chennel-{temp_chennel}","x":[],"y":[],"yerr":[]} # no average here
    if ax is None:
        fig, axDR = plt.subplots(1,1,figsize=(12,9))
        axDR.set_xlabel("Time past (min)", fontsize=26)
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
        axDR.set_ylabel(y_name,fontsize=26)
        axDR.set_ylim(min(dr_temp_array)-std(dr_temp_array), max(dr_temp_array)+std(dr_temp_array))
        axDR.xaxis.set_tick_params(labelsize=26)
        axDR.yaxis.set_tick_params(labelsize=26)
        dr_info["x"] = dr_time_array.tolist()
        dr_info["y"] = dr_temp_array.tolist()
        dr_info["yerr"] = zeros(dr_temp_array.shape[0]).tolist()
    return axDR, dr_info

def scat_DR_avg_temp(need_log_info:dict,DR_log_folder_path:str="",ax:plt.Axes=None,temp_chennel:int=6)->tuple[plt.Axes, plt.Axes]:
    """ 
    `need_log_info` follows the form: {temp:{"keep_time_min", "start_date", "start_time", "avg_min_from_the_end"}}
    """
    button:bool = 1
    dr_info = {"exp":f"Chennel-{temp_chennel}"}
    if ax is None :
        plt.close()
        fig, axDR = plt.subplots(1,1,figsize=(12,9))
        axDR.set_xlabel("Radiator temp. (K)", fontsize=26)
        x_axis = []
        x_axis += [int(temp.split("K")[0]) for temp in need_log_info]
   

    dr_temp_avg = []
    dr_temp_std = []
    x_axis = []
    for temperature in need_log_info:
        exp_keep_time_min:int = need_log_info[temperature]["keep_time_min"]
        start_date:str = need_log_info[temperature]["start_date"]
        start_time:str = need_log_info[temperature]["start_time"]
        avg_min_from_the_end:int = need_log_info[temperature]["avg_min_from_the_end"]

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
        axDR.set_ylabel(y_name,fontsize=26)
        axDR.set_ylim(min(dr_temp_avg)-std(dr_temp_avg), max(dr_temp_avg)+std(dr_temp_avg))
        axDR.xaxis.set_tick_params(labelsize=26) 
        axDR.yaxis.set_tick_params(labelsize=26)    
    return axDR, dr_info


def plot_stable_temp_dep(temp_folder_names:list, exp_type:str, slice_min_from_the_end:list=[0])->tuple[plt.Axes, dict, str, dict]:
    """exp_type should be one of ['T1', 'T2', 'effT']"""
    keep_time = {}
    x_axis = []
    exp_value = {"avgs":[],"stds":[]}
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
        exp_value["avgs"].append(mean(array(every_value_dict[f"{exp_type}s"][histo_idx_this_min:])))
        exp_value["stds"].append( std(array(every_value_dict[f"{exp_type}s"][histo_idx_this_min:])))

    fig, ax = plt.subplots(1,1,figsize=(15,10))   
    plt.grid()  
    ax: plt.Axes
    ax.yaxis.label.set_color("red")
    ax.errorbar(x_axis,exp_value["avgs"],yerr=exp_value["stds"],fmt="o-",c='red')
    ax.set_xlabel("Radiator temp. (K)",fontsize=26)
    
    if exp_type in ["T1", "T2"]:
        ax.set_ylabel(f"{exp_type} (us)",fontsize=26)
    else:
        ax.set_ylabel("Effective temp (mK)",fontsize=26)
    
    ax.xaxis.set_tick_params(labelsize=26)
    ax.yaxis.set_tick_params(labelsize=26)


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
        if ("T1" in key_names_b4) and ("T1_sd" in key_names_b4) and (which.lower() == 't1'):
            ax.scatter([4],ref_dict_b4["T1"], s=150, marker='X')
            ax.errorbar([4],ref_dict_b4["T1"],yerr=ref_dict_b4["T1_sd"],fmt="X",label=f"{which.upper()}-Before turn ON radiator")
            x.append(4)
            y.append(ref_dict_b4["T1"])
            yerr.append(ref_dict_b4["T1_sd"])
        if  ("T2" in key_names_b4) and ("T2_sd" in key_names_b4) and (which.lower() == 't2'):
            ax.scatter([4],ref_dict_b4["T2"], s=150, marker='X')
            ax.errorbar([4],ref_dict_b4["T2"],yerr=ref_dict_b4["T2_sd"],fmt="X",label=f"{which.upper()}-Before turn ON radiator")
            x.append(4)
            y.append(ref_dict_b4["T2"])
            yerr.append(ref_dict_b4["T2_sd"])
        if ("effT" in key_names_b4) and ("effT_sd" in key_names_b4) and (which.lower() == 'efft'):
            ax.scatter([4],ref_dict_b4["effT"], s=150, marker='X')
            ax.errorbar([4],ref_dict_b4["effT"],yerr=ref_dict_b4["effT_sd"],fmt="X",label=f"{which.upper()}-Before turn ON radiator")
            x.append(4)
            y.append(ref_dict_b4["effT"])
            yerr.append(ref_dict_b4["effT_sd"])
    if ref_dict_recover != {}:
        key_names_rc = list(ref_dict_recover.keys())
        if ("T1" in key_names_rc) and ("T1_sd" in key_names_rc) and (which.lower() == 't1'):
            ax.scatter([4],ref_dict_recover["T1"], s=150, marker='D')
            ax.errorbar([4],ref_dict_recover["T1"],yerr=ref_dict_recover["T1_sd"],fmt="D",label=f"{which.upper()}-LONG enough After turn OFF radiator")
            x.append(4)
            y.append(ref_dict_recover["T1"])
            yerr.append(ref_dict_recover["T1_sd"])
        if  ("T2" in key_names_rc) and ("T2_sd" in key_names_rc) and (which.lower() == 't2'):
            ax.scatter([4],ref_dict_recover["T2"], s=150, marker='D')
            ax.errorbar([4],ref_dict_recover["T2"],yerr=ref_dict_recover["T2_sd"],fmt="D",label=f"{which.upper()}-LONG enough After turn OFF radiator")
            x.append(4)
            y.append(ref_dict_recover["T2"])
            yerr.append(ref_dict_recover["T2_sd"])
        if ("effT" in key_names_rc) and ("effT_sd" in key_names_rc) and (which.lower() == 'efft'):
            ax.scatter([4],ref_dict_recover["effT"], s=150, marker='D')
            ax.errorbar([4],ref_dict_recover["effT"],yerr=ref_dict_recover["effT_sd"],fmt="D",label=f"{which.upper()}-LONG enough After turn OFF radiator")
            x.append(4)
            y.append(ref_dict_recover["effT"])
            yerr.append(ref_dict_recover["effT_sd"])
    handles, labels = ax.get_legend_handles_labels()
    return ax, handles, labels, {"exp":f"ref-{which}","x":x,"y":y,"yerr":yerr}
    


#                         ==============
#                         =            =
# ========================= user layer ================================
#                         =            =
#                         ==============

def time_trend_artist(tempera_folder:str, exp_catas:list, time_past_sec_array:ndarray, ref_before:dict, ref_recove:dict, DR_time_info:dict, log_folder:str, DRtemp_che:int=6):
    start_date = DR_time_info["start_date"]
    start_time = DR_time_info["start_time"]
    temperature = os.path.split(tempera_folder)[-1]
    if temperature[:2] == "re":
        action = "OFF"
    else:
        action = "ON"
    for exp in exp_catas:
        ax, trend_info = plot_time_behavior_sep(tempera_folder,exp,time_axis=time_past_sec_array,radiator_act=action)
        ax, ref_info = plot_ref_onAxes(ax,ref_before,ref_recove,exp)
        handles, labels = ax.get_legend_handles_labels()
        plt.grid(axis='both')
        ax, dr_info = plot_DR_tempera(start_date,start_time,DRtemp_che,log_folder,ax,time_axis=time_past_sec_array)
        new_folder_path = timelabelfolder_creator(tempera_folder,f"{target_q}_{temperature}TIME_{exp}")
        pic_values_saver(new_folder_path,'time',trend_info, ref_info, dr_info)
        plt.legend(handles, labels, fontsize=30,bbox_to_anchor=(1,1.2))
        plt.tight_layout()
        plt.savefig(os.path.join(new_folder_path,f"{target_q}_{exp}_timeMonitor.png"))
        plt.close()

def temp_depen_artist(temperature_list:list, log_info_dict:dict, exp_catas:list, ref_before:dict={}, ref_recove:dict={}, log_folder="", tempera_che:int=6):
    for exp in exp_catas:
        ax, keep_time_info, conditional_folder_path, axes_info = plot_stable_temp_dep(temperature_list,exp,[log_info_dict[temp]["avg_min_from_the_end"] for temp in temperature_list])
        for tempera in log_info_dict:
            log_info_dict[tempera]["keep_time_min"] = keep_time_info[tempera]
        ax, handles, labels, ref_info = scat_ref_temp_dep(ax,ref_before,ref_recove,exp)
        
        ax, dr_info = scat_DR_avg_temp(log_info_dict,DR_log_folder_path=log_folder,ax=ax,temp_chennel=tempera_che)
        new_folder_path = timelabelfolder_creator(conditional_folder_path,f"{target_q}_TEMP_{exp}")
        pic_values_saver(new_folder_path,'temp',axes_info, ref_info, dr_info)
        plt.legend(handles, labels, fontsize=30,bbox_to_anchor=(1,1.15))
        plt.tight_layout()
        plt.savefig(os.path.join(new_folder_path,f"{target_q}_{exp}_tempDependence.png"))
        plt.close()
        # plt.show()
    
# without shielding 
"""
# after recovering
effT_recov = 73.1
std_effT_recov = 13.4
T1_recov = 10.6
std_T1_recov = 0.6

# before turn on
effT_b4 = 67.01
std_effT_b4 = 0.39
T1_b4 = 11.1
std_T1_b4 = 1
"""




if __name__ == '__main__':
    # // If you wanna plot DR temperature, "start_date" and "start_time" in log_info_dict and also log_folder ALL should be given !
    # *********** Manully settings ***********
    analysis:bool = 0
    plot_time_trend:bool = 1
    always_plot_timeTrend:bool = 1
    plot_temp_dependence: bool = 1
    save_ref:bool = 0

    target_q = 'q0'
    exp_catas = ["effT", "T1"]
    sameple_folder = "Radiator_wisconsinQ1"
    conditional_folder = "Radiator_WS"

    # ?  log_info_dict at least should have the temperature with its 'avg_min_from_the_end'
    # // If there is only one temperature in `log_info_dict`, it will plot time trend. Or u can set 'always_plot_timeTrend = True` to enforce it plot no matter how many temperatures are there.
    # log_info_dict = {"20K-2":{"start_date":"", "start_time":"", "avg_min_from_the_end":60},
    #                  "30K-2":{"start_date":"", "start_time":"", "avg_min_from_the_end":60},
    #                  "40K-2":{"start_date":"", "start_time":"", "avg_min_from_the_end":60},
    #                  "60K-2":{"start_date":"", "start_time":"", "avg_min_from_the_end":60}
    #                 }
    log_info_dict = {#"10K":{"start_date":"2024-05-13", "start_time":"17:25", "avg_min_from_the_end":60},
                    #  "20K":{"start_date":"2024-05-14", "start_time":"10:45", "avg_min_from_the_end":60},
                    #  "40K":{"start_date":"2024-05-14", "start_time":"16:15", "avg_min_from_the_end":60},
                    #  "60K":{"start_date":"2024-05-15", "start_time":"09:15", "avg_min_from_the_end":60},
                     "re0K":{"start_date":"2024-05-15", "start_time":"15:43", "avg_min_from_the_end":60}
                    }
    
    # ? For references.
    # ? (1) If dict is empty = {}, it won't plot that reference 
    # ? (2) avg and SD should be all given but which exp is optional, 
        # * {"T1", "T1_sd"}, {"effT", "effT_sd"} are all aceptable. 
        # ! {"T1"}, {"T1","effT_sd"} are forbidden !
    ref_before = {"T1":35.9,"T1_sd":2,"effT":65.2,"effT_sd":1.8} # WS
    ref_recove = {"T1":29.2,"T1_sd":3.2,"effT":56.3,"effT_sd":3.1}
    # ref_before = {"T1":11.1,"T1_sd":1,"effT":67.01,"effT_sd":0.39} # WOS
    # ref_recove = {"T1":10.6,"T1_sd":0.6,"effT":73.1,"effT_sd":13.4}


    # ? If this folder is "", it won't plot the MXC temperature. 
    log_folder = "/Users/ratiswu/Downloads/DR_temperature_log"




    # ! Don't touch below !
    temperature_list = list(log_info_dict.keys())
    if len(temperature_list) == 1 or always_plot_timeTrend :
        for tempera in temperature_list:
            tempera_folder = os.path.join(meas_raw_dir,sameple_folder,conditional_folder,tempera)
            if analysis:
                main_analysis(target_q, tempera_folder)
            time_past_sec_array = get_time_axis(target_q,tempera_folder)
            # plot time trend with a given temperature (time as x-axis, in min)
            if plot_time_trend:
                time_trend_artist(tempera_folder, exp_catas, time_past_sec_array, ref_before, ref_recove, log_info_dict[tempera], log_folder, DRtemp_che=6)
    else:
        if plot_temp_dependence:
            temp_depen_artist(temperature_list, log_info_dict, exp_catas, ref_before, ref_recove, log_folder, tempera_che=6)


    # save reference 
    if save_ref:
        if conditional_folder == "":
            from datetime import datetime as d
            ref_name = str(d.now().date())
        else:
            ref_name = conditional_folder

        ref_path = os.path.join(meas_raw_dir,sameple_folder,"references.json")
        if os.path.exists(ref_path):
            old_ref = {}
            with open(ref_path) as record_file: 
                old_ref = json.load(record_file)
            
            old_ref[ref_name] = {"ref_before":ref_before,"ref_recove":ref_recove}
            with open(ref_path,'w') as record_file: 
                json.dump(old_ref,record_file)
        else:
            new_ref = {str(ref_name):{"ref_before":ref_before,"ref_recove":ref_recove}}
            with open(ref_path,'w') as record_file: 
                json.dump(new_ref,record_file)

                
                


            
        
        
        
   