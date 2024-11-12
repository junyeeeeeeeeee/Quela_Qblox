"""
This program is only for analyzing a series of radiation tests like 0K, 20K 40K 60K and re0K with/without shielding. This didn't support comparison between with and without shielding
"""
import os, sys, time, json
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from pandas import DataFrame as DF
from xarray import Dataset, open_dataset # Dataset.to_dict(SS_ds)
from numpy import array, ndarray, mean, std, round, arange
# from Modularize.support.Pulse_schedule_library import hist_plot
import matplotlib.pyplot as plt
from qblox_drive_AS.support.Path_Book import meas_raw_dir
from qblox_drive_AS.analysis.Radiator.RadiatorSetAna import collect_allSets_inTempera, get_time_axis, timelabelfolder_creator, find_nearest, exp_item_translator, ax_set_y_label
from qblox_drive_AS.analysis.Radiator.RadiatorSetAna import exp_items, label_font_size, tick_num_size, fig_size, legend_font_size, units

exp_colors = {"T2":"#1E90FF",
              "T1":"#D2691E",  
            "effT":"#FFA500",
      "thermalPop":"#BDB76B", 
          "gamma1":"#FF0000",
          "gamma2":"#0000FF", 
        "gammaPhi":"#9400D3"}

# +++++++++++++++++++++++ Functions +++++++++++++++++++++++++++++++++++++++++++++++
def get_conditional_folders(sample_folder_name:str)->list:
    sample_path = os.path.join(meas_raw_dir,sample_folder_name)
    conditional_folder_paths = []
    for name in os.listdir(sample_path):
        if os.path.isdir(os.path.join(sample_path,name)):
            if len(name.split("_")) >= 2:
                if name.split("_")[1] not in list(exp_items.values()):
                    conditional_folder_paths.append(os.path.join(sample_path,name))
            else:
                conditional_folder_paths.append(os.path.join(sample_path,name))
    return conditional_folder_paths

def search_tempera_folder(sample_folder_name:str, tempera:str)->list:
    temp_folder_paths = []
    for conditional_folder in get_conditional_folders(sample_folder_name):
        temp_folder_paths += [os.path.join(conditional_folder,name) for name in os.listdir(conditional_folder) if (os.path.isdir(os.path.join(conditional_folder,name)) and name == tempera)]
    
    return temp_folder_paths

def save_picvalues_Vcompa(folder_path:str,axes_infos:list, exp_type:str, mode:str):
    """ 
    mode indicates the axis in this pic is time with the mode='time', or temperature with the mode='temp'.\n
    exp_type refer to the exp_items in RadiatorSetAna.py
    """
    to_csv_dict = {}
    if exp_type in list(exp_items.values()):
        y_unit = units[exp_type] 
    else:
        y_unit = ""
    x_unit = "min" if mode.lower() in ["time"] else "K"
    for line_info in axes_infos:
        to_csv_dict[f"{line_info['condition']}_x_({x_unit})"] = line_info["x"]
        to_csv_dict[f"{line_info['condition']}_y_({y_unit})"] = line_info["y"]
        to_csv_dict[f"{line_info['condition']}_yerr_({y_unit})"] = line_info["yerr"]


    df = DF.from_dict(to_csv_dict, orient='index').transpose()
    DF.to_csv(df, os.path.join(folder_path,f"{exp_type}_dataValues.csv"))
    

# +++++++++++++++++++++++ plot function inside +++++++++++++++++++++++++++++++++++++++++++++++
def time_monitoring_compa_figure_sameTemp(temperature_folder_list:list, target_q:str, mode:str, radiator_act:str="ON",savefig:bool=True):
    """
    Compare one of T1/T2/eff_T trend along time axis. Up to 6 trends.\n
    arg `mode` assign which value to plot : refer to exp_items in RadiatorSetAna.py
    """
    if mode not in list(exp_items.values()): raise KeyError(f"Un-supported mode={mode} was given!")
    if len(temperature_folder_list) > 6: raise ValueError("The number of given temperature folder should be equal or less than 6")
    colors = ["#1E90FF","#FF6347","#DAA520","#3CB371","#EE82EE","#000000"]
    handles = []
    labels  = []
    
    fig, ax = plt.subplots(1,1,figsize=fig_size)
    a_set_time = 7 # mins
    ax:plt.Axes
    plt.grid()
    temp_name = os.path.split(temperature_folder_list[0])[-1]
    compa_folder = os.path.split(os.path.split(temperature_folder_list[0])[0])[0] # "Modularize/Meas_raw/Radiator_wisconsinQ1" with "Modularize/Meas_raw/Radiator_wisconsinQ1/Radiator_WOS/re0K-1" is the given temperature_folder
    pic_info =[]
    for temperature_folder in temperature_folder_list:
        n = len(pic_info)
        time_axis = get_time_axis(target_q, temperature_folder)
        condition_name = os.path.split(os.path.split(temperature_folder)[0])[-1]  # "Radiator_WOS" with "Modularize/Meas_raw/Radiator_wisconsinQ1/Radiator_WOS/re0K-1" is the given temperature_folder
        tempera_into = collect_allSets_inTempera(temperature_folder)
        time_axis_min = arange(1,len(tempera_into[mode]["avg"]))*a_set_time if time_axis.shape[0]==0 else round(time_axis/60,1)  
        line_info = {"condition":condition_name,"x":time_axis_min.tolist()}
        
        ax.errorbar(time_axis_min,tempera_into[mode]["avg"],yerr=tempera_into[mode]["std"],fmt="o-",color=colors[n],label=condition_name)
        line_info["y"], line_info["yerr"] = tempera_into[mode]["avg"], tempera_into[mode]["std"]

        hs, ls = ax.get_legend_handles_labels()
        labels += [ls[-1]]
        handles += [hs[-1]] 
        pic_info.append(line_info)


    ax = ax_set_y_label(ax,mode,label_font_size)
    
    ax.set_xlabel(f"time after radiator {radiator_act} (min)", fontsize=label_font_size)
    ax.xaxis.set_tick_params(labelsize=tick_num_size)
    ax.yaxis.set_tick_params(labelsize=tick_num_size)
    plt.legend(handles, labels, fontsize=legend_font_size,bbox_to_anchor=(1,1.2))
    plt.tight_layout()
    new_folder_path = timelabelfolder_creator(compa_folder,f"{target_q}_{mode}_{temp_name}_TIME")
    save_picvalues_Vcompa(new_folder_path,pic_info,mode,"time")
    if savefig:
        plt.savefig(os.path.join(new_folder_path,f"{target_q}_{temp_name}_{mode}_comparison.png"))
        plt.close()
    else:
        plt.show()



def plot_conditional_temp_dep(sample_folder_name:str,temp_folder_names:list,target_q:str,exp_type:str, slice_min_from_the_end:list=[0])->tuple[plt.Axes, dict, str, dict]:
    """ exp_type should be one of exp_items in RadiatorSetAna.py """
    if exp_type not in list(exp_items.values()): raise KeyError(f"Un-supported exp_type={exp_type} was given!")
    
    rec = []
    conditionla_folders = get_conditional_folders(sample_folder_name)
    for conditional_folder in conditionla_folders:
        condition = os.path.split(conditional_folder)[-1]
        values = {"condition":condition,"x":[],"y":[],"yerr":[]}
        for idx, temperature in enumerate(temp_folder_names):
            every_value_dict = {}
            tempera_folder = os.path.join(conditional_folder,temperature)
            if os.path.exists(tempera_folder):
                time_past_min_array = get_time_axis(target_q,tempera_folder)/60
                
                set_number = time_past_min_array.shape[0]
                json_file = os.path.join(tempera_folder,"results","jsons","every_values.json")
                with open(json_file) as JJ:
                    every_value_dict = json.load(JJ)
                    values["x"].append(int(temperature.split("K")[0]))

                slice_from_this_min = int(time_past_min_array[-1]-slice_min_from_the_end[idx])
                histo_count_in_set = int(len(every_value_dict[f"{exp_type}s"])/set_number) # data number of a exp in every_value_dict = set_number * histo_count_in_set
                this_min_idx_in_set, _ = find_nearest(time_past_min_array, slice_from_this_min) # this index is equal to the set index
                histo_idx_this_min = (this_min_idx_in_set+1)*histo_count_in_set 
                target_y = array(every_value_dict[f"{exp_type}s"][histo_idx_this_min:])         
                values["y"].append(mean(target_y[target_y != 0]))
                values["yerr"].append(std(target_y[target_y != 0]))
        if len(values["x"]) != 0:
            rec.append(values)

    fig, ax = plt.subplots(1,1,figsize=fig_size)   
    plt.grid()  
    ax: plt.Axes
    for condition_values in rec:
        ax.errorbar(condition_values["x"],condition_values["y"],yerr=condition_values["yerr"],label=condition_values["condition"],fmt='o-')
    ax.set_xlabel("Radiator temp. (K)",fontsize=label_font_size)
    
    ax = ax_set_y_label(ax,exp_type,label_font_size)
    
    ax.xaxis.set_tick_params(labelsize=tick_num_size)
    ax.yaxis.set_tick_params(labelsize=tick_num_size)
    hs, ls = ax.get_legend_handles_labels()
    plt.legend(hs, ls, fontsize=legend_font_size,bbox_to_anchor=(1,1.2))
    plt.tight_layout()
    new_folder_path = timelabelfolder_creator(os.path.join(meas_raw_dir,sample_folder_name),f"{target_q}_{exp_type}_TEMP")
    save_picvalues_Vcompa(new_folder_path,rec,exp_type,"temp")
    plt.savefig(os.path.join(new_folder_path,f"{target_q}_{exp_type}_TEMPcomparison.png"))
    plt.close()

# // TOTest
def gamma_compa_temp_dep(sample_folder_name:str, temp_folder_names:dict,slice_min_from_the_end:list,target_q:str,subtract_ref:bool=True,mode:str="gamma1"):
    """ mode: please check `exp_items` in Modularize.analysis.RadiatorSetAna """
    
    ref_dict = {}
    with open(os.path.join(meas_raw_dir,sample_folder_name,f"{target_q}_references.json")) as J:
        file_dict = json.load(J)
        for condition in file_dict:
            ref_dict[condition] = {}
            ref_dict[condition][mode] = {"avg":0}
            if subtract_ref:
                if mode in ["gamma1", "gamma2", "gammaPhi"]:
                    try:
                        ref_dict[condition][mode]["avg"] = file_dict[condition]["ref_before"][mode]
                    except:
                        print(f"no ref called {mode} in reference.json !")
    
    rec = []
    conditionla_folders = [get_conditional_folders(sample_folder_name)[-1]]
    print(conditionla_folders)
    for conditional_folder in conditionla_folders:
        condition = os.path.split(conditional_folder)[-1]
        values = {"condition":condition,"x":[],"y":[],"yerr":[]}
        for idx, temperature in enumerate(temp_folder_names):
            every_value_dict = {}
            tempera_folder = os.path.join(conditional_folder,temperature)
            if os.path.exists(tempera_folder):
                time_past_min_array = get_time_axis(target_q,tempera_folder)/60
                
                set_number = time_past_min_array.shape[0]
                json_file = os.path.join(tempera_folder,"results","jsons","every_values.json")
                with open(json_file) as JJ:
                    every_value_dict = json.load(JJ)
                    values["x"].append(int(temperature.split("K")[0]))

                slice_from_this_min = int(time_past_min_array[-1]-slice_min_from_the_end[idx])
                histo_count_in_set = int(len(every_value_dict[f"{mode}s"])/set_number) # data number of a exp in every_value_dict = set_number * histo_count_in_set
                this_min_idx_in_set, _ = find_nearest(time_past_min_array, slice_from_this_min) # this index is equal to the set index
                histo_idx_this_min = (this_min_idx_in_set+1)*histo_count_in_set  
                
                target_y = array(every_value_dict[f"{mode}s"][histo_idx_this_min:]) 

                values["y"].append(mean(target_y[target_y != 0]- ref_dict[condition][mode]["avg"]))
                values["yerr"].append(std(target_y[target_y != 0]- ref_dict[condition][mode]["avg"])) 

        if len(values["x"]) != 0:
            rec.append(values)

    fig, ax = plt.subplots(1,1,figsize=fig_size)   
    plt.grid()  
    ax: plt.Axes
    y_for_lim = []
    for condition_values in rec:
        ax.errorbar(condition_values["x"],condition_values["y"],yerr=condition_values["yerr"],label=condition_values["condition"],fmt='o-')
        y_for_lim += condition_values["y"]
    ax.set_xlabel("Radiator temp. (K)",fontsize=label_font_size)
    ax.set_ylim(min(array(y_for_lim))-std(array(y_for_lim)),max(array(y_for_lim))+std(array(y_for_lim)))
    
    if subtract_ref:
        if mode == "gamma1":
            ax.set_ylabel("$\Gamma_{IR}$ in $\Gamma_{1}$ (MHz)",fontsize=label_font_size)
        elif mode == "gamma2":
            ax.set_ylabel("$\Gamma_{IR}$ in $\Gamma_{2}$ (MHz)",fontsize=label_font_size)
        else:
            ax.set_ylabel("$\Gamma_{IR}$ in $\Gamma_{\phi}$ (MHz)",fontsize=label_font_size)
        file_name = f"Delta{mode}_IR"
    else:
        if mode == "gamma1":
            ax.set_ylabel("$\Gamma_{1}$ (MHz)",fontsize=label_font_size)
        elif mode == "gamma2":
            ax.set_ylabel("$\Gamma_{2}$ (MHz)",fontsize=label_font_size)
        else:
            ax.set_ylabel("$\Gamma_{\phi}$ (MHz)",fontsize=label_font_size)

    file_name = mode
    
    ax.xaxis.set_tick_params(labelsize=tick_num_size)
    ax.yaxis.set_tick_params(labelsize=tick_num_size)
    hs, ls = ax.get_legend_handles_labels()
    plt.legend(hs, ls, fontsize=legend_font_size,bbox_to_anchor=(1,1.2))
    plt.tight_layout()
    # new_folder_path = timelabelfolder_creator(os.path.join(meas_raw_dir,sample_folder_name),f"{target_q}_{file_name}_TEMP")
    # save_picvalues_Vcompa(new_folder_path,rec,file_name,"temp")
    # plt.savefig(os.path.join(new_folder_path,f"{target_q}_{file_name}_TEMPcomparison.png"))
    # plt.close()
    plt.show()


#                         ==============
#                         =            =
# ========================= user layer ================================
#                         =            =
#                         ==============


def plot_time_monitor_compa(sameple_folder:str, compa_tempera:list, target_q:str, exp_to_plot:list, savefig:bool=True):
    for temp in compa_tempera:
        temperature2compa = search_tempera_folder(sameple_folder,temp)
        action = "OFF" if temp[:2] == "re" else "ON"
        for exp in exp_item_translator(exp_to_plot):
            time_monitoring_compa_figure_sameTemp(temperature2compa,target_q,exp,action, savefig)


def plot_temp_monitor_compa(sample_folder_name:str,exp_catas:list,temp_names:dict,target_q:str):
    """
    temp_names = {"10K":100}, where 100 is the minutes from the end to average.
    """
    avg_time_minutes = list(temp_names.values())
    for exp in exp_item_translator(exp_catas):
        plot_conditional_temp_dep(sample_folder_name,temp_names,target_q,exp,avg_time_minutes)


def gamma_compa_artist(sample_folder_name:str, temp_folder_names:dict, exp_catas:list, target_q:str, ref:bool=True):
    """
    temp_names = {"10K":100}, where 100 is the minutes from the end to average.
    """
    avg_time_minutes = list(temp_folder_names.values())
    for exp in exp_item_translator(exp_catas):
        gamma_compa_temp_dep(sample_folder_name, temp_folder_names,slice_min_from_the_end=avg_time_minutes,target_q=target_q,subtract_ref=ref,mode=exp)



if __name__ == '__main__':
    target_q:str = 'q0'
    sameple_folder:str = "Radiator_wisconsinQ1"
    compa_tempera:dict = {"10K":60,
                          "20K":60,
                          "30K":60,
                          "40K":60,
                          "60K":60}
    exps:list = ["4"]       # {"1":"T1","2":"T2","3":"effT","4":"gamma1","5":"gamma2","6":"thermalPop","7":"gammaPhi"}
    

    # # plot time trend (time past as x-axis)
    # plot_time_monitor_compa(sameple_folder, ["re0K"], target_q, exps)

    # # plot temperature dependence (radiator temperature as x-axis) 
    # plot_temp_monitor_compa(sameple_folder,exps,compa_tempera,target_q)

    # # plot Gamma temperature dependence (radiator temperature as x-axis) 
    gamma_compa_artist(sameple_folder,compa_tempera,exps,target_q) # test in next run (un-finished)
