from pandas import read_csv, concat
from pandas import DataFrame as DF
from numpy import array, asarray, ndarray, min, max, std
from datetime import datetime as dt
from datetime import timedelta as td
import matplotlib.pyplot as plt
import json, os

def file_name_sort(file_path_list:list):
    def get_date(path:str):
        return dt.strptime(f"20{os.path.split(path)[-1].split('.')[0].split(' ')[-1]}","%Y-%m-%d")
    
    file_path_list.sort(key=get_date)

    return file_path_list

def how_long_calc(start:dict,end:dict)->dict:
    """
    start: {"date":"2024-05-13","time":"18:06"}.\n
    end:   {name:{"date":"2024-05-13","time":"23:06"}, ...}
    """
    start_dt = dt.strptime(f'{start["date"]} {start["time"]}:00',"%Y-%m-%d %H:%M:%S")
    total_min = {}
    for name in end:
        end_dt = dt.strptime(f'{end[name]["date"]} {end[name]["time"]}:00',"%Y-%m-%d %H:%M:%S")
        total_min[name] = (end_dt-start_dt).total_seconds() / 60

    return total_min

    
def mins_2_time(mins:float):
    hr = str(mins//60) if mins//60 > 9 else f"0{mins//60}"
    min = str(mins%60) if mins%60 > 9  else f"0{mins%60}"
    return f"{hr}:{min}"

def date_2_logdate(date):
    date = str(date.date())
    year = date.split("-")[0]
    mont = date.split("-")[1]
    day  = date.split("-")[-1]
    return f"{year}-{mont}-{day}"

def find_nearest(ary, value:float):
    """ find the element  which is closest to the given target_value in the given array"""
    ary = asarray(ary)
    idx = (abs(ary - value)).argmin()
    return idx, ary[idx]

def current_day_filter(log:dict,the_day:str)->tuple[int, dict]:
    """ The day : 2024-05-14 """
    df = DF.from_dict(log)
    
    wantted_day = f"{the_day[2:].split('-')[-1]}-{the_day[2:].split('-')[1]}-{the_day[2:].split('-')[0]}"
    
    the_day_df = df[df['date'].str.contains(wantted_day)]
    original_row_idx = df.index[df['date'] == " "+wantted_day].tolist()[0]
    
    return original_row_idx, the_day_df.to_dict()


def find_closest_time(log:dict,time:str, date:str):
    """ time in hr:min,\n
        date in 2024-05-15
    """
    
    original_idx, log = current_day_filter(log,date)
    time_col = log["time"]
    dt_col = []
    for i in time_col:
        dt_col.append(dt.strptime(time_col[i],"%H:%M:%S"))

    idx, colsest_time = find_nearest(array(dt_col),dt.strptime(time,"%H:%M:%S"))
    
    
    return original_idx+idx+1, colsest_time
    

def time_calc(start:str,keep:str,date:str)->tuple[str, str]:
    st = dt.strptime(date+" "+start, '%Y-%m-%d %H:%M')
    kp = td(hours=int(keep.split(":")[0]),minutes=int(keep.split(":")[-1]))
    end_time = st+kp
    
    return str(end_time.date()), str(end_time.time())

def log_arranger(folder_path:str,chennel:int=6,mode:str="T",to_json:bool=True,enforce_refresh:bool=False)->dict:
    """
    mode 'T' for temperature, 'P' for pressure.
    """
    if mode.upper() == "T":
        data_unit = "Kelvin"
    else:
        data_unit = "Bar"


    if os.path.exists(os.path.join(folder_path,f"CH{chennel}_{mode.upper()}_SumInfo.json")) and not enforce_refresh:
        with open(os.path.join(folder_path,f"CH{chennel}_{mode.upper()}_SumInfo.json")) as JJ:
            log = json.load(JJ)
    else:
        files = file_name_sort([os.path.join(folder_path,name) for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[-1]=='log' and name.split(" ")[0]==f'CH{chennel}' and name.split(" ")[1]==mode.upper())])
        log = read_csv(files[0],names=["date","time",data_unit])
        for file in files[1:]:
            log = concat([log,read_csv(file,names=["date","time",data_unit])])
        log = log.to_dict(orient='list')
        
        if to_json:
            with open(os.path.join(folder_path,f"CH{chennel}_{mode.upper()}_SumInfo.json"),'w') as rec:
                json.dump(log,rec)
    
    return log

# def reduce_to_min(time_array:ndarray)->ndarray:
#     date_time_array = []
#     for i in time_array:
#         date_time_array.append(dt.strptime(i,"%H:%M:%S"))
    
#     min_array = []
#     for i in date_time_array:
#         min_array.append((i-date_time_array[0]).seconds)
#     return array(min_array)/60

def reduce_to_min(time_array:ndarray, date_array:ndarray)->ndarray:
    datetimes = []
    for idx, value in enumerate(time_array):
        date = f'{date_array[idx].split("-")[-1]}-{date_array[idx].split("-")[1]}-{date_array[idx].split("-")[0][1:]}'
        datetimes.append(dt.strptime(f"20{date} {value}","%Y-%m-%d %H:%M:%S"))
    past_mins = []
    for a_dt in datetimes:
        past_mins.append(round((a_dt-datetimes[0]).total_seconds() / 60,1))

    return array(past_mins)
    

def Kelvin_collector(folder_path:str,start_date:str,start_time:str,execu_time_min:int,chennel:int)-> tuple[ndarray, ndarray]:
    """
    Giving a folder contains the log files from DR, with the wantted start date, start time, execution time (min) and temperature chennel,\n
    return the time array in mins and temperature array.\n
    arg Examples:
    1) start_date: 2024-05-15\n
    2) start_time: 10:20 (in 24 hr unit)\n
    3) execu_time_min: 600 (mins)\n
    4) tempera_chennel: 6 (MXC plate) or 2 (4K plate)
    """
    
    log = log_arranger(folder_path,chennel,mode="T")
    keep_time = mins_2_time(execu_time_min)
    
    date, time = time_calc(start_time,keep_time,start_date)
    
    row_s, _ = find_closest_time(log,f"{start_time}:00",start_date)
    row_e, _ = find_closest_time(log,time,date)

    temp_array = array(list(log["Kelvin"])[row_s:row_e])
    time_array = reduce_to_min(array(list(log["time"])[row_s:row_e]),array(list(log["date"])[row_s:row_e]))

    return time_array, temp_array

 
def Pressure_collector(folder_path:str,start_date:str,start_time:str,execu_time_min:int,chennel:int)-> tuple[ndarray, ndarray]:
    """
    Giving a folder contains the log files from DR, with the wantted start date, start time, execution time (min) and pressure chennel,\n
    return the time array in mins and pressure array.\n
    arg Examples:
    1) start_date: 2024-05-15\n
    2) start_time: 10:20 (in 24 hr unit)\n
    3) execu_time_min: 600 (mins)\n
    4) pressure_chennel: 6 (MXC plate) or 2 (4K plate)
    """
    
    log = log_arranger(folder_path,chennel,mode="P")
    keep_time = mins_2_time(execu_time_min)
    
    date, time = time_calc(start_time,keep_time,start_date)
    
    row_s, _ = find_closest_time(log,f"{start_time}:00",start_date)
    row_e, _ = find_closest_time(log,time,date)

    pres_array = array(list(log["Bar"])[row_s:row_e])
    time_array = reduce_to_min(array(list(log["time"])[row_s:row_e]),array(list(log["date"])[row_s:row_e]))

    return time_array, pres_array   


def plot_temperaChens_trend(DR_log_folder_path:str,start_date:str,start_time:str,monitor_time_min:int, chennels:list=[6], pic_name:str="", label_vlines:dict={}):
    """
    label_vlines is a dict follows the form: {"label name":time(minutes)}
    """
    if len(chennels) > 6: raise ValueError("Given chennels should be equal or less than 6")
    colors = ["#1E90FF","#FF6347","#DAA520","#3CB371","#EE82EE","#000000"]
    fig, ax = plt.subplots(1,1,figsize=(18,12))
    plt.grid()
    ax:plt.Axes
    ax_mk = ax.twinx()
    ax_mk:plt.Axes
    n = 0
    temp_beyond_K = []
    handles = []
    labels  = []

    for temp_chen in chennels:
        dr_time_array, dr_temp_array = Kelvin_collector(DR_log_folder_path,start_date,start_time,monitor_time_min,temp_chen)
        if str(temp_chen) == "6" :
            dr_temp_array *= 1000
            ax_mk.plot(dr_time_array, dr_temp_array, c=colors[n],lw=3, label='MXC')
            ax_mk.set_ylabel("MXC temp. (mK)",fontsize=26)
            ax_mk.set_ylim(min(dr_temp_array)-std(dr_temp_array), max(dr_temp_array)+std(dr_temp_array))
            ax_mk.xaxis.set_tick_params(labelsize=26)
            ax_mk.yaxis.set_tick_params(labelsize=26)
            ax_mk.spines['right'].set_color(colors[n])
            ax_mk.yaxis.label.set_color(colors[n])
            hs, ls = ax_mk.get_legend_handles_labels()
            
        else:
            if str(temp_chen) == '2':
                lable='4K-plate'
            else:
                lable=f'Chennel={temp_chen}'
            ax.plot(dr_time_array, dr_temp_array, c=colors[n],lw=3, label=lable)
            temp_beyond_K += dr_temp_array.flatten().tolist()
            hs, ls = ax.get_legend_handles_labels()
    
        labels += [ls[-1]]
        handles += [hs[-1]]    
        n += 1

    
    ax.xaxis.set_tick_params(labelsize=26)
    ax.yaxis.set_tick_params(labelsize=26)
    ax.set_ylim(min(array(temp_beyond_K))-std(array(temp_beyond_K)),max(array(temp_beyond_K))+std(array(temp_beyond_K)))
    ax.set_xlabel("time past (min)",fontsize=26)
    ax.set_ylabel("Temp. (K)",fontsize=26)
    plt.legend(handles, labels, fontsize=20,bbox_to_anchor=(1,1.15))

    if label_vlines != {}:
        vline_colors = ["#1E90FF","#FF6347","#DAA520","#3CB371","#EE82EE","#000000"]
        H, L = [], []
        for name in label_vlines:
            ax.vlines([float(label_vlines[name])],ymin=min(array(temp_beyond_K))-std(array(temp_beyond_K)),ymax=max(array(temp_beyond_K))+std(array(temp_beyond_K)),linestyles="--",lw=3,label=name,colors=vline_colors[len(H)])
            h, l = ax.get_legend_handles_labels()  
            H.append(h[-1])
            L.append(l[-1])

        ax.legend(H, L, fontsize=20, ncol=int(len(L)/2),bbox_to_anchor=(0.4,1.15))
    plt.tight_layout()
    if pic_name == "":
        plt.show()
    else:
        plt.savefig(os.path.join(DR_log_folder_path,f"{pic_name}.png"))
        plt.close()


if __name__ == "__main__":
    log_folder = "/Users/ratiswu/Downloads/DR_temperature_log_WOS"
    start_date = "2024-04-25"
    start_time = "10:00" # hour:min, 24hr unit
    keep_time_mins = int(6*24*60) # mins
    temp_chennel = 2
    
    # time_array,temp_array = Kelvin_collector(log_folder,start_date,start_time,keep_time_mins,temp_chennel)
    
    plot_temperaChens_trend(log_folder,start_date,start_time,keep_time_mins,chennels=[2,6],label_vlines=how_long_calc({"date":start_date,"time":start_time},
                                                                                                                      {"10K_on":{"date":"2024-04-25","time":"16:10"},
                                                                                                                       "20K_on":{"date":"2024-04-26","time":"17:50"},
                                                                                                                       "30K_on":{"date":"2024-04-29","time":"12:15"},
                                                                                                                       "40K_on":{"date":"2024-04-29","time":"16:40"},
                                                                                                                       "60K_on":{"date":"2024-04-30","time":"11:45"},
                                                                                                                       "60K_off":{"date":"2024-04-30","time":"23:30"}}),pic_name="")
    
    
    
    




