import os, json
import matplotlib.pyplot as plt
from numpy import array, ndarray, cos, sin, deg2rad, imag, real, pi, abs
from xarray import Dataset, open_dataset
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis


    
def official_cavfit(ds_path:str, new_folder_additional_name:str=""):
    
    from quantify_core.data.handling import set_datadir
    ds_folder = os.path.split(ds_path)[0]
    new_folder_path = os.path.join(ds_folder,f"officialfit_results_{new_folder_additional_name}")
    if not os.path.isdir(new_folder_path):
        os.mkdir(new_folder_path)
    set_datadir(new_folder_path)
    rs_ds = open_dataset(ds_path)
    x = ResonatorSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
    


def timelabel_sort(file_name_list:list)->list:
    import datetime as dt
    def timelable2time(file_name:str):
        time_label = file_name.split("_")[-1].split(".")[0]
        hours = time_label.split("H")[-1].split("M")[0]
        minutes=time_label.split("M")[-1].split("S")[0]
        seconds=time_label.split("S")[-1]
        return dt.datetime.strptime(f"{hours}:{minutes}:{seconds}","%H:%m:%S")
    file_name_list.sort(key=lambda name: timelable2time(name))  
    return file_name_list

def find_nearest(ary:ndarray, value:float):
    """ find the element  which is closest to the given target_value in the given array"""
    idx = (abs(ary - value)).argmin()
    return idx, ary[idx]

def get_S21_andX(ds:Dataset)->dict:
    S21 = ds.y0 * cos(deg2rad(ds.y1)) + 1j * ds.y0 * sin(deg2rad(ds.y1))
    return {"x":array(ds.x0),"S21":S21,"mag":array(ds.y0),"pha":array(ds.y1)}


def officallyPlot_a_givenAtte_cav(CavQua_nc_folder_path:str, which_dBm:int=0):
    target_q = os.path.split(CavQua_nc_folder_path)[-1].split("_")[0]
    RT_atte_dB = os.path.split(CavQua_nc_folder_path)[-1].split("_")[-2].split("e")[-1].split("d")[0]
    quality_rec_json = [os.path.join(CavQua_nc_folder_path,name) for name in os.listdir(CavQua_nc_folder_path) if (os.path.isfile(os.path.join(CavQua_nc_folder_path,name)) and name.split(".")[0]=='Quality_results')][0]
    raw_ncs = [name for name in os.listdir(CavQua_nc_folder_path) if (os.path.isfile(os.path.join(CavQua_nc_folder_path,name)) and name.split("_")[1]=='CavitySpectro')]
    sorted_raw_ncs = timelabel_sort(raw_ncs)
    quality_rec_dict = {}
    with open(quality_rec_json) as JJ:
        quality_rec_dict = json.load(JJ)
    all_dBm = quality_rec_dict[target_q]["dBm"]
    
    idx, nearest_dBm = find_nearest(array(all_dBm), which_dBm)
    target_nc_path = os.path.join(CavQua_nc_folder_path,sorted_raw_ncs[idx])
    data = get_S21_andX(open_dataset(target_nc_path)) # {"x","S21","mag","pha"}
    official_cavfit(target_nc_path, f"power{nearest_dBm}dB")
    
def plot_quality_results_for_aQ(CavQua_nc_folder_path:str, sep_qua:bool=False):
    colors = ["#0000FF","#FF0000","#008000"]
    markers = ["o", "D", "X"]
    target_q = os.path.split(CavQua_nc_folder_path)[-1].split("_")[0]
    RT_atte_dB = os.path.split(CavQua_nc_folder_path)[-1].split("_")[-2].split("e")[-1].split("d")[0]
    quality_rec_json = [os.path.join(CavQua_nc_folder_path,name) for name in os.listdir(CavQua_nc_folder_path) if (os.path.isfile(os.path.join(CavQua_nc_folder_path,name)) and name.split(".")[0]=='Quality_results')][0]
    
    quality_rec_dict = {}
    with open(quality_rec_json) as JJ:
        quality_rec_dict = json.load(JJ)
    
    x_axis_dBm = quality_rec_dict[target_q]["dBm"]
    if not sep_qua:
        fig, ax = plt.subplots(1,1,figsize=(15,10))
        ax:plt.Axes
    n = 0
    for qualities_name in quality_rec_dict[target_q]:
        if qualities_name != "dBm" and qualities_name.split("_")[-1] != "sd":
            if not sep_qua:
                ax.errorbar(x_axis_dBm,quality_rec_dict[target_q][qualities_name],yerr=quality_rec_dict[target_q][f"{qualities_name}_sd"],fmt=markers[n],c=colors[n],label=qualities_name)
                ax.scatter(x_axis_dBm,quality_rec_dict[target_q][qualities_name],s=60,c=colors[n],marker=markers[n])
               
            else:
                fig, ax = plt.subplots(1,1,figsize=(15,10))
                ax:plt.Axes
                ax.errorbar(x_axis_dBm,quality_rec_dict[target_q][qualities_name],yerr=quality_rec_dict[target_q][f"{qualities_name}_sd"],fmt=markers[n],c=colors[n])
                ax.scatter(x_axis_dBm,quality_rec_dict[target_q][qualities_name],s=60,c=colors[n],marker=markers[n])
                ax.set_xlabel("Input power (dBm)",fontsize=26)
                ax.set_ylabel(qualities_name,fontsize=26)
                ax.xaxis.set_tick_params(labelsize=26)
                ax.yaxis.set_tick_params(labelsize=26)
                ax.set_yscale("log")
                plt.grid()
                plt.title(f"{qualities_name} for {target_q}, RTatte= {RT_atte_dB} dB",fontsize=30)
                plt.tight_layout()
                plt.savefig(os.path.join(CavQua_nc_folder_path, f"{target_q}_{qualities_name}_RTatte{RT_atte_dB}dB.png"))
                plt.close()
            n += 1 

    if not sep_qua:
        ax.set_yscale("log")
        ax.set_xlabel("Input power (dBm)",fontsize=26)
        ax.set_ylabel("Qualities",fontsize=26)
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        plt.grid()
        plt.legend(fontsize=30,ncol=3)
        plt.title(f"Qualities for {target_q}, RTatte= {RT_atte_dB} dB",fontsize=30)
        plt.tight_layout()
        plt.savefig(os.path.join(CavQua_nc_folder_path, f"{target_q}_qualities_RTatte{RT_atte_dB}dB.png"))
        plt.close()


def plot_allCav_quality():
    pass


if __name__ == "__main__":
    CavQua_nc_folder_path = "Modularize/Meas_raw/2024_5_30/q2_CavityQuality_RTatte0dB_H11M45S55"

    # # plot quality factors
    plot_quality_results_for_aQ(CavQua_nc_folder_path, sep_qua=False)

    # # plot cavity and its fitting
    # officallyPlot_a_givenAtte_cav(CavQua_nc_folder_path,0)

