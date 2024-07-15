import os, json
import matplotlib.pyplot as plt
from numpy import array, ndarray, cos, sin, deg2rad, imag, real, pi, abs
from xarray import Dataset, open_dataset

def dBm2photons():
    pass

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

def plot_all_cav(result_js_path:str,specific_qs:list=[],allCavInOneFig:bool=False):
    """
    plot all the cavity along dBm axis. plot item includes 'qi', 'qc' and 'ql'.
    """
    parent = os.path.split(result_js_path)[0]
    results = {}
    with open(os.path.join(result_js_path)) as J:
        results = json.load(J)
    ro_attes = list(results.keys())
    if len(specific_qs) == 0:
        specific_qs = list(results[ro_attes[0]]["output_dBm"].keys())
    
    # collects values
    collections = {}
    for q in specific_qs:
        collections[q] = {'Qi_dia_corr':[], 'Qi_dia_corr_err':[], 'Qc_dia_corr':[], 'absQc_err':[], 'Ql':[], 'Ql_err':[], "output_dBm":[]}
        for atte in ro_attes:
            for item in collections[q]:
                collections[q][item].append(results[atte][q][item])
    
    # plotting
    plot_item = {"Qi_dia_corr":"Internal Q (Qi)", "Ql":"Loaded Q (Ql)", "Qc_dia_corr":"Coupling Q (Qc)"}
    if allCavInOneFig:
        for topic in plot_item:
            fig, ax = plt.subplots(dpi=720)
            ax:plt.Axes
            for q in collections.keys():
                if f"{topic}_err" in collections[q]:
                    ax.errorbar(x=collections[q]["output_dBm"], y=collections[q][topic], yerr=collections[q][f"{topic}_err"],label=q)
                else:
                    ax.errorbar(x=collections[q]["output_dBm"], y=collections[q][topic], yerr=collections[q]["absQc_err"],label=q)
            ax.set_ylabel(plot_item[topic],fontsize=26)
            ax.set_xlabel("Power (dBm)",fontsize=26)
            ax.xaxis.set_tick_params(labelsize=26)
            ax.yaxis.set_tick_params(labelsize=26)
            plt.savefig(os.path.join(parent, f"AllCav_{topic.split("_")[0]}.png"))
            plt.close()
    else:
        for topic in plot_item:
            for q in collections.keys():
                fig, ax = plt.subplots(dpi=720)
                ax:plt.Axes
                if f"{topic}_err" in collections[q]:
                    ax.errorbar(x=collections[q]["output_dBm"], y=collections[q][topic], yerr=collections[q][f"{topic}_err"],label=q)
                else:
                    ax.errorbar(x=collections[q]["output_dBm"], y=collections[q][topic], yerr=collections[q]["absQc_err"],label=q)
                ax.set_ylabel(plot_item[topic],fontsize=26)
                ax.set_xlabel("Power (dBm)",fontsize=26)
                ax.xaxis.set_tick_params(labelsize=26)
                ax.yaxis.set_tick_params(labelsize=26)
                plt.savefig(os.path.join(parent, f"{q}_{topic.split("_")[0]}.png"))
                plt.close()

    
    





if __name__ == "__main__":
    CavQua_nc_folder_path = "Modularize/Meas_raw/2024_5_30/q2_CavityQuality_RTatte0dB_H11M45S55"

    # # plot quality factors
    plot_quality_results_for_aQ(CavQua_nc_folder_path, sep_qua=False)

    # # plot cavity and its fitting
    # officallyPlot_a_givenAtte_cav(CavQua_nc_folder_path,0)

