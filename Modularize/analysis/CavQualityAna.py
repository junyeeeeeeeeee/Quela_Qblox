import os, json
import matplotlib.pyplot as plt
from numpy import array, ndarray, cos, sin, deg2rad, imag, real, pi, abs, sqrt
from xarray import Dataset, open_dataset
import pandas as pd
from qcat.analysis.resonator.photon_dep.res_data import PhotonDepResonator
def dBm2photons():
    pass

def timelabel_sort(file_name_list:list)->list:
    import datetime as dt
    def timelable2time(file_name:str):
        time_label = file_name.split("_")[-1].split(".")[0]
        hours = time_label.split("H")[-1].split("M")[0] if len(time_label.split("H")[-1].split("M")[0]) == 2 else f'0{time_label.split("H")[-1].split("M")[0]}'
        minutes=time_label.split("M")[-1].split("S")[0] if len(time_label.split("M")[-1].split("S")[0]) == 2 else f'0{time_label.split("M")[-1].split("S")[0]}'
        seconds=time_label.split("S")[-1] if len(time_label.split("S")[-1]) == 2 else f'0{time_label.split("S")[-1]}'
        return dt.datetime.strptime(f"{hours}:{minutes}:{seconds}","%H:%M:%S")
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
            plt.savefig(os.path.join(parent, f"AllCav_{topic.split('_')[0]}.png"))
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
                plt.savefig(os.path.join(parent, f"{q}_{topic.split('_')[0]}.png"))
                plt.close()

    
def cav_photonDepAna_bridge(folder_path:str):
    other_info = {}
    info_file = [os.path.join(folder_path,name) for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[0]=='Additional_info')][0] 
    with open(info_file) as J:
            other_info = json.load(J)
    ncs = [os.path.join(folder_path,name) for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split("_")[1]=='CavitySpectro')]
    ncs = timelabel_sort(ncs)
    power = other_info["SA_dBm"]
    RT_atte = other_info["RT_atte_dB"]
    ro_elements = other_info["ro_elements"]
    applied_attes = other_info["applied_atte"]

    for q_idx, qubit in enumerate(list(ro_elements.keys())):
        result_folder = os.path.join(folder_path,f"{qubit}_cavResults")
        resonator = PhotonDepResonator(qubit)
        if not os.path.exists(result_folder):
            os.mkdir(result_folder)
        for atte_idx, nc_file in enumerate(ncs):
            ds = open_dataset(nc_file)
            
            S21 = array(ds[f"y{2*q_idx}"] * cos(deg2rad(ds[f"y{2*q_idx+1}"])) + 1j * ds[f"y{2*q_idx}"] * sin(deg2rad(ds[f"y{2*q_idx+1}"])))
            resonator.import_array(array(ro_elements[qubit]), array(S21), float(power[qubit])-float(RT_atte)-float(applied_attes[atte_idx]))
        result = resonator.refined_analysis( result_folder )
        plot_qualities(result_folder)


def plot_qualities(result_folder:str):
    plot_items = ["Qi_dia_corr_fqc", "Qc_dia_corr", "Ql"]
    errors = ["Qi_dia_corr_err", "absQc_err", "Ql_err"]
    x = "photons"
    qubit = os.path.split(result_folder)[-1].split("_")[0]
    csvs = [os.path.join(result_folder,name) for name in os.listdir(result_folder) if (os.path.isfile(os.path.join(result_folder,name)) and name.split(".")[-1]=='csv')]
    for csv in csvs:
        
        fit_name = os.path.split(csv)[-1].split("_")[0]
        csv_dict = pd.read_csv(csv).to_dict()
        if fit_name != "free":
            for plot_idx, plot_topic in enumerate(plot_items):
                x_axis = array(list(csv_dict[x].values()))
                y = array(list(csv_dict[plot_topic].values()))
                yerr = array(list(csv_dict[errors[plot_idx]].values()))

                fig, ax = plt.subplots(figsize=(8,5),dpi=150)
                ax:plt.Axes
                ax.errorbar(x_axis,y,yerr,fmt=".")
                ax.scatter(x_axis,y,s=50)
                ax.set_yscale("log")
                ax.set_xscale("log")
                ax.set_xlabel("Photon number",fontsize=26)
                ax.set_title(f"{qubit}_{fit_name}_{plot_topic[:2]}",fontsize=26)
                ax.xaxis.set_tick_params(labelsize=26)
                ax.yaxis.set_tick_params(labelsize=26)
                plt.grid()
                plt.tight_layout()
                plt.savefig(os.path.join(result_folder,f"{qubit}_{plot_topic[:2]}.png"))
                plt.close()


def a_stte_cavFit(nc_file:str, ro_elements:dict):
    from Modularize.m2_CavitySpec import multiplexing_CS_ana
    ds = open_dataset(nc_file)
    multiplexing_CS_ana(None, ds, ro_elements)

def plot_raw_amp(folder_path:str):
    other_info = {}
    info_file = [os.path.join(folder_path,name) for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[0]=='Additional_info')][0] 
    with open(info_file) as J:
            other_info = json.load(J)
    ncs = [os.path.join(folder_path,name) for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split("_")[1]=='CavitySpectro')]
    ncs = timelabel_sort(ncs)[-2:]
    power = other_info["SA_dBm"]
    RT_atte = other_info["RT_atte_dB"]
    ro_elements = other_info["ro_elements"]
    applied_attes = other_info["applied_atte"][-2:]
    for q_idx, qubit in enumerate(list(ro_elements.keys())):
        fig, ax = plt.subplots(figsize=(8,5),dpi=150)
        ax:plt.Axes
        for idx, nc in enumerate(ncs):
            ds = open_dataset(nc)
            S21 = array(ds[f"y{2*q_idx}"] * cos(deg2rad(ds[f"y{2*q_idx+1}"])) + 1j * ds[f"y{2*q_idx}"] * sin(deg2rad(ds[f"y{2*q_idx+1}"])))
            ax.plot(array(ro_elements[qubit]),sqrt(real(S21)**2+imag(S21)**2),label=f"{round(float(power[qubit])-float(applied_attes[idx]),1)} dBm")
            ax.set_xlabel("Freq. (Hz)",fontsize=26)
            ax.set_title(f"{qubit}_RTatte = {int(RT_atte)}dB",fontsize=26)
            ax.xaxis.set_tick_params(labelsize=26)
            ax.yaxis.set_tick_params(labelsize=26)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path,f"{qubit}_rawCav_special.png"))
        plt.close()


if __name__ == "__main__":
    # CavQua_nc_folder_path = "Modularize/Meas_raw/2024_7_17/Multiplex_CavityQuality_RTatte120dB_H21M35S32"
    # cav_photonDepAna_bridge(CavQua_nc_folder_path)

    # # plot cavity and its fitting
    # officallyPlot_a_givenAtte_cav(CavQua_nc_folder_path,0)

    plot_raw_amp("Modularize/Meas_raw/2024_7_18/Multiplex_CavityQuality_RTatte120dB_H18M22S10")
    
