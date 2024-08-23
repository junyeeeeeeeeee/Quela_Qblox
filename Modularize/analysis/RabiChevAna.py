import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import matplotlib.pyplot as plt
from xarray import open_dataset
from numpy import array, linspace
from Modularize.support import QDmanager
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array
from datetime import datetime

def time_label_sort(nc_file_name:str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0],"H%HM%MS%S")


def plot_chevron(QD_agent:QDmanager, target_q:str, freq_spn_Hz:float, data_folder:str, y_is_detune:bool=True):

    ncs = sorted([name for name in os.listdir(data_folder) if (os.path.isfile(os.path.join(data_folder,name)) and name.split(".")[-1] == "nc")],key=lambda name:time_label_sort(name))
    fres_adj = linspace(-freq_spn_Hz,freq_spn_Hz,len(ncs))
    
    z = []
    for idx, name in enumerate(ncs):
        file = os.path.join(data_folder, name)
        ds = open_dataset(file)
        x = array(ds.x0)
        I,Q = dataset_to_array(dataset=ds,dims=1)
        z.append(IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[target_q][0],ref_Q=QD_agent.refIQ[target_q][-1]))

    if max(abs(x)) < 0.01:
        x_text = "Driving pulse length (ns)"
        x *= 1e9 
    else:
        x_text = "Driving pulse amplitude (V)"

    if y_is_detune:
        y = fres_adj*1e-6
        y_text = "Driving pulse detuning (MHz)"
    else:
        y = (fres_adj+QD_agent.quantum_device.get_element(target_q).clock_freqs.f01())*1e-9
        y_text = "Driving pulse frequency (GHz)"

    # +QD_agent.quantum_device.get_element(target_q).clock_freqs.f01())*1e-9
    fig, ax = plt.subplots(figsize=(12,7))
    ax:plt.Axes
    c = ax.pcolormesh(x, y, z, cmap="RdBu")
    ax.set_ylabel(y_text,fontsize=20)
    ax.set_xlabel(x_text,fontsize=20)
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    fig.colorbar(c, ax=ax, label="Contrast (V)")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__": 
    QD_file = "Modularize/QD_backup/2024_8_23/DR4#81_SumInfo.pkl"
    data_folder = 'Modularize/Meas_raw/2024_8_23/RabiChevron_q0_H0M42S22'
    target_q = 'q0'
    freq_span = 30e6


    QD_agent = QDmanager(QD_file)
    QD_agent.QD_loader()

    plot_chevron(QD_agent,target_q,freq_span,data_folder,True)

