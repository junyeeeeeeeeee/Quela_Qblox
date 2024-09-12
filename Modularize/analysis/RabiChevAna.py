import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import matplotlib.pyplot as plt
from xarray import open_dataset
from numpy import array, linspace, mean
from Modularize.support import QDmanager
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array
from datetime import datetime

def time_label_sort(nc_file_name:str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0],"H%HM%MS%S")


def plot_chevron(QD_agent:QDmanager, target_q:str, nc_path:str, y_is_detune:bool=True):

    ds = open_dataset(nc_path)
    I,Q = dataset_to_array(dataset=ds,dims=2)
    z = IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[target_q][0],ref_Q=QD_agent.refIQ[target_q][-1]).transpose() 
    x = array(ds.x0).reshape(z.shape)[0]
    freqs = array(ds.x1).reshape(z.shape).transpose()[0]

    if max(abs(x)) < 0.01:
        title = ' Time Rabi Chevron'
        x_text = "Driving pulse length (ns)"
        x *= 1e9 
    else:
        x_text = "Driving pulse amplitude (V)"
        title = ' Power Rabi Chevron'

    if y_is_detune:
        y = (freqs - mean(freqs))*1e-6
        y_text = "Driving pulse detuning (MHz)"
    else:
        y = freqs*1e-9
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
    plt.title(title,fontsize=20)
    plt.tight_layout()
    plt.show()

def plot_fringe(QD_agent:QDmanager, target_q:str, nc_path:str, y_is_detune:bool=True):
    ds = open_dataset(nc_path)
    I,Q = dataset_to_array(dataset=ds,dims=2)
    z = IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[target_q][0],ref_Q=QD_agent.refIQ[target_q][-1]).transpose() 
    x = array(ds.x0).reshape(z.shape)[0]
    freqs = array(ds.x1).reshape(z.shape).transpose()[0]

    if y_is_detune:
        y = (freqs - mean(freqs))*1e-6
        y_text = "Driving pulse detuning (MHz)"
    else:
        y = freqs*1e-9
        y_text = "Driving pulse frequency (GHz)"
    
    title = 'Ramsey Fringe'
    x_text = "Free evolution time (us)"
    x *= 1e6

    fig, ax = plt.subplots(figsize=(12,7))
    ax:plt.Axes
    c = ax.pcolormesh(x, y, z, cmap="RdBu")
    ax.set_ylabel(y_text,fontsize=20)
    ax.set_xlabel(x_text,fontsize=20)
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    fig.colorbar(c, ax=ax, label="Contrast (V)")
    plt.title(title,fontsize=20)
    plt.tight_layout()
    plt.show()



if __name__ == "__main__": 
    QD_file = "Modularize/QD_backup/2024_9_11/DR4#81_SumInfo.pkl"
    data_path = 'Modularize/Meas_raw/2024_9_12/DR4q4_RabiChevron_H12M45S34.nc'
    target_q = 'q4'



    QD_agent = QDmanager(QD_file)
    QD_agent.QD_loader()

    plot_chevron(QD_agent,target_q,data_path,True)

