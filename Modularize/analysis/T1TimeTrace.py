import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import Dataset, open_dataset
from Modularize.support import QDmanager,Data_manager
import matplotlib.pyplot as plt
from numpy import array 
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis
from datetime import datetime 
def time_label_sort(nc_file_name:str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0],"H%HM%MS%S")

folder = "Modularize/Meas_raw/spin_t2"
QD_file_path = "Modularize/QD_backup/2024_8_29/DR4#81_SumInfo.pkl"
qs = ['q1']

files = sorted([name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split(".")[-1] == "nc")],key=lambda name:time_label_sort(name))
QD_agent = QDmanager(QD_file_path)
QD_agent.QD_loader()
time = []
t1 = []

for idx, file in enumerate(files) :
    path = os.path.join(folder,file)
    nc = open_dataset(path)
    samples = array(nc.x0)
    I,Q= dataset_to_array(dataset=nc,dims=1)
    data= IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[qs[0]][0],ref_Q=QD_agent.refIQ[qs[0]][-1])
    try:
        data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=25e-6)
        if data_fit.attrs['T1_fit']*1e6 > 10:
            t1.append(data_fit.attrs['T1_fit']*1e6)
        time.append(24.5*(idx+1))
    except:
        plt.plot(data)
        plt.show()

Data_manager().save_histo_pic(QD_agent,{str(qs[0]):t1},qs[0],mode="t2")