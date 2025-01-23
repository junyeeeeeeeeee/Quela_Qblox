from xarray import open_dataset, Dataset, DataArray
import matplotlib.pyplot as plt
from numpy import array, linspace, arange, moveaxis, empty_like
from qblox_drive_AS.SOP.RabiOsci import sort_elements_2_multiples_of
from qcat.visualization.readout_fidelity import plot_readout_fidelity
nc_path:str = "qblox_drive_AS/Meas_raw/20250122/H17M37S50/T1_20250122173836.nc"
ds = open_dataset(nc_path)

shot = 500
repeat = 1
time_pts = 100 
time_samples = sort_elements_2_multiples_of(linspace(0,400e-6,100)*1e9,4)*1e-9
target_q = ["q0","q1","q2","q3"]
print(ds.coords)
print(ds.data_vars)
print("Attributes:\n")
print(ds.attrs)

print("\n\n")
print(array(ds.data_vars['y1']).reshape(repeat,shot,time_pts))

dict_ = {}
for q_idx, q in enumerate(target_q):
    i_data = array(ds[f'y{2*q_idx}']).reshape(repeat,shot,time_pts)
    q_data = array(ds[f'y{2*q_idx+1}']).reshape(repeat,shot,time_pts)
    dict_[q] = (["mixer","prepared_state","repeat","index","time_idx"],array([[i_data],[q_data]]))
    time_values = list(time_samples)*2*repeat*shot
    dict_[f"{q}_x"] = (["mixer","prepared_state","repeat","index","time_idx"],array(time_values).reshape(2,1,repeat,shot,time_pts))

dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":arange(repeat),"prepared_state":array([1]),"index":arange(shot),"time_idx":arange(time_pts)})
dataset.attrs = ds.attrs
dataset.attrs["method"] = "oneshot"
dataset.to_netcdf("qblox_drive_AS/Meas_raw/20250122/H17M37S50/T1_OS_20250122173836.nc")  

from qblox_drive_AS.support.QDmanager import QDmanager

ds = open_dataset("qblox_drive_AS/Meas_raw/20250122/H17M37S50/T1_OS_20250122173836.nc")
md_ds = open_dataset("qblox_drive_AS/Meas_raw/20250122/H16M33S50/SingleShot_20250122163431.nc")
# QD_agent = QDmanager("qblox_drive_AS/QD_backup/20250122/DR1_FQWJ.pkl")
QD_agent = QDmanager("qblox_drive_AS/QD_backup/20250122/DR1_FQWJ_activateMD.pkl")
QD_agent.QD_loader()
print(QD_agent.StateDiscriminator.elements)
# # check model
# for q in md_ds.data_vars:
#     QD_agent.StateDiscriminator.check_model_alive(md_ds[q]*1000,q)


p_rec = {}
for q in ds.data_vars:
    if q.split("_")[-1] != "x":
        
        raw_data = array(ds.data_vars[q])*1000 # Shape ("mixer","prepared_state","repeat","index","time_idx")
        reshaped_data = moveaxis(moveaxis(raw_data,2,0),-1,2)   # raw shape -> ("repeat","mixer","time_idx","prepared_state","index")
        
        md = QD_agent.StateDiscriminator.summon_discriminator(q)
        
        p_rec[q] = []
        for repetition in reshaped_data:
            time_proba = []
                
            da = DataArray(repetition, coords= [("mixer",array(["I","Q"])), ("time_idx",array(ds.coords["time_idx"])), ("prepared_state",array([1])), ("index",array(ds.coords["index"]))] )
            md.discriminator._import_data(da)
            md.discriminator._start_analysis()
            ans = md.discriminator._export_result()
            
            for i in ans:
                p = list(i.reshape(-1)).count(1)/len(list(i.reshape(-1)))
                time_proba.append(p)
            p_rec[q].append(time_proba)

                

from qcat.visualization.qubit_relaxation import plot_qubit_relaxation
from qcat.analysis.qubit.relaxation import qubit_relaxation_fitting
for q in p_rec:
    time_us = array(ds.data_vars[f"{q}_x"])[0][0][0][0]*1e6
    # print(time_us)
    for exci_proba in array(p_rec[q]):
        try:
            ans = qubit_relaxation_fitting(time_us,exci_proba)
            fig, ax = plt.subplots()
            ax = plot_qubit_relaxation(time_us,exci_proba, ax, ans)
            ax.set_title(f"{q} T1 = {round(ans.params['tau'].value,1)} µs" )
            ax.set_ylabel("|1> probability")
            ax.yaxis.minorticks_on()
            ax.grid()
        except:
            plt.scatter(time_us,exci_proba)
            plt.title(q)
        plt.savefig(f"{q}_OST1.png")

