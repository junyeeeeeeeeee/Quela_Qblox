import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from xarray import Dataset, open_dataset, DataArray
from qblox_drive_AS.analysis.Radiator.RadiatorSetAna import OSdata_arranger
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support.QDmanager import QDmanager
from numpy import array, moveaxis, mean, std, median, arange,ndarray, pi, sin, cos, vstack, arctan2, column_stack, hstack, empty_like
import matplotlib.pyplot as plt
from qblox_drive_AS.analysis.Radiator.RadiatorSetAna import sort_files
from qcat.analysis.state_discrimination import p01_to_Teff
from qblox_drive_AS.support import rotate_onto_Inphase, rotate_data

def a_OSdata_analPlot(nc_path:str, QD_agent:QDmanager=None, target_q:str=None, plot:bool=True, pic_path:str='', save_pic:bool=False): # 
    folder = os.path.join(os.path.split(nc_path)[0],'OS_pic')
    if not os.path.exists(folder):
        os.mkdir(folder)
    if save_pic:
        pic_save_path = os.path.join(folder,os.path.split(nc_path)[1].split(".")[0]) if pic_path == '' else pic_path
    else:
        pic_save_path = None
    SS_ds = open_dataset(nc_path)
    ss_dict = Dataset.to_dict(SS_ds)
    # print(ss_dict)
    pe_I, pe_Q = ss_dict['data_vars']['e']['data']
    pg_I, pg_Q = ss_dict['data_vars']['g']['data']
    pgI_collection = [pg_I]
    pgQ_collection = [pg_Q]
    peI_collection = [pe_I]
    peQ_collection = [pe_Q]

    OS_data = 1000*array([[pgI_collection,peI_collection],[pgQ_collection,peQ_collection]]) # can train or predict 2*2*histo_counts*shot
    tarin_data, fit_arrays = OSdata_arranger(OS_data)
    # train GMM
    gmm2d_fidelity = GMMROFidelity()
    gmm2d_fidelity._import_data(tarin_data[0])
    gmm2d_fidelity._start_analysis()
    g1d_fidelity = gmm2d_fidelity.export_G1DROFidelity()
    if QD_agent is not None:
        transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    else: 
        transi_freq = None
    
   
    _, angle = rotate_onto_Inphase(gmm2d_fidelity.centers[0],gmm2d_fidelity.centers[1])
    
    p00 = g1d_fidelity.g1d_dist[0][0][0]
    p01 = g1d_fidelity.g1d_dist[0][0][1]
    p11 = g1d_fidelity.g1d_dist[1][0][1]
    if transi_freq is not None:
        effT_mK = p01_to_Teff(p01, transi_freq)*1000
    RO_fidelity_percentage = (p00+p11)*100/2
    if plot:
        z = moveaxis(array(tarin_data[0]),0,1) # (IQ, state, shots) -> (state, IQ, shots)
        rotated_data = empty_like(z)
        for state_idx, state_data in enumerate(z):
            rotated_data[state_idx] = rotate_data(state_data,angle)
    
        da = DataArray(moveaxis(rotated_data,0,1), coords= [("mixer",["I","Q"]), ("prepared_state",[0,1]), ("index",arange(array(tarin_data[0]).shape[2]))] )
        gmm2d_fidelity._import_data(da)
        gmm2d_fidelity._start_analysis()
        g1d_fidelity = gmm2d_fidelity.export_G1DROFidelity()
        print(f"!! center:\n{gmm2d_fidelity.centers}")
        plot_readout_fidelity(da, gmm2d_fidelity, g1d_fidelity,transi_freq,pic_save_path)
        plt.close()

    return p01, effT_mK, RO_fidelity_percentage, angle

def share_model_OSana(QD_agent:QDmanager,target_q:str,folder_path:str,pic_save:bool=True):
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    files = [name for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split("_")[1].split("(")[0]=="SingleShot")]
    files = [os.path.join(folder_path,name) for name in sort_files(files)][:21]
    pop_rec, efft_rec = [], []
    if pic_save:
        pic_folder = os.path.join(folder_path,"OS_detail_pic")
        if not os.path.exists(pic_folder):
            os.mkdir(pic_folder)
    else:
        pic_folder = None

    for idx, file in enumerate(files):
        SS_ds = open_dataset(file)
        ss_dict = Dataset.to_dict(SS_ds)
        # print(ss_dict)
        pe_I, pe_Q = ss_dict['data_vars']['e']['data']
        pg_I, pg_Q = ss_dict['data_vars']['g']['data']
        pgI_collection = [pg_I]
        pgQ_collection = [pg_Q]
        peI_collection = [pe_I]
        peQ_collection = [pe_Q]

        OS_data = 1000*array([[pgI_collection,peI_collection],[pgQ_collection,peQ_collection]]) # can train or predict 2*2*histo_counts*shot
        tarin_data, fit_arrays = OSdata_arranger(OS_data)
        # train GMM
        if idx == 0:
            gmm2d_fidelity = GMMROFidelity()
            gmm2d_fidelity._import_data(tarin_data[0])
            gmm2d_fidelity._start_analysis()

        gmm2d_fidelity.discriminator._import_data( fit_arrays[0] )
        gmm2d_fidelity.discriminator._start_analysis()


    return pop_rec, efft_rec





if __name__ == "__main__":
    nc_path = 'Modularize/Meas_raw/20241029/DR2q0_SingleShot(0)_H12M46S03.nc'
    a_OSdata_analPlot(nc_path)