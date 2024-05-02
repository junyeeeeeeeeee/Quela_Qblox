from qcat.state_discrimination.discriminator import train_model # type: ignore
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from xarray import Dataset, open_dataset
from Modularize.analysis.RadiatorSetAna import OSdata_arranger
import os
from Modularize.support.QDmanager import QDmanager
from numpy import array, moveaxis, mean, std
import matplotlib.pyplot as plt

def a_OSdata_analPlot(QD_agent:QDmanager,target_q:str,nc_path:str, plot:bool=True, pic_path:str=''):
    folder = os.path.join(os.path.split(nc_path)[0],'pic')
    pic_save_path = os.path.join(folder,f'{os.path.split(nc_path)[1].split(".")[0]}.png') if pic_path == '' else pic_path
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
    dist_model = train_model(tarin_data[0])
    dist_model.relabel_model(array([QD_agent.refIQ[target_q]]).transpose())
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    
    # predict all collection to calculate eff_T for every exp_idx
    
    analysis_data = fit_arrays[0] #your (2,2,N) data to analysis

    new_data = moveaxis( analysis_data ,1,0)
    p0_pop = dist_model.get_state_population(new_data[0].transpose())
    p1_pop = dist_model.get_state_population(new_data[1].transpose())
    
    fig , effT_mK, snr_dB = plot_readout_fidelity(analysis_data, transi_freq, output=pic_save_path,plot=plot)
    
    plt.close()
    return effT_mK, snr_dB

if __name__ == "__main__":
    pass