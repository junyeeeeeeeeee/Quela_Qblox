import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qcat.analysis.state_discrimination.discriminator import train_GMModel # type: ignore
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from xarray import Dataset, open_dataset
from Modularize.analysis.Radiator.RadiatorSetAna import OSdata_arranger
from Modularize.support.UserFriend import *
from Modularize.support.QDmanager import QDmanager
from numpy import array, moveaxis, mean, std, median, arange
import matplotlib.pyplot as plt
from Modularize.analysis.Radiator.RadiatorSetAna import sort_files

def a_OSdata_analPlot(QD_agent:QDmanager,target_q:str,nc_path:str, plot:bool=True, pic_path:str='', save_pic:bool=False):
    folder = os.path.join(os.path.split(nc_path)[0],'OS_pic')
    if not os.path.exists(folder):
        os.mkdir(folder)
    pic_save_path = os.path.join(folder,os.path.split(nc_path)[1].split(".")[0]) if pic_path == '' else pic_path
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
    dist_model = train_GMModel(tarin_data[0])
    dist_model.relabel_model(array([QD_agent.refIQ[target_q]]).transpose())
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    
    # predict all collection to calculate eff_T for every exp_idx
    
    analysis_data = fit_arrays[0] #your (2,2,N) data to analysis

    new_data = moveaxis( analysis_data ,1,0)
    p0_pop = dist_model.get_state_population(new_data[0].transpose())
    p1_pop = dist_model.get_state_population(new_data[1].transpose())
    if save_pic:
        fig , p01,effT_mK, snr_dB = plot_readout_fidelity(analysis_data, transi_freq, output=pic_save_path,plot=plot,detail_output=True)
    else:
        fig , p01, effT_mK, snr_dB = plot_readout_fidelity(analysis_data, transi_freq, plot=plot, detail_output=True)
    
    plt.close()
    return p01, effT_mK, snr_dB

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
            dist_model = train_GMModel(tarin_data[0])
            dist_model.relabel_model(array([QD_agent.refIQ[target_q]]).transpose())

        analysis_data = fit_arrays[0] #your (2,2,N) data to analysis
        new_data = moveaxis( analysis_data ,1,0)
        p0_pop = dist_model.get_state_population(new_data[0].transpose())
        p1_pop = dist_model.get_state_population(new_data[1].transpose())
        fig , p01, effT_mK, snr_dB = plot_readout_fidelity(analysis_data, transi_freq, output=os.path.join(pic_folder,f"OS_({idx})") if pic_folder is not None else pic_folder, plot=False)
        pop_rec.append(p01*100)
        efft_rec.append(effT_mK)


    return pop_rec, efft_rec
if __name__ == "__main__":
    QD_path = "Modularize/QD_backup/2024_8_1/DR4#81_SumInfo.pkl"
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    pops, effts = share_model_OSana(QD_agent,'q0','Modularize/Meas_raw/2024_8_1',False)
    times = arange(21)*48/60
    plt.plot(times, pops)
    plt.ylim(0,1.5)
    plt.xlabel("time after MXC warm up (min)")
    plt.ylabel("Thermal pops (%)")
    plt.show()

    # folder = 'Modularize/Meas_raw/2024_5_16'
    # files = [os.path.join(folder,name) for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split("_")[1].split("(")[0]=="SingleShot")]
    
    # effts = []
    # x = []
    # for nc in files:
    #     try:
    #         efft, _ = a_OSdata_analPlot(QD_agent,"q0",nc,save_pic=True,plot=False)
    #         x.append(int(nc.split("(")[-1].split(")")[0]))
    #         effts.append(efft)
    #     except:
    #         warning_print(f"idx = {os.path.split(nc)[-1].split('_')[1].split('(')[-1].split(')')[0]} can't plot")
  
    # plt.scatter(x,effts)
    # plt.xlabel("exp index",fontsize=20)
    # plt.ylabel("Effective temp. (mK)",fontsize=20) 
    # plt.xticks(fontsize=20)
    # plt.yticks(fontsize=20)
    # plt.hlines(median(array(effts)),colors="red",xmin=0,xmax=99,label=f'median={round(median(array(effts)),1)}')
    # plt.hlines(median(array(effts))+std(array(effts)),colors="red",linestyles="--",xmin=0,xmax=99)
    # plt.hlines(median(array(effts))-std(array(effts)),colors="red",linestyles="--",xmin=0,xmax=99)
    # plt.title("Effective temp raw data distribution",fontsize=20)
    # plt.legend()
    # plt.show()
    
    # print(std(array(effts)))