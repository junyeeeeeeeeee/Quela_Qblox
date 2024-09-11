from numpy import array, arctan2, unwrap, diff, median, ndarray, max, min, average, log10
from numpy.linalg import norm 
import pandas as pd
import matplotlib.pyplot as plt

def plot_S21_fromVNAcsv(csv_file:str,S_chennel:str='S21'):
    dict = pd.read_csv(csv_file,sep=";").to_dict() # ['freq[Hz]', 're:Trc1_S21', 'im:Trc1_S21']
    freq = array(list(dict["freq[Hz]"].values()))*1e-9
    amp_dB =20*log10(abs(array(list(dict[f"re:Trc1_{S_chennel.upper()}"].values())) + 1j*array(list(dict[f"im:Trc1_{S_chennel.upper()}"].values()))))
    pha = array([median(diff(unwrap(arctan2(array(list(dict[f"im:Trc1_{S_chennel.upper()}"].values())), array(list(dict[f"re:Trc1_{S_chennel.upper()}"].values()))))))]+list(diff(unwrap(arctan2(array(list(dict[f"im:Trc1_{S_chennel.upper()}"].values())), array(list(dict[f"re:Trc1_{S_chennel.upper()}"].values())))))))

    fig, ax = plt.subplots(2,1,figsize=(13,8))
    ax0:plt.Axes = ax[0]
    ax0.plot(freq,amp_dB,'red')
    ax0.xaxis.minorticks_on()
    ax0.xaxis.set_tick_params(labelsize=18)
    ax0.yaxis.set_tick_params(labelsize=18)
    ax0.set_ylabel("Amplitude (dB)",fontsize=20)
    ax0.set_xlabel("Frequency (GHz)",fontsize=20)
    ax0.grid()
    ax1:plt.Axes = ax[1]
    ax1.plot(freq,pha,'blue')
    ax1.xaxis.minorticks_on()
    ax1.xaxis.set_tick_params(labelsize=18)
    ax1.yaxis.set_tick_params(labelsize=18)
    ax1.set_ylabel("Phase difference (a.u.)",fontsize=20)
    ax1.set_xlabel("Frequency (GHz)",fontsize=20)
    ax1.grid()
    plt.tight_layout()
    plt.show()


def plot_S21_fromPYQUMcsv(csv_file:str):
    x = pd.read_csv(csv_file).to_dict('list')
    fig, ax = plt.subplots(2,1,figsize=(13,8))
    ax0:plt.Axes=ax[0]
    ax0.plot(x['<b>frequency(GHz)</b>'],x['Amplitude'],'red')
    ax0.set_xlabel("frequency (GHz)",fontsize=20)
    ax0.set_ylabel("Amplitude (dB)",fontsize=20)
    ax0.xaxis.set_tick_params(labelsize=18)
    ax0.yaxis.set_tick_params(labelsize=18)
    ax0.minorticks_on()
    ax0.grid()
    ax1:plt.Axes=ax[1]
    ax1.plot(x['<b>frequency(GHz)</b>'],x['UPhase'],'blue')
    ax1.set_xlabel("frequency (GHz)",fontsize=20)
    ax1.set_ylabel("Phase difference (a.u.)",fontsize=20)
    ax1.xaxis.set_tick_params(labelsize=18)
    ax1.yaxis.set_tick_params(labelsize=18)
    ax1.minorticks_on()
    ax1.grid()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_S21_fromVNAcsv("/Users/ratiswu/Desktop/WJ_S21.csv","S43")
    plot_S21_fromPYQUMcsv("/Users/ratiswu/Downloads/4 to 8.csv")
