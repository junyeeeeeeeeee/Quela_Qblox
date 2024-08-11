from numpy import array, arctan2, unwrap, diff, median, ndarray, max, min, average, log10
from numpy.linalg import norm 
import pandas as pd
import matplotlib.pyplot as plt

def L2(ary:ndarray):
    return ary / norm(ary)

def normalize(ary:ndarray):
    ary_max = max(ary)
    ary_min = min(ary)
    
    return (ary-ary_min)/(ary_max-ary_min)


TJ_dict = pd.read_csv("/Users/ratiswu/Desktop/TJ_S21.csv",sep=";").to_dict() # ['freq[Hz]', 're:Trc1_S21', 'im:Trc1_S21']
WJ_dict = pd.read_csv("/Users/ratiswu/Desktop/WJ_S21.csv",sep=";").to_dict() # ['freq[Hz]', 're:Trc1_S43', 'im:Trc1_S43']
WJ_freq = array(list(WJ_dict["freq[Hz]"].values()))*1e-9
WJ_amp =20*log10(abs(array(list(WJ_dict["re:Trc1_S43"].values())) + 1j*array(list(WJ_dict["im:Trc1_S43"].values()))))
WJ_pha = array([median(diff(unwrap(arctan2(array(list(WJ_dict["im:Trc1_S43"].values())), array(list(WJ_dict["re:Trc1_S43"].values()))))))]+list(diff(unwrap(arctan2(array(list(WJ_dict["im:Trc1_S43"].values())), array(list(WJ_dict["re:Trc1_S43"].values())))))))
TJ_amp = 20*log10(abs(array(list(TJ_dict["re:Trc1_S21"].values())) + 1j*array(list(TJ_dict["im:Trc1_S21"].values()))))
TJ_pha = array([0]+list(diff(unwrap(arctan2(array(list(TJ_dict["im:Trc1_S21"].values())), array(list(TJ_dict["re:Trc1_S21"].values())))))))

fig, ax = plt.subplots(2,1,figsize=(13,8),sharex=True)
ax0:plt.Axes = ax[0]
# ax0Pha:plt.Axes = ax0.twinx()
# ax0Pha.spines['right'].set_color("blue")
# ax0Pha.yaxis.label.set_color("blue")
ax1:plt.Axes = ax[1]
# ax1Pha:plt.Axes = ax1.twinx()
# ax1Pha.spines['right'].set_color("blue")
# ax1Pha.yaxis.label.set_color("blue")

ax0.set_title("WJ",fontsize=20)
ax0.plot(WJ_freq,normalize(WJ_pha),'blue',label='phase')
# ax0Pha.plot(WJ_freq,normalize(WJ_pha),'blue',label='phase')
ax1.plot(WJ_freq,normalize(TJ_pha),'blue',label='phase')
# ax1Pha.plot(WJ_freq,normalize(TJ_pha),'blue',label='phase')
ax0.grid()
ax0.legend(fontsize=16)
# ax0Pha.legend(loc='lower left',fontsize=16)
ax1.set_title("TJ",fontsize=20)
ax1.grid()
ax1.legend(fontsize=16)
# ax1Pha.legend(loc='lower left',fontsize=16)
ax1.set_xlabel("Frequency (GHz)",fontsize=20)

for ax in [ax0, ax1]: #[ax0, ax0Pha, ax1, ax1Pha]
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=20)
plt.tight_layout()
plt.savefig("/Users/ratiswu/Desktop/WTJ_S21_pha.png")
plt.close()