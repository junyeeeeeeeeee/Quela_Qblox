import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os 
# def e_fn(x, a, tau, c):
#     return a*np.exp(-x/tau) + c

# file = "/Users/ratiswu/Downloads/A/ScalinQ/q0_re0KTIME_effT_H22M49S19/effT_dataValues.csv"
# df = pd.read_csv(file, sep=',')
# x = np.array(df["exp_x_(min)"])[:100]
# y = np.array(df["exp_y_(mK)"])[:100]

# p, e = curve_fit(e_fn, x, y)
# fig = plt.figure(figsize=(9,7))
# plt.scatter(x,y, label='data')
# plt.plot(x,e_fn(x,*p),c='red', label='fit')
# plt.xlabel("time past (min)", fontsize=26)
# plt.xticks(fontsize=20)
# plt.ylabel("Eff. Temp. (mK)", fontsize=26)
# plt.yticks(fontsize=20)
# plt.title("Temp monitor after off 80mK heater", fontsize=20)
# plt.legend(fontsize=20)
# plt.grid()
# plt.savefig(os.path.join(os.path.split(file)[0], "fit.png"))
# plt.close()
# print(p)

from numpy import arange, ndarray, array

n = 13
t = 172
spin= 3

gap = (t) // n + ((t // n) %(4*(spin+1))) # multiple by 8 ns
samples = arange((4*(spin+1)),t,gap)

def modify_time_point(ary:ndarray,factor:int):
    x = []
    for i in ary:
        if i % factor == 0:
            if i not in x :
                x.append(i)
            else:
                pass
        else:
            multiples = i // factor
            multiples_of_factor = factor*(multiples+1)
            if multiples_of_factor not in x:
                x.append(multiples_of_factor)
            else:
                pass
    return array(x)

good_samples = modify_time_point(samples, spin+1)
print(samples)
print(good_samples)
        


    