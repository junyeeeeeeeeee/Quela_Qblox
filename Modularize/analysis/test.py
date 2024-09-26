import matplotlib.pyplot as plt
import pandas as pd
from numpy import array
import os 

driving_times = array([0,1,2,3,4])
freqs = array([4718.1, 4716.5, 4718.6, 4723.5, 4723.5])
BWs = array([11.8, 11.9, 11.2, 8.2, 8.8])/2

plt.plot(driving_times, freqs, c='blue')
plt.scatter(driving_times,freqs,c='blue')
plt.xlabel("Driving time in log10 base (µs)")
plt.ylabel("Transition freq. (MHz)")
plt.title("Z-pulse 2tone with 100µs reset time")
plt.grid()
plt.tight_layout()
plt.show()