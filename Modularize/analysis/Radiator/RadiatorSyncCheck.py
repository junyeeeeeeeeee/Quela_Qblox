"""
When RadiatorSet program is running and we didn't set a timer, we can use this to plot the time trend to check it in the mean time.
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from Modularize.analysis.Radiator.RadiatorSetAna import main_analysis, time_trend_artist, get_time_axis
from Modularize.support.Path_Book import meas_raw_dir

target_q:str = 'q0'
exp_catas:list = ["effT"]
sample_folder_name:str = "Radiator_wisconsinQ1"
conditional_folder_name:str = "Radiator_WS" 
log_info_dict:dict = {"10K":{}}

for temperature in log_info_dict:
    tempera_folder = os.path.join(meas_raw_dir,sample_folder_name,conditional_folder_name,temperature)
    main_analysis(target_q, tempera_folder,mode='jump')
    time_past_sec_array = get_time_axis(target_q,tempera_folder)
    # plot time trend with a given temperature (time as x-axis, in min)

    time_trend_artist(tempera_folder, target_q, exp_catas, time_past_sec_array, {}, {}, log_info_dict[temperature], "", DRtemp_che="",show=True)