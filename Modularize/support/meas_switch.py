"""
This py support you to measure at a different bias, the detail info is recorded in Notebook in QD file.
"""
import os, json, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import QDmanager
from Modularize.support.UserFriend import *


def save_this_meas_option(dr:str,last_ip:str,target_q:str,z_bias:float='sweet'):
    """
    Save the info now in QDmanager with a given z_bias. If z_bias is `'off-sweet'` or `'sweet'`, system will use the tuneaway or sweet bias in Fluxmanager.
    """
    QD_path = find_latest_QD_pkl_for_dr(dr,last_ip)
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    options = QD_agent.Notewriter.get_all_meas_options(target_q)
    for opt_idx, option in enumerate(options):
        slightly_print("=====================================")
        print(f"idx={opt_idx}: {option}")

    chosen_idx = mark_input("Please input an index (idx) to write that corresponding option or input 'n' to cancel:")
    if chosen_idx.lower() not in ['n', 'no']:
        # save
        if z_bias == 'sweet':
            eyeson_print("Z_bias is automatically picked as sweet spot bias in Fluxmanager!")
            QD_agent.keep_meas_option(target_q, QD_agent.Fluxmanager.get_sweetBiasFor(target_q), int(chosen_idx))
        elif z_bias == 'off-sweet':
            eyeson_print("Z_bias is automatically picked as tuneaway spot bias in Fluxmanager!")
            
            QD_agent.keep_meas_option(target_q, QD_agent.Fluxmanager.get_tuneawayBiasFor(target_q), int(chosen_idx))
        else:
            QD_agent.keep_meas_option(target_q, z_bias, int(chosen_idx))
        QD_agent.QD_keeper()
    else:
        print("This deal is got cancelled !")


def switch_a_meas_point(dr:str,last_ip:str,target_q:str):
    QD_path = find_latest_QD_pkl_for_dr(dr,last_ip)
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()
    
    options = QD_agent.Notewriter.get_all_meas_options(target_q)
    for opt_idx, option in enumerate(options):
        slightly_print("=====================================")
        print(f"idx={opt_idx}: {option}")
    chosen_idx = mark_input("Please input an index (idx) to choose that corresponding option or input 'n' to cancel:")
    if chosen_idx.lower() not in ['n','no']:
        eyeson_print(f"your option is : {options[int(chosen_idx)]}")
        permission = mark_input("Sure to keep doing? [y/n]")
        if permission.lower() in ['y','yes']:
            QD_agent.write_with_meas_option(target_q,chosen_idx)
            QD_agent.QD_keeper()
        else:
            print("This deal is got cancelled !")
    else:
            print("This deal is got cancelled !")
    

# NOTE: Before you switch, please save it !!
if __name__ == "__main__":

    """ fill in """
    DRandIP = {"dr":"dr4","last_ip":"81"}
    ro_elements = {
        'q1':{"mode_idx":0,"z_bias_forSave":'sweet'}
    }                                                      #["sweet","off-sweet"] or a float number
 
    """ Running and Storing"""
    for qubit in ro_elements:
        mode = ["save","switch"][ro_elements[qubit]["mode_idx"]]

        if mode == "save":
            save_this_meas_option(DRandIP["dr"],DRandIP["last_ip"],qubit,ro_elements[qubit]["z_bias_forSave"])
        else:
            switch_a_meas_point(DRandIP["dr"],DRandIP["last_ip"],qubit)

    


