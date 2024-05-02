"""
This py support you to measure at a different bias, the detail info is recorded in Notebook in QD file.
"""
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import QDmanager


def save_this_meas_option(dr:str,last_ip:str,target_q:str,z_bias:float=None):
    """
    Save the info now in QDmanager with a given z_bias. If option_idx is `'0'` or `'sweet'`, system will use the sweet bias in Fluxmanager.
    """
    QD_path = find_latest_QD_pkl_for_dr(dr,last_ip)
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    options = QD_agent.Notewriter.get_all_meas_options(target_q)
    for opt_idx, option in enumerate(options):
        print(f"idx={opt_idx}: {option}")
    chosen_idx = input("Please input an index (idx) to write that corresponding option or input 'n' to cancel:")
    if chosen_idx.lower() not in ['n', 'no']:
        # save
        if z_bias is None:
            print("Z_bias is automatically picked as sweet spot bias in Fluxmanager!")
            QD_agent.keep_meas_option(target_q, QD_agent.Fluxmanager.get_sweetBiasFor(target_q), int(chosen_idx))
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
        print(f"idx={opt_idx}: {option}")
    chosen_idx = input("Please input an index (idx) to choose that corresponding option or input 'n' to cancel:")
    if chosen_idx.lower() not in ['n','no']:
        print(f"your option is : {options[int(chosen_idx)]}")
        permission = input("Sure to keep doing? [y/n]")
        if permission.lower() in ['y','yes']:
            QD_agent.write_with_meas_option(target_q,chosen_idx)
            QD_agent.QD_keeper()
        else:
            print("This deal is got cancelled !")
    else:
            print("This deal is got cancelled !")
    
    
if __name__ == "__main__":

    """ fill in """
    DRandIP = {"dr":"dr1","last_ip":"11"}
    ro_elements = {
        'q0':{"mode_idx":1,"z_bias_if_needed":None}   #["save","switch"] # save sweet spot bias with z_bias=None
    }

    """ Running and Storing"""
    for qubit in ro_elements:
        mode = ["save","switch"][ro_elements[qubit]["mode_idx"]]

        if mode == "save":
            save_this_meas_option(DRandIP["dr"],DRandIP["last_ip"],qubit,ro_elements[qubit]["z_bias_if_needed"])
        else:
            switch_a_meas_point(DRandIP["dr"],DRandIP["last_ip"],qubit)

    


