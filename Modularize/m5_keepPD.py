import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.support import QDmanager, cds
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr

def fillin_PDans(QD_path:str,ans:dict):
    """
    Fill in the power dep answer to the quantum_device.\n
    format:\n
    `ans = {"q0":{"dressF_Hz":,"dressP":,"bareF_Hz":},...}`
    """
    QDagent = QDmanager(QD_path)
    QDagent.QD_loader()
    for q in ans:
        qubit = QDagent.quantum_device.get_element(q)
        if ans[q]["dressP"] != "": qubit.measure.pulse_amp(ans[q]["dressP"]) 
        if ans[q]["dressF_Hz"] != "": qubit.clock_freqs.readout(ans[q]["dressF_Hz"])
        if ans[q]["bareF_Hz"] != "": QDagent.Notewriter.save_bareFreq_for(target_q=q,bare_freq=ans[q]["bareF_Hz"])
        if ans[q]["ro_atte"] != "": QDagent.Notewriter.save_DigiAtte_For(atte_dB=ans[q]["ro_atte"],target_q=q,mode='ro')

    QDagent.refresh_log("PD answer stored!")
    QDagent.QD_keeper()
    return QDagent


if __name__ == "__main__":
    """ Fill in """
    DRandIP = {"dr":"dr3","last_ip":"13"}
    PDans = {
        "q0":{"dressF_Hz":5979326009.709601,"dressP":0.01,"bareF_Hz":5977463041.694341,"ro_atte":30},
        "q1":{"dressF_Hz":6092282375.162393,"dressP":0.01,"bareF_Hz":6088746679.042166,"ro_atte":30},
        "q2":{"dressF_Hz":5921075416.364415,"dressP":0.01,"bareF_Hz":5919998796.780729,"ro_atte":30},
        "q3":{"dressF_Hz":6102516437.684792,"dressP":0.01,"bareF_Hz":6099496282.311073,"ro_atte":30},
        "q4":{"dressF_Hz":6013562459.208041,"dressP":0.01,"bareF_Hz":6010520418.566977,"ro_atte":30},
    }
    # 1 = Store
    # 0 = not store
    chip_info_restore = 1
    
    
    """ Storing """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent = fillin_PDans(QD_path, PDans)
    chip_info = cds.Chip_file(QD_agent=QD_agent)
    if chip_info_restore:
        chip_info.update_PDans(result=PDans)
    
    
