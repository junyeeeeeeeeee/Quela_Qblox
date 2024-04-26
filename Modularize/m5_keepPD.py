import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
<<<<<<< HEAD

=======
>>>>>>> origin/RatisWu
from Modularize.support import QDmanager

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


if __name__ == "__main__":
    qd_path = "Modularize/QD_backup/2024_4_25/DR2#10_SumInfo.pkl"
    PDans = {"q1":{"dressF_Hz":6.014388e9,"dressP":0.01,"bareF_Hz":6.01318e9,"ro_atte":34}} # "q0":[5.259e9,0.7,5.2589e9],"q1":[5.5278e9,0.1,5.5277e9],"q2":[5.3596e9,0.1,5.3594e9],"q3":[5.6366e9,0.1,5.6365e9]
    fillin_PDans(qd_path, PDans)
    
    
