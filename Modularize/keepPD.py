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
        qubit.measure.pulse_amp(ans[q]["dressP"])
        qubit.clock_freqs.readout(ans[q]["dressF_Hz"])
        QDagent.Notewriter.save_bareFreq_for(target_q=q,bare_freq=ans[q]["bareF_Hz"])

    QDagent.refresh_log("PD answer stored!")
    QDagent.QD_keeper()


if __name__ == "__main__":
    qd_path = 'Modularize/QD_backup/2024_3_15/DR2#171_SumInfo.pkl'
    PDans = {"q2":{"dressF_Hz":5.8359e9,"dressP":0.25,"bareF_Hz":5.8353e9}} # "q0":[5.259e9,0.7,5.2589e9],"q1":[5.5278e9,0.1,5.5277e9],"q2":[5.3596e9,0.1,5.3594e9],"q3":[5.6366e9,0.1,5.6365e9]
    fillin_PDans(qd_path, PDans)
