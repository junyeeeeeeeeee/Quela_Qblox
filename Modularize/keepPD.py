from Modularize.support import QDmanager

def fillin_PDans(QD_path:str,ans:dict):
    """
    Fill in the power dep answer to the quantum_device.\n
    format:\n
    `ans = {"q0":[ROf(Hz), ROamp, bareFreq(Hz)],...}`
    """
    QDagent = QDmanager(QD_path)
    QDagent.QD_loader()
    for q in ans:
        qubit = QDagent.quantum_device.get_element(q)
        qubit.measure.pulse_amp(ans[q][1])
        qubit.clock_freqs.readout(ans[q][0])
        QDagent.Notewriter.save_bareFreq_for(target_q=q,bare_freq=ans[q][-1])

    QDagent.refresh_log("PD answer stored!")
    QDagent.QD_keeper()


if __name__ == "__main__":
    qd_path = 'Modularize/QD_backup/2024_2_27/SumInfo.pkl'
    PDans = {"q0":[5.72e9,0.15, 5.7188e9],"q1":[6.008e9,0.15,6.0073e9],"q2":[5.8385e9,0.15,5.8372e9],"q3":[6.107e9,0.2,6.1062e9]}
    fillin_PDans(qd_path, PDans)
