def fillin_PDans(QD_path:str,ans:dict):
    """
    Fill in the power dep answer to the quantum_device.\n
    format:\n
    `ans = {"q0":[ROf(Hz), ROamp],...}`
    """
    import pickle
    with open(QD_path, 'rb') as inp:
        gift = pickle.load(inp)
    quantum_device = gift["QD"]
    fluxDict = gift["Flux"]
    log = gift["Log"]
    Hcfg = gift["Hcfg"]

    for q in ans:
        qubit = quantum_device.get_element(q)
        qubit.measure.pulse_amp(ans[q][-1])
        qubit.clock_freqs.readout(ans[q][0])

    merged_file = {"QD":quantum_device,"Flux":fluxDict,"Hcfg":Hcfg,"Log":log}
    with open(QD_path, 'wb') as file:
        pickle.dump(merged_file, file)
    print("PD answer filled in successfully!")


if __name__ == "__main__":
    qd_path = ''
    PDans = {"q0":[5.72e9,0.2]}
    fillin_PDans(qd_path, PDans)
