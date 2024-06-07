import os, sys, json
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from qblox_instruments import Cluster

def wideCS(readout_module:Cluster, lo_start_freq:int, lo_stop_freq:int, num_data:int):

    nco_freq = 1e6

    num_averages = 10
    integration_length = 1024
    holdoff_length = 200
    waveform_length = integration_length + holdoff_length

    # Acquisitions
    acquisitions = {"acq": {"num_bins": 1, "index": 0}}

    # Sequence program
    seq_prog = f"""
        move    {num_averages},R0           # Average iterator.
        nop
        reset_ph
        set_awg_offs 10000, 10000          # set amplitude of signal
        nop
    loop: 
        wait     {holdoff_length}          # Wait time of flight
        acquire  0,0,{integration_length}  # Acquire data and store them in bin_n0 of acq_index.
        loop     R0,@loop                  # Run until number of average iterations is done.
        stop                               # Stop the sequencer
        """

    # Add sequence to single dictionary and write to JSON file.
    sequence = {
        "waveforms": {},
        "weights": {},
        "acquisitions": acquisitions,
        "program": seq_prog,
    }
    with open("sequence.json", "w", encoding="utf-8") as file:
        json.dump(sequence, file, indent=4)
        file.close()

    # Upload sequence
    readout_module.sequencer0.sequence("sequence.json")
    readout_module.disconnect_outputs()
    readout_module.disconnect_inputs()

    # Configure channel map
    readout_module.sequencer0.connect_sequencer("io0")

    readout_module.sequencer0.marker_ovr_en(True)
    readout_module.sequencer0.marker_ovr_value(3)  # Enables output on QRM-RF

    # Set offset in mV
    readout_module.out0_offset_path0(5.5)
    readout_module.out0_offset_path1(5.5)

    # Configure scope mode
    readout_module.scope_acq_sequencer_select(0)
    readout_module.scope_acq_trigger_mode_path0("sequencer")
    readout_module.scope_acq_trigger_mode_path1("sequencer")

    # Configure the sequencer
    readout_module.sequencer0.mod_en_awg(True)
    readout_module.sequencer0.demod_en_acq(True)
    readout_module.sequencer0.nco_freq(nco_freq)
    readout_module.sequencer0.integration_length_acq(integration_length)
    readout_module.sequencer0.sync_en(True)

    # NCO delay compensation
    readout_module.sequencer0.nco_prop_delay_comp_en(True)

    lo_sweep_range = np.linspace(lo_start_freq, lo_stop_freq, num_data)

    I_data = []
    Q_data = []

    for lo_val in lo_sweep_range:
        # Update the LO frequency.
        readout_module.out0_in0_lo_freq(lo_val)

        # Clear acquisitions
        readout_module.sequencer0.delete_acquisition_data("acq")

        readout_module.arm_sequencer(0)
        readout_module.start_sequencer()

        # Wait for the sequencer to stop with a timeout period of one minute.
        readout_module.get_acquisition_state(0, timeout=1)

        # Move acquisition data from temporary memory to acquisition list.
        readout_module.store_scope_acquisition(0, "acq")

        # Get acquisition list from instrument.
        data = readout_module.get_acquisitions(0)["acq"]

        # Store the acquisition data.
        # The result still needs to be divided by the integration length to make sure
        # the units are correct.
        I_data.append(data["acquisition"]["bins"]["integration"]["path0"][0] / integration_length)
        Q_data.append(data["acquisition"]["bins"]["integration"]["path1"][0] / integration_length)

    # Change data type
    I_data = np.asarray(I_data)
    Q_data = np.asarray(Q_data)

    amplitude = np.sqrt(I_data**2 + Q_data**2)
    phase = np.arctan2(Q_data, I_data) * 180 / np.pi

    mean_amp = np.mean(amplitude)
    print(mean_amp)

    plt.rcParams["axes.labelsize"] = 18
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16

    fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, figsize=(15, 7))

    cav_freq = scipy.signal.argrelmin(amplitude, order = round(num_data/8))
    cav_freq = cav_freq[0]
    for i in cav_freq:
        if amplitude[i] > mean_amp*2/3:        # Discard the local minimum which is caused by noise
            cav_freq = np.delete(cav_freq, np.where(cav_freq == i))

    ax1.plot((lo_sweep_range + nco_freq) / 1e9, amplitude, color="#00839F", linewidth=2)
    ax1.plot(((cav_freq / num_data)*(lo_stop_freq - lo_start_freq) + lo_start_freq + nco_freq) /1e9, amplitude[cav_freq], "x")
    print(((cav_freq / num_data)*(lo_stop_freq - lo_start_freq) + lo_start_freq + nco_freq) /1e9)
    ax1.set_ylabel("Amplitude (V)")

    ax2.plot((lo_sweep_range + nco_freq) / 1e9, np.diff(np.unwrap(phase,180,period=360),append=np.unwrap(phase,180,period=360)[-1]), color="#00839F", linewidth=2)
    ax2.set_ylabel("Phase ($\circ$)")
    ax2.set_xlabel("Frequency (GHz)")
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, QRM_nco_init
    from Modularize.support.UI_Window import init_meas_window
    
    """ Fill in """
    QD_path, dr, ip, mode = "", "dr3", "192.168.1.13","n" #init_meas_window()
    qrmRF_slot_idx:int  = 6
    lo_start_freq:float = 5.4  * 1e9
    lo_stop_freq:float = 6.1   * 1e9
    num_data:int = 2100


    """ Preparations """
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,
                                                        dr_loc=dr,
                                                        cluster_ip=ip,
                                                        mode=mode,
                                                        qubit_number=5)
    # Set the system attenuations
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=0)
    QRM_nco_init(cluster)
    # Readout select
    readout_module = cluster.modules[qrmRF_slot_idx-1]


    """ Running """
    wideCS(readout_module=readout_module, lo_start_freq=lo_start_freq, lo_stop_freq=lo_stop_freq, num_data=num_data)


    """ Close """
    shut_down(cluster,Fctrl)
    
