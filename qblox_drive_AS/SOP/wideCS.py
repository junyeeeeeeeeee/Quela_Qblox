import os, sys, json
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
import matplotlib.pyplot as plt
from xarray import Dataset, open_dataset
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
    readout_module.sequencers[0].sequence("sequence.json")
    readout_module.disconnect_outputs()
    readout_module.disconnect_inputs()

    # Configure channel map
    readout_module.sequencers[0].connect_sequencer("io0")

    readout_module.sequencers[0].marker_ovr_en(True)
    readout_module.sequencers[0].marker_ovr_value(3)  # Enables output on QRM-RF

    # Set offset in mV
    readout_module.out0_offset_path0(5.5)
    readout_module.out0_offset_path1(5.5)

    # Configure scope mode
    readout_module.scope_acq_sequencer_select(0)
    readout_module.scope_acq_trigger_mode_path0("sequencer")
    readout_module.scope_acq_trigger_mode_path1("sequencer")

    # Configure the sequencer
    readout_module.sequencers[0].mod_en_awg(True)
    readout_module.sequencers[0].demod_en_acq(True)
    readout_module.sequencers[0].nco_freq(nco_freq)
    readout_module.sequencers[0].integration_length_acq(integration_length)
    readout_module.sequencers[0].sync_en(True)

    # NCO delay compensation
    readout_module.sequencers[0].nco_prop_delay_comp_en(True)

    lo_sweep_range = np.linspace(lo_start_freq, lo_stop_freq, num_data)

    I_data = []
    Q_data = []

    for lo_val in lo_sweep_range:
        # Update the LO frequency.
        readout_module.out0_in0_lo_freq(lo_val)

        # Clear acquisitions
        
        readout_module.sequencers[0].delete_acquisition_data("acq")

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

    
    dataset = Dataset({"data":(["mixer","freq"],np.array([I_data,Q_data]))},coords={"mixer":np.array(["I","Q"]),"freq":lo_sweep_range + nco_freq})

    return dataset

def plot_S21(dataset:Dataset,save_path:str=None):
    I_data = dataset.data_vars["data"][0]
    Q_data = dataset.data_vars["data"][1]
    f = dataset.coords["freq"]
    
    amplitude = np.sqrt(I_data**2 + Q_data**2)
    phase = np.arctan2(Q_data, I_data) * 180 / np.pi

    plt.rcParams["axes.labelsize"] = 18
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(15, 7))

    ax1:plt.Axes = axs[0]
    ax1.grid()
    ax1.plot(f / 1e9, amplitude, color="#00839F", linewidth=2)
    ax1.set_ylabel("Amplitude (V)")
    ax2:plt.Axes = axs[1]
    ax2.plot(f / 1e9, np.diff(np.unwrap(phase,180,period=360),append=np.unwrap(phase,180,period=360)[-1]), color="#00839F", linewidth=2)
    ax2.set_ylabel("Phase ($\circ$)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.grid()
    fig.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    ds = open_dataset(r"C:\GitHub\Quela_Qblox\qblox_drive_AS\Meas_raw\20241114\BroadBandCS_20241114215922.nc")
    plot_S21(ds)


    
