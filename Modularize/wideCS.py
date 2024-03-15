import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
import matplotlib as plt
import scipy.signal
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import Data_manager, QDmanager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.Pulse_schedule_library import One_tone_sche, pulse_preview
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis

def plot_spectrum(start_freq, stop_freq, num_data, freq_sweep_range, I_data, Q_data):
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

    ax1.plot(freq_sweep_range / 1e9, amplitude, color="#00839F", linewidth=2)
    ax1.plot(((cav_freq / num_data)*(stop_freq - start_freq) + start_freq) /1e9, amplitude[cav_freq], "x")
    print(((cav_freq / num_data)*(stop_freq - start_freq) + start_freq) /1e9)
    ax1.set_ylabel("Amplitude (V)")

    ax2.plot(freq_sweep_range / 1e9, phase, color="#00839F", linewidth=2)
    ax2.set_ylabel("Phase ($\circ$)")
    ax2.set_xlabel("Frequency (GHz)")
    fig.tight_layout()
    plt.show()

# Should be simplfied!!!

def select_module_widget(
    device, select_all=False, select_qrm_type: bool = True, select_rf_type: bool = False
):
    '''
    Create a widget to select modules of a certain type

    default is to show only QRM baseband

    Args:
        devices : Cluster we are currently using
        select_all (bool): ignore filters and show all modules
        select_qrm_type (bool): filter QRM/QCM
        select_rf_type (bool): filter RF/baseband
    '''
    options = [[None, None]]

    for module in device.modules:
        if module.present():
            if select_all or (
                module.is_qrm_type == select_qrm_type and module.is_rf_type == select_rf_type
            ):
                options.append(
                    [
                        f"{device.name} "
                        f"{module.short_name} "
                        f"({module.module_type}{'_RF' if module.is_rf_type else ''})",
                        module,
                    ]
                )
    widget = widgets.Dropdown(options=options)
    display(widget)

    return widget


if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down
    from numpy import NaN
    import json
    import Modularize.chip_data_store as cds
    import Modularize.UIwindow as UW
    
    # Variables
    chip_info_restore = True

    # Create or Load chip information
    chip_info = cds.Chip_file()
    
    # Reload the QuantumDevice or build up a new one
    QD_path, dr, ip, mode, vpn = UW.init_meas_window()
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,
                                                        dr_loc=dr,
                                                        cluster_ip=ip,
                                                        mode=mode,
                                                        vpn=vpn)

    # Set the system attenuations
    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),ro_out_att=0)
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

    # Initial value
    num_averages = 10
    integration_length = 1024
    holdoff_length = 200
    waveform_length = integration_length + holdoff_length

    lo_start_freq = 5.2e9
    lo_stop_freq = 5.7e9
    num_data = 501
    nco_freq = 1e6

    # Acquisitions
    acquisitions = {"acq": {"num_bins": 1, "index": 0}}

    # Readout select
    select_readout_module = select_module_widget(cluster, select_qrm_type=True, select_rf_type=True)
    readout_module = select_readout_module.value

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

    lo_data_0 = []
    lo_data_1 = []

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
        lo_data_0.append(data["acquisition"]["bins"]["integration"]["path0"][0] / integration_length)
        lo_data_1.append(data["acquisition"]["bins"]["integration"]["path1"][0] / integration_length)

    # Change data type
    lo_data_0 = np.asarray(lo_data_0)
    lo_data_1 = np.asarray(lo_data_1)

    plot_spectrum(lo_start_freq + nco_freq, lo_stop_freq + nco_freq, num_data, lo_sweep_range + nco_freq, lo_data_0, lo_data_1)

    shut_down(cluster,Fctrl)
    
