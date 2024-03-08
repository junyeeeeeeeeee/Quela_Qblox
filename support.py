import warnings
from pathlib import Path
import matplotlib.pyplot as plt
import plotly.io as pio
import ipywidgets as widgets
import numpy as np
import quantify_core.data.handling as dh
from IPython.display import display
from qblox_instruments import Cluster, ClusterType, PlugAndPlay
from qcodes import Instrument
from qcodes.parameters import ManualParameter
from quantify_core.analysis.single_qubit_timedomain import (
    RabiAnalysis,
    RamseyAnalysis,
    T1Analysis,
)
from quantify_core.measurement.control import MeasurementControl
from quantify_core.visualization.pyqt_plotmon import PlotMonitor_pyqt as PlotMonitor
from quantify_scheduler.device_under_test.quantum_device import QuantumDevice
from quantify_scheduler.device_under_test.transmon_element import BasicTransmonElement
from quantify_scheduler.gettables import ScheduleGettable
from quantify_scheduler.instrument_coordinator import InstrumentCoordinator
from quantify_scheduler.instrument_coordinator.components.qblox import ClusterComponent
from quantify_scheduler.operations.gate_library import Measure, Reset
from quantify_scheduler.operations.pulse_library import SetClockFrequency, SquarePulse
from quantify_scheduler.resources import ClockResource
from quantify_scheduler.schedules import heterodyne_spec_sched_nco, rabi_sched, t1_sched
from quantify_scheduler.schedules.timedomain_schedules import ramsey_sched
from quantify_scheduler.schedules.schedule import Schedule
from quantify_scheduler.backends.graph_compilation import SerialCompiler
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis, QubitSpectroscopyAnalysis
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from utils.tutorial_analysis_classes import (
    QubitFluxSpectroscopyAnalysis,
    ResonatorFluxSpectroscopyAnalysis,
)
from utils.tutorial_utils import (
    set_drive_attenuation,
    set_readout_attenuation,
    show_args,
    show_drive_args,
    show_readout_args,
)

# close all instruments
def shut_down(cluster:Cluster,flux_map:dict):
    '''
        Disconnect all the instruments.
    '''
    reset_offset(flux_map)
    cluster.reset() 
    Instrument.close_all() 

# connect to clusters
def connect_clusters():
    with PlugAndPlay() as p:            # Scan for available devices and display
        device_list = p.list_devices()  # Get info of all devices

    names = {dev_id: dev_info["description"]["name"] for dev_id, dev_info in device_list.items()}
    ip_addresses = {dev_id: dev_info["identity"]["ip"] for dev_id, dev_info in device_list.items()}

    connect_options = widgets.Dropdown(         # Create widget for names and ip addresses *** Should be change into other interface in the future
        options=[(f"{names[dev_id]} @{ip_addresses[dev_id]}", dev_id) for dev_id in device_list],
        description="Select Device",
    )
    display(connect_options)
    return connect_options, ip_addresses

# Create quantum device with given q number
def create_quantum_device(HARDWARE_CONFIG:dict,num_qubits : int) -> QuantumDevice:
    """
        Create the QuantumDevice with the given qubit number and the hardware config.
    """
    quantum_device = QuantumDevice("academia_sinica_device")
    quantum_device.hardware_config(HARDWARE_CONFIG)
    
    # store references
    quantum_device._device_elements = list()

    for i in range(num_qubits):
        qubit = BasicTransmonElement(f"q{i}")
        qubit.measure.acq_channel(i)
        quantum_device.add_element(qubit)
        quantum_device._device_elements.append(qubit)
    
    return quantum_device

# Configure_measurement_control_loop
def configure_measurement_control_loop(
    device: QuantumDevice, cluster: Cluster, live_plotting: bool = False
) -> (MeasurementControl, InstrumentCoordinator):
    # Close QCoDeS instruments with conflicting names
    for name in [
        "PlotMonitor",
        "meas_ctrl",
        "ic",
        "ic_generic",
        f"ic_{cluster.name}",
    ] + [f"ic_{module.name}" for module in cluster.modules]:
        try:
            Instrument.find_instrument(name).close()
        except KeyError as kerr:
            pass

    meas_ctrl = MeasurementControl("meas_ctrl")
    ic = InstrumentCoordinator("ic")
    ic.timeout(60*60)

    # Add cluster to instrument coordinator
    ic_cluster = ClusterComponent(cluster)
    ic.add_component(ic_cluster)

    if live_plotting:
        # Associate plot monitor with measurement controller
        plotmon = PlotMonitor("PlotMonitor")
        meas_ctrl.instr_plotmon(plotmon.name)

    # Associate measurement controller and instrument coordinator with the quantum device
    device.instr_measurement_control(meas_ctrl.name)
    device.instr_instrument_coordinator(ic.name)

    return (meas_ctrl, ic)

def qubits_meas(quantum_device:QuantumDevice, qubit_want_to_meas:list)->list:
    ro_elements = []
    for q in range(5):
        if qubit_want_to_meas[q] == 1:
            ro_elements.append(quantum_device._device_elements[q].name)
        else:
            pass
    return ro_elements

def pulse_preview(quantum_device:QuantumDevice,sche_func:Schedule, sche_kwargs:dict, **kwargs):
    
    pio.renderers.default='browser'
    device_compiler = SerialCompiler("Device compiler", quantum_device)
    comp_sched = device_compiler.compile(
        sche_func(**sche_kwargs)
    )
    sche_func(**sche_kwargs).plot_circuit_diagram()
    comp_sched.plot_pulse_diagram(plot_backend="plotly",**kwargs).show()

# def set attenuation for all qubit
def set_atte_for(quantum_device:QuantumDevice,atte_value:int,mode:str,target_q:list=['q1']):
    """
        Set the attenuations for RO/XY by the given mode and atte. values.\n
        atte_value: integer multiple of 2,\n
        mode: 'ro' or 'xy',\n
        target_q: ['q1']
    """
    # Check atte value
    if atte_value%2 != 0:
        raise ValueError(f"atte_value={atte_value} is not the multiple of 2!")
    # set atte.
    if mode.lower() == 'ro':
        for q_name in target_q:
            set_readout_attenuation(quantum_device, quantum_device.get_element(q_name), out_att=atte_value, in_att=0)
    elif mode.lower() == 'xy':
        for q_name in target_q:
            set_drive_attenuation(quantum_device, quantum_device.get_element(q_name), out_att=atte_value)
    else:
        raise KeyError (f"The mode='{mode.lower()}' is not 'ro' or 'xy'!")
    
def two_tone_spec_sched_nco(
    qubit_name: str,
    spec_pulse_amp: float,
    spec_pulse_duration: float,
    spec_pulse_port: str,
    spec_pulse_clock: str,
    spec_pulse_frequencies: np.ndarray,
    repetitions: int = 1,
) -> Schedule:
    """
    Generate a batched schedule for performing fast two-tone spectroscopy using the
    `SetClockFrequency` operation for doing an NCO sweep.

    Parameters
    ----------
    spec_pulse_amp
        Amplitude of the spectroscopy pulse in Volt.
    spec_pulse_duration
        Duration of the spectroscopy pulse in seconds.
    spec_pulse_port
        Location on the device where the spectroscopy pulse should be applied.
    spec_pulse_clock
        Reference clock used to track the spectroscopy frequency.
    spec_pulse_frequencies
        Sample frequencies for the spectroscopy pulse in Hertz.
    repetitions
        The amount of times the Schedule will be repeated.
    """
    sched = Schedule("two-tone", repetitions)
    sched.add_resource(ClockResource(name=spec_pulse_clock, freq=spec_pulse_frequencies.flat[0]))

    for acq_idx, spec_pulse_freq in enumerate(spec_pulse_frequencies):
        sched.add(Reset(qubit_name))
        sched.add(SetClockFrequency(clock=spec_pulse_clock, clock_freq_new=spec_pulse_freq))
        sched.add(
            SquarePulse(
                duration=spec_pulse_duration,
                amp=spec_pulse_amp,
                port=spec_pulse_port,
                clock=spec_pulse_clock,
            )
        )
        sched.add(Measure(qubit_name, acq_index=acq_idx))

    return sched
    
def CSresults_alignPlot(quantum_device:QuantumDevice, results:dict):
    item_num = len(list(results.keys()))
    fig, ax = plt.subplots(1,item_num,figsize=plt.figaspect(1/item_num), sharey = False)
    if item_num == 1:
        for idx, q in enumerate(list(results.keys())):
            fr = results[q].run().quantities_of_interest["fr"].nominal_value
            dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax)
            ax.axvline(fr, color = "red", ls = "--")
            ax.set_title(f"{q} resonator")
    else:  
        for idx, q in enumerate(list(results.keys())):
            fr = results[q].run().quantities_of_interest["fr"].nominal_value
            dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax[idx])
            ax[idx].axvline(fr, color = "red", ls = "--")
            ax[idx].set_title(f"{q} resonator")
    
    fig.suptitle(f"Resonator spectroscopy, {quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()
    
def PDresults_alignPlot(quantum_device:QuantumDevice, results:dict, show_mode:str='pha'):
    item_num = len(list(results.keys()))
    fig, ax = plt.subplots(1,item_num,figsize=plt.figaspect(1/item_num), sharey = False)
    if item_num == 1:
        for idx, q in enumerate(list(results.keys())):
            if show_mode == 'pha':
                dh.to_gridded_dataset(results[q].dataset).y1.plot(ax = ax)
            else:
                dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax)
            ax.axhline(quantum_device.get_element(q).clock_freqs.readout(), color = "red", ls = "--")
            ax.set_title(f"{q} resonator")
    else:  
        for idx, q in enumerate(list(results.keys())):
            if show_mode == 'pha':
                dh.to_gridded_dataset(results[q].dataset).y1.plot(ax = ax[idx])
            else:
                dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax[idx])
            ax[idx].axhline(quantum_device.get_element(q).clock_freqs.readout(), color = "red", ls = "--")
            ax[idx].set_title(f"{q} resonator")
        
    fig.suptitle(f"Resonator Dispersive, {quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()

def FD_results_alignPlot(quantum_device:QuantumDevice, results:dict, show_mode:str='pha'):
    item_num = len(list(results.keys()))
    fig, ax = plt.subplots(1,item_num,figsize=plt.figaspect(1/item_num), sharey = False)
    if item_num == 1:
        for idx, q in enumerate(list(results.keys())):
            dressed_f = results[q].quantities_of_interest["freq_0"]
            offset = results[q].quantities_of_interest["offset_0"].nominal_value
            if show_mode == 'pha':
                dh.to_gridded_dataset(results[q].dataset).y1.plot(ax = ax)
            else:
                dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax)
            ax.axhline(dressed_f, color = "red", ls = "--")
            ax.axvline(offset, color = "red", ls = "--")
            ax.set_title(f"{q} resonator")
    else: 
        for idx, q in enumerate(list(results.keys())):
            dressed_f = results[q].quantities_of_interest["freq_0"]
            offset = results[q].quantities_of_interest["offset_0"].nominal_value
            if show_mode == 'pha':
                dh.to_gridded_dataset(results[q].dataset).y1.plot(ax = ax[idx])
            else:
                dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax[idx])
            ax[idx].axhline(dressed_f, color = "red", ls = "--")
            ax[idx].axvline(offset, color = "red", ls = "--")
            ax[idx].set_title(f"{q} resonator")
        
    fig.suptitle(f"Resonator Flux dependence, {quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()
    
def QFD_results_alignPlot(quantum_device:QuantumDevice, results:dict, show_mode:str='pha'):
    item_num = len(list(results.keys()))
    fig, ax = plt.subplots(1,item_num,figsize=plt.figaspect(1/item_num), sharey = False)
    if item_num == 1:
        for idx, q in enumerate(list(results.keys())):
            if show_mode == 'pha':
                dh.to_gridded_dataset(results[q].dataset).y1.plot(ax = ax)
            else:
                dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax)
            ax.set_title(f"{q} resonator")
    else: 
        for idx, q in enumerate(list(results.keys())):
            if show_mode == 'pha':
                dh.to_gridded_dataset(results[q].dataset).y1.plot(ax = ax[idx])
            else:
                dh.to_gridded_dataset(results[q].dataset).y0.plot(ax = ax[idx])
            ax[idx].set_title(f"{q} resonator")
        
    fig.suptitle(f"Qubit Flux dependence, {quantum_device.cfg_sched_repetitions()} repetitions")
    fig.tight_layout()
    plt.show()

def reset_offset(flux_callable_map:dict):
    for i in flux_callable_map:
        flux_callable_map[i](0.0)


        