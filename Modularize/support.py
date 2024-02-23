from quantify_scheduler.json_utils import SchedulerJSONDecoder
from typing import Tuple
import warnings
from pathlib import Path
from typing import List, Union, Literal
import matplotlib.pyplot as plt
import ipywidgets as widgets
import numpy as np
from numpy.random import randint
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
from quantify_scheduler.operations.pulse_library import SetClockFrequency, SquarePulse, DRAGPulse
from quantify_scheduler.resources import ClockResource
from quantify_scheduler.schedules import heterodyne_spec_sched_nco, rabi_sched, t1_sched
from quantify_scheduler.schedules.timedomain_schedules import ramsey_sched
from quantify_scheduler.schedules.schedule import Schedule
from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis, QubitSpectroscopyAnalysis
from quantify_core.analysis.base_analysis import Basic2DAnalysis
from utils.tutorial_analysis_classes import (
    QubitFluxSpectroscopyAnalysis,
    ResonatorFluxSpectroscopyAnalysis,
)
from quantify_scheduler.operations.gate_library import X90, Measure, Reset, Rxy, X, Y
from utils.tutorial_utils import (
    set_drive_attenuation,
    set_readout_attenuation,
    show_args,
    show_drive_args,
    show_readout_args,
)

# TODO: Test this class to store the bias info
# dict to save the filux bias information
class FluxBiasDict():
    """
    This class helps to memorize the flux bias. 
    """
    def __init__(self,qb_number:int):
        self.__bias_dict = {}
        self.q_num = qb_number
        self.init_bias()
    def get_bias_dict(self)->dict:
        """
        Return the dictionary contain the bias info.
        """
        return self.__bias_dict
    def init_bias(self):
        """
        Initialize the dict when you create this object.
        """
        for i in range(self.q_num):
            self.__bias_dict[f"q{i}"] = {}
            for bias_position in ["SweetSpot","TuneAway"]:
                self.__bias_dict[f"q{i}"][bias_position] = 0.0
    def save_sweetspotBias_for(self,target_q:str='q0',bias:float=0.0):
        """
            Set the sweet spot bias for the given qubit, target_q label starts from 0.\n
            Ex. target_q = 'q0'
        """
        self.__bias_dict[target_q]["SweetSpot"] = bias
    def save_tuneawayBias_for(self,target_q:str='q0',bias:float=0.0):
        """
            Set the tune away bias for the given qubit, target_q label starts from 0.\n
            Ex. target_q = 'q0'
        """
        self.__bias_dict[target_q]["TuneAway"] = bias
    def activate_from_dict(self,old_bias_dict:dict):
        """
        Activate the dict which super from the old record.
        """
        for q_old in old_bias_dict:
            for bias_position in old_bias_dict[q_old]:
                self.__bias_dict[q_old][bias_position] = old_bias_dict[q_old][bias_position]
    def get_sweetBiasFor(self,target_q:str):
        """
        Return the sweetSpot bias for the target qubit.
        """
        return self.__bias_dict[target_q]["SweetSpot"]
    def get_tuneawayBiasFor(self,target_q:str):
        """
        Return the tuneAway bias for the target qubit.
        """
        return self.__bias_dict[target_q]["TuneAway"]
    

class QDmanager():
    def __init__(self,QD_path:str=''):
        self.path = QD_path
        self.refIQ = {}
        self.Hcfg = {}
        self.Log = "" 
    
    def memo_refIQ(self,ref_dict:dict):
        """
        Memorize the reference IQ according to the given ref_dict, which the key named in "q0"..., and the value is composed in list by IQ values.\n
        Ex. ref_dict={"q0":[0.03,-0.004],...} 
        """
        for q in ref_dict:
            self.refIQ[q] = ref_dict[q]
    
    def refresh_log(self,message:str):
        """
        Leave the message for this file.
        """
        self.Log = message

    def QD_loader(self):
        """
        Load the QuantumDevice, Bias config, hardware config and Flux control callable dict from a given json file path contain the serialized QD.
        """
        import pickle
        with open(self.path, 'rb') as inp:
            gift = pickle.load(inp)
        # string and int
        self.Log = gift["Log"]
        self.q_num = len(list(gift["Flux"].keys()))
        # class    
        self.Fluxmanager :FluxBiasDict = FluxBiasDict(qb_number=self.q_num)
        self.Fluxmanager.activate_from_dict(gift["Flux"])
        self.quantum_device :QuantumDevice = gift["QD"]
        # dict
        self.Hcfg = gift["Hcfg"]
        self.refIQ = gift["refIQ"]
   

        self.quantum_device.hardware_config(self.Hcfg)
        print("Old friends loaded!")
    
    def QD_keeper(self, special_path:str=''):
        """
        Save the merged dictionary to a json file with the given path. \n
        Ex. merged_file = {"QD":self.quantum_device,"Flux":self.Fluxmanager.get_bias_dict(),"Hcfg":Hcfg,"refIQ":self.refIQ,"Log":self.Log}
        """
        import pickle
        import os
        if self.path == '':
            from Modularize.path_book import qdevice_backup_dir
            qd_folder = build_folder_today(qdevice_backup_dir)
            self.path = os.path.join(qd_folder,"SumInfo.pkl")
        Hcfg = self.quantum_device.generate_hardware_config()
        merged_file = {"QD":self.quantum_device,"Flux":self.Fluxmanager.get_bias_dict(),"Hcfg":Hcfg,"refIQ":self.refIQ,"Log":self.Log}
        
        with open(self.path if special_path == '' else special_path, 'wb') as file:
            pickle.dump(merged_file, file)
            print(f'Summarized info had successfully saved to the given path!')
    
    def build_new_QD(self,qubit_number:int,Hcfg:dict):
        print("Building up a new quantum device system....")
        self.q_num = qubit_number
        self.Hcfg = Hcfg
        self.quantum_device = QuantumDevice("academia_sinica_device")
        self.quantum_device.hardware_config(self.Hcfg)
        
        # store references
        self.quantum_device._device_elements = list()

        for i in range(qubit_number):
            qubit = BasicTransmonElement(f"q{i}")
            qubit.measure.acq_channel(i)
            self.quantum_device.add_element(qubit)
            self.quantum_device._device_elements.append(qubit)

        self.Fluxmanager :FluxBiasDict = FluxBiasDict(self.q_num)


# initialize a measurement
def init_meas(QuantumDevice_path:str='',qubit_number:int=5, mode:str='new')->Tuple[QDmanager, Cluster, MeasurementControl, InstrumentCoordinator]:
    """
    Initialize a measurement by the following 2 cases:\n
    ### Case 1: QD_path isn't given, create a new QD accordingly.\n
    ### Case 2: QD_path is given, load the QD with that given path.\n
    args:\n
    mode: 'new'/'n' or 'load'/'l'. 'new' need a self defined hardware config. 'load' load the given path. 
    """
    
    import quantify_core.data.handling as dh
    meas_datadir = '.data'
    dh.set_datadir(meas_datadir)
    if mode.lower() in ['new', 'n']:
        from Experiment_setup import hardware_cfg
        cfg, pth = hardware_cfg, ''
    elif mode.lower() in ['load', 'l']:
        cfg, pth = {}, QuantumDevice_path 
    else:
        raise KeyError("The given mode can not be recognized!")
    
    Qmanager = QDmanager(pth)
    if pth == '':
        Qmanager.build_new_QD(qubit_number,cfg)
        Qmanager.refresh_log("new-born!")
    else:
        Qmanager.QD_loader()

    # Connect to the Qblox cluster
    connect, ip = connect_clusters()
    cluster = Cluster(name = "cluster0", identifier = ip.get(connect.value))
    meas_ctrl, ic = configure_measurement_control_loop(Qmanager.quantum_device, cluster)
    
    return Qmanager, cluster, meas_ctrl, ic

# Configure_measurement_control_loop
def configure_measurement_control_loop(
    device: QuantumDevice, cluster: Cluster, live_plotting: bool = False
) ->Tuple[MeasurementControl,InstrumentCoordinator]:
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

# close all instruments
def shut_down(cluster:Cluster,flux_map:dict):
    '''
        Disconnect all the instruments.
    '''
    reset_offset(flux_map)
    cluster.reset() 
    Instrument.close_all() 
    print("All instr are closed and zeroed all flux bias!")

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
    
def CSresults_alignPlot(quantum_device:QuantumDevice, results:dict):
    item_num = len(list(results.keys()))
    fig, ax = plt.subplots(1,item_num,figsize=plt.figaspect(1/item_num), sharey = False)
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

def reset_offset(flux_callable_map:dict):
    for i in flux_callable_map:
        flux_callable_map[i](0.0)





# TODO: ToTest ZZinteractions
def ZZinteractions_sched(
    times: Union[np.ndarray, float],
    ctrl_qubit: str,
    meas_qubit: str,
    artificial_detuning: float = 0,
    repetitions: int = 1,
) -> Schedule:
    r"""
    Generate a schedule for performing a Ramsey experiment to measure the
    dephasing time :math:`T_2^{\star}`.

    Schedule sequence
        .. centered:: Reset -- pi/2 -- Idle(tau) -- pi/2 -- Measure

    See section III.B.2. of :cite:t:`krantz_quantum_2019` for an explanation of the Bloch-Redfield
    model of decoherence and the Ramsey experiment.

    Parameters
    ----------
    times
        an array of wait times tau between the pi/2 pulses.
    artificial_detuning
        frequency in Hz of the software emulated, or ``artificial`` qubit detuning, which is
        implemented by changing the phase of the second pi/2 (recovery) pulse. The
        artificial detuning changes the observed frequency of the Ramsey oscillation,
        which can be useful to distinguish a slow oscillation due to a small physical
        detuning from the decay of the dephasing noise.
    qubit
        the name of the qubit e.g., :code:`"q0"` to perform the Ramsey experiment on.
    repetitions
        The amount of times the Schedule will be repeated.

    Returns
    -------
    :
        An experiment schedule.

    """
    # ensure times is an iterable when passing floats.
    times = np.asarray(times)
    times = times.reshape(times.shape or (1,))

    schedule = Schedule("Ramsey", repetitions)

    if isinstance(times, float):
        times = [times]

    for i, tau in enumerate(times):
        schedule.add(Reset(meas_qubit,ctrl_qubit), label=f"Reset {i}")
        schedule.add(X(ctrl_qubit))
        schedule.add(X90(meas_qubit))

        # the phase of the second pi/2 phase progresses to propagate
        recovery_phase = np.rad2deg(2 * np.pi * artificial_detuning * tau)
        schedule.add(
            Rxy(theta=90, phi=recovery_phase, qubit=meas_qubit), ref_pt="start", rel_time=tau
        )
        schedule.add(Measure(meas_qubit, acq_index=i), label=f"Measurement {i}")
    return schedule      

# TODO: RB
def SQRB_schedule(
    qubit:str,
    gate_num:int,
    repetitions: int=1,
) -> Schedule:

    sched = Schedule('RB',repetitions)
    sched.add(Reset(qubit))
    for idx in range(gate_num):
        pass
    sched.add(Measure(qubit))
    
    return sched

# generate time label for netCDF file name
def get_time_now()->str:
    """
    Since we save the Xarray into netCDF, we use the current time to encode the file name.\n
    Ex: 19:23:34 return H19M23S34 
    """
    import datetime
    current_time = datetime.datetime.now()
    return f"H{current_time.hour}M{current_time.minute}S{current_time.second}"

# build the folder for the data today
def build_folder_today(parent_path:str):
    """
    Build up and return the folder named by the current date in the parent path.\n
    Ex. parent_path='D:/Examples/'
    """
    import os
    import datetime
    current_time = datetime.datetime.now()
    folder = f"{current_time.year}_{current_time.month}_{current_time.day}"
    new_folder = os.path.join(parent_path, folder) 
    if not os.path.isdir(new_folder):
        os.mkdir(new_folder) 
        print(f"Folder {current_time.year}_{current_time.month}_{current_time.day} had been created!")
    else:
        print(f"Folder {current_time.year}_{current_time.month}_{current_time.day} exist!")
    return new_folder



# def add log message into the sumInfo
def leave_LogMSG(MSG:str,sumInfo_path:str):
    """
    Leave the log message in the sumInfo with the given path.
    """
    import pickle
    with open(sumInfo_path, 'rb') as inp:
        gift = pickle.load(inp)
    gift["Log"] = MSG
    with open(sumInfo_path, 'wb') as file:
        pickle.dump(gift, file)
        print("Log message had been added!")


# set attenuations
def init_system_atte(quantum_device:QuantumDevice,qb_list:list,ro_out_att:int=20,xy_out_att:int=10):
    """
    Attenuation setting includes XY and RO. We don't change it once we set it.
    """
    # atte. setting
    set_atte_for(quantum_device,ro_out_att,'ro',qb_list)
    set_atte_for(quantum_device,xy_out_att,'xy',qb_list) 