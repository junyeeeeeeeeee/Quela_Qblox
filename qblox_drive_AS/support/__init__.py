import pickle, os
from typing import Callable
from qblox_drive_AS.Configs.ClusterAddress_rec import ip_register, port_register
from qcodes.instrument import find_or_create_instrument
from typing import Tuple
import ipywidgets as widgets
from IPython.display import display
from qblox_instruments import Cluster, PlugAndPlay, ClusterType
# from qblox_instruments.qcodes_drivers.qcm_qrm import QcmQrm
from qcodes import Instrument
from quantify_core.measurement.control import MeasurementControl
from quantify_core.visualization.pyqt_plotmon import PlotMonitor_pyqt as PlotMonitor
from quantify_scheduler.device_under_test.quantum_device import QuantumDevice
from quantify_scheduler.instrument_coordinator import InstrumentCoordinator
from quantify_scheduler.instrument_coordinator.components.qblox import ClusterComponent
from utils.tutorial_utils import (
    set_drive_attenuation,
    set_readout_attenuation,
)

from qblox_drive_AS.support.QDmanager import QDmanager, Data_manager
import qblox_drive_AS.support.UI_Window as uw
import qblox_drive_AS.support.Chip_Data_Store as cds
from numpy import ndarray
from numpy import asarray, real, vstack, array
from numpy import arctan2, pi, cos, sin
from qblox_drive_AS.support.UserFriend import *


def multiples_of_x(raw_number:float, x:float):
    multiples = int(raw_number//x) + 1
    return x * multiples
    

def find_nearest(ary:ndarray, value:float):
    """ find the element  which is closest to the given target_value in the given array"""
    ary = asarray(ary)
    idx = (abs(ary - value)).argmin()
    return float(ary[idx])

# initialize a measurement
def init_meas(QuantumDevice_path:str)->Tuple[QDmanager, Cluster, MeasurementControl, InstrumentCoordinator, dict]:
    """
    Initialize a measurement by the following 2 cases:\n
    ### Case 1: QD_path isn't given, create a new QD accordingly.\n
    ### Case 2: QD_path is given, load the QD with that given path.\n
    args:\n
    mode: 'new'/'n' or 'load'/'l'. 'new' need a self defined hardware config. 'load' load the given path. 
    """
    from qblox_drive_AS.support.UserFriend import warning_print
    import quantify_core.data.handling as dh
    meas_datadir = '.data'
    dh.set_datadir(meas_datadir)


    cfg, pth = {}, QuantumDevice_path 
    dr_loc = get_dr_loca(QuantumDevice_path)
    cluster_ip = ip_register[dr_loc.lower()]

    
    if cluster_ip in list(port_register.keys()):
        # try maximum 3 connections to prevent connect timeout error 
        try:
            cluster = Cluster(name = f"cluster{dr_loc.lower()}",identifier = f"qum.phys.sinica.edu.tw", port=int(port_register[cluster_ip]))
        except:
            try:
                warning_print("First cluster connection trying")
                cluster = Cluster(name = f"cluster{dr_loc.lower()}",identifier = f"qum.phys.sinica.edu.tw", port=int(port_register[cluster_ip]))
            except:
                warning_print("Second cluster connection trying")
                cluster = Cluster(name = f"cluster{dr_loc.lower()}",identifier = f"qum.phys.sinica.edu.tw", port=int(port_register[cluster_ip]))          
    else:
        try:
            warning_print("cluster IP connection trying")
            cluster = Cluster(name = f"cluster{dr_loc.lower()}", identifier = cluster_ip)
        except:
            raise KeyError("Check your cluster ip had been log into Experiment_setup.py with its connected DR, and also is its ip-port")
    
    ip = ip_register[dr_loc.lower()]
    
    # enable_QCMRF_LO(cluster) # for v0.6 firmware
    QRM_nco_init(cluster)
    Qmanager = QDmanager(pth)
    Qmanager.QD_loader()
    bias_controller = Qmanager.activate_str_Fctrl(cluster)

    meas_ctrl, ic = configure_measurement_control_loop(Qmanager.quantum_device, cluster)
    reset_offset(bias_controller)
    cluster.reset()

    return Qmanager, cluster, meas_ctrl, ic, bias_controller

def get_ip_specifier(QD_path:str):
    specifier = QD_path.split("#")[-1].split("_")[0]
    return specifier 

def get_dr_loca(QD_path:str):
    loc = os.path.split(QD_path)[-1].split("#")[0]
    return loc

# Configure_measurement_control_loop
""""# for firmware v0.6.2
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
"""
# for firmware v0.7.0
def configure_measurement_control_loop(
    device: QuantumDevice, cluster: Cluster, live_plotting: bool = False
    ) ->Tuple[MeasurementControl,InstrumentCoordinator]:
    meas_ctrl = find_or_create_instrument(MeasurementControl, recreate=True, name="meas_ctrl")
    ic = find_or_create_instrument(InstrumentCoordinator, recreate=True, name="ic")
    ic.timeout(60*60*120) # 120 hr maximum
    # Add cluster to instrument coordinator
    ic_cluster = ClusterComponent(cluster)
    ic.add_component(ic_cluster)

    if live_plotting:
        # Associate plot monitor with measurement controller
        plotmon = find_or_create_instrument(PlotMonitor, recreate=False, name="PlotMonitor")
        meas_ctrl.instr_plotmon(plotmon.name)

    # Associate measurement controller and instrument coordinator with the quantum device
    device.instr_measurement_control(meas_ctrl.name)
    device.instr_instrument_coordinator(ic.name)

    return (meas_ctrl, ic)



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

# Multi-clusters online version connect_clusters()
def connect_clusters_withinMulti(dr_loc:str,ip:str='192.168.1.10'):
    """
    This function is only for who doesn't use jupyter notebook to connect cluster.
    args: \n
    ip: 192.168.1.10\n
    So far the ip for Qblox cluster is named with 192.168.1.170 and 192.168.1.171
    """

    permissions = {}
    with PlugAndPlay() as p:            # Scan for available devices and display
        device_list = p.list_devices()
    for devi, info in device_list.items():
        permissions[info["identity"]["ip"]] = info["identity"]["ser"]
    if ip in permissions:
        print(f"{ip} is available to connect to!")
        return ip, permissions[ip]
    else:
        raise KeyError(f"{ip} is NOT available now!")
 
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


def reset_offset(flux_callable_map:dict):
    for i in flux_callable_map:
        flux_callable_map[i](0.0)


# def add log message into the sumInfo
def leave_LogMSG(MSG:str,sumInfo_path:str):
    """
    Leave the log message in the sumInfo with the given path.
    """
    with open(sumInfo_path, 'rb') as inp:
        gift = pickle.load(inp)
    gift["Log"] = MSG
    with open(sumInfo_path, 'wb') as file:
        pickle.dump(gift, file)
        print("Log message had been added!")


# set attenuations
def init_system_atte(quantum_device:QuantumDevice,qb_list:list,ro_out_att:int=20,xy_out_att:int=20):
    """
    Attenuation setting includes XY and RO. We don't change it once we set it.
    """
    # atte. setting
    set_atte_for(quantum_device,ro_out_att,'ro',qb_list)
    set_atte_for(quantum_device,xy_out_att,'xy',qb_list) 


# LO debug
def get_connected_modules(cluster: Cluster, filter_fn: Callable):
    def checked_filter_fn(mod: ClusterType) -> bool:
        if filter_fn is not None:
            return filter_fn(mod)
        return True

    return {
        mod.slot_idx: mod for mod in cluster.modules if mod.present() and checked_filter_fn(mod)
    }


def enable_QCMRF_LO(cluster):
    QCM_RFs = get_connected_modules(cluster,lambda mod: mod.is_qcm_type and mod.is_rf_type)
    for slot_idx in QCM_RFs:
        QCM_RFs[slot_idx].out0_lo_en(True)
        QCM_RFs[slot_idx].out1_lo_en(True)

    
    print(f"QCM_RF: {list(QCM_RFs.keys())} had been enabled the LO source successfully!")


def QRM_nco_init(cluster):
    QRM_RFs = get_connected_modules(cluster,lambda mod: mod.is_qrm_type and mod.is_rf_type)
    for slot_idx in QRM_RFs:
        for i in range(6):
            getattr(QRM_RFs[slot_idx], f"sequencer{i}").nco_prop_delay_comp_en(True)      
            getattr(QRM_RFs[slot_idx], f"sequencer{i}").nco_prop_delay_comp(50)
    
    print(f" NCO in QRM_RF: {list(QRM_RFs.keys())} had initialized NCO successfully!")


def advise_where_fq(QD:QDmanager,target_q:str,guess_g_Hz:float=48e6):
    fb = QD.Notewriter.get_bareFreqFor(target_q)
    fd = QD.quantum_device.get_element(target_q).clock_freqs.readout()
    g = guess_g_Hz
    x = fd-fb
    fq_Hz = fb - (g**2)/x
    return fq_Hz


def check_QD_info(QD_agent:QDmanager,target_q:str):
    from utils.tutorial_utils import show_readout_args, show_drive_args
    qubit = QD_agent.quantum_device.get_element(target_q)
    show_readout_args(qubit)
    show_drive_args(qubit)

def coupler_zctrl(Fctrl:dict,cp_elements:dict)->dict:
    """
    control coupler Z bias.
    ------------------------------
    # * Args:\n
    cp_elements follows the form: `{ coupler_name:bias (V)}`, like `{"c0":0.2}`
    ------------------------------
    # * Example:\n
    coupler_zctrl(dr2,cluster,cp_elements={"c0":0.2})
    ------------------------------
    """
    
    for cp in cp_elements:
        slightly_print(f"{cp} biased @ {round(cp_elements[cp],2)} V")
        Fctrl[cp](cp_elements[cp])
    
    return Fctrl

def compose_para_for_multiplexing(QD_agent:QDmanager,ro_elements,mode:str)->dict:
    """
    Get the dict about the required values for all qubits in quantum_device.
    The required value can be assigned by the arg `mode`.
    ------
    ### Args:\n
    * ro_elements: a dict with the keyname in qubit name. ex: {"q0":[ ],"q1":[ ],...}\n
    * mode:\n 
        'r1' for RO-amp, 'r2' for acq-delay, 'r3' for RO-duration, 'r4' for integration time.\n
        'd1' for xy-amp,                     'd3' for xy-duration,
    ----
    ### Returns:\n
    A dict with the same keyname as the `ro_elements`, and also with the value about the required mode.  
    """
    if type(ro_elements) is dict:
        qubits = list(ro_elements.keys())
    elif type(ro_elements) is list:
        qubits:list = ro_elements
    else:
        raise ValueError(f"The type of ro_elements should be list or dict but `{type(ro_elements)}` was recieved!")
    ans = {}
    for qubit in qubits:
        match mode.lower():
            case "r1":
                ans[qubit] = QD_agent.quantum_device.get_element(qubit).measure.pulse_amp()
            case "r2":
                ans[qubit] = QD_agent.quantum_device.get_element(qubit).measure.acq_delay()
            case "r3":
                ans[qubit] = QD_agent.quantum_device.get_element(qubit).measure.pulse_duration()
            case "r4":
                ans[qubit] = QD_agent.quantum_device.get_element(qubit).measure.integration_time()
            case 'd1':
                ans[qubit] = QD_agent.quantum_device.get_element(qubit).rxy.amp180()
            case 'd3':
                ans[qubit] = QD_agent.quantum_device.get_element(qubit).rxy.duration()
            case _:
                raise KeyError(f"Un-supported mode = {mode} was given !")
    
    return ans

def rotation_matrix(angle)->ndarray:
    rotate_matrix = array([
        [cos(angle), -sin(angle)],
        [sin(angle),  cos(angle)]
    ])
    return rotate_matrix

def rotate_onto_Inphase(point_1:ndarray, point_2:ndarray, angle_degree:float=None)->tuple[ndarray,float]:
    """
        Give 2 points, rotate them makes them lie one the x-axis and return the rotate angle in degree.\n
        If you also give the rotate angle in degree, we rotate them according to this angle.
    """
    vector_21 = point_2 - point_1
    
    if angle_degree is None:
        angle = -arctan2(vector_21[1], vector_21[0])  # Angle to align with x-axis
    else:
        angle = -angle_degree*pi/180
 

    # Apply rotation to both translated points
    P1_rotated = rotation_matrix(angle) @ point_1 
    P2_rotated = rotation_matrix(angle) @ point_2 

    return vstack((P1_rotated,P2_rotated)), -angle*180/pi


def rotate_data(data:ndarray, angle_degree:float)->ndarray:
    """ data shape (IQ, shots)"""

    angle = -angle_degree*pi/180
    data_rotated = rotation_matrix(angle) @ data

    return data_rotated




if __name__ == "__main__":
    a = {"a1":10,"a2":12}
    b = {}
    print({**a,**b})