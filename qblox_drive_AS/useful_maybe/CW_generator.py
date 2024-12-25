from __future__ import annotations

import contextlib
import json
from typing import TYPE_CHECKING, Callable

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import interact

from qblox_instruments import Cluster, ClusterType


def get_connected_modules(cluster: Cluster, filter_fn: Callable | None = None):
    def checked_filter_fn(mod: ClusterType) -> bool:
        if filter_fn is not None:
            return filter_fn(mod)
        return True

    return {
        mod.slot_idx: mod for mod in cluster.modules if mod.present() and checked_filter_fn(mod)
    }


def get_cluster_modules(cluster_ip:str, cluster_name:str, slot_idx:int, port_idx:int, ro_atte:int=0):
    # Close the chosen QCodes instrument to prevent name clash
    with contextlib.suppress(KeyError):
        Cluster.find_instrument(cluster_name).close()
    cluster = Cluster(
        name=cluster_name,
        identifier=cluster_ip,
        dummy_cfg={
            2: ClusterType.CLUSTER_QCM,
            4: ClusterType.CLUSTER_QCM,
            8: ClusterType.CLUSTER_QRM_RF,
            12: ClusterType.CLUSTER_QCM_RF,
            14: ClusterType.CLUSTER_QCM_RF,
            16: ClusterType.CLUSTER_QCM_RF,
        }
        if cluster_ip is None
        else None,
        )
    connected_module = get_connected_modules(cluster)[slot_idx]
    cluster.reset()
    if port_idx == 0:
        connected_module.out0_att(ro_atte)
    else:
        connected_module.out1_att(ro_atte)
    return cluster, connected_module


def define_sequencer(amp:float, connected_module):
    acquisitions = {"acq": {"num_bins": 1, "index": 0}}
    # Sequence program.
    seq_prog = """
        wait_sync 4

    loop: play    0,0,1200
        jmp     @loop
    """
    waveforms = {"dc": {"data": [amp for i in range(0, 1500)], "index": 0}}

    # Add sequence to single dictionary and write to JSON file.
    sequence = {
        "waveforms": waveforms,
        "weights": {},
        "acquisitions": acquisitions,
        "program": seq_prog,
    }
    with open("sequence.json", "w", encoding="utf-8") as file:
        json.dump(sequence, file, indent=4)
        file.close()

    connected_module.sequencer0.sequence("sequence.json")

def turn_on_sequencer(connected_module, port_idx:int, RF_freq:float, IF_freq:float):
    if connected_module.is_qrm_type:
        connected_module.out0_in0_lo_freq(RF_freq-IF_freq)
    else:
        if port_idx == 0:
            connected_module.out0_lo_freq(RF_freq-IF_freq) # if port =1, .out1_lo_freq
        else:
            connected_module.out1_lo_freq(RF_freq-IF_freq)

    connected_module.sequencer0.marker_ovr_en(True)
    connected_module.sequencer0.marker_ovr_value(3)  # Enables output on QRM-RF

    # Configure the sequencer
    connected_module.sequencer0.mod_en_awg(True)
    connected_module.sequencer0.nco_freq(IF_freq)
    connected_module.sequencer0.sync_en(True)

    connected_module.arm_sequencer(0)
    connected_module.start_sequencer(0)

def turn_off_sequencer(cluster, connected_module):
    # Stop sequencer.
    connected_module.stop_sequencer()
    cluster.reset()


def CW_executor(cluster:Cluster, slot_idx:int, port_idx:int, ro_atte:float, amp:float, RF_freq:float, LO_freq:float, cluster_ip:str="", cluster_name:str=""):
    IF_freq = RF_freq - LO_freq
    if cluster is None:
        cluster, connected_module = get_cluster_modules(cluster_ip, cluster_name, slot_idx, port_idx, ro_atte)
    else:
        connected_module = get_connected_modules(cluster)[slot_idx]
        cluster.reset()
        if port_idx == 0:
            connected_module.out0_att(ro_atte)
        else:
            connected_module.out1_att(ro_atte)
    define_sequencer(amp, connected_module)
    turn_on_sequencer(connected_module, port_idx, RF_freq, IF_freq)

    return cluster, connected_module



if __name__ == "__main__":

    """ Fill in """
    cluster_ip:str = "192.168.1.11"
    cluster_name:str = "cluster11"
    slot_idx:int = 6
    out_voltage:float = 0.3
    atte:int = 0 # dB, multiple of 2
    port_idx:int = 0 # 0 for qrm always
    want_freq = 6.1e9 # LO + IF = RF
    LO_freq = 6e9 # inside your experiment_setup.py


    """ execute """
    cluster, connected_module = CW_executor(None, slot_idx, port_idx, atte, out_voltage, want_freq, LO_freq, cluster_ip, cluster_name)
    stop = input("input 'n' to stop sequencer: ")
    if stop.lower() in ["n", "no"]:
        turn_off_sequencer(cluster, connected_module)
        cluster.close()

    
