import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from numpy import arange, linspace, logspace
from quantify_scheduler.backends.graph_compilation import SerialCompiler
### Note: The qubits are labeled with the physical position on the chip from the left, and starts from 0.

#%%
# Please register the dr and its corresponding cluster ip here first!
ip_register = {
    "dr1":"192.168.1.11",
    "dr2":"192.168.1.10",
    "dr3":"192.168.1.13",
    "dr4":"192.168.1.81",
    "drke":"192.168.1.242"
} # all keys in lower
port_register = {
    "192.168.1.10":"5010",
    "192.168.1.11":"5011",
    "192.168.1.12":"5012",
    "192.168.1.13":"5013",
    "192.168.1.81":"5081",
    "192.168.1.242":"5242"
}


#%%
# Hardware settings
Hcfg_dr1 = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    f"clusterdr1": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        f"clusterdr1_module4": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q0:mw",
                        "clock": "q0.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q1:mw",
                        "clock": "q1.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            }
        },
        f"clusterdr1_module8": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q2:mw",
                        "clock": "q2.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q3:mw",
                        "clock": "q3.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            }
        },
        f"clusterdr1_module10": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q4:mw",
                        "clock": "q4.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q5:mw",
                        "clock": "q5.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            }
        },
        # ============ FLUX ============#
        f"clusterdr1_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]},
        },
        # ============ READOUT ============#
        f"clusterdr1_module6": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 5.81e9,       # *** Should be set as a parameter later on
                "portclock_configs": [
                    {
                        "port": "q:res",
                        "clock": "q0.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q1.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q2.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q3.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q4.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q5.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                ],
            },
        },
    },
}

Hcfg_dr2 = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    f"clusterdr2": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        f"clusterdr2_module12": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q0:mw",
                        "clock": "q0.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q1:mw",
                        "clock": "q1.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        f"clusterdr2_module14": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q2:mw",
                        "clock": "q2.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q3:mw",
                        "clock": "q3.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        f"clusterdr2_module16": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q4:mw",
                        "clock": "q4.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q5:mw",
                        "clock": "q5.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        # ============ FLUX ============#
        f"clusterdr2_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q0:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_3": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
        },
        f"clusterdr2_module4": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]}
        },
        # ============ READOUT ============#
        f"clusterdr2_module8": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 6.17e9,       # *** Should be set as a parameter later on
                "portclock_configs": [
                    {
                        "port": "q0:res",
                        "clock": "q0.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q1:res",
                        "clock": "q1.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q2:res",
                        "clock": "q2.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q3:res",
                        "clock": "q3.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q4:res",
                        "clock": "q4.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                ],
            },
        },
    },
}

Hcfg_dr3 = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    "clusterdr3": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        "clusterdr3_module8": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 3e9,
                "portclock_configs": [
                    {
                        "port": "q0:mw",
                        "clock": "q0.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 3e9,
                "portclock_configs": [
                    {
                        "port": "q1:mw",
                        "clock": "q1.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        "clusterdr3_module12": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 3e9,
                "portclock_configs": [
                    {
                        "port": "q2:mw",
                        "clock": "q2.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 3e9,
                "portclock_configs": [
                    {
                        "port": "q3:mw",
                        "clock": "q3.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        "clusterdr3_module16": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q4:mw",
                        "clock": "q4.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        # ============ FLUX ============#
        "clusterdr3_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q0:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_3": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
        },
        "clusterdr3_module4": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]},
        },
        # "clusterdr3_module6": {
        #     "instrument_type": "QCM",
        #     "real_output_0": {"portclock_configs": [{"port": "q5:fl", "clock": "cl0.baseband"}]},
        #     "real_output_1": {"portclock_configs": [{"port": "q6:fl", "clock": "cl0.baseband"}]},
        #     "real_output_2": {"portclock_configs": [{"port": "q7:fl", "clock": "cl0.baseband"}]},
        #     "real_output_3": {"portclock_configs": [{"port": "q8:fl", "clock": "cl0.baseband"}]},

        # },
        # ============ READOUT ============#
        "clusterdr3_module18": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 6.06e9,       # *** Should be set as a parameter later on
                "portclock_configs": [
                    {
                        "port": "q:res",
                        "clock": "q0.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q1.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q2.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q3.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q4.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                ],
            },
        },
    },
}

Hcfg_dr4 = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    f"clusterdr4": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        f"clusterdr4_module10": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q0:mw",
                        "clock": "q0.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q1:mw",
                        "clock": "q1.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            }
        },
        f"clusterdr4_module12": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q2:mw",
                        "clock": "q2.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q3:mw",
                        "clock": "q3.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        f"clusterdr4_module14": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q4:mw",
                        "clock": "q4.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q5:mw",
                        "clock": "q5.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },
        

        # ============ FLUX ============#
        "clusterdr4_module4": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]},
        },
        # ============ READOUT ============#
        f"clusterdr4_module18": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq":5.6e9,       # *** Should be set as a parameter later on
                "portclock_configs": [
                    {
                        "port": "q:res",
                        "clock": "q0.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q1.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q2.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q3.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q4.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q5.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                ],
            },
        },
    },
}

Hcfg_drke = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    f"clusterdrke": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        f"clusterdrke_module4": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q0:mw",
                        "clock": "q0.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q1:mw",
                        "clock": "q1.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            }
        },
        f"clusterdrke_module8": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q2:mw",
                        "clock": "q2.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
            "complex_output_1": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q3:mw",
                        "clock": "q3.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            }
        },
        f"clusterdrke_module10": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4e9,
                "portclock_configs": [
                    {
                        "port": "q4:mw",
                        "clock": "q4.01",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    }
                ],
            },
        },

        # ============ FLUX ============#
        f"clusterdrke_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
        },
        # ============ READOUT ============#
        f"clusterdrke_module6": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq":5.8e9,       # *** Should be set as a parameter later on
                "portclock_configs": [
                    {
                        "port": "q:res",
                        "clock": "q0.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q1.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q2.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q3.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                    {
                        "port": "q:res",
                        "clock": "q4.ro",
                        "mixer_amp_ratio": 1.0,
                        "mixer_phase_error_deg": 0.0,
                    },
                ],
            },
        },
    },
}

def get_FluxController(cluster, ip:str)->dict:
    which_dr = ''
    for dr in ip_register:
        if ip_register[dr] == ip:
            which_dr = dr
    if which_dr == '':
        raise ValueError("Please register the dr location with it's cluster ip in Experiment_setup.py in support folder! Or check the given ip_label is correct!")
    
    Fctrl = {}
    if which_dr.lower() == 'dr2':
        Fctrl: callable = {
            "q0":cluster.module2.out0_offset,
            "q1":cluster.module2.out1_offset,
            "q2":cluster.module2.out2_offset,
            "q3":cluster.module2.out3_offset,
            "q4":cluster.module4.out0_offset
        }
    elif which_dr.lower() == 'dr3':
        Fctrl: callable = {
            "q0":cluster.module2.out0_offset,
            "q1":cluster.module2.out1_offset,
            "q2":cluster.module2.out2_offset,
            "q3":cluster.module2.out3_offset,
            "q4":cluster.module4.out0_offset,
        }   
    elif which_dr.lower() == 'dr1':
        Fctrl: callable = {
            "q2":cluster.module2.out0_offset,
            "q3":cluster.module2.out1_offset,
            "q4":cluster.module2.out2_offset,
        }
    elif which_dr.lower() == 'dr4':
        Fctrl: callable = {
            "q2":cluster.module4.out0_offset,
            "q3":cluster.module4.out1_offset,
            "q4":cluster.module4.out2_offset,
        }
    elif which_dr.lower() == 'drke':
        Fctrl: callable = {
            "q2":cluster.module2.out0_offset,
            "q3":cluster.module2.out1_offset,
        }

    else:
        raise KeyError ("please input ip label like '170' or '171'!")
    return Fctrl

def get_CouplerController(cluster, ip:str)->dict:
    which_dr = ''
    for dr in ip_register:
        if ip_register[dr] == ip:
            which_dr = dr
    if which_dr == '':
        raise ValueError("Please register the dr location with it's cluster ip in Experiment_setup.py in support folder! Or check the given ip_label is correct!")
    
    Cctrl = {}
    if which_dr.lower() == 'dr3':
        Cctrl = {
            "c0":cluster.module6.out0_offset,
            "c1":cluster.module6.out1_offset,
            "c2":cluster.module6.out2_offset,
            "c3":cluster.module6.out3_offset
        }
    elif which_dr.lower() == 'dr4':
        Cctrl = {
        }
    return Cctrl

# # Cluster registerations
# clusters_online = {
#     "192.168.1.170":{"loc":'DR3',"ser":'00015_2247_002'},
#     "192.168.1.171":{"loc":'DR2',"ser":'00015_2406_014'}
# }

# Hcfg map
hcfg_map = {"dr2":Hcfg_dr2,'dr1':Hcfg_dr1,'dr3':Hcfg_dr3,'dr4':Hcfg_dr4,'drke':Hcfg_drke} # all keys in lower


