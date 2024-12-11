""" Thsi script helps you build a new QD_file """
import os, sys, rich
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support import QDmanager
from quantify_scheduler.backends.qblox_backend import QbloxHardwareCompilationConfig



cluster_IP:str = "192.168.1.11"
dr_name:str = "dr1"
qubit_number_onChip:int = 2
coupler_number_onChip:int = 0
chip_name:str = "5Q4C"
chip_type:str = "5Q4C"
old_QD_path:str = "qblox_drive_AS/QD_backup/20241211/DR1#11_SumInfo.pkl" # set the path in string When you want to update the Hcfg. Otherwise, set it None


# Hcfg = {
#     "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
#     f"cluster{dr_name}": {
#         "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
#         "ref": "internal",  # Use shared clock reference of the cluster
#         "instrument_type": "Cluster",
#         # ============ DRIVE ============#
#         f"cluster{dr_name}_module4": {
#             "instrument_type": "QCM_RF",
#             "complex_output_0": {
#                 "output_att": 0,
#                 "dc_mixer_offset_I": 0.0,
#                 "dc_mixer_offset_Q": 0.0,
#                 "lo_freq": 3e9,
#                 "portclock_configs": [
#                     {
#                         "port": "q0:mw",
#                         "clock": "q0.01",
#                         "mixer_amp_ratio": 1.0,
#                         "mixer_phase_error_deg": 0.0,
#                     }
#                 ],
#             },
#             "complex_output_1": {
#                 "output_att": 0,
#                 "dc_mixer_offset_I": 0.0,
#                 "dc_mixer_offset_Q": 0.0,
#                 "lo_freq": 3e9,
#                 "portclock_configs": [
#                     {
#                         "port": "q1:mw",
#                         "clock": "q1.01",
#                         "mixer_amp_ratio": 1.0,
#                         "mixer_phase_error_deg": 0.0,
#                     }
#                 ],
#             },
#         },

#         # ============ FLUX ============#
#         f"cluster{dr_name}_module2": {
#             "instrument_type": "QCM",
#             "real_output_0": {"portclock_configs": [{"port": "q0:fl", "clock": "cl0.baseband"}]},
#             "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
#             "real_output_2": {"portclock_configs": [{"port": "c0:fl", "clock": "cl0.baseband"}]},
#             "real_output_3": {"portclock_configs": [{"port": "c1:fl", "clock": "cl0.baseband"}]},
#         },
        
#         # ============ READOUT ============#
#         f"cluster{dr_name}_module6": {
#             "instrument_type": "QRM_RF",
#             "complex_output_0": {
#                 "output_att": 0,
#                 "input_att": 0,
#                 "dc_mixer_offset_I": 0.0,
#                 "dc_mixer_offset_Q": 0.0,
#                 "lo_freq": 4.62e9,       # *** Should be set as a parameter later on
#                 "portclock_configs": [
#                     {
#                         "port": "q:res",
#                         "clock": "q0.ro",
#                         "mixer_amp_ratio": 1.0,
#                         "mixer_phase_error_deg": 0.0,
#                     },
#                     {
#                         "port": "q:res",
#                         "clock": "q1.ro",
#                         "mixer_amp_ratio": 1.0,
#                         "mixer_phase_error_deg": 0.0,
#                     },
#                 ],
#             },
#         },
#     },
# }
# Hcfg = QbloxHardwareCompilationConfig.model_validate(Hcfg)
# rich.print(Hcfg)
Hcfg = {
    "config_type": "quantify_scheduler.backends.qblox_backend.QbloxHardwareCompilationConfig",
    "hardware_description": {
        f"cluster{dr_name}": {
            "instrument_type": "Cluster",
            "modules": {
                "4": {
                    "instrument_type": "QCM_RF"
                },
                "2": {
                    "instrument_type": "QCM"
                },
                "6": {
                    "instrument_type": "QRM_RF"
                }
            },
            "sequence_to_file": False,
            "ref": "internal"
        }
    },
    "hardware_options": {
        "output_att": {
            "q0:mw-q0.01": 0,
            "q1:mw-q1.01": 0,
            "q:res-q0.ro": 22,
            "q:res-q1.ro": 22
        },
        "mixer_corrections": {
            "q0:mw-q0.01": {
                "auto_lo_cal": "on_lo_interm_freq_change",
                "auto_sideband_cal": "on_interm_freq_change",
                # "dc_offset_i": 0.0,
                # "dc_offset_q": 0.0,
                # "amp_ratio": 1.0,
                # "phase_error": 0.0
            },
            "q1:mw-q1.01": {
                "auto_lo_cal": "on_lo_interm_freq_change",
                "auto_sideband_cal": "on_interm_freq_change",
                # "dc_offset_i": 0.0,
                # "dc_offset_q": 0.0,
                # "amp_ratio": 1.0,
                # "phase_error": 0.0
            },
            "q:res-q0.ro": {
                "auto_lo_cal": "on_lo_interm_freq_change",
                "auto_sideband_cal": "on_interm_freq_change",
                # "dc_offset_i": 0.0,
                # "dc_offset_q": 0.0,
                # "amp_ratio": 1.0,
                # "phase_error": 0.0
            },
            "q:res-q1.ro": {
                "auto_lo_cal": "on_lo_interm_freq_change",
                "auto_sideband_cal": "on_interm_freq_change",
                # "dc_offset_i": 0.0,
                # "dc_offset_q": 0.0,
                # "amp_ratio": 1.0,
                # "phase_error": 0.0
            }
        },
        "modulation_frequencies": {
            "q0:mw-q0.01": {
                "interm_freq": 80000000.0
            },
            "q1:mw-q1.01": {
                "interm_freq": 80000000.0
            },
            "q:res-q0.ro": {
                "lo_freq": 4.62e9
            },
            "q:res-q1.ro": {
                "lo_freq": 4.62e9
            }
        }
    },
    "connectivity": {
        "graph": [
            [
                f"cluster{dr_name}.module4.complex_output_0",
                "q0:mw"
            ],
            [
                f"cluster{dr_name}.module4.complex_output_1",
                "q1:mw"
            ],
            [
                f"cluster{dr_name}.module2.real_output_0",
                "q0:fl"
            ],
            [
                f"cluster{dr_name}.module2.real_output_1",
                "q1:fl"
            ],
            [
                f"cluster{dr_name}.module6.complex_output_0",
                "q:res"
            ],
            [
                f"cluster{dr_name}.module6.complex_output_0",
                "q:res"
            ]
        ]
    }
}


if old_QD_path is None or str(old_QD_path) == "" :
    QD_agent = QDmanager()
    QD_agent.build_new_QD(qubit_number_onChip,coupler_number_onChip,Hcfg,cluster_IP,dr_name,chip_name,chip_type)
else:
    QD_agent = QDmanager(old_QD_path)
    QD_agent.QD_loader(new_Hcfg=Hcfg)
    
QD_agent.QD_keeper()


    


