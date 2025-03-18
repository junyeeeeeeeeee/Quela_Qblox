""" Thsi script helps you build a new QD_file """
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support import QDmanager


cluster_IP:str = "192.168.1.10"
dr_name:str = "dr2"
qubit_number_onChip:int = 2
coupler_number_onChip:int = 0
chip_name:str = "firmware_update_test"
chip_type:str = "5Q4C"
old_QD_path:str = "" # set the path in string When you want to update the Hcfg. Otherwise, set it None


Hcfg = {
    "config_type": "quantify_scheduler.backends.qblox_backend.QbloxHardwareCompilationConfig",
    "hardware_description": {
        f"cluster{dr_name}": {
            "instrument_type": "Cluster",
            "modules": {
                "2": {
                    "instrument_type": "QCM"
                },
                "12": {
                    "instrument_type": "QCM_RF"
                },
                "8": {
                    "instrument_type": "QRM_RF"
                },
                "14": {
                    "instrument_type": "QCM_RF"
                },
            },
            "sequence_to_file": False,
            "ref": "internal"
        }
    },
    "hardware_options": {
        "output_att": {
            "q0:mw-q0.01": 0,
            "q1:mw-q1.01": 0,
            "q:res-q0.ro": 0,
            "q:res-q1.ro": 0
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
                "lo_freq": 4e9
            },
            "q1:mw-q1.01": {
                "lo_freq": 4e9
            },
            "q:res-q0.ro": {
                "lo_freq": 6.03e9
            },
            "q:res-q1.ro": {
                "lo_freq": 6.03e9
            }
        }
    },
    "connectivity": {
        "graph": [
            [
                f"cluster{dr_name}.module12.complex_output_0",
                "q0:mw"
            ],
            [
                f"cluster{dr_name}.module12.complex_output_1",
                "q1:mw"
            ],
            [
                f"cluster{dr_name}.module2.real_output_2",
                "q0:fl"
            ],
            [
                f"cluster{dr_name}.module2.real_output_3",
                "q1:fl"
            ],
            [
                f"cluster{dr_name}.module8.complex_output_0",
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


    


