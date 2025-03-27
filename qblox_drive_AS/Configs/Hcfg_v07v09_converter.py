import rich
from quantify_scheduler.backends.qblox_backend import QbloxHardwareCompilationConfig

dr_name = 'dr2'
Hcfg = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    f"cluster{dr_name}": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        f"cluster{dr_name}_module4": {
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

        # ============ FLUX ============#
        f"cluster{dr_name}_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q0:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "c0:fl", "clock": "cl0.baseband"}]},
            "real_output_3": {"portclock_configs": [{"port": "c1:fl", "clock": "cl0.baseband"}]},
        },
        
        # ============ READOUT ============#
        f"cluster{dr_name}_module6": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 4.62e9,       # *** Should be set as a parameter later on
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
                ],
            },
        },
    },
}
Hcfg = QbloxHardwareCompilationConfig.model_validate(Hcfg)
rich.print(Hcfg)
serialized_config = Hcfg.model_dump_json(exclude_unset=True)

deserialized_config = QbloxHardwareCompilationConfig.model_validate_json(serialized_config)
rich.print(deserialized_config)