import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from numpy import arange, linspace, logspace
from quantify_scheduler.backends.graph_compilation import SerialCompiler
### Note: The qubits are labeled with the physical position on the chip from the left, and starts from 0.

#%%



#%%
# Hardware settings
Hcfg_170 = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    "cluster170": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        "cluster170_module4": {
            "instrument_type": "QCM_RF",
            "complex_output_0": {
                "output_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 3e9,
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
        "cluster170_module6": {
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
        # "cluster0_module12": {
        #     "instrument_type": "QCM_RF",
        #     "complex_output_0": {
        #         "output_att": 0,
        #         "dc_mixer_offset_I": 0.0,
        #         "dc_mixer_offset_Q": 0.0,
        #         "lo_freq": 4e9,
        #         "portclock_configs": [
        #             {
        #                 "port": "q4:mw",
        #                 "clock": "q4.01",
        #                 "mixer_amp_ratio": 1.0,
        #                 "mixer_phase_error_deg": 0.0,
        #             }
        #         ],
        #     },
        # },
        # ============ FLUX ============#
        "cluster170_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_3": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
        },
        # "cluster170_module10": {
        #     "instrument_type": "QCM",
        #     "real_output_0": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]},
        # },
        # ============ READOUT ============#
        "cluster170_module8": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 5.48e9,       # *** Should be set as a parameter later on
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


Hcfg_171 = {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    "cluster171": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        "cluster171_module12": {
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
        "cluster171_module14": {
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
        "cluster171_module16": {
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
        "cluster171_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q0:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_3": {"portclock_configs": [{"port": "q3:fl", "clock": "cl0.baseband"}]},
        },
        "cluster171_module4": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q4:fl", "clock": "cl0.baseband"}]},
        },
        # ============ READOUT ============#
        "cluster171_module8": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 5.86e9,       # *** Should be set as a parameter later on
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

def get_FluxController(cluster, ip_label:str):
    if ip_label == '170':
        Fctrl: callable = {
            "q4":cluster.module2.out0_offset,
            "q1":cluster.module2.out1_offset,
            "q2":cluster.module2.out2_offset,
            "q3":cluster.module2.out3_offset,
            # "q4":cluster.module10.out0_offset
        }
    elif ip_label == '171':
        Fctrl: callable = {
            "q0":cluster.module2.out0_offset,
            "q1":cluster.module2.out1_offset,
            "q2":cluster.module2.out2_offset,
            "q3":cluster.module2.out3_offset,
            "q4":cluster.module4.out0_offset
        }
    else:
        raise KeyError ("please input ip label like '170' or '171'!")
    return Fctrl

# Cluster registerations
clusters_online = {
    "192.168.1.170":{"loc":'DR3',"ser":'00015_2247_002'},
    "192.168.1.171":{"loc":'DR2',"ser":'00015_2406_014'}
}

# Hcfg map
hcfg_map = {"170":Hcfg_170,"171":Hcfg_171}


'''
if __name__ == '__main__':
    meas_datadir = '.data'
    dh.set_datadir(meas_datadir)
    # Connect to Cluster
    warnings.simplefilter("ignore")
    connect, ip = connect_clusters()
    Instrument.close_all()              # Close all existing QCoDeS Instrument instances
    dev_id = connect.value   
    cluster = Cluster(name = "cluster0",identifier = f"qum.phys.sinica.edu.tw", port=5025)
    print(f"{connect.label} connected")

#     # Reset the cluster
#     cluster.reset()        
#     print(cluster.get_system_state())

#     # quantum device, measctrler, coordinator set up 
#     quantum_device = create_quantum_device(hardware_cfg, num_qubits=5)


#     meas_ctrl, instrument_coordinator = configure_measurement_control_loop(quantum_device, cluster)

#     # Flux bias setting
#     flux_settable_map: callable = {
#         "q1":cluster.module2.out0_offset,
#         "q2":cluster.module2.out1_offset,
#         "q3":cluster.module2.out2_offset,
#         "q4":cluster.module2.out3_offset,
#         "q5":cluster.module10.out0_offset
#     }
#     # default the offset in circuit
#     for i in flux_settable_map:
#         flux_settable_map[i](0.00)


#     # Activate NCO delay compensation
#     for i in range(6):
#         getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
#         getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)


        


#     # Attenuation setting. We don't change it once we set it.
#     ro_out_att = 20
#     xy_out_att = 10
#     # atte. setting
#     set_atte_for(quantum_device,ro_out_att,'ro',list(flux_settable_map.keys()))
#     set_atte_for(quantum_device,xy_out_att,'xy',list(flux_settable_map.keys())) 

#     #%%

        
#     ro_element = list(sorted({'q1','q2',}))#list(flux_settable_map.keys())

#     #=========Resonant frequnecy guess=========#

#     ro_bare=dict(
#         q1 = 5.720 * 1e9,
#         q2 = 5.8 * 1e9,
#         q3 = 5.9 * 1e9,
#         q4 = 6 * 1e9,
#         q5 = 6.1 * 1e9,
#     )


#     f01_guess=dict(
#         q1 = 4.123 * 1e9,
#         q2 = 4.232816338702908 * 1e9,
#         q3 = 3.8418045414828643 * 1e9,
#         q4 = 4.022 * 1e9,
#         q5 = 2.5738611635902258 * 1e9,
#     )

#     #=========initial parameters=========#

#     bias = dict(
#         q1 = 0,
#         q2 = 0,
#         q3 = 0,
#         q4 = 0,
#         q5 = 0,
#     )



#     pi_amp=dict(
#         q1 = 0,
#         q2 = 0,
#         q3 = 0,
#         q4 = 0,
#         q5 = 0,
#     )

#     f01=dict(
#         q1 = 4.111 * 1e9,
#         q2 = 4.232816338702908 * 1e9,
#         q3 = 3.8418045414828643 * 1e9,
#         q4 = 4.022 * 1e9,
#         q5 = 2.5738611635902258 * 1e9,
#     )

#     R_F=dict(
#         q1 = 5.720 * 1e9,
#         q2 = 5.8 * 1e9,
#         q3 = 5.9 * 1e9,
#         q4 = 6 * 1e9,
#         q5 = 6.1 * 1e9,
#     )

#     R_amp=dict(
#         q1 = 0.1,
#         q2 = 0.1,
#         q3 = 0.1,
#         q4 = 0.1,
#         q5 = 0.1,
#     )
#     R_duration=dict(
#         q1 = 2* 1e-6,
#         q2 = 2* 1e-6,
#         q3 = 2* 1e-6,
#         q4 = 2* 1e-6,
#         q5 = 2* 1e-6,
#     )

#     R_integration=dict(
#         q1 = 2* 1e-6,
#         q2 = 2* 1e-6,
#         q3 = 2* 1e-6,
#         q4 = 2* 1e-6,
#         q5 = 2* 1e-6,
#     )

#     I_ref=0
#     Q_ref=0

#     R_inte_delay= 0
#     Reset_time= 150e-6

#     def Experi_info(q):
#         return dict(qubit=q,
#                     Z_bias_offset='%E' %bias[q],
#                     f01='%E' %f01[q],
#                     pi_amp='%E' %pi_amp[q],
#                     ReadoutF='%E' %R_F[q],
#                     Readout_amp='%E' %R_amp[q],
#                     Readout_Du='%E' %R_duration[q],
#                     R_integration_Du='%E' %R_integration[q],
#                     R_inte_delay='%E' %R_inte_delay,
#                     Reset_time= '%E' %Reset_time,
#                     I_ref='%E' %I_ref,
#                     Q_ref='%E' %Q_ref)

        

#     executor = True
#     if executor:
#         for q in ro_element:
#             qubit = quantum_device.get_element(q)
#             qubit.clock_freqs.f01(f01_guess[q])
#             qubit.clock_freqs.readout(ro_bare[q])
#             qubit.reset.duration(Reset_time)
#             #flux_settable_map[q](bias[q])
#             show_args(Experi_info(q),title='')

#     #flux_settable_map['q1']() flux check

#     shut_down(cluster,flux_settable_map)

'''

