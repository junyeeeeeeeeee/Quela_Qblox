
import sys
sys.path.append("/Users/LT/Desktop/Wei_en_program/Qblox")
from support import *
from Save import *  
from tqdm import tqdm
from numpy import arange, linspace, logspace
from quantify_scheduler.backends.graph_compilation import SerialCompiler
#%%
meas_datadir = '.data'
dh.set_datadir(meas_datadir)
# Connect to Cluster
warnings.simplefilter("ignore")
connect, ip = connect_clusters()
Instrument.close_all()              # Close all existing QCoDeS Instrument instances
dev_id = connect.value   

cluster = Cluster(name = "cluster0", identifier = ip.get(dev_id))
#cluster = Cluster(name = "cluster0",identifier = f"qum.phys.sinica.edu.tw", port=5025)
#cluster = Cluster(name = "cluster171",identifier = f"qum.phys.sinica.edu.tw", port=5171)
print(f"{connect.label} connected")

# Reset the cluster
cluster.reset()        
print(cluster.get_system_state())


#%%
# Hardware settings
hardware_cfg= {
    "backend": "quantify_scheduler.backends.qblox_backend.hardware_compile",
    "cluster0": {
        "sequence_to_file": False,  # Boolean flag which dumps waveforms and program dict to JSON file
        "ref": "internal",  # Use shared clock reference of the cluster
        "instrument_type": "Cluster",
        # ============ DRIVE ============#
        "cluster0_module4": {
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
        "cluster0_module6": {
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
        # ============ FLUX ============#
        "cluster0_module2": {
            "instrument_type": "QCM",
            "real_output_0": {"portclock_configs": [{"port": "q0:fl", "clock": "cl0.baseband"}]},
            "real_output_1": {"portclock_configs": [{"port": "q1:fl", "clock": "cl0.baseband"}]},
            "real_output_2": {"portclock_configs": [{"port": "q2:fl", "clock": "cl0.baseband"}]},
            "real_output_3": {"portclock_configs": [{"port": "qc0:fl", "clock": "cl0.baseband"}]},
        },
        # ============ READOUT ============#
        "cluster0_module8": {
            "instrument_type": "QRM_RF",
            "complex_output_0": {
                "output_att": 0,
                "input_att": 0,
                "dc_mixer_offset_I": 0.0,
                "dc_mixer_offset_Q": 0.0,
                "lo_freq": 6e9,       # *** Should be set as a parameter later on
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
                ],
            },
        },
    },
}
#%%
# quantum device, measctrler, coordinator set up 
quantum_device = create_quantum_device(hardware_cfg, num_qubits=4)


meas_ctrl, instrument_coordinator = configure_measurement_control_loop(quantum_device, cluster)


# Flux bias setting
flux_settable_map:callable = {
            "q0":cluster.module2.out0_offset,
            "q1":cluster.module2.out1_offset,
            "q2":cluster.module2.out2_offset,
            "qc0":cluster.module2.out3_offset,
        }

# default the offset in circuit
for i in flux_settable_map:
    flux_settable_map[i](0.00)


# Activate NCO delay compensation
for i in range(6):
    getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
    getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)

#%%
    


# Attenuation setting. We don't change it once we set it.
ro_out_att = 40
xy_out_att = 0
# atte. setting
set_atte_for(quantum_device,ro_out_att,'ro',['q0', 'q1', 'q2','q3'])
set_atte_for(quantum_device,xy_out_att,'xy',['q0', 'q1', 'q2','q3']) 
#LO debug fixing 
enable_QCMRF_LO(cluster)


#%%

    
elements = list(sorted({'q0','q1','q2'}))#list(flux_settable_map.keys())

#=========Resonant frequnecy guess=========#

ro_bare=dict(
    q0 = 5.762 * 1e9,
    q1 = 5.8 * 1e9,
    q2 = 5.9 * 1e9,
    q3 = 6 * 1e9,
    q4 = 6.1 * 1e9,
)


f01_guess=dict(
    q0 = 4.3 * 1e9,
    q1 = 4.232816338702908 * 1e9,
    q2 = 3.8418045414828643 * 1e9,
    q3 = 4.022 * 1e9,
    q4 = 2.5738611635902258 * 1e9,
)

#=========initial parameters=========#

bias = dict(
    q0 = 0,
    q1 = 0,
    q2 = 0,
    q3 = 0,
    q4 = 0,
)



pi_amp=dict(
    q0 = 0.1,
    q1 = 0,
    q2 = 0,
    q3 = 0,
    q4 = 0,
)

f01=dict(
    q0 = 4.3 * 1e9,
    q1 = 4.232816338702908 * 1e9,
    q2 = 3.8418045414828643 * 1e9,
    q3 = 4.022 * 1e9,
    q4 = 2.5738611635902258 * 1e9,
)

R_F=dict(
    q0 = 5.7 * 1e9,
    q1 = 5.8 * 1e9,
    q2 = 5.9 * 1e9,
    q3 = 6 * 1e9,
    q4 = 6.1 * 1e9,
)

R_amp=dict(
    q0 = 0.02,
    q1 = 0.1,
    q2 = 0.2,
    q3 = 0.2,
    q4 = 0.2,
)
R_duration=dict(
    q0 = 5* 1e-6,
    q1 = 5* 1e-6,
    q2 = 5* 1e-6,
    q3 = 5* 1e-6,
    q4 = 5* 1e-6,
)

R_integration=dict(
    q0 = 4* 1e-6,
    q1 = 4* 1e-6,
    q2 = 4* 1e-6,
    q3 = 4* 1e-6,
    q4 = 4* 1e-6,
)

I_ref=dict(
    q0 = 0,
    q1 = 0,
    q2 = 0,
    q3 = 0,
    q4 = 0,
)
Q_ref=dict(
    q0 = 0,
    q1 = 0,
    q2 = 0,
    q3 = 0,
    q4 = 0,
)

R_inte_delay= 200e-9
Reset_time= 400e-6



def Experi_info(q,glob):
    return dict(qubit=q,
                Z_bias_offset='%E' %glob['bias'][q],
                f01='%E' %glob['f01'][q],
                pi_amp='%E' %glob['pi_amp'][q],
                ReadoutF='%E' %glob['R_F'][q],
                Readout_amp='%E' %glob['R_amp'][q],
                Readout_Du='%E' %glob['R_duration'][q],
                R_integration_Du='%E' %glob['R_integration'][q],
                R_inte_delay='%E' %glob['R_inte_delay'],
                Reset_time= '%E' %glob['Reset_time'],
                I_ref='%E' %glob['I_ref'][q],
                Q_ref='%E' %glob['Q_ref'][q])


def manual_update(glob,q,para_key,update_value):
    if para_key== 'R_inte_delay' or para_key=='Reset_time':
        glob[para_key] =update_value
    else:   
        glob[para_key][q] =update_value
    qubit = quantum_device.get_element(q)
    qubit.clock_freqs.f01(glob['f01'][q])
    qubit.clock_freqs.readout(glob['R_F'][q])
    qubit.reset.duration(glob['Reset_time'])
    flux_settable_map[q](glob['bias'][q])


executor = True
if executor: # initial setting
    for q in elements:
        qubit = quantum_device.get_element(q)
        qubit.clock_freqs.f01(f01_guess[q])
        qubit.clock_freqs.readout(ro_bare[q])
        qubit.reset.duration(Reset_time)
        show_args(Experi_info(q,globals()),title='')
        
        
Save_filepath= "/Users/LT/Desktop/QC data/MMD_202403/Q1/"

#save_exper_info(Save_filepath+"initial_parameters_save",globals())


#  manual_update(globals(),q,para_key,update_value)
#  show_args(Experi_info(q),title='Current parameters')

#load_data()


