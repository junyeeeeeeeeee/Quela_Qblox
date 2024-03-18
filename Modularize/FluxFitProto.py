from numpy import array
from Modularize.FluxQubit import Zgate_two_tone_spec
from Modularize.support import QDmanager, Data_manager
from Modularize.QuFluxFit import get_arrays_from_netCDF, sortAndDecora, plot_HeatScat, calc_Gcoef_inFbFqFd, calc_fq_g_excluded


if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from numpy import pi, absolute
    
    # Reload the QuantumDevice or build up a new one
    QD_path = 'Modularize/QD_backup/2024_3_14/DR2#171_SumInfo.pkl'
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    
    # Set system attenuation
    # init_system_atte(QDmanager.quantum_device,list(Fctrl.keys()))
    for i in range(6):
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp_en(True)
        getattr(cluster.module8, f"sequencer{i}").nco_prop_delay_comp(50)
    
    execute = True
    for qb in ["q1"]:
        # for i in Fctrl:
        #     if i != qb:
        #         tuneaway = QDmanager.Fluxmanager.get_tuneawayBiasFor(i)
        #         if abs(tuneaway) <= 0.3:
        #             Fctrl[i](tuneaway)
        #         else:
        #             raise ValueError(f"tuneaway bias wrong! = {tuneaway}")

        center = QD_agent.Fluxmanager.get_sweetBiasFor(target_q=qb)
        half_period = QD_agent.Fluxmanager.get_PeriodFor(target_q=qb)/8
        window_shifter = 0
        results, origin_f01, path = Zgate_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=center-half_period+window_shifter,Z_amp_end=center+half_period+window_shifter,q=qb,run=execute)
        reset_offset(Fctrl)
        print('Flux qubit done!')
        shut_down(cluster,Fctrl)


