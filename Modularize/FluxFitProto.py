from numpy import array
from Modularize.FluxQubit import Zgate_two_tone_spec
from Modularize.support import QDmanager, Data_manager
from Modularize.QuFluxFit import get_arrays_from_netCDF, sortAndDecora, plot_HeatScat, calc_Gcoef_inFbFqFd, calc_fq_g_excluded


if __name__ == "__main__":
    from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
    from numpy import pi, absolute
    
    
    

