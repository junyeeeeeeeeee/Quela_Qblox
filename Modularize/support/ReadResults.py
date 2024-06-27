import xarray as xr
import quantify_core.data.handling as dh
from Modularize.support.QDmanager import QDmanager
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from numpy import sqrt, array, moveaxis
from xarray import open_dataset
import matplotlib.pyplot as plt
from Modularize.support.Pulse_schedule_library import dataset_to_array
# from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
# from quantify_core.analysis.base_analysis import Basic2DAnalysis



def plot_QbFlux(Qmanager:QDmanager, nc_path:str, target_q:str):
    ref = Qmanager.refIQ[target_q]
    # plot flux-qubit 
    f,z,i,q = convert_netCDF_2_arrays(nc_path)
    amp = array(sqrt((i-array(ref)[0])**2+(q-array(ref)[1])**2)).transpose()
    fig, ax = plt.subplots()
    c = ax.pcolormesh(z, f, amp, cmap='RdBu')
    fig.colorbar(c, ax=ax)
    plt.show()


# from quantify_scheduler.helpers.collections import find_port_clock_path
# QD_agent = QDmanager('Modularize/QD_backup/2024_3_19/DR2#171_SumInfo.pkl')
# QD_agent.QD_loader()
# print(QD_agent.Fluxmanager.get_sweetBiasFor('q2'))
# print(QD_agent.Fluxmanager.get_PeriodFor('q2'))
# qd = QD_agent.quantum_device
# hcfg = QD_agent.Hcfg
# qubit = qd.get_element('q4')
# output_path = find_port_clock_path(
#         hcfg, port=qubit.ports.readout(), clock=qubit.name + ".ro"
#     )

# cluster_key, module_key, output_key, _, _ = tuple(output_path)
# readout_module = hcfg[cluster_key][module_key]
# print(readout_module[output_key]["output_att"])


# meas_datadir = 'tempt'
# dh.set_datadir(meas_datadir)
# print(ds.attrs["tuid"])
# ana = Basic2DAnalysis(tuid=ds.attrs["tuid"], dataset=ds).run()
# print(ana.quantities_of_interest)



# from Modularize.support import QDmanager
# from numpy import array
# QD_path = 'Modularize/QD_backup/2024_3_31/DR2#171_SumInfo.pkl'
# QD_agent = QDmanager(QD_path)
# QD_agent.QD_loader()
# for i in ["q0","q1","q2","q3","q4"]:
#     qu = QD_agent.quantum_device.get_element(i)
#     rof = qu.clock_freqs.readout()
#     xyf = qu.clock_freqs.f01()
#     xyl = qu.rxy.amp180()
#     T1 = QD_agent.Notewriter.get_T1For(target_q=i)
#     T2 = QD_agent.Notewriter.get_T2For(target_q=i)
#     print(f"***{i}:")
#     print(f"rof : {rof}")
#     print(f"xyf : {xyf}")
#     print(f"xyl : {xyl}")
#     print(f"T1 : {T1}")
#     print(f"T2 : {T2}")
#     print("===========================\n")
# def aprx_fq(disper_MHz:float,bareF_MHz:float,g_MHz:float=45.0):
#     return bareF_MHz-(g_MHz**2)/disper_MHz


# print("aprx fq=",aprx_fq(x,bare))

if __name__ == '__main__':
    # QD_agent = QDmanager('Modularize/QD_backup/2024_6_11/DR1SCA#11_SumInfo.pkl')
    # QD_agent.QD_loader()
    # qubit = QD_agent.quantum_device.get_element("q0")
    # qubit.measure.acq_delay(280e-9)
    # qubit.measure.pulse_duration(1e-6)
    # qubit.measure.integration_time(500e-9)
    # QD_agent.QD_keeper()


    def zgateT1_Qblox2QM_adapter():
        """
        trnaslate the raw data form into the shape (IQ, Z, evoTime)sssssss
        """
        ds = open_dataset("Modularize/Meas_raw/2024_6_27/DR1SCAq0_zT1(1)_H14M8S35.nc")
        i, Q = dataset_to_array(dataset=ds,dims=2)
        new = moveaxis(i,0,-1)
        print(new.shape)