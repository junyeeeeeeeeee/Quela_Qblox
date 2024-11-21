"""
Use the results from m1 and a light attenuation (10 ~ 16 is recommended) to find the BARE cavity frequency.\n
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from xarray import Dataset
import matplotlib.pyplot as plt
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from numpy import array, arange, real, imag, arctan2,column_stack
from quantify_core.measurement.control import MeasurementControl
from qcat.analysis.resonator.photon_dep.res_data import ResonatorData
from qblox_drive_AS.support import Data_manager, QDmanager, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import One_tone_multi_sche, pulse_preview


def Cavity_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:dict,n_avg:int=10,run:bool=True)->Dataset:
    """
        Doing the multiplexing cavity search according to the arg `ro_elements`\n
        Please fill up the initial value about measure for qubit in QuantumDevice first, like: amp, duration, integration_time and acqusition_delay!\n
        ----
        ### Args:\n
        * ro_elements: {'q0': data point array, 'q1':[], ...}\n
        ----
        ## Warning:\n
        The sweep frequency data-point for each cavity in ro_elements is better to set equally.
    """
    datapoint_idx = arange(0,len(list(list(ro_elements.values())[0])))

    quantum_device = QD_agent.quantum_device
    sche_func = One_tone_multi_sche
    for q in ro_elements:
        quantum_device.get_element(q).clock_freqs.readout(NaN) # avoid cluster clock warning
    

    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True


    spec_sched_kwargs = dict(   
        frequencies=ro_elements,
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'r1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'r3'),
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'r4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'r2'),
        powerDep=False,
    )
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(ro_elements.keys())),
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(datapoint_idx)
        
        rs_ds = meas_ctrl.run("One-tone")
        dict_ = {}
        for q_idx, q in enumerate(list(ro_elements.keys())):
            i_data = array(rs_ds[f'y{2*q_idx}'])
            q_data = array(rs_ds[f'y{2*q_idx+1}'])
            dict_[q] = (["mixer","freq"],array([i_data,q_data]))
            dict_[f'{q}_freq'] = (["mixer","freq"],array([ro_elements[q],ro_elements[q]]))
        
        DS = Dataset(dict_,coords={"mixer":array(["I","Q"]),"freq":datapoint_idx})
      
    else:
        n_s = 2
        preview_para = {}
        for q in ro_elements:
            preview_para[q] = ro_elements[q][:n_s]
        
        spec_sched_kwargs['frequencies']= preview_para
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        DS = {}

    return DS

def QD_RO_init(QD_agent:QDmanager, ro_elements:dict):
    for target_q in list(ro_elements.keys()):
        qubit = QD_agent.quantum_device.get_element(target_q)
        qubit.measure.pulse_duration(100e-6)
        qubit.measure.integration_time(100e-6)


def multiplexing_CS_ana(QD_agent:QDmanager, ds:Dataset, save_pic_folder:str=None)->dict:
    """
    # Return\n
    A dict sorted by q_name with its fit results.\n
    Ex. {'q0':{..}, ...}\n
    ----------------------------
    # fit results key names: \n
    ['Qi_dia_corr', 'Qi_no_corr', 'absQc', 'Qc_dia_corr', 'Ql', 'fr', 'theta0', 'phi0', 'phi0_err', 'Ql_err', 'absQc_err', 'fr_err', 'chi_square', 'Qi_no_corr_err', 'Qi_dia_corr_err', 'A', 'alpha', 'delay', 'input_power']
    """
    fit_results = {}

    for idx, q in enumerate(ds.data_vars):
        if str(q).split("_")[-1] != "freq":
            S21 = array(ds[q])[0] + array(ds[q])[1]*1j
            freq = array(ds[f"{q}_freq"])[0][5:]
            res_er = ResonatorData(freq=freq,zdata=array(S21)[5:])
            result, data2plot, fit2plot = res_er.fit()
            fig, ax = plt.subplots(2,2,figsize=(12,12))
            ax0:plt.Axes = ax[0][0] 
            ax0.grid()       
            ax0.plot(freq,result['A']*abs(data2plot))
            ax0.plot(freq,result['A']*abs(fit2plot),c="red",label='fitting')
            ax0.vlines(result['fr'],result['A']*min(data2plot),result['A']*max(data2plot),linestyles="--")
            ax0.set_title(f"{q} cavity @ {round(float(result['fr'])*1e-9,5)} GHz")
            ax0.legend()
            ax1:plt.Axes = ax[0][1] 
            ax1.grid()       
            ax1.plot(freq,arctan2(imag(data2plot),real(data2plot)))
            ax1.plot(freq,arctan2(imag(fit2plot),real(fit2plot)),c="red",label='fitting')
            ax1.set_title("Phase")
            ax1.legend()
            ax2:plt.Axes = ax[1][0]  
            ax2.grid()      
            ax2.scatter(real(array(S21)[1:]),imag(array(S21)[1:]),label='data')
            ax2.set_title("S21 raw data")
            ax2.legend()
            ax3:plt.Axes = ax[1][1]  
            ax3.grid()      
            ax3.scatter(result['A']*real(data2plot),result['A']*imag(data2plot),label='data')
            ax3.scatter(result['A']*real(fit2plot),result['A']*imag(fit2plot),label='fit',c='red',s=10)
            ax3.set_title("S21 after fit")
            ax3.legend()
            plt.tight_layout()
            if save_pic_folder is not None:
                Data_manager().save_multiplex_pics(QD_agent,q,'cs',fig,save_pic_folder)
                plt.close()
            else:
                plt.show()
                
            fit_results[q] = result

    return fit_results


def CS_ana(QD_agent:QDmanager, cs_ds:Dataset, pic_save_folder:str=None, keep_bare:bool=True):
    Quality_values = ["Qi_dia_corr", "Qc_dia_corr", "Ql"]
    Quality_errors = ["Qi_dia_corr_err", "absQc_err", "Ql_err"]
    CS_results = multiplexing_CS_ana(QD_agent,cs_ds, Data_manager().get_today_picFolder() if pic_save_folder is None else pic_save_folder)
    for qubit in CS_results:
        qu = QD_agent.quantum_device.get_element(qubit)
        qu.clock_freqs.readout(float(CS_results[qubit]['fr']))
        if keep_bare:
            QD_agent.Notewriter.save_bareFreq_for(target_q=qubit,bare_freq=CS_results[qubit]['fr'])
        print(f"{qubit}:")
        print("Res @ ",round(qu.clock_freqs.readout()*1e-9,4)," GHz")
        for Qua_idx, Qua in enumerate(Quality_values):
            print(f"{Qua[:2]} = {round(float(CS_results[qubit][Qua])/1000,2)} 土 {round(float(CS_results[qubit][Quality_errors[Qua_idx]])/1000,2)} k")




    
    

