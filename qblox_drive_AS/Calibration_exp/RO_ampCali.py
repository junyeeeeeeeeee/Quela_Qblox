import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from numpy import array, arange, moveaxis
from qblox_drive_AS.support import QDmanager, Data_manager, compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import multi_ROL_Cali_sche, pulse_preview

def rolCali(QD_agent:QDmanager,meas_ctrl:MeasurementControl,rol_coef_samples:dict,n_avg:int=500,run:bool=True)->Dataset:
    sche_func= multi_ROL_Cali_sche
    
    for q in rol_coef_samples:
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        rol_data_idx = arange(rol_coef_samples[q].shape[0])
    
    
    option_rol = ManualParameter(name="ROamp", unit="V", label="Amplitude")
    option_rol.batched = True

    def state_dep_sched(ini_state:str):
        sched_kwargs = dict(
            R_amp_coefs=rol_coef_samples,
            ini_state=ini_state,
            pi_amp=compose_para_for_multiplexing(QD_agent,rol_coef_samples,'d1'),
            pi_dura=compose_para_for_multiplexing(QD_agent,rol_coef_samples,'d3'),
            R_duration=compose_para_for_multiplexing(QD_agent,rol_coef_samples,'r3'),
            R_amp = compose_para_for_multiplexing(QD_agent,rol_coef_samples,'r1'),
            R_integration=compose_para_for_multiplexing(QD_agent,rol_coef_samples,'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,rol_coef_samples,'r2'),
            )
        
        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func,
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(rol_coef_samples.keys())),
            )
            
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            meas_ctrl.gettables(gettable)
            meas_ctrl.settables(option_rol)
            meas_ctrl.setpoints(rol_data_idx)
            
            ds = meas_ctrl.run("Rof-Calibrate")
            dict_ = {}
            for q_idx, q in enumerate(rol_coef_samples):
                i_data = array(ds[f'y{2*q_idx}'])
                q_data = array(ds[f'y{2*q_idx+1}'])
                dict_[q] = array([i_data.tolist(),q_data.tolist()]) # shape (mixer, rof_idx)
                dict_[f"{q}_rol"] = array([array(rol_coef_samples[q])]*2) # shape (mixer, rof_idx)
            
            return dict_
    
        else:
            preview_para = {}
            for q in rol_coef_samples:
                preview_para[q] = rol_coef_samples[q][:2]
            sched_kwargs['R_amp_coefs']= preview_para
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

            return {}
    
    slightly_print("Running |0>")
    data_dict_g = state_dep_sched('g')
    slightly_print("Running |1>")
    data_dict_e = state_dep_sched('e')
    rofcali_ds = ""
    if run:
        dict_ = {}
        for var in data_dict_g :
            if var.split("_")[-1] != "rol":
                state_data = []
                rof_data = [array(data_dict_g[f"{var}_rol"]).tolist()]*2
                for item in [array(data_dict_g[var]).tolist(), array(data_dict_e[var]).tolist()]: # shape (state, mixer, rof)
                    state_data.append(item)
                dict_[var] = (["mixer","state","rol"],moveaxis(array(state_data),0,1)) 
                dict_[f"{var}_rol"] = (["mixer","state","rol"],moveaxis(array(rof_data),0,1))
        
        rofcali_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"state":["g","e"],"rol":rol_data_idx})
        rofcali_ds.attrs["execution_time"] = Data_manager().get_time_now()

    
    return rofcali_ds

