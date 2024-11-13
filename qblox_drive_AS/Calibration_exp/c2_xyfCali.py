import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from numpy import array
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager

from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from qblox_drive_AS.SOP.m12_T2 import T2_waiter,ramsey_executor
from qblox_drive_AS.analysis.raw_data_demolisher import T2_dataReducer
from qblox_drive_AS.analysis.Multiplexing_analysis import Multiplex_analyzer

def XYF_calibration_waiter(ro_element:dict, initialized_detune:float=20e6,evoT=0.5e-6, actual_detune:dict=None):
    if actual_detune is None:
        for q in ro_element:
            ro_element[q]["evo_T"] = evoT
            ro_element[q]["detune"] = initialized_detune
    else:
        for q in ro_element:
            ro_element[q]["evo_T"] = evoT
            ro_element[q]["detune"] = actual_detune[q]

    return ro_element


if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    DRandIP = {"dr":"dr2","last_ip":"10"}
    ro_elements = {
        "q0":{"xy_IF":250e6},
        "q1":{"xy_IF":250e6},
    }
    couplers = []


    """ Iteration (Do NOT touch!)"""
    ana_target = ["c2","m12"]
    avg_n = 1000
    initializing_detune = 10e6
    evoT = 0.5e-6

    
    """ Preparation """
    actual_detune = {}
    for idx, exp in enumerate(ana_target):
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        for q in ro_elements:
            init_system_atte(QD_agent.quantum_device,list([q]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(q,'xy'))
        if idx == 0:
            new_ro_elements = T2_waiter(QD_agent,XYF_calibration_waiter(ro_elements,initializing_detune,evoT),200)
        else:
            new_ro_elements = T2_waiter(QD_agent,XYF_calibration_waiter(ro_elements,None,20e-6,actual_detune),100)

        """ Running """
        nc_path = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,new_ro_elements,repeat=1,ith=0,run=execution,avg_n=avg_n,specific_folder='')


        """ Analysis """
        ds = T2_dataReducer(nc_path)
        for var in ds.data_vars:
            if var.split("_")[-1] != 'x':
                time_data = array(ds[f"{var}_x"])[0][0]
                ANA = Multiplex_analyzer(exp)
                ANA._import_data(ds[var],var_dimension=2,refIQ=QD_agent.refIQ[var],fq_Hz=QD_agent.quantum_device.get_element(var).clock_freqs.f01()+actual_detune[var] if idx == 1 else None)
                ANA._import_2nddata(time_data)
                ANA._start_analysis()

                fit_pic_folder = Data_manager().get_today_picFolder()
                ANA._export_result(fit_pic_folder)
                if idx == 0: 
                    actual_detune[var] = initializing_detune-ANA.fit_packs['freq']
                    highlight_print(f"{var}: actual detune = {round(actual_detune[var]*1e-6,4)} MHz")
                else:
                    highlight_print(f"Calibrated {var} detune = {round(ANA.fit_packs['freq']*1e-6,4)} MHz")


        """ Store """
        if idx == 1:
            permission = mark_input(f"Update these freqs for both qubits ? [y/n] or an update target qubit name ").lower()
            match permission[0]:
                case 'y':
                    for q in actual_detune:
                        QD_agent.quantum_device.get_element(q).clock_freqs.f01(QD_agent.quantum_device.get_element(q).clock_freqs.f01()+actual_detune[q])
                    QD_agent.QD_keeper()
                case "q":
                    if permission in list(actual_detune.keys()):
                        QD_agent.quantum_device.get_element(permission).clock_freqs.f01(QD_agent.quantum_device.get_element(permission).clock_freqs.f01()+actual_detune[permission])
                        QD_agent.QD_keeper()
                case _:
                    print("Do nothing ~")


        """ Close """
        print('XYF calibration done!')
        shut_down(cluster,Fctrl)

    





            

        
        
            
                
           