import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from numpy import arange
from qblox_drive_AS.SOP.m13_T1  import T1_executor
from qblox_drive_AS.SOP.m12_T2  import ramsey_executor, modify_time_point
from qblox_drive_AS.SOP.m14_SingleShot import SS_executor
from qblox_drive_AS.support import init_meas, shut_down, init_system_atte, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import set_LO_frequency
from qblox_drive_AS.support.QDmanager import QDmanager, Data_manager

def timeMonitor_waiter(QD_agent:QDmanager,ro_element:dict, shots:int, time_pts:int=100)->tuple[dict, dict, dict]:
    """
    ro_element = {"q0":{"T2detune":0.5e6,"freeTime":{"T1":120e-6,"T2":20e-6}, "oneshot":True, "xy_IF":250e6, "spin_num":0},... } for example. \n
    All the keys must be provided.
    """
    T1_elements, T2_elements, OS_elements = {}, {}, {}
    for q in ro_element:
        # T2 settings
        if ro_element[q]["freeTime"]["T2"] != 0:
            if T2_elements == {}:
                T2_elements = {"xyf":{},"spin_num":{},"time_samples":{}}
            # detune settings
            new_freq = QD_agent.quantum_device.get_element(q).clock_freqs.f01() + (ro_element[q]['T2detune'] if 'T2detune' in list(ro_element[q].keys()) else 0)
            T2_elements["xyf"][q] = new_freq
            # echo settings
            T2_elements["spin_num"][q] = ro_element[q]['spin_num'] if 'spin_num' in list(ro_element[q].keys()) else 0
            # evolution time settings
            gap = (ro_element[q]["freeTime"]["T2"])*1e9 // time_pts + (((ro_element[q]["freeTime"]["T2"])*1e9 // time_pts) %(4*(T2_elements["spin_num"][q]+1))) # multiple by 8 ns
            time_samples = modify_time_point(arange(0,ro_element[q]["freeTime"]["T2"],gap*1e-9), T2_elements["spin_num"][q]*8e-9 if T2_elements["spin_num"][q]>=1 else 4e-9)
            T2_elements["time_samples"][q] = time_samples
        # T1 settings
        if ro_element[q]["freeTime"]["T1"] != 0:
            if T1_elements == {}:
                T1_elements = {"time_samples":{}}
            # evolution time settings 
            gap = (ro_element[q]["freeTime"]["T1"]*1e9 // time_pts) + ((ro_element[q]["freeTime"]["T1"]*1e9 // time_pts)%4)
            T1_elements["time_samples"][q] = arange(0,ro_element[q]["freeTime"]["T1"],gap*1e-9)
        # OneShot settings
        if shots != 0 and ro_element[q]["oneshot"] :
            OS_elements[q] = {"ro_amp_factor":1}
            
        # LO settings 
        IF_key = [name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'][0] if len([name for name in list(ro_element[q].keys()) if str(name).lower() == 'xy_if'])==1 else ""      
        IF = ro_element[q][IF_key] if IF_key != "" else 150e6
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=QD_agent.quantum_device.get_element(q).clock_freqs.f01()+IF)

    return T1_elements, T2_elements, OS_elements


if __name__ == "__main__":

    """ fill in """
    QD_path = 'Modularize/QD_backup/20241107/DR2#10_SumInfo.pkl'
    ro_elements = {
        "q0":{"T2detune":0.1e6,"freeTime":{"T1":100e-6,"T2":60e-6},"oneshot":True, "xy_IF":250e6, "spin_num":0}, # all the keies are required exist 
        "q1":{"T2detune":0.1e6,"freeTime":{"T1":100e-6,"T2":60e-6},"oneshot":True, "xy_IF":250e6, "spin_num":0}
    }
    couplers = []
    tracking_time_min = 30        # if you wanna interupt it manually, set it 'free'

    
    """ Optional paras """
    n_avg:int = 300
    time_pts:int = 100
    shot_num:int = 10000


    """ Preparations """
    start = time.time()
    cut_time = time.time()
    specific_folder = Data_manager().creat_datafolder_today(f"TimeMonitor_MultiQ_{Data_manager().get_time_now()}")
    dr = os.path.split(QD_path)[-1].split("#")[0].lower()
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
    T1_elements, T2_elements, OS_elements = timeMonitor_waiter(QD_agent,ro_elements,shot_num,time_pts)


    """ Running """
    if str(tracking_time_min).lower() == 'free':
        tracking_time_min = 500 * 24 * 60 # keep running for 500 days, waiting interupted manually

    set_idx = 0
    while (cut_time-start)/60 < tracking_time_min:
        if set_idx != 0:
            # avoid error on open quantum_device twice 
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
        for qubit in ro_elements:
            init_system_atte(QD_agent.quantum_device,[qubit],xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
        
        # T1
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        if len(list(T1_elements.keys())) != 0:
            _ = T1_executor(QD_agent,cluster,meas_ctrl,Fctrl,T1_elements,ith=set_idx,run=True,specific_folder=specific_folder,avg_times=n_avg)
        
        # T2
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        if len(list(T2_elements.keys())) != 0:
            _ = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,T2_elements,ith=set_idx,run=True,specific_folder=specific_folder,avg_n=int(1.5*n_avg))
        
        # OneShot
        Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
        if len(list(OS_elements.keys())) != 0:
            _ = SS_executor(QD_agent,cluster,Fctrl,OS_elements,execution=True,data_folder=specific_folder,exp_label=set_idx,shots=shot_num)
            

        """ Close """
        cut_time = time.time()
        shut_down(cluster)
        set_idx += 1
        
            
            


    