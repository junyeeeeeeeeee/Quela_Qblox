import os, sys, json, time
from datetime import datetime
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qblox_instruments import Cluster
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.m13_T1  import T1_executor
from Modularize.m12_T2  import ramsey_executor
from Modularize.m14_SingleShot import SS_executor
from Modularize.support import QDmanager, init_meas, shut_down, init_system_atte
from quantify_core.measurement.control import MeasurementControl

def create_special_folder(parent_dir:str,folder_idx:int):
    folder_name = f"Radiator({folder_idx})"
    new_folder_path = os.path.join(parent_dir,folder_name)
    os.mkdir(new_folder_path)
    print(f"dir '{folder_name}' had been created!")
    return new_folder_path


def radiation_test(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,freeDura:dict,T2_detu:float=0e6,histo_counts:int=10,run:bool=True,new_folder:str=''):
    # do T1
    _, mean_T1_us, std_T1_us = T1_executor(QD_agent,cluster,meas_ctrl,Fctrl,specific_qubits,freeDura=freeDura["T1"],histo_counts=histo_counts,run=run,specific_folder=new_folder)
    # do T2
    _, mean_T2_us, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,specific_qubits,artificial_detune=T2_detu,freeDura=freeDura["T2"],histo_counts=histo_counts,run=run,plot=False,specific_folder=new_folder)
    # do single shot
    for ith in range(histo_counts):
        SS_executor(QD_agent,cluster,Fctrl,specific_qubits,execution=run,data_folder=new_folder,exp_label=ith,plot=False)


if __name__ == "__main__":
    # 2 sets, 2 histo_counts, take 2.7 mins
    """ fill in """
    Temp = '0K-1'
    QD_path = 'Modularize/QD_backup/2024_5_13/DR1#11_SumInfo.pkl'
    ro_elements = {
        "q0":{"T2detune":0e6,"freeTime":{"T1":120e-6,"T2":20e-6},"histo_counts":2} # histo_counts min = 2 when for test
    }
    data_parent_dir = os.path.join(meas_raw_dir,Temp)
    tracking_time_min= 0.1 * 60

    """ Preparations """
    os.mkdir(data_parent_dir)
    exp_start_time = datetime.now()
    exp_start_time = f"{exp_start_time.strftime('%Y-%m-%d')} {exp_start_time.strftime('%H:%M')}"
    start = time.time()
    cut_time = time.time()

    """ Running """
    other_info= {}
  
    for qubit in ro_elements:
        set_idx = 0
        while (cut_time-start)/60 < tracking_time_min:
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
        
            if set_idx == 0:
                other_info[qubit]={"start_time":exp_start_time,"refIQ":QD_agent.refIQ[qubit],"time_past":[],"f01":QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()}
            
            set_folder = create_special_folder(parent_dir=data_parent_dir,folder_idx=set_idx)
            evoT = ro_elements[qubit]["freeTime"]
            ramsey_detune = ro_elements[qubit]["T2detune"]
            histo_count = ro_elements[qubit]["histo_counts"]

            radiation_test(QD_agent, cluster, meas_ctrl, Fctrl, qubit, freeDura=evoT, T2_detu=ramsey_detune, histo_counts=histo_count, new_folder=set_folder)
            cut_time = time.time()
            other_info[qubit]["time_past"].append(cut_time-start)
            """ Close """
            shut_down(cluster,Fctrl)
            set_idx += 1



    """ Storing """
    with open(f"{data_parent_dir}/otherInfo.json","w") as record_file:
        json.dump(other_info,record_file)


    end = time.time()
    print(f"{set_idx}*{ro_elements['q0']['histo_counts']} Cost time: {round((end-start)/60,1)} mins")

    