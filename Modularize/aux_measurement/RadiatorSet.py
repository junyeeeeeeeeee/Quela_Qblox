import os, sys, json, time
from datetime import datetime
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_instruments import Cluster
from Modularize.support.Path_Book import meas_raw_dir
from Modularize.m13_T1  import T1_executor
from Modularize.m12_T2  import ramsey_executor
from Modularize.m14_SingleShot import SS_executor
from Modularize.support import QDmanager, init_meas, shut_down, init_system_atte, coupler_zctrl
from quantify_core.measurement.control import MeasurementControl

def create_set_folder(parent_dir:str,folder_idx:int):
    folder_name = f"Radiator({folder_idx})"
    new_folder_path = os.path.join(parent_dir,folder_name)
    os.mkdir(new_folder_path)
    print(f"dir '{folder_name}' had been created!")
    return new_folder_path

def create_temperature_folder(temperature:str,within_specific_path:str="")->str:
    if within_specific_path != "":
        temp_folder_path = os.path.join(within_specific_path,temperature)
    else:
        temp_folder_path = os.path.join(meas_raw_dir,temperature)
    os.mkdir(temp_folder_path)

    return temp_folder_path



def radiation_test(QD_agent:QDmanager,cluster:Cluster,meas_ctrl:MeasurementControl,Fctrl:dict,specific_qubits:str,freeDura:dict,T2_detu:float=0e6,histo_counts:int=10,run:bool=True,new_folder:str=''):
    # do T1
    _, mean_T1_us, std_T1_us = T1_executor(QD_agent,cluster,meas_ctrl,Fctrl,specific_qubits,freeDura=freeDura["T1"],histo_counts=histo_counts,run=run,specific_folder=new_folder)
    # do T2
    _, mean_T2_us, _, average_actual_detune = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,specific_qubits,artificial_detune=T2_detu,freeDura=freeDura["T2"],histo_counts=histo_counts,run=run,plot=False,specific_folder=new_folder)
    # do single shot
    for ith in range(histo_counts):
        SS_executor(QD_agent,cluster,Fctrl,specific_qubits,execution=run,data_folder=new_folder,exp_label=ith,plot=False)

def time_monitor(monitoring_info:dict, other_info:dict, qubit:str, data_parent_dir:str, start_time):
    from Modularize.analysis.Radiator.RadiatorSetAna import a_set_analysis,live_time_monitoring_plot
    monitoring_info = a_set_analysis(set_folder,monitoring_info,other_info[qubit]["refIQ"],other_info[qubit]["f01"])
    monitoring_info["x_minutes"].append(round((time.time()-start_time)/60,1))
    live_time_monitoring_plot(monitoring_info,data_parent_dir)
    return monitoring_info

if __name__ == "__main__":
    # 2 sets, 2 histo_counts, take 2.7 mins
    """ fill in """
    Temp = '0K'                        # avoid named start with 're', this name is only for radiator off. If it's for reference please use 4K as temperature name 
    QD_path = 'Modularize/QD_backup/2024_5_13/DR1#11_SumInfo.pkl'
    ro_elements = {
        "q0":{"T2detune":0e6,"freeTime":{"T1":120e-6,"T2":20e-6},"histo_counts":2} # histo_counts min = 2 when for test
    }
    couplers = ['c0','c1']
    tracking_time_min = "free"         # if you wanna interupt it manually, set 'free'


    """ Preparations """
    data_parent_dir = create_temperature_folder(Temp)
    exp_start_time = datetime.now()
    exp_start_time = f"{exp_start_time.strftime('%Y-%m-%d')} {exp_start_time.strftime('%H:%M')}"
    start = time.time()
    cut_time = time.time()
    dr = os.path.split(QD_path)[-1].split("#")[0].lower()

    """ Running """
    other_info= {}
    if str(tracking_time_min).lower() == 'free':
        tracking_time_min = 500 * 24 * 60 # keep running for 500 days, waiting interupted manually
  
    for qubit in ro_elements:
        set_idx = 0
        monitoring_info = {}
        while (cut_time-start)/60 < tracking_time_min:
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
            Cctrl = coupler_zctrl(dr,cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            if set_idx == 0:
                other_info[qubit]={"start_time":exp_start_time,"refIQ":QD_agent.refIQ[qubit],"time_past":[],"f01":QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()}
            else:
                with open(os.path.join(data_parent_dir,"otherInfo.json")) as JJ:
                    other_info = json.load(JJ)
            
            set_folder = create_set_folder(parent_dir=data_parent_dir,folder_idx=set_idx)
            evoT = ro_elements[qubit]["freeTime"]
            ramsey_detune = ro_elements[qubit]["T2detune"]
            histo_count = ro_elements[qubit]["histo_counts"]

            radiation_test(QD_agent, cluster, meas_ctrl, Fctrl, qubit, freeDura=evoT, T2_detu=ramsey_detune, histo_counts=histo_count, new_folder=set_folder)
            # time_monitor(monitoring_info, other_info, qubit, data_parent_dir, start)
            cut_time = time.time()
            other_info[qubit]["time_past"].append(cut_time-start)
            """ Close """
            shut_down(cluster,Fctrl,Cctrl)
            set_idx += 1

            """ Storing """
            with open(os.path.join(data_parent_dir,"otherInfo.json"),"w") as record_file:
                json.dump(other_info,record_file)

    end = time.time()
    print(f"{set_idx}*{ro_elements['q0']['histo_counts']} Cost time: {round((end-start)/60,1)} mins")

    