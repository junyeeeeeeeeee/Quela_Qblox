import os, sys, json, time
from datetime import datetime
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support.Path_Book import meas_raw_dir
from qblox_drive_AS.SOP.m13_T1  import T1_executor
from qblox_drive_AS.SOP.m12_T2  import ramsey_executor
from qblox_drive_AS.SOP.m14_SingleShot import SS_executor
from qblox_drive_AS.support import init_meas, shut_down, init_system_atte, coupler_zctrl

def create_set_folder(parent_dir:str,folder_idx:int):
    folder_name = f"Radiator({folder_idx})"
    new_folder_path = os.path.join(parent_dir,folder_name)
    if not os.path.exists(new_folder_path):
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


def time_monitor(monitoring_info:dict, other_info:dict, qubit:str, data_parent_dir:str, start_time):
    from qblox_drive_AS.analysis.Radiator.RadiatorSetAna import a_set_analysis,live_time_monitoring_plot
    monitoring_info = a_set_analysis(set_folder,monitoring_info,other_info[qubit]["refIQ"],other_info[qubit]["f01"])
    monitoring_info["x_minutes"].append(round((time.time()-start_time)/60,1))
    live_time_monitoring_plot(monitoring_info,data_parent_dir)
    return monitoring_info

if __name__ == "__main__":
    # 2 sets, 2 histo_counts, take 2.7 mins
    """ fill in """
    Temp = '4K'                        # avoid named start with 're', this name is only for radiator off. If it's for reference please use 4K as temperature name 
    special_parent_dir = "Modularize/Meas_raw/Radiator/ScalinQ_Q4"
    QD_path = 'Modularize/QD_backup/2024_6_11/DR1SCA#11_SumInfo.pkl'
    ro_elements = {
        "q0":{"T2detune":0.2e6,"freeTime":{"T1":120e-6,"T2":20e-6},"histo_counts":10} # histo_counts min = 2 when for test
    }
    couplers = ['c0']
    tracking_time_min = "free"         # if you wanna interupt it manually, set 'free'

    """ Optional paras """
    doing_exp = ["T1","T2","OS"]


    """ Preparations """
    data_parent_dir = create_temperature_folder(Temp,within_specific_path=special_parent_dir)
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
        while (cut_time-start)/60 < tracking_time_min:
            set_folder = create_set_folder(parent_dir=data_parent_dir,folder_idx=set_idx)
            for exp_idx, exp in enumerate(doing_exp):
                for ith_histo in range(ro_elements[qubit]["histo_counts"]):
                    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
                    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
                    Cctrl = coupler_zctrl(dr,cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
                    if set_idx == 0 and exp_idx == 0 and ith_histo == 0:
                        other_info[qubit]={"start_time":exp_start_time,"refIQ":QD_agent.refIQ[qubit],"time_past":[],"f01":QD_agent.quantum_device.get_element(qubit).clock_freqs.f01()}
                    
                    if exp == "T1":
                        _ = T1_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,freeDura=ro_elements[qubit]["freeTime"]["T1"],ith=ith_histo,run=True,specific_folder=set_folder)
                    elif exp == "T2":
                        _ = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=ro_elements[qubit]["T2detune"],freeDura=ro_elements[qubit]["freeTime"]["T2"],ith=ith_histo,run=True,specific_folder=set_folder)
                    elif exp == "OS":
                        SS_executor(QD_agent,cluster,Fctrl,qubit,execution=True,data_folder=set_folder,exp_label=ith_histo,plot=False)
                    else:
                        print(f"*** Can't support this exp called '{exp}' in Radiator test set !")
                    
                    """ Close """
                    shut_down(cluster,Fctrl,Cctrl)
            
            
            cut_time = time.time()
            other_info[qubit]["time_past"].append(cut_time-start)
            set_idx += 1

            """ Storing """
            with open(os.path.join(data_parent_dir,"otherInfo.json"),"w") as record_file:
                json.dump(other_info,record_file)

    end = time.time()
    print(f"{set_idx}*{ro_elements['q0']['histo_counts']} Cost time: {round((end-start)/60,1)} mins")

    