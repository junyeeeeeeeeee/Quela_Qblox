import os, sys, json, time
from datetime import datetime
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support.Path_Book import meas_raw_dir
from qblox_drive_AS.SingleReadout.m13_T1  import T1_executor
from qblox_drive_AS.SingleReadout.m12_T2  import ramsey_executor
from qblox_drive_AS.SingleReadout.m14_SingleShot import SS_executor
from qblox_drive_AS.support import init_meas, shut_down, init_system_atte, coupler_zctrl


if __name__ == "__main__":

    """ fill in """
    T1_folder_path = 'Modularize\Meas_raw\T1_timeDep'
    T2_folder_path = 'Modularize\Meas_raw\T2_timeDep'
    OS_folder_path = ''

    QD_path = 'Modularize\QD_backup\2024_9_25\DR4#81_SumInfo.pkl'
    ro_elements = {
        "q4":{"T2detune":0.5e6,"freeTime":{"T1":120e-6,"T2":20e-6}} # histo_counts min = 2 when for test
    }
    couplers = []
    tracking_time_min = "free"         # if you wanna interupt it manually, set 'free'

    """ Optional paras """
    doing_exp = {"T1":True,"T2":True,"OS":False}
    XY_IF:float = 250e6
    n_avg:int = 300
    shots = 5e3


    """ Preparations """
    
    exp_start_time = datetime.now()
    exp_start_time = f"{exp_start_time.strftime('%Y-%m-%d')} {exp_start_time.strftime('%H:%M')}"
    start = time.time()
    cut_time = time.time()
    dr = os.path.split(QD_path)[-1].split("#")[0].lower()

    """ Running """
    if str(tracking_time_min).lower() == 'free':
        tracking_time_min = 500 * 24 * 60 # keep running for 500 days, waiting interupted manually

    time_recs = {'T1_times_rec':[],'T2_times_rec':[],'OS_times_rec':[]}
    paths = [T1_folder_path, T2_folder_path, OS_folder_path]
    for qubit in ro_elements:
        set_idx = 0
        while (cut_time-start)/60 < tracking_time_min:
            # build time record json in folder
            if set_idx == 0:
                for idx, folder_path in enumerate(paths):
                    if folder_path != '':
                        with open(os.path.join(folder_path,"timeInfo.json"),"w") as record_file:
                            json.dump({list(time_recs.keys())[idx]:time_recs[list(time_recs.keys())[idx]]},record_file)
            
            for idx, exp in enumerate(list(doing_exp.keys())):
                if doing_exp[exp]:
                    # load the time record
                    if set_idx != 0 and paths[idx] != '':
                        with open(os.path.join(paths[idx],"timeInfo.json")) as recorded_file:
                            time_recs[list(time_recs.keys())[idx]] = json.load(recorded_file)[list(time_recs.keys())[idx]]

                    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)
                    init_system_atte(QD_agent.quantum_device,list(Fctrl.keys()),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
                    Fctrl = coupler_zctrl(Fctrl,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
                    Fctrl['c3'](0.13)

                    if exp == "T1" and T1_folder_path != '' and doing_exp[exp]:
                        _ = T1_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,freeDura=ro_elements[qubit]["freeTime"]["T1"],ith=set_idx,run=True,specific_folder=T1_folder_path,avg_times=n_avg,IF=XY_IF)
                        
                    elif exp == "T2" and T2_folder_path != '' and doing_exp[exp]:
                        _ = ramsey_executor(QD_agent,cluster,meas_ctrl,Fctrl,qubit,artificial_detune=ro_elements[qubit]["T2detune"],freeDura=ro_elements[qubit]["freeTime"]["T2"],ith=set_idx,run=True,specific_folder=T2_folder_path,avg_n=int(1.5*n_avg),second_phase='y',IF=XY_IF)
                        
                    elif exp == "OS" and OS_folder_path != '' and doing_exp[exp]:
                        SS_executor(QD_agent,cluster,Fctrl,qubit,execution=True,data_folder=OS_folder_path,exp_label=set_idx,plot=False,IF=XY_IF,shots=shots)
                        
                    else:
                        print(f"*** Can't support this exp called '{exp}' in Radiator test set !")
                    now = time.time()
                    time_recs[list(time_recs.keys())[idx]].append((now-start)/60)
                    # save the new time record json
                    if paths[idx] != '':
                        with open(os.path.join(paths[idx],"timeInfo.json"),"w") as recorded_file:
                                json.dump({list(time_recs.keys())[idx]:time_recs[list(time_recs.keys())[idx]]},recorded_file)
                    """ Close """
                    shut_down(cluster,Fctrl)
            set_idx += 1
            
            
            


    