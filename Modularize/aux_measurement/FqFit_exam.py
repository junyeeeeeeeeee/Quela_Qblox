from numpy import array
from qblox_instruments import Cluster
from Modularize.support import QDmanager
from Modularize.m8_Cnti2Tone import Two_tone_spec
from Modularize.m7_RefIQ import Single_shot_ref_spec
from quantify_core.measurement.control import MeasurementControl
from Modularize.support import init_meas, init_system_atte, shut_down
from Modularize.support.Pulse_schedule_library import twotone_comp_plot
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr




def fqExam_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,cluster:Cluster,Fctrl:dict,bias:float,specific_qubits:str,xyf_guess:float,xyAmp_guess:list=[],xyf_span:float=500e6,xy_if:float=100e6,run:bool=True):
    
    if run:
        init_system_atte(QD_agent.quantum_device,list([specific_qubits]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'xy'))
        original_rof = QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout()
        this_rof = QD_agent.Fluxmanager.sin_for_cav(specific_qubits,array([bias]))[0]
        QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(this_rof)
        guess_fq = [xyf_guess-500e6, xyf_guess, xyf_guess+500e6]

        if len(xyAmp_guess) == 0:
            xyAmp_guess = [0, QD_agent.Notewriter.get_2tone_piampFor(specific_qubits)]

        Fctrl[specific_qubits](bias) 
        ref_results = Single_shot_ref_spec(QD_agent,q=specific_qubits,want_state='g',shots=10000)
        Fctrl[specific_qubits](0.0)
        
        I_ref, Q_ref= ref_results[specific_qubits]['fit_pack'][0],ref_results[specific_qubits]['fit_pack'][1]
        ref = [I_ref,Q_ref]

        cluster.reset()
        for XYF in guess_fq:
            ori_data = []
            for XYL in xyAmp_guess:
                # cluster.reset() # *** important
                Fctrl[specific_qubits](bias) 
                QS_results = Two_tone_spec(QD_agent,meas_ctrl,xyamp=XYL,IF=xy_if,f01_guess=XYF,q=specific_qubits,xyf_span_Hz=xyf_span,points=50,n_avg=500,run=True,ref_IQ=ref) # 
                Fctrl[specific_qubits](0.0)
                cluster.reset()
                if XYL != 0:
                    twotone_comp_plot(QS_results[specific_qubits], ori_data, True)
                else:
                    twotone_comp_plot(QS_results[specific_qubits], ori_data, False)
                    ori_data = QS_results[specific_qubits].data_vars['data']

        QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(original_rof)    
        return QS_results[specific_qubits]
                
    else:
        qu = specific_qubits
        QS_results = Two_tone_spec(QD_agent,meas_ctrl,xyamp=0.1,IF=xy_if,f01_guess=4e9,q=qu,xyf_span_Hz=xyf_span,points=50,n_avg=500,run=False,ref_IQ=QD_agent.refIQ[qu])

        return 0


# test record fig with _o.png  (threshold -= 0.5 )
# q0 (3.79) : set 3.6 get 3.59 , set 3.4 get 3.74, set 3.0 get 3.4
# q1 (4.25) : set 4.1 get 4.14 , set 3.8 get 4.0 , set 3.5 get 3.87
# q2 (3.23) : set 3.1 get 2.97 , set 2.7 get 2.57, 
# q3 (3.57) : set 3.4 get 3.25 , set 3.1 get 3.37, set 2.8 get 3.23
# q4 (4.76) : set 4.6 get 4.63 , set 4.3 get 4.71, set 4.0 get 4.49

# test record fig with _oo.png (fitting error < 1/4)
# q0 (3.79) : set 3.6 get 3.59 , set 3.4 get 3.74, set 3.0 get 3.4
# q1 (4.25) : set 4.1 get 4.16 , set 3.8 get 4.02 , set 3.5 get 3.89
# q2 (3.23) : set 3.1 get 3.12 , set 2.7 get 2.57, 
# q3 (3.57) : set 3.4 get 3.24 , set 3.1 get 3.09, set 2.8 get 3.23
# q4 (4.76) : set 4.6 get 4.59 , set 4.3 get 4.71, set 4.0 get 4.49


# test on 5Q4C 5/15
# q0 (4.49) : set 4.3 get 4.41, set 4.0 get 4.29, set 3.7 get 4.1

### 
# assign = [4.5, 4.45, 4.4, 4.35, 4.3, 4.25, 4.2, 4.15, 4.1, 4.05, 4.0, 3.95, 3.9]
# z = [-0.12363310112050828, -0.10091577892001294, -0.08899541749401903, -0.0797748633555727, -0.07198633775825916, -0.06512428584980591, -0.058925813740038364, -0.0532334880965163, -0.04794412618388684, -0.04298576453405387, -0.03830583786473896, -0.033864559428041244, -0.029630949672850578]
# xyf = [4.5077, 4.4852, 4.4634, 4.4396, 4.4165, 4.3948, 4.3691, 4.35, 4.331, 4.311, 4.2894, 4.2672, 4.2484]



if __name__ == "__main__":

    """ Fill in """
    DRandIP = {"dr":"dr3","last_ip":"13"}
    desired_fq = 4.0e9
    target_q = 'q0'
    x_amp = []
    execution = True


    """ Preparations """
    
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    desired_bias = QD_agent.Fluxmanager.get_biasWithFq_from(target_q,desired_fq)
    print(f"fq={desired_fq*1e-9}GHz at bias={desired_bias}V")

    """ Running """
    if desired_bias != 'n':
        tt_results = fqExam_executor(QD_agent,meas_ctrl,cluster,Fctrl,bias=desired_bias,specific_qubits=target_q,xyf_guess=desired_fq,xyAmp_guess=x_amp,run=execution)


    """ Storing """



    """ Close """
    print('Fq Fit examination done!')
    shut_down(cluster,Fctrl)
    
    # assign = [4.5, 4.45, 4.4, 4.35, 4.3, 4.25, 4.2, 4.15, 4.1, 4.05, 4.0, 3.95, 3.9]
    # z = [-0.12363310112050828, -0.10091577892001294, -0.08899541749401903, -0.0797748633555727, -0.07198633775825916, -0.06512428584980591, -0.058925813740038364, -0.0532334880965163, -0.04794412618388684, -0.04298576453405387, -0.03830583786473896, -0.033864559428041244, -0.029630949672850578]
    # xyf = [4.5077, 4.4852, 4.4634, 4.4396, 4.4165, 4.3948, 4.3691, 4.35, 4.331, 4.311, 4.2894, 4.2672, 4.2484]


    
    # import json
    # import matplotlib.pyplot as plt

    # d = {}

    # with open("measANDfit_points.json") as dd:
    #     d = json.load(dd)
    
    # # d["exam"] = {"z":z,"xyf":xyf,"assign":assign}
    # # with open("measANDfit_points.json", "w") as record_file:
    # #     json.dump(d,record_file)


    # plt.rc('font',size=30)
    # plt.plot(d["fit"]["z"],d["fit"]["xyf"],c='orange',label='fitting',lw=2.5)
    # plt.scatter(d["meas"]["z"],d["meas"]["xyf"],s=120*1.5,label='Z2tone data')
    # plt.vlines(x=d["exam"]["z"],ymin=d["exam"]["assign"],ymax=d["exam"]["xyf"],linestyles='--',colors='black')
    # plt.scatter(d["exam"]["z"],d["exam"]["assign"],c='red',marker="*",s=120*2,label='wantted fq')
    # plt.scatter(d["exam"]["z"],d["exam"]["xyf"],c='pink',s=120*1.5,label='2tone get')
        
    # plt.legend()
    # plt.grid()
    # plt.xlabel("Flux (V)")
    # plt.ylabel("XY Frequency (GHz)")
    
    # plt.show()
    # print(d["meas"]["z"])
 