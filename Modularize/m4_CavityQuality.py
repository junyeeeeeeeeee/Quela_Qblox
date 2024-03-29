from Modularize.m2_CavitySpec import Cavity_spec
from Modularize.support import Data_manager, QDmanager
from quantify_core.measurement.control import MeasurementControl
from Modularize.support import init_meas, init_system_atte, shut_down

def qualityFit_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_amp:float,specific_qubits:str,ro_span_Hz:float=10e6,run:bool=True,f_shifter:float=0):
    rof = {str(specific_qubits):QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout()+f_shifter}
    init_system_atte(QD_agent.quantum_device,[specific_qubits],ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(specific_qubits,'ro'))
    if run:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=True)[specific_qubits]
    else:
        qb_CSresults = Cavity_spec(QD_agent,meas_ctrl,ro_bare_guess=rof,ro_amp=ro_amp,q=specific_qubits,ro_span_Hz=ro_span_Hz,run=False)[specific_qubits]
    
    QD_agent.quantum_device.get_element(specific_qubits).clock_freqs.readout(rof[specific_qubits]-f_shifter)

    return qb_CSresults

def show_quality_for(CS_results:dict,target_q:str)->dict:
    qi = round(CS_results[target_q].quantities_of_interest['Qi'].nominal_value*1e-3,2) # in k
    ql = round(CS_results[target_q].quantities_of_interest['Ql'].nominal_value*1e-3,2)
    qc = round(CS_results[target_q].quantities_of_interest['Qc'].nominal_value*1e-3,2)

    print(f"{target_q}: Qi= {qi}k, Qc= {qc}k, Ql= {ql}k")
    return {"QI":qi,"QC":qc,"QL":ql}

if __name__ == "__main__":

    """ Fill in """
    execution = True
    QD_path = 'Modularize/QD_backup/2024_3_29/DR2#171_SumInfo.pkl'
    ro_elements = {
        "q4":{"ro_amp":0.7}
    }
    freq_shift = -3e6


    """ Preparations """ 
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')


    """ Running """
    CS_results = {}
    for qubit in ro_elements:
        CS_results[qubit] = qualityFit_executor(QD_agent=QD_agent,meas_ctrl=meas_ctrl,specific_qubits=qubit,ro_amp=ro_elements[qubit]["ro_amp"],run = execution, f_shifter=freq_shift)
        cluster.reset()
        print(f"{qubit}: Cavity @ {round(CS_results[qubit].quantities_of_interest['fr'].nominal_value*1e-9,5)} GHz")
        _ = show_quality_for(CS_results,qubit)

    """ Storing (future) """


    """ Close """
    print('Cavity quality fit done!')
    shut_down(cluster,Fctrl)
