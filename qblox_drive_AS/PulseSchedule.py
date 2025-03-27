from numpy import array 
from typing import Optional, Literal
from quantify_scheduler.backends.qblox.operations.gate_library import ConditionalReset
from quantify_scheduler.operations.acquisition_library import SSBIntegrationComplex,Trace, ThresholdedAcquisition
from quantify_scheduler.operations.pulse_library import SquarePulse, DRAGPulse
from quantify_scheduler.enums import BinMode
from qblox_drive_AS.support.Pulse_schedule_library import Schedule


electrical_delay = 280e-9

def X_pi_p(sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
    amp= pi_amp[q]
    delay_c= -pi_Du-freeDu
    return sche.add(DRAGPulse(G_amp=amp, D_amp=0, duration= pi_Du, phase=0, port=q+":mw", clock=f"{q}.01",sigma=pi_Du/4),rel_time=delay_c,ref_op=ref_pulse_sche,ref_pt=ref_point)
       

def Readout(sche,q,R_amp,R_duration,powerDep=False):
    if powerDep is True:
        amp= R_amp
        Du= R_duration[q]
    else:    
        amp= R_amp[q]
        Du= R_duration[q]

    return sche.add(SquarePulse(duration=Du,amp=amp,port="q:res",clock=q+".ro",t0=4e-9))

def Multi_Readout(sche,q,ref_pulse_sche,R_amp,R_duration,powerDep=False,):
    if powerDep is True:
        amp= R_amp
        Du= R_duration[q]
    else:    
        amp= R_amp[q]
        Du= R_duration[q]

    return sche.add(SquarePulse(duration=Du,amp=amp,port="q:res",clock=q+".ro",t0=4e-9),ref_pt="start",ref_op=ref_pulse_sche,)

    
def Integration(sche,q,R_inte_delay:float,R_inte_duration,ref_pulse_sche,acq_index,acq_channel:int=0,single_shot:bool=False,get_trace:bool=False,trace_recordlength:float=5*1e-6,
            RO_protocol:Optional[
            Literal[
                "SSBIntegrationComplex",
                "ThresholdedAcquisition",
            ]
        ] = "SSBIntegrationComplex",):
    if single_shot== False:     
        bin_mode=BinMode.AVERAGE
    else: bin_mode=BinMode.APPEND
    # Trace acquisition does not support APPEND bin mode !!!
    if get_trace==False:
        if RO_protocol == "SSBIntegrationComplex":
            return sche.add(SSBIntegrationComplex(
                duration=R_inte_duration[q]-4e-9,
                port="q:res",
                clock=q+".ro",
                acq_index=acq_index,
                acq_channel=acq_channel,
                bin_mode=bin_mode,
                ),rel_time=R_inte_delay
                ,ref_op=ref_pulse_sche,ref_pt="start")
        else:
            return sche.add(ThresholdedAcquisition(
                duration=R_inte_duration[q]-4e-9,
                port="q:res",
                clock=q+".ro",
                acq_index=acq_index,
                acq_channel=acq_channel,
                bin_mode=bin_mode,
                ),rel_time=R_inte_delay
                ,ref_op=ref_pulse_sche,ref_pt="start")
    else:  
        return sche.add(Trace(
                duration=trace_recordlength,
                port="q:res",
                clock=q+".ro",
                acq_index=acq_index,
                acq_channel=acq_channel,
                bin_mode=BinMode.AVERAGE,
                ),rel_time=R_inte_delay
                ,ref_op=ref_pulse_sche,ref_pt="start")


from numpy import linspace
free_evolution_time = {"q0":linspace(0,300e-6,100), "q1":linspace(0,300e-6,100)}
pi_amp = {"q0":0.2, "q1":0.2}
pi_dura = {"q0":200e-9, "q1":200e-9}
R_amp = {"q0":0.2, "q1":0.2}
R_duration = {"q0":10e-6, "q1":10e-6}
R_integration = {"q0":10e-6, "q1":10e-6}
R_inte_delay = {"q0":0, "q1":0}

def PulseSchedule(
        freeduration:dict,
        pi_amp: dict,
        pi_dura:dict,
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        repetitions:int=1,
        singleshot:bool=False 
        )->Schedule:

        qubits2read = list(freeduration.keys())
        sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
        sched = Schedule("T1", repetitions=repetitions)

        for acq_idx in range(sameple_idx):
            for qubit_idx, q in enumerate(qubits2read):
                freeDu = freeduration[q][acq_idx]
    
                sched.add(
                    ConditionalReset(q, acq_index=acq_idx, acq_channel=qubit_idx),
                    label=f"Reset {q} {acq_idx}",
                )
            
                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
            
            
                X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu+electrical_delay)
                
                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=singleshot,get_trace=False,trace_recordlength=0,RO_protocol="ThresholdedAcquisition")
        
        return sched