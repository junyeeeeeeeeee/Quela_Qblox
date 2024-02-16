from numpy import array, linspace

from Modulize.Pulse_schedule_library import One_tone_sche, pulse_preview
from utils.tutorial_analysis_classes import ResonatorFluxSpectroscopyAnalysis
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from quantify_scheduler.gettables import ScheduleGettable
from support import QuantumDevice
from quantify_core.measurement.control import MeasurementControl

n_s=2 
def Cavity_spec(qubit_spec:dict,quantum_device:QuantumDevice,meas_ctrl:MeasurementControl,ro_bare_guess:dict,ro_span_Hz:int=5e6,n_avg:int=1000,points:int=200,run:bool=True,q:str='q1',Experi_info:dict={}):
    """
        Do the cavity search by the given QuantumDevice with a given target qubit q
    """
    sche_func = One_tone_sche
        
    analysis_result = {}
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    #TODO: fill the variables with quit_spec (dynamic spec)
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp=R_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        powerDep=False,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     )
    
    if run:
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables(freq)
        meas_ctrl.setpoints(ro_f_samples)
        
        
        
        rs_ds = meas_ctrl.run("One-tone")
        rs_ds
        analysis_result[q] = ResonatorFluxSpectroscopyAnalysis(tuid=rs_ds.attrs["tuid"], dataset=rs_ds).run()
        print(f"{q} Cavity:")
        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        sweep_para= array(ro_f_samples[:n_s])
        spec_sched_kwargs['frequencies']= sweep_para.reshape(sweep_para.shape or (1,))
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        

        show_args(exp_kwargs, title="One_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
    return analysis_result