""" switch waveform in the future """
from quantify_scheduler.operations.pulse_library import DRAGPulse, GaussPulse
from quantify_scheduler.schedules.schedule import Schedule


XY_waveform:str = 'gauss'
RO_waveform:str = "gaussian_edge"
s_factor = 4                     # sigma = puse duration / s_factor
half_pi_ratio:float = 0.5             # pi/2 pulse amp is pi-pulse amp * half_pi_ratio, should be less than 1



class Waveformer():
    """ 
    2024/8/16 start building up for xy-control \n
    ### kwargs name:\n
    1. `drag_ratio`: Conventionally it will be 0.5* anharmonicity ~ 0.5*0.22*pi
    2. `duraOVERsigma`: The default is 4 which means duration/sigma = 4
    3. `halfPI_ratio`: The default is 0.5 which means pi/2 pulse amp is 0.5* pi pulse amp.
    4. `waveform`: It should be in ['gauss', 'DRAG']. Now 'gauss' is the default waveform.
    
    """
    def __init__(self,q_num:int=5,c_num:int=0,**kwargs):
        self.__log = {}
        if "drag_ratio" in list(kwargs.keys()):
            drag_ratio = kwargs["drag_ratio"]
        else:
            drag_ratio = 1
        
        if "duraOVERsigma" in list(kwargs.keys()):
            duraOVERsigma = kwargs["duraOVERsigma"]
        else:
            duraOVERsigma = 4
        
        if "halfPI_ratio" in list(kwargs.keys()):
            halfPI_ratio = kwargs["halfPI_ratio"]
        else:
            halfPI_ratio = 0.5

        if "waveform" in list(kwargs.keys()):
            waveform = kwargs["waveform"]
        else:
            waveform = "gauss"

        for q_idx in range(q_num):
            self.__log[f"q{q_idx}"] = {"waveform":waveform,"duraOVERsigma":duraOVERsigma,"drag_ratio":drag_ratio,"halfPI_ratio":halfPI_ratio} 

    
    def XY_waveform_controller(self, pulse_sche:Schedule, amp:float, duration:float, phase:float, q:str, rel_time:float, ref_op:Schedule, ref_pt:str="start"):
        """
        waveform switch\n
        sigma_factor determines the relations between pulse sigma and pulse duration is `sigma = duration/sigma_factor`.
        ### * kwargs: *
        1. if use DRAG, there should ne a argument named 'drag_amp'. It specifies the derivative part amplitude in DRAG as default is the same as amp.
        """
        match self.__log[q]["waveform"].lower():
            case 'drag':
                diff_amp = float(self.__log[q]["drag_ratio"]) * amp
                return pulse_sche.add(DRAGPulse(G_amp=amp, D_amp=diff_amp, duration= duration, phase=phase, port=q+":mw", clock=q+".01",sigma=duration/float(self.__log[q]["duraOVERsigma"])),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
            case 'gauss':
                return pulse_sche.add(GaussPulse(G_amp=amp, phase=phase,duration=duration, port=q+":mw", clock=q+".01",sigma=duration/float(self.__log[q]["duraOVERsigma"])),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
            case _:
                pass

    
    def get_waveform_log(self):
        return self.__log
    
    def activate(self, old_log:dict):
        self.__log = old_log
    