""" switch waveform in the future """
from quantify_scheduler.operations.pulse_library import DRAGPulse, GaussPulse, SquarePulse
from quantify_scheduler.schedules.schedule import Schedule


XY_waveform:str = 'gauss'
RO_waveform:str = "gaussian_edge"
s_factor = 4                     # sigma = puse duration / s_factor
half_pi_ratio:float = 0.5             # pi/2 pulse amp is pi-pulse amp * half_pi_ratio, should be less than 1



class Waveformer():
    """ 
    2024/8/16 start building up for xy-control \n
    #### Arg:\n
    * q_num: int, how many qubit.\n
    * c_num: int, how many coupler.
    
    
    """
    def __init__(self,q_num:int,c_num:int,**kwargs):
        self.__xylog = {}
        self.__zlog = {}
        self.xy_waveform_type:list = ['gauss', 'drag']
        self.z_waveform_type: list = ['square', 'gauss']

        for q_idx in range(q_num):
            self.__xylog[f"q{q_idx}"] = {"waveform":'drag',"duraOVERsigma":4,"drag_ratio":0,"halfPI_ratio":0.5} 
            self.__zlog[f"q{q_idx}"]  = {"waveform":'square',"duraOVERsigma":4} 
        
        for c_idx in range(c_num):
            self.__zlog[f"c{c_idx}"]  = {"waveform":'square',"duraOVERsigma":4} 

    def XY_waveform_controller(self, pulse_sche:Schedule, amp:float, duration:float, phase:float, q:str, rel_time:float, ref_op:Schedule, ref_pt:str="start"):
        """
        waveform switch\n
        sigma_factor determines the relations between pulse sigma and pulse duration is `sigma = duration/sigma_factor`.
        ### * kwargs: *
        1. if use DRAG, there should ne a argument named 'drag_amp'. It specifies the derivative part amplitude in DRAG as default is the same as amp.
        """
        match str(self.__xylog[q]["waveform"]).lower():
            case 'drag':
                diff_amp = float(self.__xylog[q]["drag_ratio"]) * amp
                return pulse_sche.add(DRAGPulse(G_amp=amp, D_amp=diff_amp, duration= duration, phase=phase, port=q+":mw", clock=q+".01",sigma=duration/float(self.__xylog[q]["duraOVERsigma"])),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
            case 'gauss':
                return pulse_sche.add(GaussPulse(G_amp=amp, phase=phase,duration=duration, port=q+":mw", clock=q+".01",sigma=duration/float(self.__xylog[q]["duraOVERsigma"])),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
            case _:
                pass
    
    def Z_waveform_controller(self, pulse_sche:Schedule, amp:float, duration:float, q:str, rel_time:float, ref_op:Schedule, ref_pt:str="start", phase:float=0):
        match str(self.__zlog[q]["waveform"]).lower():
            case 'square':
                return pulse_sche.add(SquarePulse(duration= duration,amp=amp, port=q+":fl", clock="cl0.baseband"),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
            case 'gauss':
                return pulse_sche.add(GaussPulse(G_amp=amp, phase=phase,duration=duration, port=q+":fl", clock="cl0.baseband",sigma=duration/float(self.__zlog[q]["duraOVERsigma"])),rel_time=rel_time,ref_op=ref_op,ref_pt=ref_pt)
            case _:
                pass

    def set_waveform_for(self,target_q:str, waveform:str, mode:str='xy'):
        """ 
        set waveform for `target_q` for the mode gate. 
        #### Args:\n
        * waveform: str,\n
            * xy-gate use ['gauss', 'drag']\n
            * z-gate use ['square', 'gauss']\n
        * mode: str, what gate ? 'xy' or 'z'.
        """
        match mode.lower():
            case 'xy':
                if waveform.lower() in self.xy_waveform_type:
                    self.__xylog[target_q]['waveform'] = waveform
                else:
                    raise ValueError(f"XY-gate Waveform must be in {self.xy_waveform_type}, but '{waveform}' was given.")
            case 'z':
                if waveform.lower() in self.z_waveform_type:
                    self.__zlog[target_q]['waveform'] = waveform
                else:
                    raise ValueError(f"Z-gate Waveform must be in {self.z_waveform_type}, but '{waveform}' was given.")
    
    def set_duraOVERsigma_for(self,target_q:str,duraOVERsigma:float):
        self.__xylog[target_q]["duraOVERsigma"] = duraOVERsigma

    def set_dragRatio_for(self,target_q:str,drag_ratio:float):
        self.__xylog[target_q]["drag_ratio"] = drag_ratio

    def set_halfPIratio_for(self,target_q:str,halfPI_ratio:float):
        self.__xylog[target_q]["halfPI_ratio"] = halfPI_ratio

    def get_XYwaveform_log(self):
        return self.__xylog
    def get_Zwaveform_log(self):
        return self.__zlog
    
    def activate(self, old_log:dict, mode:str):

        match mode.lower():
            case 'xy':
                self.__xylog = old_log
            case 'z':
                self.__zlog = old_log
            


    def X_theta(self,sche,amp,Du,q,ref_pulse_sche,freeDu):

        if Du!=0:
            delay_c= -Du-freeDu
            return self.XY_waveform_controller(sche,amp,Du,0,q,delay_c,ref_pulse_sche,ref_pt='start')
        else: pass

    def Y_theta(self,sche,amp,Du,q,ref_pulse_sche,freeDu):
        if Du!=0:
            delay_c= -Du-freeDu
            return self.XY_waveform_controller(sche,amp,Du,90,q,delay_c,ref_pulse_sche,ref_pt='start')
        else: pass

    def X_pi_2_p(self,sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu):
        amp= pi_amp[q]*half_pi_ratio
        delay_c= -pi_Du-freeDu
        return self.XY_waveform_controller(sche,amp,pi_Du,0,q,delay_c,ref_pulse_sche,ref_pt='start')

    def Y_pi_2_p(self,sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu):
        amp= pi_amp[q]*half_pi_ratio
        delay_c= -pi_Du-freeDu
        return self.XY_waveform_controller(sche,amp,pi_Du,90,q,delay_c,ref_pulse_sche,ref_pt='start')

    def X_pi_p(self,sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
        amp= pi_amp[q]
        delay_c= -pi_Du-freeDu
        return self.XY_waveform_controller(sche,amp,pi_Du,0,q,delay_c,ref_pulse_sche,ref_pt=ref_point)

    def Y_pi_p(self,sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
        amp= pi_amp[q]
        delay_c= -pi_Du-freeDu
        return self.XY_waveform_controller(sche,amp,pi_Du,90,q,delay_c,ref_pulse_sche,ref_pt=ref_point)

    def X_12pi_p(self,sche,pi_amp,q,pi_Du:float,ref_pulse_sche,freeDu, ref_point:str="start"):
        amp= pi_amp[q]
        delay_c= -pi_Du-freeDu
        return self.XY_waveform_controller(sche,amp,pi_Du,0,q,delay_c,ref_pulse_sche,ref_pt=ref_point)

    def Z(self,sche,Z_amp,Du,q,ref_pulse_sche,freeDu,ref_position='start'):
        if Du!=0:
            delay_z= -Du-freeDu
            return self.Z_waveform_controller(sche,Z_amp,Du,q,delay_z,ref_pulse_sche,ref_pt=ref_position)

    def Zc(self,sche,Z_amp,Du,q,ref_pulse_sche,freeDu):
        if Du!=0:
            delay_z= -Du-freeDu
            return self.Z_waveform_controller(sche,Z_amp,Du,q,delay_z,ref_pulse_sche,ref_pt="start")
        
