import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support import QDmanager
from qblox_drive_AS.support.UserFriend import *
from quantify_scheduler.helpers.collections import find_port_clock_path
from qblox_drive_AS.support.Pulse_schedule_library import set_LO_frequency
from qblox_drive_AS.support.QDmanager import find_path_by_port

def set_integration_time(QD_agent:QDmanager, inte_time_s:dict=None):
    if inte_time_s is not None:
        if len(list(inte_time_s.keys())) != 0:
            for q in inte_time_s:
                QD_agent.quantum_device.get_element(q).measure.integration_time(inte_time_s[q])
                QD_agent.quantum_device.get_element(q).measure.pulse_duration(inte_time_s[q])

def set_drag_coef(QD_agent:QDmanager,Coefs:dict={}):
    if Coefs is not None:
        if len(list(Coefs.keys())) != 0:
            for q in Coefs:
                QD_agent.Waveformer.set_dragRatio_for(q, Coefs[q])

def setGlobally_reset_time(QD_agent:QDmanager, reset_time_s:float=None):
    if reset_time_s is not None:
        for q in QD_agent.quantum_device.elements():
            QD_agent.quantum_device.get_element(q).reset.duration(reset_time_s)

def set_drivin_IF(QD_agent:QDmanager, driving_IF_Hz:dict=None):
    if driving_IF_Hz is not None:
        if len(list(driving_IF_Hz.keys())) != 0:
            for q in driving_IF_Hz:
                QD_agent.Notewriter.save_xyIF_for(q,driving_IF_Hz[q])

def setGlobally_driving_atte(QD_agent:QDmanager,xy_atte:int=None):
    """ Atte must be multiples of 2 """
    if xy_atte is not None:
        slightly_print(f"Setting all qubits xy-attenuation at {xy_atte} dB")
        for q in QD_agent.quantum_device.elements():
            QD_agent.Notewriter.save_DigiAtte_For(xy_atte,q,'xy')

def set_ROF(QD_agent:QDmanager, ROFs:dict={}):
    if ROFs is not None:
        if len(list(ROFs.keys())) != 0:
            for q in ROFs:
                QD_agent.quantum_device.get_element(q).clock_freqs.readout(ROFs[q])


def set_sweet_bias(QD_agent:QDmanager, offsets:dict={}):
    if offsets is not None:
        if len(list(offsets.keys())) != 0:
            for q in offsets:
                QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=q,bias=offsets[q])

def set_roLOfreq(QD_agent:QDmanager,LO_Hz:float,target_q:str='q0'):
    """ ## *Warning*: 
        Set the LO for those qubits who shares the same readout module with the `target_q`.
    """
    if LO_Hz is not None:
        slightly_print(f"Set {find_port_clock_path(QD_agent.quantum_device.hardware_config(),'q:res',f'{target_q}.ro')[1]} RO-LO = {round(LO_Hz*1e-9,2)} GHz")
        set_LO_frequency(QD_agent.quantum_device,target_q,'readout',LO_Hz)

def set_roAtte(QD_agent:QDmanager,ro_atte:int, target_q:str='q0'):
    """ ## *Warning*: 
        * Set the readout attenuation for those qubits who shares the same readout module with the `target_q`.\n
        ## Args:\n
        ro_atte (int): multiple of 2.
    """
    if ro_atte is not None:
        slightly_print(f"Set {find_port_clock_path(QD_agent.quantum_device.hardware_config(),'q:res',f'{target_q}.ro')[1]} RO-attenuation = {int(ro_atte)} dB")
        who_onTheSameBoat:dict = find_path_by_port(QD_agent.quantum_device.hardware_config(),"q:res")
        for q in who_onTheSameBoat:
            if who_onTheSameBoat[q][0] == who_onTheSameBoat[target_q][0]:
                print(f"set {q} ...")
                QD_agent.Notewriter.save_DigiAtte_For(ro_atte,q,'ro')

def set_sweet_g(QD_agent:QDmanager,g_dict:dict=None):
    """ save g is the `g_dict` for q_key, g unit in Hz"""
    if g_dict is not None:
        if len(list(g_dict.keys())) != 0:
            for q in g_dict:
                QD_agent.Notewriter.save_sweetG_for(g_dict[q],q)

def set_ROamp_by_coef(QD_agent:QDmanager, roAmp_coef_dict:dict=None):
    if roAmp_coef_dict is not None:
        if len(list(roAmp_coef_dict.keys())) != 0:
            for q in roAmp_coef_dict:
                QD_agent.quantum_device.get_element(q).measure.pulse_amp( QD_agent.quantum_device.get_element(q).measure.pulse_amp()*float(roAmp_coef_dict[q]))

def update_coupler_bias(QD_agent:QDmanager,cp_elements:dict):
    """
    Update the idle bias in Fluxmanager for couplers.\n
    --------------------------
    ### Args:\n
    cp_elements = {"c0":0.2}
    """
    if cp_elements is not None:
        if len(list(cp_elements.keys())) != 0:
            for cp in cp_elements:
                slightly_print(f"Set coupler {cp} at {cp_elements[cp]} V.")
                QD_agent.Fluxmanager.save_idleBias_for(cp, cp_elements[cp])

if __name__ == "__main__":

    QD_path = "qblox_drive_AS/QD_backup/20241128/DR1#11_SumInfo.pkl"
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()

    """ Set RO amp by a coef. """
    set_ROamp_by_coef(QD_agent, roAmp_coef_dict={}) # roAmp_coef_dict = {"q0":0.93, "q1":0.96, ...}, set None or {} to bypass 

    """ Set RO freq """
    set_ROF(QD_agent, ROFs={})                      # ROFs = {"q0":6.0554e9, .....}
    
    """ Set Integration time """ 
    set_integration_time(QD_agent, inte_time_s={}) # inte_time_s = {"q0":1e-6, "q1":0.75e-6, ...}, set None or {} to bypass 

    """ Set reset time (All qubits global) """
    setGlobally_reset_time(QD_agent, reset_time_s=None)      # reset_time_s = 250e-6, all the qubit in the quantum_device will share the same value

    """ Set driving attenuation (All qubits blobal) """     # xy_atte = 10 recommended, all qubits are shared with a same value. Must be multiples of 2.
    setGlobally_driving_atte(QD_agent, xy_atte=None)

    """ Set driving IF """
    set_drivin_IF(QD_agent, driving_IF_Hz={})   # driving_IF_Hz = {"q0":-150e6, "q1":-100e6, ...}, set None or {} to bypass  !!! Always be negative !!!

    """ Set RO-LO, RO-atte ( target_q QRM-RF modlue global) """
    set_roLOfreq(QD_agent, LO_Hz=None, target_q='q0') # LO is global in the same QRM-RF module, set None to bypass 
    set_roAtte(QD_agent, ro_atte=None, target_q='q0') # RO-attenuation is global in the same QRM-RF module, set None to bypass 

    """ Set coupler bias """
    update_coupler_bias(QD_agent,cp_elements={})  # cp_elements = {"c0":0.1, "c2":0.05, ...}, set None or {} to bypass 
    
    """ Set grq for sweet spot """
    set_sweet_g(QD_agent,g_dict={})    # g_dict = {"q0":45e6, "q1":90e6, ...}, set None or {} to bypass 
    
    """ Set sweet spot bias """
    set_sweet_bias(QD_agent, offsets={})          # offsets = {"q0": 0.08, ....}

    """ Set Darg coef """
    set_drag_coef(QD_agent, Coefs={})


    QD_agent.QD_keeper() 