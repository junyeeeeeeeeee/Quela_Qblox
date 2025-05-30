import os, rich
from qblox_drive_AS.support.QDmanager import QDmanager, BasicTransmonElement
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support.StatifyContainer import Statifier
from qblox_drive_AS.support import set_LO_frequency
from typing import Literal, Optional

class QD_modifier():
    def __init__(self,QD_path:str):
        self.to_modifiy_item = []
        self.QD_path = QD_path
        self.QD_agent = QDmanager(QD_path)
        self.QD_agent.QD_loader()
        self.export(target_q='all')
        
    
    def comment(self, message:str=None):
        if message is not None:
            if message != "":
                self.QD_agent.refresh_log(str(message))
                self.to_modifiy_item.append("Comments")

    def init_12_settings(self, switch:Optional[Literal["ON", "OFF"]]="OFF"):
        if switch == "ON":
            qubits = list(self.QD_agent.quantum_device.elements())
            HCfg:dict = self.QD_agent.quantum_device.hardware_config()

            for q in qubits:
                # 12 Driving attenuation
                HCfg['hardware_options']['output_att'][f"{q}:mw-{q}.12"] = 0
                # 12 LO settings (The same LO as 01-driving)
                HCfg['hardware_options']["modulation_frequencies"][f"{q}:mw-{q}.12"] = {"lo_freq":self.QD_agent.quantum_device.get_element(q).clock_freqs.f01()-self.QD_agent.Notewriter.get_xyIFFor(q)}

            self.to_modifiy_item.append("12_initializations")


    def active_reset_switch(self, switch:Optional[Literal["ON", "OFF"]]="OFF"):
        if switch in ["ON", "OFF"]:
            status = True if switch == "ON" else False
            if self.QD_agent.activeReset != status:
                print("Current ActiveReset status: ", self.QD_agent.activeReset)
                self.QD_agent.activeReset = False # status
                self.to_modifiy_item.append("ActiveReset ON" if status else "ActiveReset OFF")
    
    def reset_discriminator(self, switch:Optional[Literal["ON", "OFF"]]="OFF"):
        if switch == "ON":
            self.QD_agent.StateDiscriminator = Statifier()
            self.to_modifiy_item.append("Discriminators")

    def reset_rotation_angle(self, target_qs:list):
        if len(target_qs) > 0:
            for q in target_qs:
                self.QD_agent.rotate_angle[q] = [0]
            self.to_modifiy_item.append("Rotation_angle")
    
    def set_XYF(self,xyfs_Hz:dict={}):
        if xyfs_Hz is not None:
            if len(list(xyfs_Hz.keys())) != 0:
                for q in xyfs_Hz:
                    Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    Xmon.clock_freqs.f01(xyfs_Hz[q])
                self.to_modifiy_item.append("XYF")

    def set_RamseyT2detuing(self, detunes:dict={}):
        if detunes is not None:
            if len(list(detunes.keys())) != 0:
                for q in detunes:
                    self.QD_agent.Notewriter.save_artiT2Detune_for(q,detunes[q])
                self.to_modifiy_item.append("RamseyT2_detuning")
    
    def set_f12ArtifDetuing(self, detunes:dict={}):
        if detunes is not None:
            if len(list(detunes.keys())) != 0:
                for q in detunes:
                    self.QD_agent.Notewriter.save_12artiDetune_for(q,detunes[q])
                self.to_modifiy_item.append("f12_artificial_detuning")

    def set_integration_time(self, inte_time_s:dict=None):
        if inte_time_s is not None:
            if len(list(inte_time_s.keys())) != 0:
                for q in inte_time_s:
                    Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    Xmon.measure.integration_time(inte_time_s[q])
                    Xmon.measure.pulse_duration(inte_time_s[q]+Xmon.measure.acq_delay())

                self.to_modifiy_item.append("inte_time")

    def set_EF_duration(self, EF_driving_duras:dict={}):
        if EF_driving_duras is not None:
            if len(list(EF_driving_duras.keys())) != 0:
                for q in EF_driving_duras:
                    self.QD_agent.Notewriter.save_12duration_for(q, EF_driving_duras[q])
                self.to_modifiy_item.append("EF-driving-duration")

    def set_drag_coef(self, Coefs:dict={}):
        if Coefs is not None:
            if len(list(Coefs.keys())) != 0:
                for q in Coefs:
                    self.QD_agent.Waveformer.set_dragRatio_for(q, Coefs[q])
                self.to_modifiy_item.append("DRAG_coef")

    def setGlobally_reset_time(self, reset_time_s:float=None):
        
        if reset_time_s is not None:
            for q in self.QD_agent.quantum_device.elements():
                Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                Xmon.reset.duration(reset_time_s)
            self.to_modifiy_item.append("reset_time")

    def set_drivin_IF(self, driving_IF_Hz:dict=None):
        if driving_IF_Hz is not None:
            if len(list(driving_IF_Hz.keys())) != 0:
                for q in driving_IF_Hz:
                    self.QD_agent.Notewriter.save_xyIF_for(q,driving_IF_Hz[q])
                self.to_modifiy_item.append("Driving_IF")

    def setGlobally_driving_atte(self,xy_atte:int=None):
        """ Atte must be multiples of 2 """
        if xy_atte is not None:
            slightly_print(f"Setting all qubits xy-attenuation at {xy_atte} dB")
            for q in self.QD_agent.quantum_device.elements():
                self.QD_agent.Notewriter.save_DigiAtte_For(xy_atte,q,'xy')
            self.to_modifiy_item.append("Driving_Atte")

    def set_ROF(self, ROFs:dict={}):
        if ROFs is not None:
            if len(list(ROFs.keys())) != 0:
                for q in ROFs:
                    Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    Xmon.clock_freqs.readout(ROFs[q])
                self.to_modifiy_item.append("ROF")

    def set_sweet_bias(self, offsets:dict={}):
        if offsets is not None:
            if len(list(offsets.keys())) != 0:
                for q in offsets:
                    self.QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=q,bias=offsets[q])
                self.to_modifiy_item.append("Sweet_Offset")

    def set_roLOfreq(self,LO_Hz:float,target_qs:list=['q0']):
        """ ## *Warning*: 
            Set the LO for those qubits who shares the same readout module with the `target_q`.
        """
        if LO_Hz is not None:
            slightly_print(f"Set RO-LO = {round(LO_Hz*1e-9,2)} GHz")
            for q in target_qs:
                set_LO_frequency(self.QD_agent.quantum_device,q,'readout',LO_Hz)
            self.to_modifiy_item.append("RO_LO")

    def set_roAtte(self, ro_atte:int, target_qs:list=['q0']):
        """ ## *Warning*: 
            * Set the readout attenuation for those qubits who shares the same readout module with the `target_q`.\n
            ## Args:\n
            ro_atte (int): multiple of 2.
        """
        if ro_atte is not None:
            slightly_print(f"Set RO-attenuation = {int(ro_atte)} dB")
            for q in target_qs:
                print(f"set {q} ...")
                self.QD_agent.Notewriter.save_DigiAtte_For(ro_atte,q,'ro')
            self.to_modifiy_item.append("RO_Atte")

    def set_sweet_g(self, g_dict:dict=None):
        """ save g is the `g_dict` for q_key, g unit in Hz"""
        if g_dict is not None:
            if len(list(g_dict.keys())) != 0:
                for q in g_dict:
                    self.QD_agent.Notewriter.save_sweetG_for(g_dict[q],q)
                self.to_modifiy_item.append("Sweet_g")

    def set_ROamp_by_coef(self,  roAmp_coef_dict:dict=None):
        if roAmp_coef_dict is not None:
            if len(list(roAmp_coef_dict.keys())) != 0:
                for q in roAmp_coef_dict:
                    Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    Xmon.measure.pulse_amp( Xmon.measure.pulse_amp()*float(roAmp_coef_dict[q]))
                self.to_modifiy_item.append("RO_amp")

    def set_XY_amp(self, pi_amps:dict={}):
        if pi_amps is not None:
            if len(list(pi_amps.keys())) != 0:
                for q in pi_amps:
                    Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    Xmon.rxy.amp180(pi_amps[q]) #rxy.duration
                self.to_modifiy_item.append("PI-amp")
    
    def set_XY_duration(self, pi_duras:dict={}):
        if pi_duras is not None:
            if len(list(pi_duras.keys())) != 0:
                for q in pi_duras:
                    Xmon:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    Xmon.rxy.duration(pi_duras[q]) #rxy.duration
                self.to_modifiy_item.append("PI-duration")

    def update_coupler_bias(self,cp_elements:dict):
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
                    self.QD_agent.Fluxmanager.save_idleBias_for(cp, cp_elements[cp])
                self.to_modifiy_item.append("coupler_bias")

    def save_modifications(self):
        if len(self.to_modifiy_item) != 0:
            highlight_print(f"These items will be modified: {self.to_modifiy_item}")
            permission = mark_input("Is it correct ? [y/n] :")
            
            if permission.lower() in ["y", "yes"]:
                self.export(target_q='all')   
                self.QD_agent.QD_keeper() 

            else:
                slightly_print("Your modifications got denied !")
        else:
            slightly_print("Nothing changed ~ ")
    
    def show_hcfg(self, switch:Optional[Literal["ON", "OFF"]]="OFF"):
        if switch == "ON":
            Hcfg = self.QD_agent.quantum_device.hardware_config()
            rich.print(Hcfg)

    def export(self, target_q:list=[]):
        """ Give a list contains qubit names you want to read info like: ['q0','q1'] or you can set target_q = 'all' to print all qubits."""
        if isinstance(target_q,str):
            if target_q.lower() == 'all':
                qs = self.QD_agent.quantum_device.elements()
        else:
            qs = target_q

        if len(qs) != 0:
            with open(os.path.join("qblox_drive_AS","Configs","QD_info.toml"), "w") as file:
                file.write(f"QD_file: {self.QD_path}\n")
                file.write(f"Chip name: {self.QD_agent.chip_name}\n")
                file.write(f"Comments: {self.QD_agent.Log if self.QD_agent.Log != '' else '---'}\n\n")
                file.write(f"RO-atte = {self.QD_agent.Notewriter.get_DigiAtteFor(qs[0],'ro')} dB\n")
                file.write(f"XY-atte = {self.QD_agent.Notewriter.get_DigiAtteFor(qs[0],'xy')} dB\n")
                if not self.QD_agent.activeReset:
                    file.write(f"Active Reset OFF, Reset time = {round(self.QD_agent.quantum_device.get_element(qs[0]).reset.duration()*1e6)} us\n\n")
                else:
                    file.write(f"Active Reset ON\n\n")
                for q in qs:
                    file.write(f'[{q}]\n')  
                    qubit:BasicTransmonElement = self.QD_agent.quantum_device.get_element(q)
                    file.write(f"    bare    = {self.QD_agent.Notewriter.get_bareFreqFor(q)*1e-9} GHz\n")
                    file.write(f"    ROF     = {qubit.clock_freqs.readout()*1e-9} GHz\n")
                    file.write(f"    RO-amp  = {round(qubit.measure.pulse_amp(),3)} V\n")
                    file.write(f"    ROT     = {round(qubit.measure.integration_time()*1e6,2)} us\n")
                    file.write(f"    f01     = {qubit.clock_freqs.f01()*1e-9} GHz\n")
                    file.write(f"    GE-amp  = {qubit.rxy.amp180()} V\n")
                    file.write(f"    GE-dura = {round(qubit.rxy.duration()*1e9,0)} ns\n")
                    file.write(f"    EF-dura = {round(self.QD_agent.Notewriter.get_12durationFor(q)*1e9,0)} ns\n")
                    file.write(f"    Motzoi  = {qubit.rxy.motzoi()}\n")
                    file.write(f"    x       = {(qubit.clock_freqs.readout()-self.QD_agent.Notewriter.get_bareFreqFor(q))*1e-6} MHz\n")
                    file.write(f"    g       = {self.QD_agent.Notewriter.get_sweetGFor(q)*1e-6} MHz\n")
                    file.write(f"    Ec      = {self.QD_agent.Notewriter.get_EcFor(q)*1e-6} MHz\n")
                    file.write("\n")

if __name__ == "__main__":




    QD_path = "qblox_drive_AS/QD_backup/20250418/DR1#11_SumInfo.pkl"


    QMaster = QD_modifier(QD_path)

    ### Readout
    """ Active Reset switch """
    QMaster.active_reset_switch(switch="OFF")   # On -> turn on the active reset. Otherwise, OFF for turning it off. 

    """ Reset Rotation angle """
    QMaster.reset_rotation_angle(target_qs=["q1","q3"])    # target_qs = ['q0', 'q1', ...]

    """ Set RO amp by a coef. """
    QMaster.set_ROamp_by_coef(roAmp_coef_dict={}) # roAmp_coef_dict = {"q0":0.93, "q1":0.96, ...}, set None or {} to bypass 

    """ Set RO freq """
    QMaster.set_ROF(ROFs={})                      # ROFs = {"q0":6.0554e9, .....}

    """ Set RO-LO, RO-atte ( target_q QRM-RF modlue global) """
    QMaster.set_roLOfreq(LO_Hz=None, target_qs=["q0","q1","q2","q3"]) # LO is global in the same QRM-RF module, set None to bypass 
    QMaster.set_roAtte(ro_atte=None, target_qs=["q0","q1","q2","q3"]) # RO-attenuation is global in the same QRM-RF module, set None to bypass 
    
    """ Set Integration time """ 
    QMaster.set_integration_time(inte_time_s={}) # inte_time_s = {"q0":1e-6, "q1":0.75e-6, ...}, set None or {} to bypass 

    """ Set reset time (All qubits global) """
    QMaster.setGlobally_reset_time(reset_time_s=None)      # reset_time_s = 250e-6, all the qubit in the quantum_device will share the same value

    ### Driving 
    """ Set XY Frequency """
    QMaster.set_XYF(xyfs_Hz = {})                         # xyfs_Hz = {"q0":4e9, "q1":4.5e9}, ** unit: Hz

    """ Set pi pulse amplitude """
    QMaster.set_XY_amp(pi_amps = {})                      # pi_amps = {"q0":0.2, "q1":0.15, ...}

    """ Set pi pulse duration """
    QMaster.set_XY_duration(pi_duras = {})                # pi_duras = {"q0":40e-9, "q1":48e-9, ...}

    """ Set driving attenuation (All qubits blobal) """     # xy_atte = 10 recommended, all qubits are shared with a same value. Must be multiples of 2.
    QMaster.setGlobally_driving_atte(xy_atte=None)

    """ Set driving IF """
    QMaster.set_drivin_IF(driving_IF_Hz={})   # driving_IF_Hz = {"q0":-150e6, "q1":-100e6, ...}, set None or {} to bypass  !!! Always be negative !!!

    """ Set Darg coef """
    QMaster.set_drag_coef(Coefs={})    # Coefs = {"q0": -0.5, ....}

    """ Set T2 use detuing, unit: Hz """
    QMaster.set_RamseyT2detuing(detunes={})   # detunes = {"q0": -0.5e6, ....}

    ### z and others
    """ Set coupler bias """
    QMaster.update_coupler_bias(cp_elements={})  # cp_elements = {"c0":0.1, "c2":0.05, ...}, set None or {} to bypass 
    
    """ Set grq for sweet spot """
    QMaster.set_sweet_g(g_dict={})    # g_dict = {"q0":45e6, "q1":90e6, ...}, set None or {} to bypass 
    
    """ Set sweet spot bias """
    QMaster.set_sweet_bias(offsets={})          # offsets = {"q0": 0.08, ....}


    """ Some comments for this QD file """
    QMaster.comment(message="")              # message is a str, like "Hello QD comments !", etc.


    QMaster.show_hcfg(switch="OFF")   # You can see the HCFG in the terminal if you switch on

    """ Reset the whole Statifier """
    QMaster.reset_discriminator(switch="OFF")  # trun on the switch if you wanna initialize it


    ### 12 state settings
    """ 1-2 state hardware setting initialize """
    QMaster.init_12_settings(switch="OFF")     # hardware config settings for 12 transition clock   

    """ 1-2 driving pulse duration """
    QMaster.set_EF_duration(EF_driving_duras= {})   # EF_driving_duras = {"q0":40e-9, ...}. Keep it empty if you don't wanna change it.

    """ 1-2 state frequency artificial detuning """
    QMaster.set_f12ArtifDetuing(detunes = {})    # artificial detuning for f12, unit in Hz. Keep it empty if you don't wanna change it.
    
    QMaster.save_modifications()