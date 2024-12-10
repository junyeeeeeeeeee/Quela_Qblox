class Notebook():
    def __init__(self,q_number:str):
        self.__InfoDict = {}
        self.q_num = q_number
        self.cata = ["bareF","T1","T2","T2*","CoefInG","sweetG","digitalAtte","realAtte","meas_options","2tone_piamp","xy_if", "arti_detu_ramsey"] # all in float

        self.init_dict()

    def init_dict(self):
        for i in range(self.q_num):
            self.__InfoDict[f"q{i}"] = {}
            for cata in self.cata:
                if cata in ["bareF","T1","T2","T2*","CoefInG","sweetG","2tone_piamp","arti_detu_ramsey"]: # float type
                    self.__InfoDict[f"q{i}"][cata] = 0.0
                elif cata in ["meas_options"]: # list type
                    self.__InfoDict[f"q{i}"][cata] = []  
                elif cata == 'xy_if':
                    self.__InfoDict[f"q{i}"][cata] = -150e6
                else: # dict type
                    self.__InfoDict[f"q{i}"][cata] = {}
                    self.__InfoDict[f"q{i}"][cata]['xy'] = 0
                    self.__InfoDict[f"q{i}"][cata]['ro'] = 16 if cata == 'digitalAtte' else 0


    ## About get
    def get_notebook(self,target_q:str=''):
        if target_q != '':
            return self.__InfoDict[target_q]
        else:
            return self.__InfoDict
        
    def save_artiT2Detune_for(self, target_q:str, arti_detu:float):
        """ Save artificial detuning only used on Ramsey"""
        self.__InfoDict[target_q]["arti_detu_ramsey"] = arti_detu
    def get_artiT2DetuneFor(self, target_q:str):
        """ get artificial detuning only used on Ramsey"""
        return self.__InfoDict[target_q]["arti_detu_ramsey"]
        
    # For driving IF freq
    def save_xyIF_for(self,target_q:str,driving_if:float):
        """ Always save the minus IF freq """
        self.__InfoDict[target_q]["xy_if"] = -abs(driving_if)
    def get_xyIFFor(self, target_q:str):
        """ Get the driving IF in Hz for `target_q`, and it's negative."""
        return self.__InfoDict[target_q]["xy_if"]

    ## About write and save 
    def save_2tone_piamp_for(self,target_q:str,pi_amp:float):
        self.__InfoDict[target_q]["2tone_piamp"] = pi_amp  
    def get_2tone_piampFor(self,target_q:str):
        return self.__InfoDict[target_q]["2tone_piamp"]

    # For bare cavity freq
    def save_bareFreq_for(self,bare_freq:float,target_q:str="q1"):
        self.__InfoDict[target_q]["bareF"] = bare_freq
    def get_bareFreqFor(self,target_q:str="q1"):
        return self.__InfoDict[target_q]["bareF"]
    
    # For T1
    def save_T1_for(self,T1:float,target_q:str="q1"):
        self.__InfoDict[target_q]["T1"] = T1
    def get_T1For(self,target_q:str="q1"):
        return self.__InfoDict[target_q]["T1"]
    # For T2
    def save_T2_for(self,T2:float,target_q:str="q1"):
        self.__InfoDict[target_q]["T2"] = T2
    def get_T2For(self,target_q:str="q1"):
        return self.__InfoDict[target_q]["T2"]
    # For echo T2
    def save_echoT2_for(self,T2e:float,target_q:str="q1"):
        self.__InfoDict[target_q]["T2*"] = T2e
    def get_echoT2For(self,target_q:str="q1"):
        return self.__InfoDict[target_q]["T2*"]

    # For coef A in g formula: g(MHz) = coefA*sqrt(fb*fq)/1000, fb and fq in GHz
    def save_CoefInG_for(self,A:float,target_q:str='q0'):
        """
        A is right the coef in the g formula: g = A*sqrt(fb*fq)/1000, \nwhich g will be in MHz, fb and fq are in GHz.      
        """
        self.__InfoDict[target_q]["CoefInG"] = A
    def get_CoefInGFor(self,target_q:str):
        return self.__InfoDict[target_q]["CoefInG"]
    
    # note the g at sweet spot in Hz
    def save_sweetG_for(self,g_Hz:float,target_q:str='q0'):
        self.__InfoDict[target_q]["sweetG"] = g_Hz
    def get_sweetGFor(self,target_q:str):
        """ Return the g_rq in Hz for target_q"""
        return self.__InfoDict[target_q]["sweetG"]
    
    # for attenuation
    def save_DigiAtte_For(self,atte_dB:int,target_q:str='q0',mode:str='xy'):
        """
        Save digital attenuation for target q according to the mode = 'xy' or 'ro'. 
        """
        if mode.lower() in ['xy', 'ro']:
            self.__InfoDict[target_q]["digitalAtte"][mode] = atte_dB
        else:
            raise KeyError(f"Wrong given mode={mode}. Expecting 'xy' or 'ro'.")
    def save_RealAtte_For(self,atte_dB:int,target_q:str='q0',mode:str='xy'):
        """
        Save real attenuation for target q according to the mode = 'xy' or 'ro'. 
        """
        if mode.lower() in ['xy', 'ro']:
            self.__InfoDict[target_q]["realAtte"][mode] = atte_dB
        else:
            raise KeyError(f"Wrong given mode={mode}. Expecting 'xy' or 'ro'.")
    def modify_DigiAtte_For(self,add_atte_dB:int,target_q:str, mode:str):
        """
        For a target_q, directly add the given add_atte_dB to its saved DigiAtte.\n
        add_atte_dB should be a multiple of 2.\n
        mode : 'ro' for readout atte, 'xy' for control atte.
        """
        if mode.lower() in ['xy', 'ro']:
            self.__InfoDict[target_q]["digitalAtte"][mode] = self.__InfoDict[target_q]["digitalAtte"][mode] + add_atte_dB
        else:
            raise KeyError(f"Wrong given mode={mode}. Expecting 'xy' or 'ro'.")


    def get_DigiAtteFor(self,target_q:str,mode:str):
        """
        Return digital attenuation for target q according to the mode = 'xy' or 'ro'. 
        """
        if mode.lower() in ['xy', 'ro']:
            return self.__InfoDict[target_q]["digitalAtte"][mode]
        else:
            raise KeyError(f"Wrong given mode={mode}. Expecting 'xy' or 'ro'.")
    def get_RealAtteFor(self,target_q:str,mode:str):
        """
        Return digital attenuation for target q according to the mode = 'xy' or 'ro'. 
        """
        if mode.lower() in ['xy', 'ro']:
            return self.__InfoDict[target_q]["realAtte"][mode]
        else:
            raise KeyError(f"Wrong given mode={mode}. Expecting 'xy' or 'ro'.")

    def create_meas_options(self,target_q:str):
        meas_elements = {"rof":0,"rop":0,"f01":0,"bias":0,"refIQ":[],"pi_amp":0,"2tone_pi_amp":0,"pi_dura":0,"ro_atte":0,"XY_waveform":{},"Z_waveform":{}}
        self.__InfoDict[target_q]["meas_options"].append(meas_elements)
    def get_all_meas_options(self,target_q:str)->list:
        return self.__InfoDict[target_q]["meas_options"]
    def write_meas_options(self,options_inclu_targetQ:dict,ToModified_index:int=-1):
        """
        1) options_inclu_targetQ={"q0":{"f01","rof","rop","pi_amp","2tone_pi_amp","pi_dura","refIQ","bias","ro_atte"}}.\n
        2) ToModified_index is the index of the self.__InfoDict[target_q]["meas_options"] list, default is do append it. !Don't use negative index!
        """
        if ToModified_index == -1:
            for target_q in options_inclu_targetQ:
                self.__InfoDict[target_q]["meas_options"].append(options_inclu_targetQ[target_q])
        else:
            for target_q in options_inclu_targetQ:
                if len(self.__InfoDict[target_q]["meas_options"]) > ToModified_index:
                    self.__InfoDict[target_q]["meas_options"][ToModified_index] = options_inclu_targetQ[target_q]
                else:
                    print("warning the given index > len(write_meas_options) now we have !")

    # Activate notebook
    def activate_from_dict(self,old_notebook:dict):
        for qu in old_notebook:
            for cata in self.cata:
                try:
                    self.__InfoDict[qu][cata] = old_notebook[qu][cata]
                except:
                    if cata in self.cata:
                        print(f"Old notebook didn't exist cata named '{cata}', initialize it in new notebook.")
                        if cata in ["bareF","T1","T2","T2*","CoefInG","sweetG","2tone_piamp", "arti_detu_ramsey"]:
                            self.__InfoDict[qu][cata] = 0.0
                        elif cata in ["meas_options"]:
                            self.__InfoDict[qu][cata] = []
                        elif cata == 'xy_if':
                            self.__InfoDict[qu][cata] = -150e6
                        else:
                            self.__InfoDict[qu][cata] = {}
                            self.__InfoDict[qu][cata]['xy'] = 0
                            self.__InfoDict[qu][cata]['ro'] = 0
                    else:
                        raise KeyError(f"Given cata='{cata}' can not be recognized!")

