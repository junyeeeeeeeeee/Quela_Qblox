class Notebook():
    def __init__(self,q_number:str):
        self.__dict = {}
        self.q_num = q_number
        self.cata = ["bareF","T1","T2","CoefInG","sweetG"] # all in float

        self.init_dict()

    def init_dict(self):
        for i in range(self.q_num):
            self.__dict[f"q{i}"] = {}
            for float_cata in self.cata:
                self.__dict[f"q{i}"][float_cata] = 0.0
    ## About get
    def get_notebook(self,target_q:str=''):
        if target_q != '':
            return self.__dict[target_q]
        else:
            return self.__dict
    ## About write and save   
    # For bare cavity freq
    def save_bareFreq_for(self,bare_freq:float,target_q:str="q1"):
        self.__dict[target_q]["bareF"] = bare_freq
    def get_bareFreqFor(self,target_q:str="q1"):
        return self.__dict[target_q]["bareF"]
    # For T1
    def save_T1_for(self,T1:float,target_q:str="q1"):
        self.__dict[target_q]["T1"] = T1
    def get_T1For(self,target_q:str="q1"):
        return self.__dict[target_q]["T1"]
    # For T2
    def save_T2_for(self,T2:float,target_q:str="q1"):
        self.__dict[target_q]["T2"] = T2
    def get_T2For(self,target_q:str="q1"):
        return self.__dict[target_q]["T2"]
    # For coef A in g formula: g(MHz) = coefA*sqrt(fb*fq)/1000, fb and fq in GHz
    def save_CoefInG_for(self,A:float,target_q:str='q0'):
        """
        A is right the coef in the g formula: g = A*sqrt(fb*fq)/1000, \nwhich g will be in MHz, fb and fq are in GHz.      
        """
        self.__dict[target_q]["CoefInG"] = A
    def get_CoefInGFor(self,target_q:str):
        return self.__dict[target_q]["CoefInG"]
    # note the g at sweet spot in Hz
    def save_sweetG_for(self,g_Hz:float,target_q:str='q0'):
        self.__dict[target_q]["sweetG"] = g_Hz
    def get_sweetGFor(self,target_q:str):
        return self.__dict[target_q]["sweetG"]
    # Activate notebook
    def activate_from_dict(self,old_notebook:dict):
        for qu in old_notebook:
            for cata in self.cata:
                try:
                    self.__dict[qu][cata] = old_notebook[qu][cata]
                except:
                    if cata in self.cata:
                        print(f"Old notebook didn't exist cata named '{cata}', initialize it in new notebook.")
                        self.__dict[qu][cata] = 0.0
                    else:
                        raise KeyError(f"Given cata='{cata}' can not be recognized!")

