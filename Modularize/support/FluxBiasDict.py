from numpy import ndarray, sin, sqrt, cos

# TODO: Test this class to store the bias info
# dict to save the filux bias information
def quadratic(x,a,b,c):
    return a*(x**2)+b*x+c

class FluxBiasDict():
    """
    This class helps to memorize the flux bias. 
    """
    def __init__(self,qb_number:int):
        self.__bias_dict = {}
        self.q_num = qb_number
        self.float_cata = ["SweetSpot","TuneAway","Period"]
        self.list_cata = ["cavFitParas","qubFitParas"] # For the fitting functions access
        self.dict_cata = {"qubFitData":["bias_data","XYF_data"]} # For XYF-Flux fitting data storing

        self.init_bias()

    def sin_for_cav(self,target_q:str,bias_ary:ndarray):
        """
        Return the ROF curve data fit by the Quantify's sin model about the target_q with the bias array. 
        """
        def Sin(x,amp,w,phs,offset):
            return float(amp)*sin(float(w)*x+float(phs))+float(offset)
        return Sin(bias_ary,*self.__bias_dict[target_q]["cavFitParas"])
    
    def fqEqn_for_qub(self,target_q:str,bias_ary:ndarray):
        """
        Return the XYF curve data fit by quadratic about target_q with the bias array.
        """
        def FqEqn(x,a,b,Ec,coefA,d):
            """
            a ~ period, b ~ offset, 
            """
            return sqrt(8*coefA*Ec*sqrt(cos(a*(x-b))**2+d**2*sin(a*(x-b))**2))-Ec
        return FqEqn(bias_ary,*self.__bias_dict[target_q]["qubFitParas"])


    def get_bias_dict(self)->dict:
        """
        Return the dictionary contain the bias info.
        """
        return self.__bias_dict
    
    def init_bias(self):
        """
        Initialize the dict when you create this object.
        """
        for i in range(self.q_num):
            self.__bias_dict[f"q{i}"] = {}
            for bias_position in self.float_cata:
                self.__bias_dict[f"q{i}"][bias_position] = 0.0
            
            for para_cata in self.list_cata:
                self.__bias_dict[f"q{i}"][para_cata] = []
            
            for dict_cata in self.dict_cata:
                self.__bias_dict[f"q{i}"][dict_cata] = {}
                for subcata in self.dict_cata[dict_cata]: 
                    self.__bias_dict[f"q{i}"][dict_cata][subcata] = []


    def save_sweetspotBias_for(self,target_q:str='q0',bias:float=0.0):
        """
            Set the sweet spot bias for the given qubit, target_q label starts from 0.\n
            Ex. target_q = 'q0'
        """
        self.__bias_dict[target_q]["SweetSpot"] = bias

    def save_tuneawayBias_for(self,mode:str,target_q:str='q0',bias:float=0.0):
        """
            Set the tune away bias for the given qubit, target_q label starts from 0.\n
            Ex. target_q = 'q0'\n
            mode = 'manual' | 'auto'. `'manual'` get the given bias in args, `'auto'` calculated by the stored offset and period. 
        """
        if mode == 'manual':
            self.__bias_dict[target_q]["TuneAway"] = bias
        elif mode == 'auto':
            offset = self.get_sweetBiasFor(target_q)
            T = self.get_PeriodFor(target_q)
            tuneaway = offset + T/2 if abs(offset + T/2) < abs(offset - T/2) else offset - T/2
            self.__bias_dict[target_q]["TuneAway"] = tuneaway
        else:
            raise KeyError("Wrong mode!")
        
    def save_period_for(self,target_q:str,period:float):
        """
        Save the period for the target_q with that given period.
        """
        self.__bias_dict[target_q]["Period"] = period

    def save_cavFittingParas_for(self,target_q:str,amp:float,f:float,phi:float,offset:float):
        """
        Save the fitting fuction parameters for flux dep cavity, includes amplitude, frequency, phase and offset for a sin in the form:\n
        `A*sin(fx+phi)+offset`
        """
        self.__bias_dict[target_q]["cavFitParas"] = [amp,f,phi,offset]

    def save_qubFittingParas_for(self,target_q:str,a:float,b:float,Ec:float,Ej_sum:float,d:float):
        """
        Save the fitting function parameters for z-gate 2tone qubit which are a, b, Ec, Ej_sum, d
        """
        self.__bias_dict[target_q]["qubFitParas"] = [a, b, Ec, Ej_sum, d]
    
    def activate_from_dict(self,old_bias_dict:dict):
        """
        Activate the dict which super from the old record.
        """
        for q_old in old_bias_dict:
            for catas in old_bias_dict[q_old]:
                self.__bias_dict[q_old][catas] = old_bias_dict[q_old][catas]

    def get_sweetBiasFor(self,target_q:str):
        """
        Return the sweetSpot bias for the target qubit.
        """
        return self.__bias_dict[target_q]["SweetSpot"]
    
    def get_tuneawayBiasFor(self,target_q:str):
        """
        Return the tuneAway bias for the target qubit.
        """
        return self.__bias_dict[target_q]["TuneAway"]
    
    def get_PeriodFor(self,target_q:str):
        """
        Return the period for the target_q.
        """
        return self.__bias_dict[target_q]["Period"]