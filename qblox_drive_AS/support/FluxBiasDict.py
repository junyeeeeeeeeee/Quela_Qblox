from numpy import array, ndarray, sin, sqrt, cos, pi, real
from qblox_drive_AS.support.UserFriend import *

# TODO: Test this class to store the bias info
# dict to save the filux bias information
def quadratic(x,a,b,c):
    return a*(x**2)+b*x+c

class FluxBiasDict():
    """
    This class helps to memorize the flux bias. 
    """
    def __init__(self,qb_number:int, cp_number:int):
        self.__bias_dict = {}
        self.q_num = qb_number
        self.c_num = cp_number
        self.float_cata = ["SweetSpot","TuneAway","Period"]
        self.list_cata = ["cavFitParas","qubFitParas"] # For the fitting functions access
        self.dict_cata = {"qubFitData":["bias_data","XYF_data"]} # For XYF-Flux fitting data storing
        self.bool_cata = ["offSweetSpot"]

        self.cp_float_cata = ["idle", "operation", "max"]

        self.init_bias()

    def press_offsweetspot_button(self,target_q:str,offSweetSpot:bool):
        self.__bias_dict[target_q]["offSweetSpot"] = offSweetSpot
    def get_offsweetspot_button(self,target_q:str)->bool:
        return self.__bias_dict[target_q]["offSweetSpot"]


    def sin_for_cav(self,target_q:str,bias_ary:ndarray):
        """
        Return the ROF curve data fit by the Quantify's sin model about the target_q with the bias array. 
        """
        
        def Sin(x,amp,w,phs,offset):
            return sin(float(w)*x+float(phs))*float(amp)+float(offset)
        
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
            
            for bool_cate in self.bool_cata:
                self.__bias_dict[f"q{i}"][bool_cate] = False #Fasle
        
        for i in range(self.c_num):
            self.__bias_dict[f"c{i}"] = {}
            for bias_position in self.cp_float_cata:
                self.__bias_dict[f"c{i}"][bias_position] = 0.0

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
            if q_old[0] == "q": # qubit 
                for catas in old_bias_dict[q_old]:
                    try:
                        self.__bias_dict[q_old][catas] = old_bias_dict[q_old][catas]
                    except:
                        if catas in self.float_cata:
                            self.__bias_dict[q_old][catas] = 0.0
                        elif catas in self.list_cata:
                            self.__bias_dict[q_old][catas] = []
                        elif catas in self.dict_cata:
                            self.__bias_dict[q_old][catas] = {}
                            for subcata in self.dict_cata[q_old]: 
                                self.__bias_dict[q_old][catas][subcata] = []
                        else:
                            self.__bias_dict[q_old][catas] = 0 #Fasle
            else: # coupler 
                for catas in old_bias_dict[q_old]:
                    try:
                        self.__bias_dict[q_old][catas] = old_bias_dict[q_old][catas]
                    except:
                        self.__bias_dict[q_old][catas] = 0.0

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
    
    def check_offset_and_correctFor(self,target_q:str,new_offset:float):
        """
        If we had once fit the fq vs bias even the offset changed, we can use the offset differences to correct the fitting para 'b'. 
        ## *** This function should be executed before `save_sweetspotBias_for()`
        """
        f = 2*pi/self.get_PeriodFor(target_q)
        old_offset = self.get_sweetBiasFor(target_q)
        offset_diff = new_offset - old_offset
        if len(self.__bias_dict[target_q]["qubFitParas"]) != 0:
            self.__bias_dict[target_q]["qubFitParas"][1] += offset_diff/f  # correct the offset in fitting paras with the offset difference
            print("You don't have to fit the transition freq again! Correct it by the offset differences!")


    def get_biasWithFq_from(self,target_q:str,target_fq_Hz:float, flux_guard:float=0.4):
        """
        After we fit the tarnsition freq vs bias, we can get the bias according to the given `target_fq_Hz` for the `target_q`.\n
        ### The given `target_fq_Hz` should unit in Hz.\n
        Return the bias unit in V.
        """
        import sympy as sp
        from qblox_drive_AS.support import find_nearest
        z = sp.Symbol('z',real=True)
        if self.__bias_dict[target_q]["qubFitParas"] == []:
            raise ValueError("You have NOT fit the transition frequency with bias!")

        a,b,Ec,Ej_sum,d = self.__bias_dict[target_q]["qubFitParas"][0], self.__bias_dict[target_q]["qubFitParas"][1], self.__bias_dict[target_q]["qubFitParas"][2], self.__bias_dict[target_q]["qubFitParas"][3], self.__bias_dict[target_q]["qubFitParas"][4]
        to_solve = sp.sqrt(8*Ej_sum*Ec*sp.sqrt(sp.cos(a*(z-b))**2+d**2*sp.sin(a*(z-b))**2))-Ec - target_fq_Hz*1e-9

        candidators = array(sp.solvers.solve(to_solve, z))
        
        if candidators.shape[0] == 0:
            answer = 'n'
            fq_max = self.fqEqn_for_qub(target_q,array([self.get_sweetBiasFor(target_q)]))[0]
            warning_print(f"Can NOT find a bias makes the fq @ {target_fq_Hz*1e-9} GHz !")
            warning_print(f"The max fq about this qubit = {fq_max} GHz")
        else:
            answer = find_nearest(candidators, 0) if abs(find_nearest(candidators, 0)) < flux_guard else flux_guard

        return answer

    def get_proper_zbiasFor(self,target_q:str)->float:
        """
        According to the `offsweetspot_button`, it returns the z-bias for the given target_q.\n
        The only way to press this button is in `Modularize/support/meas_switch.py`
        """
        if self.get_offsweetspot_button(target_q):
            eyeson_print(f"Now in tune away bias = {self.get_tuneawayBiasFor(target_q)} V")
            return self.get_tuneawayBiasFor(target_q)
        else:
            eyeson_print(f"Now in sweet spot bias = {self.get_sweetBiasFor(target_q)} V ")
            return self.get_sweetBiasFor(target_q)
        

    def save_idleBias_for(self,target_q:str,idle_bias:float):
        """ Save the idle bias for a qubit """
        self.__bias_dict[target_q]["idle"] = idle_bias

    def get_idleBiasFor(self, target_q:str):
        """ Get the idle bias for a qubit """
        return self.__bias_dict[target_q]["idle"]

    def save_operationBias_for(self, target_q:str, operation_bias:float):
        """ Save the operation bias for a qubit """
        self.__bias_dict[target_q]["operation"] = operation_bias

    def get_operationBiasFor(self, target_q:str):
        """ Get the operation bias for a qubit """
        return self.__bias_dict[target_q]["operation"]

    def save_CPmaxBias_for(self, target_q:str, max_freq_bias:float):
        """ Save the max freq bias for a qubit """
        self.__bias_dict[target_q]["max"] = max_freq_bias

    def get_CPmaxBiasFor(self, target_q:str):
        """ Get the max freq bias for a qubit """
        return self.__bias_dict[target_q]["max"]
    
    def build_Cctrl_instructions(self, cp_names:list, bias_where:str='idle')->dict:
        """
        #### Build up a dict for Cctrl control by the given coupler names in list and the bias you want.
        --------------
        ### Args:\n
        `cp_names` in list, like ["c0", "c2"].\n
        `bias_where` should be "idle", "operation" or "max".\n
        --------------
        ### Returns:\n
        Dict, {"c0":0.2, "c2":0.08} for example.
        """
        if bias_where.lower() in ['idle','i','id']:
            cata = 'idle'
        elif  bias_where.lower() in ['operation','o','op']:
            cata = 'operation'
        elif  bias_where.lower() in ['max','m']:
            cata = 'max'
        else:
            raise ValueError(f"Unsupported bias was given = {bias_where}!")

        dict_for_Cctrl = {}
        for cp_name in cp_names:
            dict_for_Cctrl[cp_name] = self.__bias_dict[cp_name][cata]
        
        return dict_for_Cctrl




if __name__ == "__main__":
    pass