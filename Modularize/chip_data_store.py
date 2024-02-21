import os
import json
import time
import support as sup

'''
Unsolved problem: 
1. QD_file name should be set before running this program
2. QB_name needs to be settle up 

Questions:
1. 看起來CavitySpec還差了ro_out_att和xy_out_att?
'''

# This file should be setup after finishing Experiment setup, before any experiment start.

class Chip_file():
    
    def __init__(self,chip_name:str = "new_chip", chip_type:str = "5Q", QD_path:str = "", ro_out_att:int = 20, xy_out_att:int = 10):
        '''
        chip_name: The name of the chip, needs to be the same as in fabricaation, using the date of big current exposure + # of chip, eg. "20240206#8"
        chip_type: The type of the chip, eg. "5Q", "5Q4C"
        QD_path: The file name of Quantum Device in QD_backup
        ro_out_att: attenuation of ro line
        xy_out_att: attenuation of xy line
        '''
        
        self.name = chip_name            # chip 名稱
        self.file = self.name+".json"    # chip 檔案名稱
        self.type = chip_type            # chip 類別
        self.file_path = os.getcwd()+"\chip_information"  # chip 檔案位置
        self.file_name = os.path.join(self.file_path, self.file)    # chip 檔案完整位置+名稱
        self.path_today = self.file_path+"\Timestamp"+'\\'+time.strftime('%Y%m%d',time.localtime(time.time()))             # chip 當日檔案位置
        self.file_today = os.path.join(self.path_today, time.strftime('%Y%m%d',time.localtime(time.time()))+'_'+self.file) # chip 當日檔案完整位置+名稱
        self.QD_path = QD_path
        self.ro_out_att = ro_out_att
        self.xy_out_att = xy_out_att
        self.__chip_dict = {}
        self.init_chip_file()

    def init_chip_file(self):
        """
        生成一份chip file，輸入type和chip name後，會檢查這一份chip有沒有存在，
        若有則直接讀取後更新每天的timestamp，
        若無則會根據chip type複製一份檔案並更新在chip_information和每天的timestamp
        
        type: "5Q", ...
        name: arbitrary name
        
        """
        if self.type == "5Q":
            blank_file = "blank_chip_information.json"
            
        else:
            raise ValueError("We don't have this chip type, are you live in parallel universe?")
        
        file_exist = False
        
        
        # check the total chip information
        for root, dirs, files in os.walk(self.file_path):
            for il in files:
                if il == self.file:
                    #print(f'{il} exist')
                    file_exist = True
                    
        if file_exist == False:
            with open(os.path.join(self.file_path, blank_file), "r") as blank, open(self.file_name, 'w') as new:
                new.write(blank.read())
                
            with open(self.file_name, 'r') as rd:
                self.__chip_dict = json.load(rd)
            
            self.__chip_dict["basic_information"]["chip_name"] = self.name
            self.__chip_dict["basic_information"]["chip_type"] = self.type
            self.__chip_dict["basic_information"]["chip_file"] = self.file_name
            self.__chip_dict["basic_information"]["ro_att"] = self.ro_out_att
            self.__chip_dict["basic_information"]["xy_att"] = self.xy_out_att
            self.__chip_dict["basic_information"]["create_time"] = time.strftime('%Y%m%d',time.localtime(time.time()))
            
            with open(self.file_name, 'w') as up:
                json.dump(self.__chip_dict, up, indent=4)
            
        else:
            with open(self.file_name, 'r') as qu:
                self.__chip_dict = json.load(qu)
                
        # check today's chip imformation

        if os.path.isdir(self.path_today): 
            with open(self.file_today, 'w') as up:
                json.dump(self.__chip_dict, up, indent=4)
            # if the demo_folder2 directory is  
            # not present then create it. 
        else:
            os.makedirs(self.path_today)
            with open(self.file_today, 'w') as up:
                json.dump(self.__chip_dict, up, indent=4) 
        
    def get_chip_dict(self) -> dict:
        return self.__chip_dict
    
    def update_Cavity_spec_bare(self, QB_name:dict = dict(q0=None, q1=None, q2=None, q3=None, q4=None,), result:dict = {}) -> None:
        
        '''
        QB_name: the names of qubits. eg. dict(q0, q1, q2, q3, q4)
        result: CS_result
        '''
        
        cluster, quantum_device, meas_ctrl, ic, FluxRecorder, Hcfg = sup.init_meas(QuantumDevice_path=self.QD_path,mode='l')
        for qb in QB_name:
            qubit = quantum_device.get_element(qb)
            self.__chip_dict["1Q_information"][qb]["oper"]["readout"]["dura_time"] = qubit.measure.pulse_duration()
            self.__chip_dict["1Q_information"][qb]["oper"]["readout"]["acq_delay"] = qubit.measure.acq_delay()
            self.__chip_dict["1Q_information"][qb]["oper"]["readout"]["integ_time"] = qubit.measure.integration_time()
            self.__chip_dict["1Q_information"][qb]["init"]["wait"]["time"] = qubit.reset.duration()
            self.__chip_dict["1Q_information"][qb]["char"]["bare"]["bare_freq"] = result[qb].quantities_of_interest["fr"].nominal_value
            self.__chip_dict["1Q_information"][qb]["char"]["bare"]["Qi"] = result[qb].quantities_of_interest["Qi"].nominal_value
            self.__chip_dict["1Q_information"][qb]["char"]["bare"]["Qc"] = result[qb].quantities_of_interest["Qc"].nominal_value
            
        with open(self.file_name, 'w') as up:
            json.dump(self.__chip_dict, up, indent=4)
        with open(self.file_today, 'w') as up:
            json.dump(self.__chip_dict, up, indent=4) 
        print("updated!")     
    
         