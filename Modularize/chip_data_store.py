import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import json
import time
import Modularize.support as sup
import tkinter as tk

# This file should be setup after finishing Experiment setup, before any experiment start.

'''
How to use:
1.  After finishing Experiment Setup or After any measurement,
    if you are sure that the operation parameters are right, then run "00_update_chip_info.py".
    Store the most recent .pkl file into the chip.
2.  (Built-in) After a measurement, the result will be stored into chip file.
'''

class Chip_file():
    
    def __init__(self,chip_name:str = "new_chip"):
        '''
        chip_name:  The name of the chip, needs to be the same as in fabricaation,
                    using the date of big current exposure + # of chip, eg. "20240206#8"
        '''
        
        self.name = chip_name           # chip 名稱
        self.file = self.name+".json"   # chip 檔案名稱
        self.file_path = os.getcwd()+"\chip_information"  # chip 檔案位置
        self.file_name = os.path.join(self.file_path, self.file)    # chip 檔案完整位置+名稱
        self.path_today = self.file_path+"\Timestamp"+'\\'+time.strftime('%Y%m%d',time.localtime(time.time()))             # chip 當日檔案位置
        self.file_today = os.path.join(self.path_today, time.strftime('%Y%m%d',time.localtime(time.time()))+'_'+self.file) # chip 當日檔案完整位置+名稱
        self.__chip_dict = {}
        self.init_chip_file()
        
    def init_chip_file(self) -> None:
        '''
        檢查此名稱的chip有沒有存在，
        若有，則直接讀取後更新每天的timestamp，
        若無，則會根據chip type複製一份檔案並更新在chip_information和每天的timestamp。
        
        chip_type: The type of the chip, eg. "5Q", "5Q4C"
        ro_out_att: attenuation of ro line
        xy_out_att: attenuation of xy line
        
        '''
        
        # check the total chip information
        file_exist = False
        for root, dirs, files in os.walk(self.file_path):
            for il in files:
                if il == self.file:
                    print(f'Chip "{self.name}" exist, loading...')
                    file_exist = True
                    
        if file_exist == False:
            print("File doesn't exist, creating new file...")
            # 手動輸入
            '''
            setfile = tk.Tk()
            setfile.title('Basic information window')
            chip_type_la = tk.Label(setfile, text = "chip type")
            chip_type_la.pack()
            chip_type_en = tk.Entry(setfile, width=25)
            chip_type_en.pack()
            ro_out_att_la = tk.Label(setfile, text = "RO output attenuation")
            ro_out_att_la.pack()
            ro_out_att_en = tk.Entry(setfile, width=25)
            ro_out_att_en.pack()
            xy_out_att_la = tk.Label(setfile, text = "RO output attenuation")
            xy_out_att_la.pack()
            xy_out_att_en = tk.Entry(setfile, width=25)
            xy_out_att_en.pack()
            check_button = tk.Button(setfile, text="Enter")
            check_button.pack()
            setfile.mainloop()
            
            chip_type = chip_type_en.get()
            ro_out_att = ro_out_att_en.get()
            xy_out_att = xy_out_att_en.get()
            '''
            self.create_chip_file(chip_type, ro_out_att, xy_out_att)
            
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

    def create_chip_file(self, chip_type:str = "5Q", ro_out_att:int = 20, xy_out_att:int = 10) -> None:
        """
        僅使用於init_chip_file，當該名稱的chip不存在時將會執行此函式
        """
        
        if chip_type == "5Q":
            blank_file = "blank_chip_information.json"
        else:
            raise ValueError("We don't have this chip type, are you live in parallel universe?")
        
        with open(os.path.join(self.file_path, blank_file), "r") as blank, open(self.file_name, 'w') as new:
            new.write(blank.read())
                
        with open(self.file_name, 'r') as rd:
            self.__chip_dict = json.load(rd)
        
        self.__chip_dict["basic_information"]["chip_name"] = self.name
        self.__chip_dict["basic_information"]["chip_type"] = chip_type
        self.__chip_dict["basic_information"]["chip_file"] = self.file_name
        self.__chip_dict["basic_information"]["ro_att"] = ro_out_att
        self.__chip_dict["basic_information"]["xy_att"] = xy_out_att
        self.__chip_dict["basic_information"]["create_time"] = time.strftime('%Y%m%d',time.localtime(time.time()))
        
        with open(self.file_name, 'w') as up:
            json.dump(self.__chip_dict, up, indent=4)
        
    def get_chip_dict(self) -> dict:
        return self.__chip_dict
    
    def update_QD(self, QB_name:dict = {'q0':None, 'q1':None, 'q2':None, 'q3':None, 'q4':None}, QD_path = '') -> None:
        
        '''
        QB_name: the names of qubits. eg. {'q0':None, 'q1':None, 'q2':None, 'q3':None, 'q4':None}
        QD_path: 
        '''
        
        if QD_path == '':
            print("No Quantum_device is being stored")
            pass
        else:
            Qmanager = sup.QDmanager(QD_path)
            Qmanager.QD_loader()
            for qb in QB_name:
                qubit = Qmanager.quantum_device.get_element(qb)
                self.__chip_dict["1Q_information"][qb]["oper"]["readout"]["dura_time"] = qubit.measure.pulse_duration()
                self.__chip_dict["1Q_information"][qb]["oper"]["readout"]["acq_delay"] = qubit.measure.acq_delay()
                self.__chip_dict["1Q_information"][qb]["oper"]["readout"]["integ_time"] = qubit.measure.integration_time()
                self.__chip_dict["1Q_information"][qb]["init"]["wait"]["time"] = qubit.reset.duration()
            with open(self.file_name, 'w') as up:
                json.dump(self.__chip_dict, up, indent=4)
            with open(self.file_today, 'w') as up:
                json.dump(self.__chip_dict, up, indent=4) 
            print("Quantum_device updated!")     
        
    def update_Cavity_spec_bare(self, QB_name:dict = {'q0':None, 'q1':None, 'q2':None, 'q3':None, 'q4':None}, result:dict = {}) -> None:
        
        for qb in QB_name:
            self.__chip_dict["1Q_information"][qb]["char"]["bare"]["bare_freq"] = result[qb].quantities_of_interest["fr"].nominal_value
            self.__chip_dict["1Q_information"][qb]["char"]["bare"]["Qi"] = result[qb].quantities_of_interest["Qi"].nominal_value
            self.__chip_dict["1Q_information"][qb]["char"]["bare"]["Qc"] = result[qb].quantities_of_interest["Qc"].nominal_value
        
        with open(self.file_name, 'w') as up:
            json.dump(self.__chip_dict, up, indent=4)
        with open(self.file_today, 'w') as up:
            json.dump(self.__chip_dict, up, indent=4) 
        print("updated!")     
         