import os
import json
import time

def init_chip_file(type:str, name:str, file):
    """
    生成一份chip file，輸入type和chip name後，會檢查這一份chip有沒有存在，
    若有則直接讀取後更新每天的timestamp，
    若無則會根據chip type複製一份檔案並更新在chip_information和每天的timestamp
    
    type: "5Q", ...
    name: arbitrary name
    
    """
    if type == "5Q":
        blank_file = "blank_chip_information.json"
        
    else:
        raise ValueError("We don't have this chip type, are you live in parallel space?")
    
    chip_file_path = os.getcwd()+"\chip_information"
    chip_Timestamp_path = chip_file_path+"\Timestamp"
    chip_path_today = chip_Timestamp_path+'\\'+time.strftime('%Y%m%d',time.localtime(time.time()))
    chip_file_today = os.path.join(chip_path_today, time.strftime('%Y%m%d',time.localtime(time.time()))+'_'+file)    
    
    file_exist = False
    file_name = os.path.join(chip_file_path, file)
    
    # check the total chip information
    for root, dirs, files in os.walk(chip_file_path):
        for il in files:
            if il == file:
                #print(f'{il} exist')
                file_exist = True
                
    if file_exist == False:
        with open(os.path.join(chip_file_path, blank_file), "r") as blank, open(file_name, 'w') as new:
            new.write(blank.read())
            
        with open(file_name, 'r') as rd:
            chip_information = json.load(rd)
        
        chip_information["basic_information"]["chip_name"] = name
        chip_information["basic_information"]["chip_type"] = type
        chip_information["basic_information"]["chip_file"] = file_name
        chip_information["basic_information"]["create_time"] = time.strftime('%Y%m%d',time.localtime(time.time()))
            
        with open(file_name, 'w') as up:
            json.dump(chip_information, up, indent=4)
        
    else:
        
        with open(file_name, 'r') as qu:
            chip_information = json.load(qu)
            #print(chip_information["basic_information"]["qubit_name"])
            
    # check today's chip imformation

    if os.path.isdir(chip_path_today): 
        with open(chip_file_today, 'w') as up:
            json.dump(chip_information, up, indent=4)
        # if the demo_folder2 directory is  
        # not present then create it. 
    else:
        os.makedirs(chip_path_today)
        with open(chip_file_today, 'w') as up:
            json.dump(chip_information, up, indent=4) 
    
            
    
    return chip_information

def update_chip_file(file, chip_information):
    
    chip_file_path = os.getcwd()+"\chip_information"
    chip_Timestamp_path = chip_file_path+"\Timestamp"
    chip_path_today = chip_Timestamp_path+'\\'+time.strftime('%Y%m%d',time.localtime(time.time()))
    chip_file_today = os.path.join(chip_path_today, time.strftime('%Y%m%d',time.localtime(time.time()))+'_'+file)  
    
    with open(os.path.join(chip_file_path, file), 'w') as up:
        json.dump(chip_information, up, indent=4)

        
    with open(chip_file_today, 'w') as up:
        json.dump(chip_information, up, indent=4) 
    
    print("updated!")      