import os
import json
import time

def load_chip_file(file_path, name, file):
    
    file_exist = False
    file_name = os.path.join(file_path, file)

    for root, dirs, files in os.walk(file_path):
        for il in files:
            if il == file:
                #print(f'{il} exist')
                file_exist = True
                
    if file_exist == False:
        with open(os.path.join(file_path, "blank_chip_information.json"), "r") as blank, open(file_name, 'w') as new:
            new.write(blank.read())
            
        with open(file_name, 'r') as rd:
            chip_information = json.load(rd)
        
        chip_information["basic_information"]["chip_name"] = name
        chip_information["basic_information"]["chip_file"] = file_name
        #chip_information["basic_information"]["create_time"] = 
            
        with open(file_name, 'w') as up:
            json.dump(chip_information, up, indent=4)
        
    else:
        with open(file_name, 'r') as qu:
            chip_information = json.load(qu)
            #print(chip_information["basic_information"]["qubit_name"])
            
    
    return chip_information

def update_chip_file(file_path, file, chip_information):
    with open(os.path.join(file_path, file), 'w') as up:
        json.dump(chip_information, up, indent=4)

            