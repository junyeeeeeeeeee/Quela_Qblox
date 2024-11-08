import os, sys, time, inspect, tomli
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from Modularize.Council.MainBrain import Exp_Encyclopedia

class Coordinator(Exp_Encyclopedia):
    def __init__(self, toml_path:str):
        self.__toml_path__ = toml_path
        super().__init__(exp_type="_")
        self.__toml_decoder__()
        

    def __toml_decoder__(self):
        # Load data from TOML file
        with open(self.__toml_path__, "rb") as file:
            toml_data = tomli.load(file)
        
        # Assign values to object attributes
        # Dynamically create attributes for each section in the TOML data
        for section, attributes in toml_data.items():
            setattr(self, section, attributes)
        
        paras = {}
        shared_paras = {}
        for attribute, value in vars(self).items():
            if not attribute.startswith("__"):
                if attribute not in self.__shared_attr__:
                    paras[attribute] = value
                else:
                    shared_paras[attribute] = value
        print("Unique paras: ",paras)
        print("\nShared paras: ",shared_paras)


if __name__ == "__main__":
    Paras = Coordinator("/Users/ratiswu/Documents/GitHub/Quela_Qblox/S1_ExpParasServey.toml")