import os, sys, time, inspect, tomli
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from Modularize.Council.MainBrain import Exp_Encyclopedia


class Maid(Exp_Encyclopedia):
    def __init__(self,exp_type:str, target_qs:list):
        self.exp_type=exp_type
        self.__exp_target_qs__ = target_qs
        super().__init__(exp_type)
        

    def __generate_toml__(self):
        # Get all attributes of the class excluding built-in ones
        attributes = [name for name, _ in inspect.getmembers(self) if not name.startswith("__")]
        # exclude_attrs = ["exp_type"]
        # attributes = [attr for attr in attributes if attr not in exclude_attrs]
        # Open a file to write the toml content
        self.__file_path__=f"{self.exp_type}_ParasAssignment.toml"
        with open(self.__file_path__, "w") as file:
            # For each target (q0 and q1), create a section with the same attributes
            file.write("# Set variables with the rule: [start, end, pts]\n\n")
            for attribute in self.__shared_attr__:
                file.write(f"{attribute} =        # type in  {type(getattr(self, attribute))}\n")
            file.write("\n")
            for target in self.__exp_target_qs__:
                file.write(f'[{target}]\n')  # Create a section for each target
                for attribute in attributes:
                    if attribute not in self.__shared_attr__:
                        file.write(f"  {attribute} =        # type in  {type(getattr(self, attribute))}\n")  # Inline comments
                        file.write("\n")  #
                file.write("\n")  # Add a blank line between sections for readability

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
        print(paras)
        print(shared_paras)
        
if __name__ == "__main__":
    # maid = Maid('S1', ["q0","q1"])
    # maid.__generate_toml__()

    Paras = Coordinator("/Users/ratiswu/Documents/GitHub/Quela_Qblox/S1_ParasAssignment.toml")
