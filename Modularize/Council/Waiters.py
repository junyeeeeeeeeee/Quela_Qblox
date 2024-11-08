import os, sys, time, inspect, tomli
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from Modularize.Council.MainBrain import Exp_Encyclopedia


class Canvasser(Exp_Encyclopedia):
    def __init__(self,exp_type:str, target_qs:list):
        self.exp_type=exp_type
        self.__exp_target_qs__ = target_qs
        super().__init__(exp_type)
        

    def __generate_ExpParas_servey__(self):
        # Get all attributes of the class excluding built-in ones
        attributes = [name for name, _ in inspect.getmembers(self) if not name.startswith("__")]
        exclude_attrs = []
        attributes = [attr for attr in attributes if attr not in exclude_attrs]
        # Open a file to write the toml content
        self.__file_path__=f"{self.exp_type}_ExpParasServey.toml"
        with open(self.__file_path__, "w") as file:
            # For each target (q0 and q1), create a section with the same attributes
            file.write("# Measurement paras configs editor\n\n")
            for attribute in self.__shared_attr__:
                file.write(f"{attribute} =        # type in  { type(getattr(self, attribute))}\n")
            file.write("\n")
            for target in self.__exp_target_qs__:
                file.write(f'[{target}]\n')  # Create a section for each target
                for attribute in attributes:
                    if attribute not in self.__shared_attr__:
                        file.write(f"  {attribute} =        # type in {type(getattr(self, attribute)) if type(getattr(self, attribute)) != list else 'list, rule: [ start, end, pts ] or [ fixed value ]'}\n")  # Inline comments
                        file.write("\n")  #
                file.write("\n")  # Add a blank line between sections for readability


        
if __name__ == "__main__":
    Survey = Canvasser('S1', ["q0","q1"])
    Survey.__generate_ExpParas_servey__()

    
