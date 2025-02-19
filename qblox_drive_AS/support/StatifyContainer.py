from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from xarray import DataArray


class tester():
    def __init__(self,msg:str):
        self.rec = msg
    def play(self):
        print(self.rec)


class Statifier():  # state + indentify + er = statifier 
    def __init__(self):
        self.__container:dict = {}

    @property
    def elements( self ):
        return self.__container

    def serialize(self, name:str, discriminator:GMMROFidelity, version:str=""):
        setattr(self, name, discriminator)
        self.__container[name] = version
    
    def summon_discriminator(self, name:str)->GMMROFidelity:
        if name not in self.__container:
            raise NameError(f"The name you gave didn't show in container, check it please:\n{self.__container}")
        return getattr(self, name)

    def check_model_alive(self, dataArray:DataArray, name:str, show_plot:bool=True):
        """ dataArray is from OneShot nc dataset and *1000 """
        md = self.summon_discriminator(name)
        md._import_data(dataArray)
        md._start_analysis()
        g1d_fidelity = md.export_G1DROFidelity()

        if show_plot:
            plot_readout_fidelity(dataArray, md, g1d_fidelity)



if __name__ == "__main__":
    
    Bob = tester("Age: 27")
    Alice = tester("Age: 25")
    Ratis = tester("Age: 26")
    rec = {"Bob":Bob, "Alice":Alice, "Ratis":Ratis}
    
    container = Statifier()
    for i in rec:
        container.serialize(i, rec[i], version="v0.1.1")
    
    print(container.elements)

    person:tester = container.summon_discriminator("Bob")
    person.play()