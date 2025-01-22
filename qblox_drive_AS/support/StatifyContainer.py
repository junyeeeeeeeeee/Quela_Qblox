from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity

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

if __name__ == "__main__":
    
    Bob = tester("B")
    Alice = tester("A")
    Ratis = tester("R")
    rec = {"Bob":Bob, "Alice":Alice, "Ratis":Ratis}
    
    container = Statifier()
    for i in rec:
        container.serialize(i, rec[i], version="v0.1.1")

    person:tester = container.summon_discriminator("bob")
    person.play()