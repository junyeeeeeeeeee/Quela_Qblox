class Exp_Encyclopedia():
    def __init__(self,exp_type:str):
        """
        exp_type should be in the following:\n
            S1. BroadBandCS,
            S2. ZoomCS,
            S3. PowerCavity,
            S4. FluxCavity,
            S5. Conti2Tone,
            S6. FluxQubit,
            S7. Rabi,
            S8. OneShot,
            S9. T2,
            S10. T1,
            ----------------
            C1. ROFcali,
            C2. XYFcali,
            C3. PIampcali,
            C4. HalfPIcali,
            ----------------
            A1. TimeMonitor,
            A2. ZgateT1
        """
        self.__exp__ = exp_type.lower()
        self.__expVarable_Reminder__()
        self.__shared_attr__ = ["avg_n", "exp_type"]

    def __expVarable_Reminder__(self):
        match self.__exp__:
            case "s1" | "broadbandcs" | "bbcs":
                self.freq_range = list([])
                self.avg_n:int = 100