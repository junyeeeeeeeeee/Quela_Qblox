""" Set the readout parameters here about (1) integration time, (2) reset time """
#// Reset time can be set in str 'auto', which will use 10*T1 if there is the T1 record in the QD_agent. Otherwise, set it as 250e-6. 
Readers:dict = {
    "q0":{
        "integration_time":1e-6, "reset_time":250e-6
    },
    "q1":{
        "integration_time":1e-6, "reset_time":250e-6
    }
}

def get_manully_integration_time(target_q:str):
    return Readers[target_q]["integration_time"]

def get_manully_reset_time(target_q:str):
    return Readers[target_q]["integration_time"] if str(Readers[target_q]["integration_time"]).lower() != 'auto' else 0
