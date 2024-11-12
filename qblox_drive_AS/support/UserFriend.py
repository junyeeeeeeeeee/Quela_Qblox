"""
This script sorts the functions to print different color to emphasize a important massage in the console.\n
Here are different level: slightly -> should look -> very important -> warning . 
"""
from colorama import init, Fore, Back, Style


# slightly important level
def slightly_print(msg:str):
    print(Fore.YELLOW + msg + Style.RESET_ALL)

# print should look level
def eyeson_print(msg:str):
    print(Back.WHITE + Fore.GREEN + Style.BRIGHT + msg + Style.RESET_ALL)

# very important level
def highlight_print(msg:str):
    print(Back.YELLOW + Fore.RED + Style.BRIGHT + f"*** {msg} ***" + Style.RESET_ALL)

# warning level
def warning_print(msg:str):
    print(Back.RED + Fore.YELLOW + Style.BRIGHT + f"!!! WARNING: {msg} !!!" + Style.RESET_ALL)

# interaction level
def mark_input(msg:str)->str:
    ans = input(Fore.RED + Style.BRIGHT + msg + Style.RESET_ALL)
    
    return ans

if __name__ == "__main__":
    highlight_print("hi")
    warning_print("hi")
    mark_input("hi")