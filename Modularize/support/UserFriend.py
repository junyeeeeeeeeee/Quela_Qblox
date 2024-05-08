"""
This script sorts the functions to print different color to emphasize a important massage in the console.\n
Here are different level: slightly -> should look -> very important -> warning . 
"""
from colorama import init, Fore, Back, Style
init(autoreset=True)

# slightly important level
def slightly_print(msg:str):
    init(autoreset=True)
    print(Fore.YELLOW + msg)

# print should look level
def eyeson_print(msg:str):
    init(autoreset=True)
    print(Back.WHITE + Fore.GREEN + Style.BRIGHT + msg)

# very important level
def highlight_print(msg:str):
    init(autoreset=True)
    print(Back.YELLOW + Fore.RED + Style.BRIGHT + f"*** {msg} ***")

# warning level
def warning_print(msg:str):
    init(autoreset=True)
    print(Back.RED + Fore.YELLOW + Style.BRIGHT + f"!!! WARNING: {msg} !!!")


if __name__ == "__main__":
    highlight_print("hi")
    warning_print("hi")