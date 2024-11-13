import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from numpy import arange, linspace, logspace
from quantify_scheduler.backends.graph_compilation import SerialCompiler
### Note: The qubits are labeled with the physical position on the chip from the left, and starts from 0.

#%%
# Please register the dr and its corresponding cluster ip here first!
ip_register = {
    "dr1":"192.168.1.11",
    "dr2":"192.168.1.10",
    "dr3":"192.168.1.13",
    "dr4":"192.168.1.81",
    "drke":"192.168.1.242"
} # all keys in lower
port_register = {
    "192.168.1.10":"5010",
    "192.168.1.11":"5011",
    "192.168.1.12":"5012",
    "192.168.1.13":"5013",
    "192.168.1.81":"5081",
    "192.168.1.242":"5242"
}
