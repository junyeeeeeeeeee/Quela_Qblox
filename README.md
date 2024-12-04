Welcome to Quala_Qblox
=============================================
### the storage area for miscellaneous and incomplete codes~
##### Now you are one of guinea pig of the Qblox users in Quela!
---------------------------------
## Update
### 2024-11-07
 1. All the SOP measurements (m), calibrations (c), ZgateT1 (auxA) and TimeMoniter (auxB) are multiplexed readout supported.
 2. Creating Readout setup for manually reset_time and integration_time setting.
 3. QD_backup folder remaned.
 4. Demolish the exp in to 3 parties, namely 'Waiter', 'Executor', and 'Analyzer'.

### 2024-10-08 
 1. pip install -r required_packages.txt (whole path of required_packages.txt in your local is required)
 2. conda install git (if `git` hasn't been installed)
 3. pip install git+https://github.com/sebastianprobst/resonator_tools
 4. pip install -e QCAT (whole path of QCAT repo in your local is required)
 * firmware == v0.7.0
 * python>=3.10

#### Single Qubit Charaterization 08-05-2024
You will need the following package with pip and the firmware:  
 * qblox-instruments == v0.12.0  
 * quantify-core == v0.7.4  
 * quantify-scheduler == v0.19.0  
 * colorama
 * firmware == v0.7.0

#### Oidle5Q V0_4
1. It's now able to measure any qubit(s) at the same time. Update to Oidle5QMeas_V0_4.ipynb and support.py
    * Measurements under two-tone can only measure one qubit in one measurement.
2. Function 'pulse_preview()' is updated into support.py
    * Figure out the problem why two-tone can't draw the measuring pulse: if the pulse is longer than 1us, then the pulse will not show.
    * Waiting for Qblox to fix.
    * After fixing the problem, clear the codes in Oidle5QMeas.
  



    
