Welcome to Quala_Qblox
=============================================
### the storage area for miscellaneous and incomplete codes~
##### Now you are one of guinea pig of the Qblox users in Quela!
---------------------------------
There are some disordered files in every branches, with some useless debris intermingled inside.  
However, there must be some effective things within them, so let me introduce them~

1. **Oidle5QMeas_V0_X.ipynb**  
This is the main program now to run Qblox
  
2. **Support.py**  
This is the library for running other jupyter notebook files  
  
3. **chip_store.py**  
This is the library for storing every chip information


--------------------------------------
## Update

#### Oidle5Q V0_4
1. It's now able to measure any qubit(s) at the same time. Update to Oidle5QMeas_V0_4.ipynb and support.py
    * Measurements under two-tone can only measure one qubit in one measurement.
2. Function 'pulse_preview()' is updated into support.py
    * Figure out the problem why two-tone can't draw the measuring pulse: if the pulse is longer than 1us, then the pulse will not show.
    * Waiting for Qblox to fix.
    * After fixing the problem, clear the codes in Oidle5QMeas.
