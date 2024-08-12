Welcome to Quala_Qblox
=============================================
### the storage area for miscellaneous and incomplete codes~
##### Now you are one of guinea pig of the Qblox users in Quela!
---------------------------------
There are some disordered files in every branches, with some useless debris intermingled inside.  
However, there must be some effective things within them, so let me introduce them~

Start from folder **Modularize**, remember doing the measurement step by step!

### 1. **m1_wideCS**  
Wide scan cavity with big readout amplitude, then you will roughly get cavityies' bare frequency. Remember them to do the next measurement!

### 2. **m2_CavitySpec**
Type in the rough bare frequency of **all qubits** into ro_bare. Also type the chip's name and type. This program will scan bare cavity accurately, and create a chip information file for you to read some important data. On top of that, it will create a QD file, remember all the information to use for the later measurement (If you run this program again, it will rewrite the QD file, so **if you start m3 or even further, then don't run this again!!**)

### 3. **m3_CavitySpec**
Just run, and it will return 2D map of the relationship between readout power and frequency, remember the dress readout power, frequency gap between bare and dress for the next program.

### 4. **m4_BDCavityFit**
Type in the dress readout power(usually 0.01 with ro_atte=30), and fill in window_shift with frequency gap. Remember, **Run bare and dress at the same time!** if you run only bare, then the next run on the same qubit will crash!

(BTW, if you are unlucky, fitting a wrong dress frequency and tune out readout frequency away, then you need to fix it manually with clock_freqs.readout(##) function, fill in ## with bare frequency.)

### 5. **m5_keepPD** 
**Merged into m4**

### 6. **m6_FluxCavSpec**
Run. After measurement, check the picture of the fitting plot. If the fitting is wierd or couplers are disturbing, type 'n' in terminal. You can try to add coupler's flux bias to tune out coupler(in {cx:##}), try the best tune out flux to get a good flux cavity fitting plot.

### 7. **m7_RefIQ**
This program will calibrate the IQ plane of the readout automaticly, so don't hesitate, just run.

### 8. **m8_Cnti2Tone**

### 9. **m9_Cnti2Tone**

### 10. **m10_Cnti2Tone**
**Not done**

### 11. **m11_Cnti2Tone**

### 12. **m12_T2**

### 13. **m13_T1**

--------------------------------------
## Update

#### 2024-08-12 
 * python==3.10
 * qblox-instruments==0.12.0  
 * quantify-core==0.7.4  
 * quantify-scheduler==0.20.0  
 * colorama
 * numpy==1.26.2
 * xarray==2023.12.0
 * QCAT from shiau109
 * scikit-learn
 * pip install git+https://github.com/sebastianprobst/resonator_tools
 * firmware == v0.7.0

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
  



    
