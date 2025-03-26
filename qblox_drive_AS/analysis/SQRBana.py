import os
from xarray import open_dataset, Dataset
from numpy import array, inf
from qblox_drive_AS.support.QDmanager import QDmanager, QuantumDevice
from qblox_drive_AS.support import rotate_data
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exp_decay(m, A, p, B):
    """Exponential decay function for RB fitting."""
    return A * (p ** m) + B

def fit_rb_decay(m_values, data):
    """
    Fits RB decay data to the exponential model A * (p^m) + B.
    
    Parameters:
    - m_values: Array of sequence lengths
    - data: Measured probabilities

    Returns:
    - params: Fitted (A, p, B)
    - covariance: Covariance matrix from curve_fit
    - fitted_curve: Computed fit values
    - fidelity: Average gate fidelity (1 + p) / 2
    """
    # Initial guesses
    A_guess = max(data) - min(data)  # Approximate amplitude
    p_guess = 0.98  # Typical decay rate close to 1
    B_guess = min(data)  # Offset

    p0 = [A_guess, p_guess, B_guess]  # Initial parameter estimates
    
    # Bounds: 
    # A > 0, 0 < p < 1 (physical decay), B unrestricted
    bounds = ([-inf, 0, -inf], [inf, 1, inf])

    # Fit the data
    params, covariance = curve_fit(exp_decay, m_values, data, p0=p0, bounds=bounds)

    # Extract parameters
    A_fit, p_fit, B_fit = params
    fitted_curve = exp_decay(m_values, A_fit, p_fit, B_fit)

    # Compute average gate fidelity
    F_cliford = 1 - ((1 - p_fit) / 2)
    F_gate = 1 - ((1 - p_fit) / 2)/1.875


    return params, covariance, fitted_curve, (F_cliford, F_gate)


def SQRB_ana(ds:Dataset, rotation_angle_degree:dict, pic_save_folder:str, quantum_device:QuantumDevice):
    g_n = ds.coords["gate_length"].data
    for var in ds.data_vars:
        tg = quantum_device.get_element(var).rxy.duration()

        avg_IQ_data = ds[var].data*1000
        Ir = rotate_data(avg_IQ_data, rotation_angle_degree[var][0])[0]

        arams, covariance, fitted_curve, (F_cliford, F_gate) = fit_rb_decay(g_n, Ir)

        plt.scatter(g_n, Ir)
        plt.plot(g_n,fitted_curve,c='red')
        plt.ylabel("I (mV)")
        plt.xlabel("Gates")
        plt.title(f"{var}, tg = {round(tg*1e9)} ns \n Gate Fidelity = {round(100*F_gate,2)} %")
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(pic_save_folder,f"{var}_SQRB.png"))
        plt.close()

        





