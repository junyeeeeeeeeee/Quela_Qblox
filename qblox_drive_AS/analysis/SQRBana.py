import os
from xarray import open_dataset, Dataset
from numpy import array, inf, moveaxis, average, std
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
        data = moveaxis(ds[var].data, 0, 1)*1000  # ["mixer","random_circuits","gate_length"] -> ["random_circuits","mixer","gate_length"]
        circuit_dep_data = []
        for circuit_idx, a_rand_circuit_data in enumerate(data):

            Ir = rotate_data(a_rand_circuit_data, rotation_angle_degree[var][0])[0]
            circuit_dep_data.append(Ir)

        avg_signal = average(array(circuit_dep_data), axis=0)

        arams, covariance, fitted_curve, (F_cliford, F_gate) = fit_rb_decay(g_n, avg_signal)

        # for a_result in circuit_dep_data:
        #     plt.scatter(g_n, a_result, c='#87CEFA',s=10)
        plt.errorbar(g_n, avg_signal,yerr=std(array(circuit_dep_data), axis=0), c='red', fmt='*')
        plt.plot(g_n,fitted_curve,c='orange')
        plt.ylabel("I (mV)")
        plt.xlabel("Gates")
        plt.title(f"{var}, tg = {round(tg*1e9)} ns \n Gate Fidelity = {round(100*F_gate,2)} %")
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(pic_save_folder,f"{var}_SQRB.png"))
        plt.close()


if __name__=="__main__":
    QD_agent = QDmanager("qblox_drive_AS/QD_backup/20250327/DR4#81_SumInfo.pkl")
    QD_agent.QD_loader()
    ds = open_dataset("qblox_drive_AS/Meas_raw/20250327/H14M08S31/SQRB_20250327141515.nc")


    SQRB_ana(ds,QD_agent.rotate_angle,"qblox_drive_AS/Meas_raw/20250327/H14M08S31",QD_agent.quantum_device)


