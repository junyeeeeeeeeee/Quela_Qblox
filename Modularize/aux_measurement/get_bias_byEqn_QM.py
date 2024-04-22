from numpy import asarray, ndarray, array
import sympy as sp


def find_nearest(ary:ndarray, near_target:float):
    """ find the element  which is closest to the given near_target in the given array"""
    ary = asarray(ary)
    idx = (abs(ary - near_target)).argmin()
    return float(ary[idx])



def get_biasWithFq_from(fitting_popts:list,target_fq_Hz:float, flux_guard:float=0.4):
    """
    After we fit the tarnsition freq vs bias, we can get the bias according to the given `target_fq_Hz` for the `target_q`.\n
    ### The given `target_fq_Hz` should unit in Hz.\n
    ### The given `fitting_popts` should follow the order: [a, b, Ec, Ej_sum, d]
    Return the bias unit in V.
    """
    
    z = sp.Symbol('z',real=True)
    
    a,b,Ec,Ej_sum,d = fitting_popts[0], fitting_popts[1], fitting_popts[2], fitting_popts[3], fitting_popts[4]
    to_solve = sp.sqrt(8*Ej_sum*Ec*sp.sqrt(sp.cos(a*(z-b))**2+d**2*sp.sin(a*(z-b))**2))-Ec - target_fq_Hz*1e-9
    candidators = array(sp.solvers.solve(to_solve, z))
    
    if candidators.shape[0] == 0:
        answer = 'n'
        print(f"Can NOT find a bias makes the fq @ {target_fq_Hz*1e-9} GHz !")
    else:
        answer = find_nearest(candidators, 0) if find_nearest(candidators, 0) < flux_guard else flux_guard

    return answer