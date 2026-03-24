import numpy as np

def N2_profile(stratification_kw, dimensional_N2 = 1):

    if stratification_kw == "constant":
        N2_function = lambda z : dimensional_N2 * np.ones_like(z)

    elif stratification_kw == "barotropic":
        N2_function = lambda z : dimensional_N2 * z 

    return N2_function
