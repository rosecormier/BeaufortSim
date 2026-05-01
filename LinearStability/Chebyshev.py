"""
Modification of Storer's code cheb.py.

I reference "Spectral Methods in MATLAB" by L. Trefethen.
"""

import numpy as np

from math import pi

class Parameters:
    
    def __init__(self, args, discretizeVertical = False):
    
        self.Nmax   = args["buoyancyfreq"]
        self.f0     = args["Coriolis"]
        self.Umax   = args["bkgdU"]
        self.sigmar = args["sigmar"]
    
        self.Lr     = args["Lr"] #Max. r in phys. space; half of comp. domain
        self.Nr     = args["Nr"] #No. (odd) of radial points in comp. domain
        self.halfNr = self.Nr // 2   
        self.Np     = args["Np"] #No. of azimuthal points (visualization only)
        self.kps    = np.arange(args["k_phi"][0], args["k_phi"][1], args["k_phi"][2])
        
        self.nmodes = args["nmodes"]
        
        self.discretizeVertical = discretizeVertical
        
        if not discretizeVertical:
            self.bkgd = args["bkgd"]
            self.kzs  = np.arange(args["k_z"][0], args["k_z"][1], args["k_z"][2])
        
        elif discretizeVertical:
        
            self.bkgd = "BG"
        
            self.stratification_kw = args["strat_shape"]
            self.sigmaz            = args["sigmaz"]

            self.Lz = args["Lz"] #Max. depth (i.e., -min(z)) in physical domain
            self.Nz = args["Nz"] #Number of computational gridpoints in z
    
def Chebyshev(N, xTransform = None):
    """
    Computes the Chebyshev differentiation matrix on N+1 points (i.e., N
     intervals).
    Returns:
      D = (N+1) x (N+1) Chebyshev differentiation matrix
      x = Chebyshev grid
    """

    if N == 0:
        D, x = 0, 1

    else:

        #Create a vector of N+1 Chebyshev-spaced components from 1 to -1
        x = np.cos(pi * np.arange(0, (N + 1)) / N)

        #Define the Chebyshev coeffs c_{ij}
        c = np.hstack([2, np.ones(N - 1), 2]) * (-1)**np.arange(0, (N + 1))
            
        X  = np.tile(x, ((N + 1), 1))
        dX = X.T - X.conj()
    
        #Initialize D with off-diag entries computed by eq. 6.5 (Trefethen)
        D = ((c[:, np.newaxis] * (1.0 / c.conj())[np.newaxis, :])
                 / (dX + (np.eye(N + 1))))
    
        #Update diagonal entries (currently = 0) using identity 6.6 (Trefethen)
        D = D - np.diag(D.sum(axis = 1))

        #Shift entries of x (note: important to do this AFTER building D)
        if xTransform is not None:
            x = xTransform(x)
        
    return D, x
    
class ChebyshevGeometry:

    def __init__(self, params):
        
        self.method = "Chebyshev"

        #Compute differentiation matrix and Chebyshev-spaced grid
        Dr, r = Chebyshev(params.Nr)
                        
        #Scale gridpoints and variable of differentiation to fit domain
        self.r, self.Dr = r * params.Lr, Dr / params.Lr

        #Second-order r-differentiation matrix
        self.Dr2 = np.matmul(self.Dr, self.Dr)
        
        if params.discretizeVertical:
        
            def zTransform():
                return lambda z : (z - 1) / 2
                
            #Compute differentiation matrix and Chebyshev-spaced grid
            Dz, z = Chebyshev(params.Nz, xTransform = zTransform())

            #Scale gridpoints and variable of differentiation to fit domain
            self.z, self.Dz = z * params.Lz, Dz / params.Lz

            #Second-order z-differentiation matrix
            self.Dz2 = np.matmul(self.Dz, self.Dz)