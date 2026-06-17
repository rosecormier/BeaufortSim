"""
Modification of Storer's code cheb.py.

I reference "Spectral Methods in MATLAB" by L. Trefethen.
"""

import numpy as np
import scipy.sparse as ssp

from math import pi

class Parameters:
    
    def __init__(self, args, discretizeRadial = False,
                 discretizeVertical = False,
                 nondimensional = False):

        self.Np     = args["Np"] #No. of azimuthal points (visualization only)
        self.kps    = np.arange(args["k_phi"][0], args["k_phi"][1],
                                args["k_phi"][2])
        self.nmodes = args["nmodes"]
                 
        if discretizeRadial:
            self.Lr     = args["Lr"] #Max. r in phys. space; half of comp. domain
            self.Nr     = args["Nr"] #No. (odd) of radial points in comp. domain
            self.halfNr = self.Nr // 2   
                 
            if not discretizeVertical:
                self.bkgd = args["bkgd"]
                self.kzs  = np.arange(args["k_z"][0], args["k_z"][1],
                                      args["k_z"][2])
        
        if discretizeVertical:
        
            self.Lz = args["Lz"] #Max. depth (i.e., -min(z)) in physical domain
            self.Nz = args["Nz"] #Number of computational gridpoints in z
        
            self.bkgd = "BG"
        
            self.stratification_kw = args["strat_shape"]
            
            if self.stratification_kw in ["doubleTanh", "doubleTanhTWB"]:
            
                #Stratification parameters
                self.g     = -9.81  #m/s^2
                self.rho_0 = 1025.5 #kg/m^3
                self.A_s   = 2.5    #kg/m^3
                self.z_s   = -40    #m
                self.C_s   = 15     #m
                self.A_d   = 1.05   #kg/m^3
                self.z_d   = -200   #m
                self.C_d   = 60     #m
                
                self.doubleTanhParams = [self.g, self.rho_0, self.A_s, self.z_s,
                                         self.C_s, self.A_d, self.z_d, self.C_d]
            
            self.sigmaz = args["sigmaz"]
            
            self.verticalBCs = args["zBCs"]
            
            #Store size of z-differential operators used in gen. eig. problem
            if self.verticalBCs == "homogeneous":
                self.DzSize = self.Nz - 1
            elif self.verticalBCs == "continuousBuoyancy":
                self.DzSize = self.Nz + 1
                
            if not discretizeRadial:
                self.rs      = np.arange(args["r"][0], args["r"][1],
                                         args["r"][2])
                self.rs_plot = args["rs_plot"]
                 
        self.f0 = args["Coriolis"]
    
        if nondimensional:
        
            #Nondimensional parameters are explicitly set
            self.Ro = args["Ro"] #Rossby number
            self.Bu = args["Bu"] #Burger number
            
            #"Dimensional" length scales are set (for compatibility) to unity
            
            if discretizeVertical:
                self.sigmaz = 1
                
            self.sigmar = 1
            
        else:
            
            #Dimensional parameters are explicitly set
            self.N2_far = args["N2_far"] #Far-field squared buoyancy frequency
            self.Umax   = args["bkgdU"]  #Background velocity scale
            self.sigmar = args["sigmar"] #Radial background length scale
            
            #Rossby number is computed
            self.Ro = self.Umax / (self.sigmar * self.f0)
            
            if discretizeVertical:
                #Far-field Burger number is also computed
                self.Bu = (self.N2_far * self.sigmaz / (self.f0 * self.sigmar))**2
                
        #Specify units of various variables
        
        units = {}
        
        if nondimensional:
            units["r"]   = "times $\sigma_r$"
            units["z"]   = "times $\sigma_z$"
            units["kz"]  = "per $\sigma_z$"
            units["psi"] = "times $f_0$ per Ro"
            units["u"]   = "times $\sigma_r f_0$ per Ro"
        else:
            units["r"]   = "m"
            units["z"]   = "m"
            units["kz"]  = "m$^{{-1}}$"
            units["psi"] = "s$^{{-1}}$"
            units["u"]   = "m/s"
            
        self.units = units
            
        #String, representing dimensionality, to be used in output filenames
        if nondimensional:
            self.dimString = "nondimensional"
        else:
            self.dimString = "dimensional"    
        
        self.discretizeRadial   = discretizeRadial
        self.discretizeVertical = discretizeVertical
        self.nondimensional     = nondimensional
    
def Chebyshev(N, xIntervalScaleFactor = 1, xIntervalShiftAmt = 0, 
              sparse = False):
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
        
        #Rescale x as dictated (note: important to do this BEFORE building D)
        x *= xIntervalScaleFactor
        
        #Define the Chebyshev coeffs c_{ij}
        c = np.hstack([2, np.ones(N - 1), 2]) * (-1)**np.arange(0, (N + 1))
            
        X  = np.tile(x, ((N + 1), 1))
        dX = X.T - X.conj()
    
        #Initialize D with off-diag entries computed by eq. 6.5 (Trefethen)
        D_ndarray = ((c[:, np.newaxis] * (1.0 / c.conj())[np.newaxis, :])
                 / (dX + (np.eye(N + 1))))
    
        #Update diagonal entries (currently = 0) using identity 6.6 (Trefethen)
        D = D_ndarray - np.diag(D_ndarray.sum(axis = 1))
        
        if sparse:
            D = ssp.csr_array(D)

        #Shift entries of x (note: important to do this AFTER building D)
        x += xIntervalShiftAmt
        
    return D, x
    
class ChebyshevGeometry:

    def __init__(self, params):
    
        self.method = "Chebyshev"
        
        if (not params.discretizeVertical or not params.discretizeRadial):
            self.sparse = False
    
        elif (params.discretizeRadial and params.discretizeVertical):
            self.sparse = True
                
        if params.discretizeVertical:
            #Compute differentiation matrix and Chebyshev-spaced grid
            Dz, z = Chebyshev(params.Nz, xIntervalScaleFactor = 0.5,
                              xIntervalShiftAmt = -0.5, sparse = self.sparse)

            #Scale gridpoints and variable of differentiation to fit domain
            self.z, self.Dz = z * params.Lz, Dz / params.Lz

            #Second-order z-differentiation matrix
            self.Dz2 = self.Dz @ self.Dz
            
            if not params.discretizeRadial:
                self.r = params.rs

        if params.discretizeRadial:
            #Compute differentiation matrix and Chebyshev-spaced grid
            Dr, r = Chebyshev(params.Nr, sparse = self.sparse)
                            
            #Scale gridpoints and variable of differentiation to fit domain
            self.r, self.Dr = r * params.Lr, Dr / params.Lr
    
            #Second-order r-differentiation matrix
            if self.sparse:
                self.Dr2 = self.Dr @ self.Dr
            elif not self.sparse:
                self.Dr2 = np.matmul(self.Dr, self.Dr)