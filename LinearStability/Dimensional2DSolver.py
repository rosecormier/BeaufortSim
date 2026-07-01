"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables should be specified in SI base units.
"""

import argparse
import numpy as np
import scipy.sparse.linalg as sspla

import sys
import time
import timeit

from BuildDiscreteOperators import *
from Chebyshev import Parameters, ChebyshevGeometry
from SaveToNetCDF import SaveToNetCDF
from VisualizationFunctions import RunVisFromSavedData

from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser()
parser.add_argument("-Nr", 
                    help = "Number (must be ODD) of r-grid points for direct computation",
                    type = int, default = 201)
parser.add_argument("-Nz", 
                    help = "Number of z-grid points for direct computation",
                    type = int, default = 20)
parser.add_argument("-Lr", 
                    help = "DIMENSIONAL radius of the physical domain (m)",
                    type = float, default = 2.5e6)
parser.add_argument("-Lz",
                    help = "DIMENSIONAL depth (> 0) of physical domain (m)",
                    type = float, default = 1e3)
parser.add_argument("--zBCs",
                    help = "Vertical boundary conditions on streamfunction",
                    choices = ("continuousBuoyancy", "homogeneous"),
                    type = str, default = "continuousBuoyancy")
parser.add_argument("--strat_shape",
                    help = "Shape of ambient buoyancy profile",
                    choices = ("constant", "TWB", "doubleTanh", "doubleTanhTWB"),
                    type = str, default = "linear")
parser.add_argument("--N2_far",
                    help = "Far-field squared buoyancy frequency (1/s^2)",
                    type = float, default = 1e-6)
parser.add_argument("-f0", "--Coriolis",
                    help = "Coriolis frequency f0 (Hz)",
                    type = float, default = 1.4e-4)
parser.add_argument("-U", "--bkgdU",
                    help = "Characteristic scale for background velocity (m/s)",
                    type = float, default = 5e-2)
parser.add_argument("--sigmar",
                    help = "Radial length scale of gyre (m)",
                    type = float, default = 2.5e5)
parser.add_argument("--sigmaz",
                    help = "Vertical length scale of gyre (m)",
                    type = float, default = 3e2)
parser.add_argument("-Np", 
                    help = "Number of phi-points for visualization", 
                    type = int, default = 50)
parser.add_argument("--k_phi", 
                    help = "Azimuthal wavenumbers; enter as -kp start stop step",
                    type = float, default = [1, 3, 1], nargs = 3)
parser.add_argument("--nmodes", 
                    help = "Number of modes of instability to be considered",
                    type = int, default = 1)
parser.add_argument("--useSaved",
                    help = "Flag to skip eig. solver and go straight to visualizing previously saved data",
                    action = "store_true", default = False)
args = vars(parser.parse_args())

def QG_Vortex_Stability():

    #Initialize parameters and set up geometry for Chebyshev solver
    params = Parameters(args, discretizeRadial = True, 
                        discretizeVertical = True)
    geom   = ChebyshevGeometry(params)
        
    if not args["useSaved"]:

        #Build discrete operators
        ComputeRecips(params, geom)
        BuildBkgdOperators(params, geom)
        BuildHorizontalLaplacian(params, geom)
         
        #Information about wavenumbers and modes
        kφs, nmodes = params.kps, params.nmodes
    
        #Initialize arrays to store results of eigen-computation
        growth = np.zeros([kφs.shape[0], nmodes])
        prop   = np.zeros([kφs.shape[0], nmodes])
        modes  = np.zeros([kφs.shape[0], ((params.halfNr + 1) * (params.Nz + 1)),
                           nmodes],
                          dtype = complex)
    
        #Solve generalized eigenvalue problem
    
        for kφ_idx in range(kφs.shape[0]):
    
            kφ = kφs[kφ_idx]
        
            print("Solving for kφ =", kφ)
    
            #Build matrices "A" and "B"
            B = BuildMatrixB(params, geom, kφ)
            A = BuildMatrixA(params, geom)
            
            #Find (nmodes) eigenvalue/vector pairs, indirectly
            eigVals, eigVecs = sspla.eigs(A, k = nmodes, M = B, which = "LI")
            #Note: "which = 'LI'" prioritizes finding eigenvalue with largest 
            # imaginary part.
            #Assuming indirectly-found generalized eigvals will occur in
            # conjugate pairs, the largest imaginary part of any eigenvalue 
            # equals the largest negative imaginary part of any eigenvalue -- 
            # i.e., largest growth rate.

            #Indexing that sorts eigvals by DESCENDING Im(c)
            indSort = np.argsort(eigVals.imag)[::-1]
    
            eigVals = eigVals[indSort]    #Sort eigvals
            eigVecs = eigVecs[:, indSort] #Sort eigvecs in the same order
            ωs      = eigVals * kφ        #Corresponding ω-values for this kφ
            
            #Growth rates and propagation speeds
            growth[kφ_idx, 0:nmodes] = ωs.imag
            prop[kφ_idx, 0:nmodes]   = ωs.real
            
            modesLen = len(modes[kφ_idx, :, 0:nmodes])
            
            if params.verticalBCs == "homogeneous":
                #Update 'modes' at interior points only
                modes[kφ_idx, 
                      ((np.mod(np.arange(modesLen), (params.Nz + 1)) != 0)
                       & (np.mod(np.arange(modesLen), 
                                 (params.Nz + 1)) != params.Nz)
                       & (np.arange(modesLen) > params.Nz)),
                      0:nmodes] = eigVecs[:, 0:nmodes]
            elif params.verticalBCs == "continuousBuoyancy":
                #Update 'modes' at all z-points but interior points only in r
                modes[kφ_idx, 
                      (np.arange(modesLen) > params.Nz), 
                      0:nmodes] = eigVecs[:, 0:nmodes]
    
        #Save results to nc file
        SaveToNetCDF(params, geom, growth, prop, modes)
        
    RunVisFromSavedData(params, geom)
        
if __name__ == '__main__': #For testing
    QG_Vortex_Stability()