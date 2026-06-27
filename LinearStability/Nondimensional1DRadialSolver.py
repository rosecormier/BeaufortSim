"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables should dimensionless unless otherwise specified.
"""

import argparse
import numpy as np
import scipy.linalg as spalg

import sys
import time
import timeit

from BuildDiscreteOperators import *
from Chebyshev import Parameters, ChebyshevGeometry
from SaveToNetCDF import SaveToNetCDF
from VisualizationLinearStability import RunVisFromSavedData

parser = argparse.ArgumentParser()
parser.add_argument("-Nr", 
                    help = "Number (must be ODD) of grid points for direct computation",
                    type = int, default = 501)
parser.add_argument("-Lr", 
                    help = "DIMENSIONLESS radius of the physical domain (i.e., as a multiple of sigma_r)",
                    type = float, default = 10)
parser.add_argument("-Ro",
                    help = "Rossby number of background flow", 
                    type = float, default = 1e-1)
parser.add_argument("-Bu",
                    help = 'Burger number of background flow',
                    type = float, default = 2.5e-3)
parser.add_argument("-f0", "--Coriolis",
                    help = "DIMENSIONAL Coriolis frequency f0 (Hz)",
                    type = float, default = 1.4e-4)
parser.add_argument("--bkgd",
                    help = "Background flow to use ('GM' or 'BG')",
                    type = str, default = "BG")
parser.add_argument("-Np", 
                    help = "Number of phi-points for visualization", 
                    type = int, default = 50)
parser.add_argument("--k_phi", 
                    help = "Azimuthal wavenumbers; enter as --k_phi start stop step",
                    type = float, default = [1, 3, 1], nargs = 3)
parser.add_argument("--k_z", 
                    help = "DIMENSIONLESS vertical wavenumbers (i.e., as a multiple of sigma_z^{-1}); enter as --k_z start stop step",
                    type = float, default = [0, 2e-1, 3e-3], nargs = 3)
parser.add_argument("--nmodes", 
                    help = "Number of modes of instability to be considered",
                    type = int, default = 1)
parser.add_argument("--useSaved",
                    help = "Flag to skip eig. solver and go straight to visualizing previously saved data",
                    action = "store_true", default = False)
args = vars(parser.parse_args())

def QG_Vortex_Stability():

    #Initialize parameters and set up geometry for Chebyshev solver
    params = Parameters(args, nondimensional = True)
    geom   = ChebyshevGeometry(params)

    if not args["useSaved"]:

        #Build discrete operators
        ComputeRecips(params, geom)
        BuildBkgdOperators(params, geom)
        BuildHorizontalLaplacian(params, geom)
         
        #Information about wavenumbers and modes
        kφs, kzs, nmodes = params.kps, params.kzs, params.nmodes
    
        #Initialize arrays to store results of eigen-computation
        growth = np.zeros([kzs.shape[0], kφs.shape[0], nmodes])
        prop   = np.zeros([kzs.shape[0], kφs.shape[0], nmodes])
        modes  = np.zeros([kzs.shape[0], kφs.shape[0], params.halfNr + 1,
                           nmodes], dtype = complex)
    
        #Solve generalized eigenvalue problem
        
        for kz_idx in range(kzs.shape[0]):
    
            kz  = kzs[kz_idx]
            kz2 = kz**2
      
            for kφ_idx in range(kφs.shape[0]):
    
                kφ = kφs[kφ_idx]
        
                print("Solving for kφ =", kφ, ", kz =", kz)
    
                #Build matrices "A" and "B"
                B = BuildMatrixB(params, geom, kφ, kz = kz)
                A = BuildMatrixA(params, geom)
    
                #Find generalized eigenspace (directly)
                
                t0 = timeit.timeit()
            
                #Compute eigvals c and eigvecs psi with direct solver
                eigVals, eigVecs = spalg.eig(A, B)
    
                solveTime = timeit.timeit() - t0 #Time for direct solver
                
                #Indexing that sorts eigvals by ASCENDING Im(c)
                indSort = np.argsort(eigVals.imag)
    
                eigVals = eigVals[indSort]    #Sort eigvals
                eigVecs = eigVecs[:, indSort] #Sort eigvecs in the same order
                ωs      = eigVals * kφ        #Corresponding ω-values
               
                growth[kz_idx, kφ_idx, :]    = -ωs[0:nmodes].imag
                prop[kz_idx, kφ_idx, :]      = ωs[0:nmodes].real
                modes[kz_idx, kφ_idx, 1:, :] = eigVecs[:, 0:nmodes]
    
        #Dimensionalize eigenvalues before saving
        growthDim = growth * params.f0 * params.Ro
        propDim   = prop * params.f0 * params.Ro
    
        #Save results to nc file
        SaveToNetCDF(params, geom, growthDim, propDim, modes)
    
    RunVisFromSavedData(params, geom)

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()