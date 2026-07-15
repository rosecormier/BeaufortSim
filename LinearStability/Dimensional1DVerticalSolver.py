"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

All variables should be specified in SI base units.
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
from VisualizationFunctions import RunVisFromSavedData

parser = argparse.ArgumentParser()
parser.add_argument("-Nz", 
                    help = "Number of grid points for direct computation",
                    type = int, default = 100)
parser.add_argument("-Lz", 
                    help = "DIMENSIONAL depth of the physical domain (m)",
                    type = float, default = 1e3)
parser.add_argument("--zBCs",
                    help = "Vertical boundary conditions on streamfunction",
                    choices = ("continuousBuoyancy", "constantStreamfunction", "homogeneous"),
                    type = str, default = "continuousBuoyancy")
parser.add_argument("--strat_shape",
                    help = "Shape of ambient squared buoyancy frequency profile",
                    choices = ("constant", "TWB", "doubleTanh", "doubleTanhTWB"),
                    type = str, default = "doubleTanhTWB")
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
                    help = "Azimuthal wavenumbers; enter as --k_phi start stop step",
                    type = float, default = [1, 5, 1], nargs = 3)
parser.add_argument("-r", 
                    help = "DIMENSIONAL r-values to solve at (m); enter as (-r start stop step) or as r-values themselves",
                    type = float, default = [1, 2.5e6, 1e4], nargs = "*")
parser.add_argument("--nmodes", 
                    help = "Number of modes of instability to be considered",
                    type = int, default = 1)
parser.add_argument("--useSaved",
                    help = "Flag to skip eig. solver and go straight to visualizing previously saved data",
                    action = "store_true", default = False)
parser.add_argument("--rs_plot",
                    help = "r-values to plot eigenmodes at (will choose discretized r-values closest to args provided)",
                    type = float,
                    default = [0, 2.5e5, 2.5e6], nargs = "*"
                   )
args = vars(parser.parse_args())

def QG_Vortex_Stability():

    #Initialize parameters and set up geometry for Chebyshev solver
    params = Parameters(args, discretizeVertical = True)
    geom   = ChebyshevGeometry(params)
        
    if not args["useSaved"]:

        #Information about wavenumbers and modes
        kφs, rs, nmodes = params.kps, params.rs, params.nmodes

        #Initialize arrays to store results of eigen-computation
        growth = np.zeros([rs.shape[0], kφs.shape[0], nmodes])
        prop   = np.zeros([rs.shape[0], kφs.shape[0], nmodes])
        modes  = np.zeros([rs.shape[0], kφs.shape[0], params.Nz + 1, nmodes],
                          dtype = complex)
    
        #Solve generalized eigenvalue problem
        
        BuildGlobalUDerivMatrices(params, geom)
        
        for r_idx in range(rs.shape[0]):
        
            #Build discrete operators
            ComputeRecips(params, geom, r = rs[r_idx])
            BuildBkgdOperators(params, geom, r_idx = r_idx)

            for kφ_idx in range(kφs.shape[0]):
    
                kφ = kφs[kφ_idx]
        
                print("Solving for kφ =", kφ, ", r =", rs[r_idx])
    
                #Build matrices "A" and "B"
                B = BuildMatrixB(params, geom, kφ, r_idx = r_idx)
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

                if r_idx == 0:
                    print(f"Eigval at r = {rs[r_idx]}: ", eigVals[0])
                
                growth[r_idx, kφ_idx, :] = -ωs[0:nmodes].imag
                prop[r_idx, kφ_idx, :]   = ωs[0:nmodes].real

                if params.verticalBCs == "homogeneous":
                    #Update 'modes' at interior points only
                    modes[r_idx, kφ_idx, 1:-1, 0:nmodes] = eigVecs[:, 0:nmodes]
                elif (params.verticalBCs == "continuousBuoyancy" 
                      or params.verticalBCs == "constantStreamfunction"):
                    #Update 'modes' at all z-points
                    modes[r_idx, kφ_idx, :, 0:nmodes] = eigVecs[:, 0:nmodes]
                
        #Save results to nc file
        SaveToNetCDF(params, geom, growth, prop, modes)
    
    RunVisFromSavedData(params, geom)

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()