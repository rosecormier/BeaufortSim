"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables are assumed to be dimensional.
"""

import argparse
import numpy as np
import scipy.linalg as spalg
#Note: I think we should aim to switch from matrix objects (e.g. built by sp.spdiags) to array objects, as the matrix functionality is now deprecated
import sys
import time
import timeit

from math import pi

from BuildLaplacian import BuildLaplacian
from BuildBkgdOperators import BuildBkgdOperators, GridInterior
from Chebyshev import Chebyshev
from SaveToNetCDF import SaveToNetCDF

#Parse command-line inputs

parser = argparse.ArgumentParser()
parser.add_argument('--NrEig',
                    help = 'Number (must be ODD) of r-gridpoints for eig computations',
                    type = int, default = 201)
parser.add_argument('--NzEig',
                    help = 'Number (must be EVEN) of z-gridpoints for eig computations',
                    type = int, default = 100)
#parser.add_argument('--Neigs', 
#                    help = 'Number of grid points for eigs computations',
#                    type = int, default = 1001)
parser.add_argument('-Lr', 
                    help = 'Radius of physical domain (km)',
                    type = float, default = 2.5e3)
parser.add_argument('-Lz',
                    help = 'Depth (> 0) of physical domain (km)',
                    type = float, default = 1.0)
parser.add_argument('--sigmar',
                    help = 'Horizontal decay scale of gyre (km)',
                    type = float, default = 250)
parser.add_argument('--sigmaz',
                    help = 'Vertical decay scale of gyre (km)',
                    type = float, default = 3e-1)
parser.add_argument('Umax',
                    help = 'Maximum background velocity (m/s)',
                    type = float, default = 3.5)
parser.add_argument('-f0', '--Coriolis',
                    help = 'Coriolis frequency (Hz)',
                    type = float, default = 1.4e-4)
parser.add_argument('--Nmax',
                    help = 'Maximum buoyancy frequency (Hz)',
                    type = float, default = 1e-2)
parser.add_argument('--strat_shape',
                    help = 'Shape of ambient squared buoyancy frequency profile',
                    type = str, default = "constant")
parser.add_argument('-p', '--PrintOutputs',
                    help = 'Flag to turn on display for each computation',
                    action = 'store_true')
parser.add_argument('-kp', '--k_phi', 
                    help = 'Azimuthal wavenumbers; enter as -kp start stop step',
                    type = float, default = [1, 10, 1], nargs = 3)
parser.add_argument('--modes', 
                    help = 'Number of modes of instability to be considered',
                    type = int, default = 1)
parser.add_argument('-Np',
                    help = 'Number of points for discretization of phi',
                    type = int, default = 50)
parser.add_argument('--z_idx',
                    help = 'Constant z-index to plot 2D slices at',
                    type = int, default = 0)
#parser.add_argument('--r_idx',
#                    help = 'Constant r-index to plot 2D slices at',
#                    type = int, default = args.Nr // 2 - args.Nr // (2 * args.Lr))
args = parser.parse_args()
               
class Parameters:
    
    σr   = args.sigmar * 1000 #Convert to m
    σz   = args.sigmaz * 1000 #Convert to m
    Umax = args.Umax
    f0   = args.Coriolis
    Nmax = args.Nmax
    bkgd = "BG"

    stratification_kw = args.strat_shape

    Lr     = args.Lr         #Max. r in physical space; half of computational domain
    Nr     = args.NrEig      #Number of computational gridpoints in r
    halfNr = args.NrEig // 2 #Number of physical r-gridpoints

    Nφ = args.Np #Number of azimuthal gridpoints; for visualization only

    Lz = args.Lz    #Max. depth (i.e., -min(z)) in physical domain
    Nz = args.NzEig #Number of computational gridpoints in z
   
    z_idx = args.z_idx
    r_idx = int(halfNr - halfNr // Lr)

    kφs    = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])        
    nmodes = args.modes

    printout = args.PrintOutputs

    def display(self):
        print(f"Ro = {self.Ro}")
        print(f"Bu = {self.Bu}")
        print(f"f0 = {self.f0}")
        print(f"Lr = {self.Lr}")
        print(f"Nr = {self.Nr}")
        print(f"halfNr = {self.halfNr}")
        print(f"Np = {self.Np}")
        print(f"Nz = {self.Nz}")
        print(f"Lz = {self.Lz}")
        print(f"kps = {self.kps}")
        print(f"nmodes = {self.nmodes}")

def zTransform():
    return lambda z : (z - 1) / 2

class Geometry:

    def __init__(self, method, params):
        
        self.method = method
        
        if method == "Chebyshev":
           
            #Compute r- and z-diff. matrices and Chebyshev-spaced grids
            Dr, r = Chebyshev(params.Nr)
            Dz, z = Chebyshev(params.Nz, xTransform = zTransform())
            
            #Scale gridpoints and variables of differentiation to fit domain
            self.r, self.Dr = r * params.Lr, Dr / params.Lr
            self.z, self.Dz = z * params.Lz, Dz / (params.Lz)
 
            self.Dr2 = np.matmul(self.Dr, self.Dr) #Second-order diff. matrix

def QG_Vortex_Stability():
 
    #Initialize parameters and set up geometry for Chebyshev solver
    paramsCh   = Parameters()
    geomCh     = Geometry("Chebyshev", paramsCh)
    geomCh.Lap = BuildLaplacian(paramsCh, geomCh, discretizeVertical = True)

    #Perform 2D discretization of grids
    GridInterior(paramsCh, geomCh, discretizeVertical = True, 
                 dimensional_σz = paramsCh.σz)

    #Discretize background-state-flow operators
    geomCh.Ψ_op, geomCh.Q_op = BuildBkgdOperators(paramsCh, geomCh, 
                                                  discretizeVertical = True,
                                                  dimensional_σz = paramsCh.σz)

    #Information about wavenumbers and modes
    kφs, nmodes = paramsCh.kφs, paramsCh.nmodes

    #Initialize arrays to store results of eigen-computations
    growthCh = np.zeros([kφs.shape[0], nmodes])
    propCh   = np.zeros([kφs.shape[0], nmodes])
    modesCh  = np.zeros([kφs.shape[0], (paramsCh.halfNr * paramsCh.Nz), 
                         nmodes], dtype = complex)

    ########################################
    # SOLVE GENERALIZED EIGENVALUE PROBLEM #
    ########################################
    
    for kφ_idx in range(0, kφs.shape[0]):

        kφ  = kφs[kφ_idx]
        kφ2 = kφ**2
    
        ##############################
        # BUILD MATRICES 'A' and 'B' #
        ##############################

        B_Ch = (geomCh.Lap - (kφ2 * geomCh.r2Recip)
                + ((1 / paramsCh.Bu) * (geomCh.Dz 
                    * (geomCh.N2Recip * geomCh.Dz)))
               )
        A_Ch = np.matmul(geomCh.Ψ_op, B_Ch) - geomCh.Q_op

        ############################
        # FIND EIGENSPACE DIRECTLY #
        ############################
            
        t0Ch = timeit.timeit()
        
        #Compute eigvals c and eigvecs psi with direct solver
        eigValsCh, eigVecsCh = spalg.eig(A_Ch, B_Ch)

        solveTimeCh = timeit.timeit() - t0Ch #Time for direct solver

        #Indexing that sorts eigvals by ASCENDING Im(c)
        indSortCh = np.argsort(eigValsCh.imag)
            
        eigValsCh = eigValsCh[indSortCh]    #Sort eigvals
        eigVecsCh = eigVecsCh[:, indSortCh] #Sort eigvecs in the same order
        ωsCh      = eigValsCh * kφ          #Corresponding omegas for this kφ
        
        print(eigValsCh[0:nmodes])

        #Store results of direct solve
        growthCh[kφ_idx, :]   = -ωsCh[0:nmodes].imag
        propCh[kφ_idx, :]     = ωsCh[0:nmodes].real
        modesCh[kφ_idx, :, :] = eigVecsCh[:, 0:nmodes]

    SaveToNetCDF(paramsCh, np.arange(0, nmodes), growthCh, propCh, modesCh, 
                 paramsCh.Lr * np.arange(paramsCh.halfNr, 0, -1), 
                 (2 * pi / paramsCh.Nφ) * np.arange(1, (paramsCh.Nφ + 1)), 
                 paramsCh.Lz * np.arange(0, paramsCh.Nz))
"""
    #################
    # VISUALIZATION #
    #################

    plt.rcParams.update({"text.usetex": True, "font.size": 12})

    nkφ = (np.ravel(kφs)).shape[0]

    if nmodes == 1:

        ############################################################
        # VISUALIZE GROWTH RATES AND PROP. SPEEDS FOR DIFFERENT kφ #
        ############################################################
        
        fig, axs = plt.subplots(nmodes, 2, figsize = (10, 5), constrained_layout = True)#, sharex = "col")
            
        #for i in range(nmodes):
                
        ax_growth = axs[0]
        ax_growth.scatter(kφs, np.ravel(growthDimCh[:, 0]), 
                           color = "mediumpurple", label = "Cheb solver")
        ax_growth.set(title = f"Growth rate; fastest-growing mode",
                          ylabel = "Growth rate ($s^{-1}$)")
                
        ax_prop = axs[1]
        ax_prop.scatter(kφs, np.ravel(propDimCh[:, 0]), 
                         color = "mediumpurple", label = "Cheb solver")
        ax_prop.set(title = f"Propagation speed; fastest-growing mode",
                        ylabel = "Azimuthal speed ($s^{-1}$)")

        ax_growth.set(xlabel = r"Azimuthal wavenumber $k_{\phi}$")
        ax_prop.set(xlabel = r"Azimuthal wavenumber $k_{\phi}$")
        axs[0].legend()
        #plt.show()
        fig.savefig(f"omega_vs_k_fastestgrowing_nondimBCgyre.png")
        plt.close(fig)
    
        ##########################################
        # PLOT SPATIAL STRUCTURES OF EIGEN-MODES #
        ##########################################

        dφ      = 2 * pi / paramsCh.Nφ
        φCoords = dφ * np.arange(1, (paramsCh.Nφ + 1))

        #Create meshgrid of r-φ-z coordinates for visualization
        rVisCh, φVisCh, zVisCh = np.meshgrid(geomCh.r[1:(paramsCh.halfNr + 1)],
                                             φCoords,
                                             geomCh.z[1:(paramsCh.Nz + 1)],
                                             indexing = "ij")

        #for jj in range(0, nmodes):
        
        for kφ_idx in range(nkφ):

            kφ = kφs[kφ_idx] #Wavenumber to plot for

            figEigVec, axsEigVec = plt.subplots(1, 2, figsize = (9, 5), 
                                                    sharey = "row")
            figψ_rφ, axsψ_rφ = plt.subplots(1, 2, figsize = (9, 5),
                                        subplot_kw = {"projection": "polar"})
            figψ_φz, axsψ_φz = plt.subplots(1, 2, figsize = (9, 5),
                                                sharey = "row")

                #figu_rφ, axsu_rφ = plt.subplots(2, 2, figsize = (8, 9),
                #                        subplot_kw = {"projection": "polar"})
                #figu_φz, axsu_φz = plt.subplots(2, 2, figsize = (9, 9), 
                #                              sharex = "col", sharey = "row")

            for i in range(2): #Remove grids (required for pcolormesh)
                axsEigVec[i].grid(False)
                axsψ_rφ[i].grid(False)
                axsψ_φz[i].grid(False)

                ####################################
                # PLOT THE EIGENVECTORS THEMSELVES #
                ####################################

                #Reshape and normalize eigenvector
            eigVecCh_C     = np.reshape(modesCh[kφ_idx, :, 0],
                                        (paramsCh.halfNr, paramsCh.Nz),
                                          order = "C")

            eigVecAmpCh    = np.sqrt(eigVecCh_C.real**2 + eigVecCh_C.imag**2)
            eigVecNormCh_C = eigVecCh_C / np.max(eigVecAmpCh)

            axsEigVec[0].pcolormesh(rVisCh[:, 0, :], 
                    zVisCh[:, 0, :], eigVecNormCh_C.real, 
                    cmap = "RdBu_r", vmin = -1, vmax = 1)
            axsEigVec[0].set(xlabel = "$\\tilde{{r}}$", 
                                 ylabel = "$\\tilde{{z}}$",
                                 title = "Real part")
            axsEigVec[1].pcolormesh(rVisCh[:, 0, :], 
                    zVisCh[:, 0, :], eigVecNormCh_C.imag, 
                    cmap = "RdBu_r", vmin = -1, vmax = 1)
            axsEigVec[1].set(xlabel = "$\\tilde{{r}}$",
                                 title = "Imaginary part")

            figEigVec.colorbar(ScalarMappable(norm =
                        Normalize(vmin = -1, vmax = 1), cmap = "RdBu_r"),
                                   ax = axsEigVec.ravel().tolist(),
                                   orientation = "horizontal", shrink = 0.8)
            figEigVec.suptitle(f"Components of normalized eigenvectors for fastest-growing mode with $k =$ {kφ}")
                #plt.show()
            figEigVec.savefig(f"eigvec_2D_fastestgrowingk{kφ}_nondimBCgyre.png")
            plt.close(figEigVec)

                ##################################
                # PLOT THE EIGEN-STREAMFUNCTIONS #
                ##################################

            r_idx, z_idx = paramsCh.r_idx, paramsCh.z_idx

                #Need to reshape and re-normalize eigenvector
                #Ideally, will figure out a way to combine this with previous
            eigVecCh_F     = np.reshape(modesCh[kφ_idx, :, 0],
                                            (paramsCh.halfNr, paramsCh.Nz),
                                            order = "F")
            eigVecNormCh_F = eigVecCh_F / np.max(eigVecAmpCh)

                #Array to hold streamfunction values
            ψ = np.zeros([paramsCh.halfNr, paramsCh.Nφ, paramsCh.Nz],
                             dtype = complex)

                #Evaluate streamfunction at gridpoints
            for φ_idx in range(paramsCh.Nφ):
                    ψ[:, φ_idx, :] = Streamfunction(eigVecNormCh_F,
                                                    k = kφ,
                                                    φ = φCoords[φ_idx])

                #Max. magnitude of streamfunction on rφ-slice
            ψ_rφ_max = np.max(np.abs(ψ[:, :, z_idx]))

            axsψ_rφ[0].pcolormesh(φVisCh[:, :, z_idx], rVisCh[:, :, z_idx],
                                      ψ[:, :, z_idx].real, cmap = "RdBu_r",
                                      vmin = -ψ_rφ_max, vmax = ψ_rφ_max)
            axsψ_rφ[0].set_title("Re[$\\hat{{\psi}}$]")
            axsψ_rφ[1].pcolormesh(φVisCh[:, :, z_idx], rVisCh[:, :, z_idx],
                                      ψ[:, :, z_idx].imag, cmap = "RdBu_r",
                                      vmin = -ψ_rφ_max, vmax = ψ_rφ_max)
            axsψ_rφ[1].set_title("Im[$\\hat{{\psi}}$]")

                #Max. magnitude of streamfunction on φz-slice
            ψ_φz_max = np.max(np.abs(ψ[r_idx, :, :]))

            axsψ_φz[0].pcolormesh(φVisCh[r_idx, :, :], zVisCh[r_idx, :, :],
                                      ψ[-1, :, :].real, cmap = "RdBu_r",
                                      vmin = -ψ_φz_max, vmax = ψ_φz_max)
            axsψ_φz[0].set(xlabel = "$\\phi$",
                               ylabel = "$\\tilde{{z}}$",
                               title = "Re[$\\hat{{\psi}}$]")
            axsψ_φz[1].pcolormesh(φVisCh[r_idx, :, :], zVisCh[r_idx, :, :],
                                      ψ[r_idx, :, :].imag, cmap = "RdBu_r",
                                      vmin = -ψ_φz_max, vmax = ψ_φz_max)
            axsψ_φz[1].set(xlabel = "$\\phi$", 
                               title = "Im[$\\hat{{\psi}}$]")

                #############################
                # PLOT THE EIGEN-VELOCITIES #
                #############################

                

                ###########################
                # FORMAT AND SAVE FIGURES #
                ###########################

            for i in range(2): #Restore grids for final version
                axsψ_rφ[i].grid(True)
                axsψ_φz[i].grid(True)
            
            figψ_rφ.subplots_adjust(hspace = 0.5, wspace = 0.75)
            figψ_rφ.suptitle(f"Cross section (at $\\tilde{{z}} =$ {geomCh.z[z_idx]:.2f}) of fastest-growing eigen-streamfunction with $k =$ {kφ}\n")
            figψ_rφ.colorbar(ScalarMappable(norm = 
                        Normalize(vmin = -ψ_rφ_max, vmax = ψ_rφ_max), 
                                                cmap = "RdBu_r"),
                                 ax = axsψ_rφ.ravel().tolist(),
                                 orientation = "horizontal", shrink = 0.8)
                #plt.show()
            figψ_rφ.savefig(f"streamfunc_z{geomCh.z[z_idx]:.2f}_fastestgrowingk{kφ}_nondimBCgyre.png")
            plt.close(figψ_rφ)

            figψ_φz.suptitle(f"Cross section (at $\\tilde{{r}} =$ {geomCh.r[r_idx]:.2f}) of fastest-growing eigen-streamfunction with $k =$ {kφ}\n")
            figψ_φz.colorbar(ScalarMappable(norm = 
                        Normalize(vmin = -ψ_φz_max, vmax = ψ_φz_max), 
                                                cmap = "RdBu_r"),
                                 ax = axsψ_φz.ravel().tolist(),
                                 orientation = "horizontal", shrink = 0.8)
                #plt.show()
            figψ_φz.savefig(f"streamfunc_r{geomCh.r[r_idx]:.2f}_fastestgrowingk{kφ}_nondimBCgyre.png")
            plt.close(figψ_φz)
"""
if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
