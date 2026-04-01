"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables should be specified in SI base units.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as spalg

import sys
import time
import timeit

from math import pi
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from BuildDiscreteOperators import *
from Chebyshev import Chebyshev
from Streamfunctions import Streamfunction, EigenvelocityFrom1DEigvec

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number (must be ODD) of grid points for eig computations',
                    type = int, default = 501)
parser.add_argument('-Lr', 
                    help = 'DIMENSIONAL radius of the physical domain (m)',
                    type = float, default = 2.5e6)
parser.add_argument('-N', '--buoyancyfreq',
                    help = 'Buoyancy frequency (Hz)',
                    type = float, default = 1e-2)
parser.add_argument('-f0', '--Coriolis',
                    help = 'Coriolis frequency f0 (Hz)',
                    type = float, default = 1.4e-4)
parser.add_argument('-U', '--bkgdU',
                    help = 'Characteristic scale for background velocity (m/s)',
                    type = float, default = 3.5)
parser.add_argument('-σr', '--sigmar',
                    help = 'Radial length scale of gyre (m)',
                    type = float, default = 2.5e5)
parser.add_argument('--bkgd',
                    help = 'Background flow to use ("GM" or "BG")',
                    type = str, default = "BG")
parser.add_argument('-Np', 
                    help = 'Number of points for discretization of phi', 
                    type = int, default = 50)
parser.add_argument('-kp', '--k_phi', 
                    help = 'Azimuthal wavenumbers; enter as -kp start stop step',
                    type = float, default = [1, 3, 1], nargs = 3)
parser.add_argument('-kz', '--k_z', 
                    help = 'DIMENSIONAL vertical wavenumbers (m^{-1}); enter as -kz start stop step',
                    type = float, default = [0, 1e-5, 2e-7], nargs = 3)
parser.add_argument('--modes', 
                    help = 'Number of modes of instability to be considered',
                    type = int, default = 1)
args = parser.parse_args()
                    
class Parameters:
    
    Nmax = args.buoyancyfreq
    f0   = args.Coriolis
    Umax = args.bkgdU
    σr   = args.sigmar
    bkgd = args.bkgd

    #Ro = Umax / (σr * f0) #Rossby number

    Lr     = args.Lr #Max. r in physical space; half of computational domain
    Nr     = args.Neig #Number (odd) of gridpoints in computational domain
    halfNr = args.Neig // 2   
    Nφ     = args.Np #Number of azimuthal gridpoints; for visualization only
    kφs    = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
    kzs    = np.arange(args.k_z[0], args.k_z[1], args.k_z[2])
    nmodes = args.modes

class Geometry:

    def __init__(self, params):
        
        self.method = "cheb"
           
        #Compute differentiation matrix and Chebyshev-spaced grid
        Dr, r = Chebyshev(params.Nr)
                        
        #Scale gridpoints and variable of differentiation to fit domain
        self.r, self.Dr = r * params.Lr, Dr / params.Lr
            
        #Second-order differentiation matrix
        self.Dr2 = np.matmul(self.Dr, self.Dr)

def QG_Vortex_Stability():

    #Initialize parameters and set up geometry for Chebyshev solver
    params = Parameters()
    geom   = Geometry(params)

    #Build discrete operators
    DiscretizeGrid(params, geom)
    BuildFullBkgdOperators(params, geom)
    BuildFullLaplacian(params, geom)
     
    #Information about wavenumbers and modes
    kφs, kzs, nmodes = params.kφs, params.kzs, params.nmodes

    #Initialize arrays to store results of eigen-computation
    growth = np.zeros([kzs.shape[0], kφs.shape[0], nmodes])
    prop   = np.zeros([kzs.shape[0], kφs.shape[0], nmodes])
    modes  = np.zeros([kzs.shape[0], kφs.shape[0], params.halfNr + 1, nmodes],
                      dtype = complex)

    #Solve generalized eigenvalue problem
    
    for kz_idx in range(0, kzs.shape[0]):

        kz  = kzs[kz_idx]
        kz2 = kz**2
  
        for kφ_idx in range(0, kφs.shape[0]):

            kφ  = kφs[kφ_idx]
            kφ2 = kφ**2
    
            print("Solving for kφ, =", kφ, ", kz =", kz)

            #Build matrices "A" and "B"
            
            B = BuildMatrixB(params, geom, kφ, kz)
            A = BuildMatrixA(params, geom)

            #Find eigenspace (directly)
            
            t0 = timeit.timeit()
        
            #Compute eigvals c and eigvecs psi with direct solver
            eigVals, eigVecs = spalg.eig(A, B)

            solveTime = timeit.timeit() - t0 #Time for direct solver
            
            #Indexing that sorts eigvals by ASCENDING Im(c)
            indSort = np.argsort(eigVals.imag)

            eigVals = eigVals[indSort]    #Sort eigvals
            eigVecs = eigVecs[:, indSort] #Sort eigvecs in the same order
            ωs      = eigVals * kφ        #Corresponding ω values for this kφ
           
            growth[kz_idx, kφ_idx, :]    = -ωs[0:nmodes].imag
            prop[kz_idx, kφ_idx, :]      = ωs[0:nmodes].real
            modes[kz_idx, kφ_idx, 1:, :] = eigVecs[1:(params.halfNr + 1),
                                                   0:nmodes]

            if (kφ == 1 and abs(kz - 2.6e-6) < 1e-7):
                print("First eigval =", eigVals[0])
                np.savetxt("BT_eigvec_k1_Nr501.out", modes[kz_idx, kφ_idx, :, 0])
            elif (kφ == 2 and kz == 0):
                print("First eigval =", eigVals[0])
                np.savetxt("BT_eigvec_k2_Nr501.out", modes[kz_idx, kφ_idx, :, 0])

    #Run visualization

    plt.rcParams.update({"text.usetex": True, "font.size": 17})
    
    nkφ, nkz = (np.ravel(kφs)).shape[0], (np.ravel(kzs)).shape[0]

    for jj in range(0, nmodes):
        
        if nkφ < 4:

            #Visualize growth rates and propagation speeds for different kφ

            fig, axes = plt.subplots(nkφ, 2, figsize = (13, 7), sharex = "col")

            for ii in range(0, nkφ):
                
                ax_growth = axes[ii, 0]
                ax_growth.plot(kzs, np.ravel(growth[:, ii, jj]), ".-",
                               color = "mediumpurple")
                ax_growth.set(title = f"Growth rate; $k_{{\phi}}$ = {kφs[ii]}",
                              ylabel = "Growth rate (s$^{{-1}}$)")
                ax_growth.grid(True)

                ax_prop = axes[ii, 1]
                ax_prop.plot(kzs, np.ravel(prop[:, ii, jj]), ".-",
                             color = "mediumpurple")
                ax_prop.set(title = 
                                f"Propagation speed; $k_{{\phi}}$ = {kφs[ii]}",
                            ylabel = "Angular velocity (s$^{{-1}}$)")
                ax_prop.grid(True)

            ax_growth.set(xlabel = r'Vertical wavenumber (m$^{{-1}}$)')
            ax_prop.set(xlabel = r'Vertical wavenumber (m$^{{-1}}$)')

            if nmodes == 1:
                fig.savefig(f"omega_vs_m_fastestgrowing_dimensionalBTgyre.png")
            elif nmodes > 1:
                fig.savefig(f"omega_vs_m_mode{jj}_dimensionalBTgyre.png")
            
            plt.close(fig)

        #Plot spatial structures of eigenmodes
        
        kz_idx, kφ_idx = 10, 0 #0, 1
        kz, kφ         = kzs[kz_idx], kφs[kφ_idx] #Wavenumbers to plot for

        eigVec     = modes[kz_idx, kφ_idx, :, jj]
        eigVecAmp  = np.sqrt(eigVec.real**2 + eigVec.imag**2)
        eigVecNorm = eigVec / max(eigVecAmp) #Normalize eigenvector
     
        fig, ax = plt.subplots(figsize = (10, 8))

        ax.plot(geom.r[0:(params.halfNr + 1)], eigVecNorm.real, "-",
                color = "mediumpurple", label = "Re[$\hat{\psi}$]")
        ax.plot(geom.r[0:(params.halfNr + 1)], eigVecNorm.imag, "--",
                color = "mediumpurple", label = "Im[$\hat{\psi}$]")
        
        ax.set(xlabel = "$r$ (m)", 
               ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
               title = f"Components of fastest-growing eigenvector for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz}")
        ax.legend()
        fig.savefig(f"eigvec_1Dstructure_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)
        
        #The following two figures are just for testing: plots of B*psi and A*psi, respectively
        
        fig, ax = plt.subplots()
        ax.plot(geom.r[1:(params.halfNr + 1)],
                np.matmul(B[:, 1:(params.halfNr + 1)], 
                          eigVecNorm[1:]).real[1:(params.halfNr + 1)],
                "-", color = "green", label = "Re[$B\hat{\psi}$]")
        ax.plot(geom.r[1:(params.halfNr + 1)], 
                np.matmul(B[:, 1:(params.halfNr + 1)],
                          eigVecNorm[1:]).imag[1:(params.halfNr + 1)],
                "--", color = "green", label = "Im[$B\hat{\psi}$]")
        ax.set(xlabel = "$r$ (m)", ylabel = "$B$ times normalized $\hat{\psi}$",
               title = f"Components of $B$ times fastest-growing eigenvector \nfor wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz}")
        ax.legend()
        fig.savefig(f"Bpsi_1Dstructure_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)
        
        fig, ax = plt.subplots()
        ax.plot(geom.r[1:(params.halfNr+1)], np.matmul(A[:, 1:(params.halfNr + 1)], eigVecNorm[1:]).real[1:(params.halfNr + 1)], "-",
                color = "blue",
                label = "Re[$A\hat{\psi}$]")
        ax.plot(geom.r[1:(params.halfNr + 1)], np.matmul(A[:, 1:(params.halfNr + 1)], eigVecNorm[1:]).imag[1:(params.halfNr + 1)],
                "--", color = "blue",
                label = "Im[$A\hat{\psi}$]")

        ax.set(xlabel = "$r$ (m)", ylabel = "$A$ times normalized $\hat{\psi}$",
               title = f"Components of $B$ times fastest-growing eigenvector \nfor wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz}")
        ax.legend()
        fig.savefig(f"Apsi_1Dstructure_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)

        #Set up to plot eigen-structures in r-φ plane
        
        dφ         = 2 * pi / params.Nφ
        φCoords    = dφ * np.arange(1, (params.Nφ + 1))
        φVis, rVis = np.meshgrid(φCoords, geom.r[0:(params.halfNr + 1)])
        
        #Array to hold streamfunction values
        ψ = np.zeros([(params.halfNr + 1), params.Nφ], dtype = complex)
        
        #Evaluate streamfunction at (r, φ)-coordinate pairs  
        for φ_idx in range(params.Nφ):
            for r_idx in range(params.halfNr + 1):
                ψ[r_idx, φ_idx] = Streamfunction(eigVecNorm[r_idx],
                                                 k = kφ, φ = φCoords[φ_idx])

        #Evaluate components of eigen-velocity
        eigVecMesh, φMesh = np.meshgrid(eigVecNorm, φCoords)
        ur, uφ            = EigenvelocityFrom1DEigvec(params, geom, eigVecNorm,
                                                      kφ, φ = φMesh)

        #Absolute maxmimum amplitudes of velocity components
        urMax = np.max(np.abs(np.sqrt(ur.real**2 + ur.imag**2)))
        uφMax = np.max(np.abs(np.sqrt(uφ.real**2 + uφ.imag**2)))

        #Plot streamfunction in r-φ plane

        fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                subplot_kw = {"projection": "polar"})

        for i in range(2):
            axs[i].grid(False) #Required for pcolormesh

        axs[0].pcolormesh(φVis, rVis, ψ.real, cmap = "RdBu_r", vmin = -1,
                          vmax = 1)
        axs[0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]")
        axs[1].pcolormesh(φVis, rVis, ψ.imag, cmap = "RdBu_r", vmin = -1,
                          vmax = 1)
        axs[1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]")

        for i in range(2):
            axs[i].grid(True) #Restore grids for final version of plot

        fig.subplots_adjust(hspace = 0.5, wspace = 0.75)
        fig.suptitle(f"Components of fastest-growing eigen-streamfunction in $r\phi$-plane\n for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz} m$^{{-1}}$\n\n")
        fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                    cmap = "RdBu_r"), 
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.8)
        fig.savefig(f"streamfunc2D_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)

        #Plot eigen-velocities in r-φ plane

        fig, axs = plt.subplots(2, 2, figsize = (8, 9),
                                subplot_kw = {"projection": "polar"})

        for i in range(2):
            for j in range(2):
                axs[i, j].grid(False) #Required for pcolormesh

        pcm_ur = axs[0, 0].pcolormesh(φVis, rVis, np.transpose(ur.real),
                                      cmap = "RdBu_r", vmin = -urMax, 
                                      vmax = urMax)
        axs[0, 0].set_title(f"Re[$u_r'(r, \phi)$]")
        axs[0, 1].pcolormesh(φVis, rVis, np.transpose(ur.imag),
                             cmap = "RdBu_r", vmin = -urMax, vmax = urMax)
        axs[0, 1].set_title(f"Im[$u_r'(r, \phi)$]")
        
        pcm_uφ = axs[1, 0].pcolormesh(φVis, rVis, np.transpose(uφ.real),
                                      cmap = "RdBu_r", vmin = -uφMax, 
                                      vmax = uφMax)
        axs[1, 0].set_title(f"Re[$u_{{\phi}}'(r,\phi)$]")
        axs[1, 1].pcolormesh(φVis, rVis, np.transpose(uφ.imag),
                             cmap = "RdBu_r", vmin = -uφMax, vmax = uφMax)
        axs[1, 1].set_title(f"Im[$u_{{\phi}}'(r,\phi)$]")

        for i in range(2):
            for j in range(2):
                axs[i, j].grid(True) #Restore grids for final version

        fig.subplots_adjust(hspace = 0.4, wspace = 0.8)
        fig.suptitle(f"Velocities derived from fastest-growing "
                     + "eigen-streamfunction \n in $r\phi$-plane "
                     + fr"for wavenumbers $k_{{\phi}} =$ {kφ}, $\tilde{{m}} =$ {kz}")
        fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], location = "right",
                     shrink = 0.6)
        fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], location = "right",
                     shrink = 0.6)
        fig.savefig(f"eigvelocities_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()