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
from Chebyshev import Parameters, ChebyshevGeometry
from Streamfunctions import Streamfunction, EigenvelocityFrom2DEigmode

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
                    type = float, default = 3e4)
parser.add_argument("--strat_shape",
                    help = "Shape of ambient squared buoyancy frequency profile",
                    type = str, default = "constant")
parser.add_argument("-N", "--buoyancyfreq",
                    help = "Maximum buoyancy frequency (Hz)",
                    type = float, default = 1e-2)
parser.add_argument("-f0", "--Coriolis",
                    help = "Coriolis frequency f0 (Hz)",
                    type = float, default = 1.4e-4)
parser.add_argument("-U", "--bkgdU",
                    help = "Characteristic scale for background velocity (m/s)",
                    type = float, default = 3.5)
parser.add_argument("--sigmar",
                    help = "Radial length scale of gyre (m)",
                    type = float, default = 2.5e5)
parser.add_argument("--sigmaz",
                    help = "Vertical length scale of gyre (m)",
                    type = float, default = 1e20)
parser.add_argument("-Np", 
                    help = "Number of phi-points for visualization", 
                    type = int, default = 50)
parser.add_argument("--k_phi", 
                    help = "Azimuthal wavenumbers; enter as -kp start stop step",
                    type = float, default = [1, 3, 1], nargs = 3)
parser.add_argument("--nmodes", 
                    help = "Number of modes of instability to be considered",
                    type = int, default = 1)
args = vars(parser.parse_args())

def QG_Vortex_Stability():

    #Initialize parameters and set up geometry for Chebyshev solver
    params = Parameters(args, discretizeVertical = True)
    geom   = ChebyshevGeometry(params)

    #Build discrete operators
    ComputeRecips(params, geom, discretizeVertical = True)
    BuildBkgdOperators(params, geom, discretizeVertical = True)
    BuildHorizontalLaplacian(params, geom, discretizeVertical = True)
     
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
        B = BuildMatrixB(params, geom, kφ, discretizeVertical = True)
        A = BuildMatrixA(params, geom, discretizeVertical = True)

        #Find generalized eigenspace (directly)
            
        t0 = timeit.timeit()
        
        #Compute eigvals c and eigvecs psi with direct solver
        eigVals, eigVecs = spalg.eig(A, B)

        solveTime = timeit.timeit() - t0 #Time for direct solver
            
        #Indexing that sorts eigvals by ASCENDING Im(c)
        indSort = np.argsort(eigVals.imag)

        eigVals = eigVals[indSort]    #Sort eigvals
        eigVecs = eigVecs[:, indSort] #Sort eigvecs in the same order
        ωs      = eigVals * kφ        #Corresponding ω-values for this kφ
        
        #Save growth rates and propagation speeds
        growth[kφ_idx, 0:nmodes] = -ωs[0:nmodes].imag
        prop[kφ_idx, 0:nmodes]   = ωs[0:nmodes].real
        
        #Save eigenvectors at interior points (they vanish at boundary points)
        modesLen = len(modes[kφ_idx, :, 0:nmodes])
        modes[kφ_idx, 
              ((np.mod(np.arange(modesLen), (params.Nz + 1)) != 0)
               & (np.mod(np.arange(modesLen), (params.Nz + 1)) != params.Nz)
               & (np.arange(modesLen) > params.Nz)),
              0:nmodes] = eigVecs[:, 0:nmodes]
    
    def plot_sigmar_polarGrid(ax, params):
        """
        Plot indication of radial gyre length scale, if it is within domain.
        """
        if params.sigmar < params.Lr:
            ax.plot(np.linspace(0, (2 * pi), params.Np), 
                    params.sigmar * np.ones(params.Np), color = "k", ls = "--")
                    
    def plot_sigmar_CartesianGrid(ax, params):
        """
        Plot indication of radial gyre length scale, if it is within domain.
        """
        if params.sigmar < params.Lr:
            ax.axvline(params.sigmar, color = "k", ls = "--")
    
    def plot_sigmaz(ax, params):
        """
        Plot indication of vertical gyre length scale, if it is within domain.
        """
        if params.sigmaz < params.Lz:
            ax.axhline(-params.sigmaz, color = "k", ls = "--")

    #Run visualization

    plt.rcParams.update({"text.usetex": True, "font.size": 17})
    
    nkφ = (np.ravel(kφs)).shape[0]

    for mode in range(nmodes):
        
        #Visualize growth rates and propagation speeds for different kφ

        fig, axs = plt.subplots(1, 2, figsize = (13, 5))

        ax_growth = axs[0]
        ax_growth.scatter(kφs, np.ravel(growth[:, mode]),
                          color = "mediumpurple")
        ax_growth.set(title = "Growth rate", xlabel = "Azimuthal wavenumber",
                      ylabel = "Growth rate (s$^{{-1}}$)")
        ax_growth.grid(True)

        ax_prop = axs[1]
        ax_prop.scatter(kφs, np.ravel(prop[:, mode]), color = "mediumpurple")
        ax_prop.set(title = "Propagation speed", 
                    xlabel = "Azimuthal wavenumber",
                    ylabel = "Angular velocity (s$^{{-1}}$)")
        ax_prop.grid(True)

        if nmodes == 1:
            fig.savefig(f"omega_vs_k_fastestgrowing_dimensional2Dgyre.png")
        elif nmodes > 1:
            fig.savefig(f"omega_vs_k_mode{jj}_dimensional2Dgyre.png")
        
        plt.close(fig)

        #Plot spatial structures of eigenmodes
        
        for kφ_idx in range(len(kφs)):
        
            kφ = kφs[kφ_idx] #Wavenumber to plot for
    
            eigVec     = modes[kφ_idx, :, mode]
            eigVecAmp  = np.sqrt(eigVec.real**2 + eigVec.imag**2)
            eigVecNorm = eigVec / max(eigVecAmp) #Normalize eigenvector
         
            eigMode_rz = np.reshape(eigVecNorm, ((params.halfNr + 1),
                                                 (params.Nz + 1))
                                   )
         
            zVis_rz, rVis_rz = np.meshgrid(geom.z, geom.r[:(params.halfNr + 1)])
         
            fig, axs = plt.subplots(1, 2, figsize = (12, 7), sharey = "row")
            
            for i in range(2):
                axs[i].grid(False) #Required for pcolormesh
                
            axs[0].pcolormesh(rVis_rz, zVis_rz, eigMode_rz.real, 
                              cmap = "RdBu_r", vmin = -1, vmax = 1)
            axs[0].set(xlabel = "$r$ [m]", ylabel = "$z$ [m]",
                       title = "Re[$\hat{\psi} (r,z)$]")
            axs[1].pcolormesh(rVis_rz, zVis_rz, eigMode_rz.imag,
                              cmap = "RdBu_r", vmin = -1, vmax = 1)
            axs[1].set(xlabel = "$r$ [m]",
                       title = "Im[$\hat{\psi} (r,z)$]")
    
            for i in range(2):
                
                #Show gyre length scales
                plot_sigmar_CartesianGrid(axs[i], params)
                plot_sigmaz(axs[i], params)
                
                axs[i].grid(True) #Restore grids for final version
                
            fig.suptitle(f"Components of fastest-growing eigenmode in $rz$-plane for wavenumber $k_{{\phi}}= {kφ}$\n\n")
            fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                        cmap = "RdBu_r"), 
                         ax = axs.ravel().tolist(), orientation = "horizontal",
                         shrink = 0.8,
                         label = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$")
            
            if nmodes == 1:
                fig.savefig(f"eigmode_structure_k{kφ}_fastestgrowing_dimensional2Dgyre.png")
            elif nmodes > 1:
                fig.savefig(f"eigmode_structure_k{kφ}_mode{jj}_dimensional2Dgyre.png")
            
            plt.close(fig)

            #Set up to plot eigen-structures in r-φ and φ-z planes
            
            dφ               = 2 * pi / params.Np
            φCoords          = dφ * np.arange(1, (params.Np + 1))
            φVis_rφ, rVis_rφ = np.meshgrid(φCoords,
                                           geom.r[:(params.halfNr + 1)])
            zVis_φz, φVis_φz = np.meshgrid(geom.z, φCoords)
            
            #Array to hold streamfunction values
            ψ = np.zeros([(params.halfNr + 1), params.Np, (params.Nz + 1)],
                         dtype = complex)
            
            #Evaluate streamfunction at (r, φ, z)-coordinate triples  
            for z_idx in range(params.Nz + 1):
                for φ_idx in range(params.Np):
                    for r_idx in range(params.halfNr + 1):
                        ψ[r_idx, φ_idx, z_idx] = Streamfunction(
                           eigMode_rz[r_idx, z_idx], k = kφ, φ = φCoords[φ_idx])
            
            #Evaluate components of eigen-velocity
            ur, uφ = EigenvelocityFrom2DEigmode(params, geom, eigVecNorm, kφ)
    
            #Absolute maxmimum amplitudes of velocity components
            urMax = np.max(np.abs(np.sqrt(ur.real**2 + ur.imag**2)))
            uφMax = np.max(np.abs(np.sqrt(uφ.real**2 + uφ.imag**2)))
            
            #Plot streamfunction in r-φ plane
            
            z_idx_plt = 7
    
            fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                    subplot_kw = {"projection": "polar"})
    
            for i in range(2):
                axs[i].grid(False) #Required for pcolormesh
    
            axs[0].pcolormesh(φVis_rφ, rVis_rφ, ψ[:, :, z_idx_plt].real,
                              cmap = "RdBu_r", vmin = -1, vmax = 1)
            axs[0].set(title = f"Re[$\hat{{\psi}}(r,z)$ exp($ik\phi$)]")
            axs[1].pcolormesh(φVis_rφ, rVis_rφ, ψ[:, :, z_idx_plt].imag,
                              cmap = "RdBu_r", vmin = -1, vmax = 1)
            axs[1].set(title = f"Im[$\hat{{\psi}}(r,z)$ exp($ik\phi$)]")
    
            for i in range(2):
                plot_sigmar_polarGrid(axs[i], params) #Show gyre length scale
                axs[i].grid(True) #Restore grids for final version of plot
    
            fig.subplots_adjust(hspace = 0.75, wspace = 0.5)
            fig.suptitle(f"Components of fastest-growing eigen-streamfunction for $k_{{\phi}}$ = {kφ} in plane $z=$ {geom.z[z_idx_plt]:.0f} m\n\n\n")
            fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                        cmap = "RdBu_r"), 
                         ax = axs.ravel().tolist(), orientation = "horizontal",
                         shrink = 0.8)
            fig.savefig(f"streamfunc_z{geom.z[z_idx_plt]:.0f}_k{kφ}_dimensional2Dgyre.png")
            plt.close(fig)
            
            #Plot streamfunction in φ-z plane
            
            fig, axs = plt.subplots(1, 2, figsize = (11, 7), sharey = "row")
            
            r_idx_plt = params.halfNr - 2
            
            for i in range(2):
                axs[i].grid(False) #Required for pcolormesh
    
            axs[0].pcolormesh(φVis_φz, zVis_φz, ψ[r_idx_plt, :, :].real,
                              cmap = "RdBu_r", vmin = -1, vmax = 1)
            axs[0].set(xlabel = "$\phi$",
                       ylabel = "$z$ [m]",
                       title = f"Re[$\hat{{\psi}}(r,z)$ exp($ik\phi$)]")
            axs[1].pcolormesh(φVis_φz, zVis_φz, ψ[r_idx_plt, :, :].imag,
                              cmap = "RdBu_r", vmin = -1, vmax = 1)
            axs[1].set(xlabel = "$\phi$",
                       title = f"Im[$\hat{{\psi}}(r,z)$ exp($ik\phi$)]")
    
            for i in range(2):
                plot_sigmaz(axs[i], params) #Show gyre length scale
                axs[i].grid(True) #Restore grids for final version of plot
                
            fig.suptitle(f"Components of fastest-growing eigen-streamfunction for $k_{{\phi}}$ = {kφ} in plane $r=$ {geom.r[r_idx_plt]:.0f} m\n\n\n")
            fig.subplots_adjust(hspace = 0.8)
            fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                        cmap = "RdBu_r"), 
                         ax = axs.ravel().tolist(), orientation = "horizontal",
                         shrink = 0.8)
            fig.savefig(f"streamfunc_r{geom.r[r_idx_plt]:.0f}_k{kφ}_dimensional2Dgyre.png")
            plt.close(fig)
            
            #Plot eigen-velocities in r-z plane
    
            fig, axs = plt.subplots(2, 2, figsize = (11, 7), sharex = "col",
                                    sharey = "row")
    
            for i in range(2):
                for j in range(2):
                    axs[i, j].grid(False) #Required for pcolormesh
                    
            ur_rz = np.reshape(ur, ((params.halfNr + 1), (params.Nz + 1)))
            uφ_rz = np.reshape(uφ, ((params.halfNr + 1), (params.Nz + 1)))
    
            pcm_ur = axs[0, 0].pcolormesh(rVis_rz, zVis_rz, ur_rz.real,
                                          cmap = "RdBu_r", vmin = -urMax, 
                                          vmax = urMax)
            axs[0, 0].set(ylabel = "$z$ [m]",
                          title = f"Re[$u_r'(r, z)$]")
            axs[0, 1].pcolormesh(rVis_rz, zVis_rz, ur_rz.imag,
                                 cmap = "RdBu_r", vmin = -urMax, vmax = urMax)
            axs[0, 1].set(title = f"Im[$u_r'(r, z)$]")
            
            pcm_uφ = axs[1, 0].pcolormesh(rVis_rz, zVis_rz, uφ_rz.real,
                                          cmap = "RdBu_r", vmin = -uφMax, 
                                          vmax = uφMax)
            axs[1, 0].set(xlabel = "$r$ [m]",
                          ylabel = "$z$ [m]",
                          title = f"Re[$u_{{\phi}}'(r,z)$]")
            axs[1, 1].pcolormesh(rVis_rz, zVis_rz, uφ_rz.imag,
                                 cmap = "RdBu_r", vmin = -uφMax, vmax = uφMax)
            axs[1, 1].set(xlabel = "$r$ [m]",
                          title = f"Im[$u_{{\phi}}'(r,z)$]")
    
            for i in range(2):
                for j in range(2):
                
                    #Show gyre length scales
                    plot_sigmar_CartesianGrid(axs[i, j], params)
                    plot_sigmaz(axs[i, j], params)
                
                    axs[i, j].grid(True) #Restore grids for final version
    
            fig.subplots_adjust(hspace = 0.2, wspace = 0.25)
            fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction for $k_{{\phi}} =$ {kφ} in $rz$-plane\n\n\n")
            fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], location = "right",
                         shrink = 0.8, label = "m/s", pad = 0.1)
            fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], location = "right",
                         shrink = 0.8, label = "m/s", pad = 0.1)
            fig.savefig(f"eigvelocities_k{kφ}_dimensional2Dgyre.png")
            plt.close(fig)
        
if __name__ == '__main__': #For testing
   QG_Vortex_Stability()