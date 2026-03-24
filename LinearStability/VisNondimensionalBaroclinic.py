#import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

#matplotlib.use("Agg")

from math import pi
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from netCDF4 import Dataset

from Streamfunctions import Streamfunction

def RunVisNondimensionalBaroclinic(Lr, Lz, Nr, Nz, Ro, Bu, f0, r_idx, z_idx, rCoords, zCoords):

    ###################
    # LOAD SAVED DATA #
    ###################

    ds = Dataset(f"./Data/data_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.nc")

    modes = ds.variables["mode"][:]
    kφs   = ds.variables["kφ"][:]
    #r     = ds.variables["r"][:]
    #r /= len(r)#this is because i saved them with the wrong scale factor; need to fix
    φ     = ds.variables["φ"][:]
    #z     = ds.variables["z"][:]
    #z /= len(z)#likewise
    
    r, z = rCoords, zCoords
    #Need to fix the way r and z are saved (im saving linearly spaced grids, not Chebyshev)

    growth_rates = ds.variables["growth_rate"][:, :]
    prop_speeds  = ds.variables["prop_speed"][:, :]

    eigVecsReal = ds.variables["eigVec"][:, :, :, :]["r"]
    eigVecsImag = ds.variables["eigVec"][:, :, :, :]["i"]

    eigStreamfns = ds.variables["eigStreamfn"][:, :, :, :]

    ds.close()

    ##################
    # VISUALIZE DATA #
    ##################

    plt.rcParams.update({"text.usetex": True, "font.size": 12})

    os.makedirs("./Graphs", exist_ok = True) #Make folder if nonexistent

    #Dimensionalize parts of eigenvalues for visualization
    growthDim = growth_rates * f0 * Ro
    propDim   = prop_speeds * f0 * Ro

    nmodes, nkφ = len(modes), len(kφs)
    
    if nmodes == 1:

        ##################################################
        # GROWTH RATES AND PROP. SPEEDS FOR DIFFERENT kφ #
        ##################################################

        fig, axs = plt.subplots(1, 2, figsize = (12, 5), 
                                constrained_layout = True)

        ax_growth = axs[0]
        ax_growth.scatter(kφs, growthDim[:, 0], color = "mediumpurple")
        ax_growth.set(title = "Growth rate of fastest-growing mode for each $k_{\phi}$",
                      ylabel = "Growth rate (s$^{-1}$)")
        
        ax_prop = axs[1]
        ax_prop.scatter(kφs, propDim[:, 0], color = "mediumpurple")
        ax_prop.set(title = "Propagation speed of fastest-growing mode for each $k_{\phi}$",
                    ylabel = "Azimuthal speed ($s^{-1}$)")

        ax_growth.set(xlabel = "Azimuthal wavenumber $k_{\phi}$")
        ax_prop.set(xlabel = "Azimuthal wavenumber $k_{\phi}$")
        #plt.show()
        fig.savefig(f"./Graphs/omega_vs_k_fastestgrowing_nondimBCgyre_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.png")
        plt.close(fig)
        
        #########################################
        # PLOT SPATIAL STRUCTURES OF EIGENMODES #
        #########################################

        #Create meshgrid of r-φ-z coordinates
        rVis, φVis, zVis = np.meshgrid(r, φ, z, indexing = "ij")

        for kφ_idx in range(nkφ):

            kφ = kφs[kφ_idx] #Wavenumber to plot for

            #Set up figures
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
            
            #Normalize real and imaginary parts of eigenvector
            eigVecReal     = eigVecsReal[kφ_idx, :, :, 0]
            eigVecImag     = eigVecsImag[kφ_idx, :, :, 0]
            eigVecMaxAmp   = np.max(np.sqrt(np.ravel(eigVecReal)**2 + np.ravel(eigVecImag)**2))
            eigVecRealNorm = eigVecReal / eigVecMaxAmp
            eigVecImagNorm = eigVecImag / eigVecMaxAmp

            axsEigVec[0].pcolormesh(rVis[:, 0, :],
                    zVis[:, 0, :], eigVecRealNorm,
                    cmap = "RdBu_r", vmin = -1, vmax = 1)
            axsEigVec[0].set(xlabel = "$\\tilde{{r}}$",
                                 ylabel = "$\\tilde{{z}}$",
                                 title = "Real part")
            axsEigVec[1].pcolormesh(rVis[:, 0, :],
                    zVis[:, 0, :], eigVecImagNorm,
                    cmap = "RdBu_r", vmin = -1, vmax = 1)
            axsEigVec[1].set(xlabel = "$\\tilde{{r}}$",
                                 title = "Imaginary part")

            figEigVec.colorbar(ScalarMappable(norm =
                        Normalize(vmin = -1, vmax = 1), cmap = "RdBu_r"),
                                   ax = axsEigVec.ravel().tolist(),
                                   orientation = "horizontal", shrink = 0.8)
            figEigVec.suptitle(f"Components of normalized eigenvectors for fastest-growing mode with $k_{{\phi}} =$ {kφ}")
            #plt.show()
            figEigVec.savefig(f"./Graphs/eigvec_2D_fastestgrowingk{kφ}_nondimBCgyre_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.png")
            plt.close(figEigVec)

            ##################################
            # PLOT THE EIGEN-STREAMFUNCTIONS #
            ##################################

            #Array to hold streamfunction values
            ψ = np.zeros([len(r), len(φ), len(z)], dtype = complex)

            #Evaluate streamfunction at gridpoints
            for φ_idx in range(len(φ)):
                ψ[:, φ_idx, :] = Streamfunction(eigVecRealNorm + 1j*eigVecImagNorm,
                                                    k = kφ,
                                                    φ = φ[φ_idx])

            #Max. magnitude of streamfunction on rφ-slice
            ψ_rφ_max = np.max(np.abs(ψ[:, :, z_idx]))

            axsψ_rφ[0].pcolormesh(φVis[:, :, z_idx], rVis[:, :, z_idx],
                    ψ[:, :, z_idx].real, cmap = "RdBu_r",
                                      vmin = -ψ_rφ_max, vmax = ψ_rφ_max)
            axsψ_rφ[0].set_title("Re[$\\hat{{\psi}}$]")
            axsψ_rφ[1].pcolormesh(φVis[:, :, z_idx], rVis[:, :, z_idx],
                    ψ[:, :, z_idx].imag, cmap = "RdBu_r",
                                      vmin = -ψ_rφ_max, vmax = ψ_rφ_max)
            axsψ_rφ[1].set_title("Im[$\\hat{{\psi}}$]")

            #Max. magnitude of streamfunction on φz-slice
            ψ_φz_max = np.max(np.abs(ψ[r_idx, :, :]))

            axsψ_φz[0].pcolormesh(φVis[r_idx, :, :], zVis[r_idx, :, :],
                    ψ[r_idx, :, :].real, cmap = "RdBu_r",
                                      vmin = -ψ_φz_max, vmax = ψ_φz_max)
            axsψ_φz[0].set(xlabel = "$\\phi$",
                               ylabel = "$\\tilde{{z}}$",
                               title = "Re[$\\hat{{\psi}}$]")
            axsψ_φz[1].pcolormesh(φVis[r_idx, :, :], zVis[r_idx, :, :],
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
            figψ_rφ.suptitle(f"Cross section (at $\\tilde{{z}} =$ {z[z_idx]:.2f}) of fastest-growing eigen-streamfunction with $k =$ {kφ}\n")
            figψ_rφ.colorbar(ScalarMappable(norm =
                        Normalize(vmin = -ψ_rφ_max, vmax = ψ_rφ_max),
                                                cmap = "RdBu_r"),
                                 ax = axsψ_rφ.ravel().tolist(),
                                 orientation = "horizontal", shrink = 0.8)
            figψ_rφ.savefig(f"./Graphs/streamfunc_z{z[z_idx]:.2f}_fastestgrowingk{kφ}_nondimBCgyre_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.png")
            plt.close(figψ_rφ)

            figψ_φz.suptitle(f"Cross section (at $\\tilde{{r}} =$ {r[r_idx]:.2f}) of fastest-growing eigen-streamfunction with $k =$ {kφ}\n")
            figψ_φz.colorbar(ScalarMappable(norm =
                        Normalize(vmin = -ψ_φz_max, vmax = ψ_φz_max),
                                                cmap = "RdBu_r"),
                                 ax = axsψ_φz.ravel().tolist(),
                                 orientation = "horizontal", shrink = 0.8)
            figψ_φz.savefig(f"./Graphs/streamfunc_r{r[r_idx]:.2f}_fastestgrowingk{kφ}_nondimBCgyre_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.png")
            plt.close(figψ_φz)

### Testing

from Chebyshev import Chebyshev
Nr, Nz = 21, 2

def zTransform():
    return lambda z : (z - 1) / 2

rCoords = Chebyshev(Nr)[1][1:Nr//2+1]
zCoords = Chebyshev(Nz, xTransform=zTransform())[1][1:Nz+1]
RunVisNondimensionalBaroclinic(10, 3.3, 21, 2, 4e-3, 2.5e-3, 1e-4, 5, 0, rCoords, zCoords)
