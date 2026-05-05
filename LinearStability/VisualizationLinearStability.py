import matplotlib.pyplot as plt
import numpy as np

from math import pi
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from netCDF4 import Dataset
from os import makedirs

#To do: functions to label units for dimensional case; normalize streamfunctions before plotting; function to run visualization immediately after eig solver (i.e., not from saved data)
        
def GetModeString(nmodes, mode):
    """
    String, representing mode, to be incorporated into filename.
    """
    if nmodes == 1:
        return "fastestgrowing"
    else:
        return f"mode{mode}"

def LoadCommonVariables(ds):
    """
    Load, from nc file, variables common to 1D- and 2D-solver outputs.
    """

    modes = ds.variables["mode"][:]
    kφs   = ds.variables["kφ"][:]
    r     = ds.variables["r"][:]
    φ     = ds.variables["φ"][:]
    
    eigModesReal     = ds.variables["eigMode"][:, :, :, :]["r"]
    eigModesImag     = ds.variables["eigMode"][:, :, :, :]["i"]
    eigStreamfnsReal = ds.variables["eigStreamfn"][:, :, :, :, :]["r"]
    eigStreamfnsImag = ds.variables["eigStreamfn"][:, :, :, :, :]["i"]
    eig_urReal       = ds.variables["eig_ur"][:, :, :, :, :]["r"]
    eig_urImag       = ds.variables["eig_ur"][:, :, :, :, :]["i"]
    eig_uφReal       = ds.variables["eig_uφ"][:, :, :, :, :]["r"]
    eig_uφImag       = ds.variables["eig_uφ"][:, :, :, :, :]["i"]
    
    return (modes, kφs, r, φ, eigModesReal, eigModesImag, eigStreamfnsReal,
            eigStreamfnsImag, eig_urReal, eig_urImag, eig_uφReal, eig_uφImag)

def LoadSavedData1D(params, geom):
    """
    Load results of 1D gen. eig. solver from nc file.
    """

    ds = Dataset(f"./Data/{params.dimString}_Lr{params.Lr:.1E}_Nr{params.Nr}_Ro{params.Ro:.1E}_BuInf_f{params.f0:.1E}.nc")
    
    #Load variables
    commonVariables = LoadCommonVariables(ds)
    kzs             = ds.variables["kz"][:]
    growthDim       = ds.variables["growth_rate"][:, :, :]
    propDim         = ds.variables["prop_speed"][:, :, :]
    
    ds.close()
    
    return commonVariables, kzs, growthDim, propDim

def LoadSavedData2D(params, geom):
    """
    Load results of 2D gen. eig. solver from nc file.
    """
    
    ds = Dataset(f"./Data/{params.dimString}_Lr{params.Lr:.1E}_Lz{params.Lz:.1E}_Nr{params.Nr}_Nz{params.Nz}_Ro{params.Ro:.1E}_Bu{params.Bu:.1E}_f{params.f0:.1E}.nc")

    #Load variables
    commonVariables = LoadCommonVariables(ds)
    z               = ds.variables["z"][:]
    growthDim       = ds.variables["growth_rate"][:, :]
    propDim         = ds.variables["prop_speed"][:, :]
    
    ds.close()
    
    return commonVariables, z, growthDim, propDim
    
def plot_sigmar_polarGrid(ax, params):
    """
    Plot indication of radial gyre length scale, if it is within domain, on 
     polar rφ-grid.
    """
    if params.sigmar < params.Lr:
        ax.plot(np.linspace(0, (2 * pi), params.Np), 
                params.sigmar * np.ones(params.Np), color = "k", ls = "--")
                    
def plot_sigmar_CartesianGrid(ax, params):
    """
    Plot indication of radial gyre length scale, if it is within domain, on 
     Cartesian grid with horizontal axis r.
    """
    if params.sigmar < params.Lr:
        ax.axvline(params.sigmar, color = "k", ls = "--")
    
def plot_sigmaz(ax, params):
    """
    Plot indication of vertical gyre length scale, if it is within domain.
    """
    if params.sigmaz < params.Lz:
        ax.axhline(-params.sigmaz, color = "k", ls = "--")

def PlotEigvals(params, nmodes, kφs, kzs, dimensionalGrowthRates,
                dimensionalPropSpeeds, setupString):
    """
    Visualize growth rates and propagation speeds for different wavenumbers.
    """
  
    for mode in range(nmodes):
    
        modeString = GetModeString(nmodes, mode)
        
        if params.discretizeVertical:
        
            dimString = params.dimString + "2D"
            xVariable = "k" 
        
            fig, (axGrowth, axProp) = plt.subplots(1, 2, figsize = (13, 5))
    
            axGrowth.grid(True)
            axGrowth.scatter(kφs, np.ravel(dimensionalGrowthRates[:, mode]),
                              color = "mediumpurple")
            axGrowth.set(title = "Growth rate", xlabel = "Azimuthal wavenumber",
                         ylabel = "Growth rate (s$^{{-1}}$)")
    
            axProp.grid(True)
            axProp.scatter(kφs, np.ravel(dimensionalPropSpeeds[:, mode]),
                           color = "mediumpurple")
            axProp.set(title = "Propagation speed", 
                       xlabel = "Azimuthal wavenumber",
                       ylabel = "Angular velocity (s$^{{-1}}$)")
            
        elif not params.discretizeVertical:
        
            dimString = params.dimString + "1D"
            xVariable = "m"

            nRows = min(len(kφs), 4)

            fig, axs = plt.subplots(nRows, 2, figsize = (13, (7 * nRows)),
                                    sharex = "col")
            
            for ii in range(len(kφs)):
                
                axGrowth = axs[ii, 0]
                axGrowth.grid(True)
                axGrowth.plot(kzs, 
                              np.ravel(dimensionalGrowthRates[:, ii, mode]),
                              ".-", color = "mediumpurple")
                axGrowth.set(title = f"Growth rate; $k_{{\phi}}$ = {kφs[ii]}",
                              ylabel = "Growth rate (s$^{{-1}}$)")

                axProp = axs[ii, 1]
                axProp.grid(True)
                axProp.plot(kzs, np.ravel(dimensionalPropSpeeds[:, ii, mode]),
                            ".-", color = "mediumpurple")
                axProp.set(title = f"Propagation speed; $k_{{\phi}}$ = {kφs[ii]}",
                           ylabel = "Angular velocity (s$^{{-1}}$)")
            
            #Set x-labels on lowest axes
            axGrowth.set(xlabel = r'Vertical wavenumber (m$^{{-1}}$)')
            axProp.set(xlabel = r'Vertical wavenumber (m$^{{-1}}$)')

        fig.savefig(f"./Graphs/omega_vs_{xVariable}_{modeString}_{dimString}gyre_{setupString}.png")
        plt.close(fig)
        
def PlotEigModeStructures(params, nmodes, kφs, kzs, r, z, eigModesReal, 
                          eigModesImag, setupString):
    """
    Visualize spatial structures of eigenmodes.
    """
    
    for mode in range(nmodes):
    
        modeString = GetModeString(nmodes, mode)
    
        for kφ_idx in range(len(kφs)):
        
            kφ = kφs[kφ_idx] #Wavenumber to plot for
    
            if params.discretizeVertical:
                
                #Normalize eigenvector components
                eigModeReal     = eigModesReal[kφ_idx, :, :, mode]
                eigModeImag     = eigModesImag[kφ_idx, :, :, mode]
                eigModeAmp      = np.sqrt(eigModeReal**2 + eigModeImag**2)
                eigModeRealNorm = eigModeReal / np.max(eigModeAmp)
                eigModeImagNorm = eigModeImag / np.max(eigModeAmp)
             
                #Reshape eigenvector to fit rz-grid
                eigModeRealNorm_rz = np.reshape(eigModeRealNorm, 
                                                (len(r), len(z)))
                eigModeImagNorm_rz = np.reshape(eigModeImagNorm,
                                                (len(r), len(z)))
             
                zMesh, rMesh = np.meshgrid(z, r)
             
                fig, axs = plt.subplots(1, 2, figsize = (12, 7), sharey = "row")
                
                for i in range(2):
                    axs[i].grid(False) #Required for pcolormesh
                    
                axs[0].pcolormesh(rMesh, zMesh, eigModeRealNorm_rz, 
                                  cmap = "RdBu_r", vmin = -1, vmax = 1)
                axs[0].set(xlabel = "$r$ (m)", ylabel = "$z$ (m)",
                           title = "Re[$\hat{\psi} (r,z)$]")
                axs[1].pcolormesh(rMesh, zMesh, eigModeImagNorm_rz,
                                  cmap = "RdBu_r", vmin = -1, vmax = 1)
                axs[1].set(xlabel = "$r$ (m)",
                           title = "Im[$\hat{\psi} (r,z)$]")
        
                for i in range(2):
                    
                    #Gyre length scales
                    plot_sigmar_CartesianGrid(axs[i], params)
                    plot_sigmaz(axs[i], params)
                    
                    axs[i].grid(True) #Restore grids for final version
                    
                fig.suptitle(f"Components of fastest-growing eigenmode in $rz$-plane for wavenumber $k_{{\phi}}= {kφ}$\n\n")
                fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, 
                                                             vmax = 1),
                                            cmap = "RdBu_r"), 
                             ax = axs.ravel().tolist(),
                             orientation = "horizontal", shrink = 0.8,
                             label = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$")
                fig.savefig(f"./Graphs/eigModeStructure_k{int(kφ)}_{modeString}_{params.dimString}2Dgyre_{setupString}.png")
                plt.close(fig)
                
            elif not params.discretizeVertical:
            
                for kz_idx in range(len(kzs)):
            
                    kz = kzs[kz_idx] #Wavenumber to plot for
                    
                    if params.nondimensional:
                        kz_units = "per $\sigma_z$"
                    else:
                        kz_units = "m$^{-1}$"
        
                    #Normalize eigenvector components
                    eigModeReal     = eigModesReal[kz_idx, kφ_idx, :, mode]
                    eigModeImag     = eigModesImag[kz_idx, kφ_idx, :, mode]
                    eigModeAmp      = np.sqrt(eigModeReal**2 + eigModeImag**2)
                    eigModeRealNorm = eigModeReal / max(eigModeAmp)
                    eigModeImagNorm = eigModeImag / max(eigModeAmp)
                          
                    fig, ax = plt.subplots(figsize = (10, 8))
        
                    ax.plot(r, eigModeRealNorm, "-", color = "mediumpurple", 
                            label = "Re[$\hat{\psi}$]")
                    ax.plot(r, eigModeImagNorm, "--", color = "mediumpurple", 
                            label = "Im[$\hat{\psi}$]")
        
                    plot_sigmar_CartesianGrid(ax, params) #Gyre length scale
                
                    ax.set(xlabel = "$r$ (m)",
                           ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
                           title = f"Components of fastest-growing eigenvector for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz} {kz_units}")
                    ax.legend()
                    fig.savefig(f"./Graphs/eigModeStructure_k{int(kφ)}_m{kz:.4E}_{modeString}_{params.dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)
                    
def PlotStreamfnsAndVelocities(params, geom, nmodes, kφs, kzs, r, φ, z, 
                               eigModesReal, eigModesImag,
                               eigStreamfnsReal, eigStreamfnsImag, eig_urReal,
                               eig_urImag, eig_uφReal, eig_uφImag, setupString):
    """
    Visualize spatial structures of eigen-streamfunctions and velocities.
    """
    
    for mode in range(nmodes):
    
        modeString = GetModeString(nmodes, mode)
        
        for kφ_idx in range(len(kφs)):
        
            kφ = kφs[kφ_idx] #Wavenumber to plot for
            
            if params.discretizeVertical:
            
                eigStreamfnReal = eigStreamfnsReal[kφ_idx, :, :, :, mode]
                eigStreamfnImag = eigStreamfnsImag[kφ_idx, :, :, :, mode]
                
                z_idx = 1
                
                φMesh, rMesh = np.meshgrid(φ, r)
                
                #Plot eigen-streamfunction on constant-z surface
                
                fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                        subplot_kw = {"projection": "polar"})
    
                for i in range(2):
                    axs[i].grid(False) #Required for pcolormesh

                axs[0].pcolormesh(φMesh, rMesh, eigStreamfnReal[:, :, z_idx],
                                  cmap = "RdBu_r", vmin = -1, vmax = 1)
                axs[0].set(title = f"Re[$\hat{{\psi}}$ exp($ik\phi$)]")
                axs[1].pcolormesh(φMesh, rMesh, eigStreamfnImag[:, :, z_idx],
                                  cmap = "RdBu_r", vmin = -1, vmax = 1)
                axs[1].set(title = f"Im[$\hat{{\psi}}$ exp($ik\phi$)]")
    
                for i in range(2):
                    plot_sigmar_polarGrid(axs[i], params) #Gyre length scale
                    axs[i].grid(True) #Restore grids for final version of plot
        
                fig.subplots_adjust(hspace = 0.75, wspace = 0.5)
                fig.suptitle(f"Components of fastest-growing eigen-streamfunction for $k_{{\phi}}$ = {kφ} in plane $z=$ {z[z_idx]:.0f} m\n\n\n")
                fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, 
                                                             vmax = 1),
                                            cmap = "RdBu_r"), 
                             ax = axs.ravel().tolist(), 
                             orientation = "horizontal", shrink = 0.8)
                fig.savefig(f"./Graphs/streamfn_z{z[z_idx]:.0f}_k{int(kφ)}_{modeString}_{params.dimString}2Dgyre_{setupString}.png")
                plt.close(fig)
                
                #Eigen-velocity components on constant-z surface
                urReal_z0 = eig_urReal[kφ_idx, :, :, z_idx, mode]
                urImag_z0 = eig_urImag[kφ_idx, :, :, z_idx, mode]
                uφReal_z0 = eig_uφReal[kφ_idx, :, :, z_idx, mode]
                uφImag_z0 = eig_uφImag[kφ_idx, :, :, z_idx, mode]
                
                #Abs. max. of eigen-velocity components on constant-z surface
                urMax_z0 = np.max(np.abs(np.sqrt(urReal_z0**2 + urImag_z0**2)))
                uφMax_z0 = np.max(np.abs(np.sqrt(uφReal_z0**2 + uφImag_z0**2)))
        
                #Plot eigen-velocity on constant-z surface
                
                fig, axs = plt.subplots(2, 2, figsize = (8, 8),
                                        subplot_kw = {"projection": "polar"})

                for i in range(2):
                    for j in range(2):
                        axs[i, j].grid(False) #Required for pcolormesh
                            
                pcm_ur = axs[0, 0].pcolormesh(φMesh, rMesh, urReal_z0,
                                      cmap = "RdBu_r", vmin = -urMax_z0, 
                                      vmax = urMax_z0)
                axs[0, 0].set_title(f"Re[$u_r'$]")
                axs[0, 1].pcolormesh(φMesh, rMesh, urImag_z0,
                                         cmap = "RdBu_r", vmin = -urMax_z0, 
                                         vmax = urMax_z0)
                axs[0, 1].set_title(f"Im[$u_r'$]")
                    
                pcm_uφ = axs[1, 0].pcolormesh(φMesh, rMesh, uφReal_z0,
                                                  cmap = "RdBu_r", 
                                                  vmin = -uφMax_z0, vmax = uφMax_z0)
                axs[1, 0].set_title(f"Re[$u_{{\phi}}'$]")
                axs[1, 1].pcolormesh(φMesh, rMesh, uφImag_z0,
                                         cmap = "RdBu_r", vmin = -uφMax_z0, 
                                         vmax = uφMax_z0)
                axs[1, 1].set_title(f"Im[$u_{{\phi}}'$]")
                    
                for i in range(2):
                    for j in range(2):
                        
                        #Gyre length scale
                        plot_sigmar_polarGrid(axs[i, j], params)
                        
                        axs[i, j].grid(True) #Restore grids for final version

                fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
                fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction \n in plane $z =$ {z[z_idx]:.1E} for $k_{{\phi}} =$ {kφ}")
                fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], 
                                 location = "right", shrink = 0.6,
                                 label = "m/s", pad = 0.1)
                fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], 
                                 location = "right", shrink = 0.6,
                                 label = "m/s", pad = 0.1)
                fig.savefig(f"./Graphs/velocities_z{z[z_idx]:.0f}_k{int(kφ)}_{modeString}_{params.dimString}2Dgyre_{setupString}.png")
                plt.close(fig)
                
                #Plot eigen-streamfunction on constant-r surface (r ~ sigma_r)
                
                zMesh, φMesh = np.meshgrid(z, φ)
                
                #Get index of r closest ot sigma_r
                r_idx = np.abs(r - params.sigmar).argmin()
                
                fig, axs = plt.subplots(1, 2, figsize = (11, 7), sharey = "row")
                
                for i in range(2):
                    axs[i].grid(False) #Required for pcolormesh
        
                axs[0].pcolormesh(φMesh, zMesh, eigStreamfnReal[r_idx, :, :],
                                  cmap = "RdBu_r", vmin = -1, vmax = 1)
                axs[0].set(xlabel = "$\phi$", ylabel = "$z$ (m)",
                           title = "Re[$\hat{{\psi}}$ exp($ik\phi$)]")
                axs[1].pcolormesh(φMesh, zMesh, eigStreamfnImag[r_idx, :, :],
                                  cmap = "RdBu_r", vmin = -1, vmax = 1)
                axs[1].set(xlabel = "$\phi$",
                           title = "Im[$\hat{{\psi}}$ exp($ik\phi$)]")
        
                for i in range(2):
                    plot_sigmaz(axs[i], params) #Gyre length scale
                    axs[i].grid(True) #Restore grids for final version of plot
                    
                fig.suptitle(f"Components of fastest-growing eigen-streamfunction for $k_{{\phi}}$ = {kφ} in plane $r=$ {r[r_idx]:.1E} m\n\n\n")
                fig.subplots_adjust(hspace = 0.8)
                fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1,
                                                             vmax = 1),
                                            cmap = "RdBu_r"), 
                             ax = axs.ravel().tolist(),
                             orientation = "horizontal", shrink = 0.8)
                fig.savefig(f"./Graphs/streamfn_r{r[r_idx]:.1E}_k{int(kφ)}_{modeString}_{params.dimString}2Dgyre_{setupString}.png")
                plt.close(fig)
                
                #Eigen-velocity components on constant-r surface
                urReal_r0 = eig_urReal[kφ_idx, r_idx, :, :, mode]
                urImag_r0 = eig_urImag[kφ_idx, r_idx, :, :, mode]
                uφReal_r0 = eig_uφReal[kφ_idx, r_idx, :, :, mode]
                uφImag_r0 = eig_uφImag[kφ_idx, r_idx, :, :, mode]
                
                #Abs. max. of eigen-velocity components on constant-r surface
                urMax_r0 = np.max(np.abs(np.sqrt(urReal_r0**2 + urImag_r0**2)))
                uφMax_r0 = np.max(np.abs(np.sqrt(uφReal_r0**2 + uφImag_r0**2)))
        
                #Plot eigen-velocity on constant-r surface
                
                fig, axs = plt.subplots(2, 2, figsize = (8, 8), sharey = "row",
                                        sharex = "col")

                for i in range(2):
                    for j in range(2):
                        axs[i, j].grid(False) #Required for pcolormesh
                            
                pcm_ur = axs[0, 0].pcolormesh(φMesh, zMesh, urReal_r0,
                                      cmap = "RdBu_r", vmin = -urMax_r0,
                                      vmax = urMax_r0)
                axs[0, 0].set(ylabel = "$z$ (m)", title = "Re[$u_r'$]")
                axs[0, 1].pcolormesh(φMesh, zMesh, urImag_r0, cmap = "RdBu_r",
                                     vmin = -urMax_r0, vmax = urMax_r0)
                axs[0, 1].set(title = "Im[$u_r'$]")
                    
                pcm_uφ = axs[1, 0].pcolormesh(φMesh, zMesh, uφReal_r0,
                                              cmap = "RdBu_r", vmin = -uφMax_r0,
                                              vmax = uφMax_r0)
                axs[1, 0].set(xlabel = "$\phi$", ylabel = "$z$ (m)", 
                              title = "Re[$u_{{\phi}}'$]")
                axs[1, 1].pcolormesh(φMesh, zMesh, uφImag_r0, cmap = "RdBu_r",
                                     vmin = -uφMax_r0, vmax = uφMax_r0)
                axs[1, 1].set(xlabel = "$\phi$", title = "Im[$u_{{\phi}}'$]")
                    
                for i in range(2):
                    for j in range(2):
                        plot_sigmaz(axs[i, j], params) #Gyre length scale
                        axs[i, j].grid(True) #Restore grids for final version

                fig.subplots_adjust(hspace = 0.3, wspace = 0.8)
                fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction \n on surface $r =$ {r[r_idx]:.1E} for $k_{{\phi}} =$ {kφ}")
                fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], 
                             location = "right", shrink = 0.6, label = "m/s",
                             pad = 0.1)
                fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], 
                             location = "right", shrink = 0.6, label = "m/s", 
                             pad = 0.1)
                fig.savefig(f"./Graphs/velocities_r{r[r_idx]:.0f}_k{int(kφ)}_{modeString}_{params.dimString}2Dgyre_{setupString}.png")
                plt.close(fig)
            
            elif not params.discretizeVertical:
            
                for kz_idx in range(len(kzs)):
                
                    kz = kzs[kz_idx] #Wavenumber to plot for
                    
                    if params.nondimensional:
                        kz_units = "per $\sigma_z$"
                    else:
                        kz_units = "m$^{-1}$"

                    eigStreamfnReal = eigStreamfnsReal[kz_idx, kφ_idx, :, :, 
                                                       mode]
                    eigStreamfnImag = eigStreamfnsImag[kz_idx, kφ_idx, :, :, 
                                                       mode]

                    #Plot eigen-streamfunction in rφ-plane
                    
                    φMesh, rMesh = np.meshgrid(φ, r)

                    fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                            subplot_kw = {"projection": "polar"}
                                           )
                                
                    for i in range(2):
                        axs[i].grid(False) #Required for pcolormesh

                    axs[0].pcolormesh(φMesh, rMesh, eigStreamfnReal,
                                      cmap = "RdBu_r", vmin = -1, vmax = 1)
                    axs[0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]")
                    axs[1].pcolormesh(φMesh, rMesh, eigStreamfnImag,
                                      cmap = "RdBu_r", vmin = -1, vmax = 1)
                    axs[1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]")
                    
                    for i in range(2):
                        plot_sigmar_polarGrid(axs[i], params) #Gyre length scale     
                        axs[i].grid(True) #Restore grids for final version of plot
            
                    fig.subplots_adjust(hspace = 0.5, wspace = 0.75)
                    fig.suptitle(f"Components of fastest-growing eigen-streamfunction in $r\phi$-plane\n for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz:.4E} {kz_units}\n\n")
                    fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, 
                                                                 vmax = 1),
                                                cmap = "RdBu_r"), 
                                 ax = axs.ravel().tolist(), 
                                 orientation = "horizontal", shrink = 0.8)
                    fig.savefig(f"./Graphs/streamfn_k{int(kφ)}_m{kz:.4E}_{modeString}_{dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)
                    
                    urReal = eig_urReal[kz_idx, kφ_idx, :, :, mode]
                    urImag = eig_urImag[kz_idx, kφ_idx, :, :, mode]
                    uφReal = eig_uφReal[kz_idx, kφ_idx, :, :, mode]
                    uφImag = eig_uφImag[kz_idx, kφ_idx, :, :, mode]
                    
                    #Absolute maxmimum of eigen-velocity components
                    urMax = np.max(np.abs(np.sqrt(urReal**2 + urImag**2)))
                    uφMax = np.max(np.abs(np.sqrt(uφReal**2 + uφImag**2)))
        
                    #Plot eigen-velocity in rφ-plane
                    
                    fig, axs = plt.subplots(2, 2, figsize = (8, 8),
                                            subplot_kw = {"projection": 
                                                          "polar"})

                    for i in range(2):
                        for j in range(2):
                            axs[i, j].grid(False) #Required for pcolormesh
                            
                    pcm_ur = axs[0, 0].pcolormesh(φMesh, rMesh, urReal,
                                      cmap = "RdBu_r", vmin = -urMax, 
                                      vmax = urMax)
                    axs[0, 0].set_title(f"Re[$u_r'(r, \phi)$]")
                    axs[0, 1].pcolormesh(φMesh, rMesh, urImag,
                                         cmap = "RdBu_r", vmin = -urMax, 
                                         vmax = urMax)
                    axs[0, 1].set_title(f"Im[$u_r'(r, \phi)$]")
                    
                    pcm_uφ = axs[1, 0].pcolormesh(φMesh, rMesh, uφReal,
                                                  cmap = "RdBu_r", 
                                                  vmin = -uφMax, vmax = uφMax)
                    axs[1, 0].set_title(f"Re[$u_{{\phi}}'(r,\phi)$]")
                    axs[1, 1].pcolormesh(φMesh, rMesh, uφImag,
                                         cmap = "RdBu_r", vmin = -uφMax, 
                                         vmax = uφMax)
                    axs[1, 1].set_title(f"Im[$u_{{\phi}}'(r,\phi)$]")
                    
                    for i in range(2):
                        for j in range(2):
                        
                            #Gyre length scale
                            plot_sigmar_polarGrid(axs[i, j], params)
                        
                            axs[i, j].grid(True) #Restore grids for final version

                    fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
                    fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction \n in $r\phi$-plane for wavenumbers $k_{{\phi}} =$ {kφ}, $m =$ {kz:.4E} {kz_units}")
                    fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], 
                                 location = "right", shrink = 0.6,
                                 label = "m/s", pad = 0.1)
                    fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], 
                                 location = "right", shrink = 0.6,
                                 label = "m/s", pad = 0.1)
                    fig.savefig(f"./Graphs/velocities_k{int(kφ)}_m{kz:.4E}_{modeString}_{params.dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)

def RunVisualization(params, geom, modes, kφs, kzs, r, φ, z, 
                     dimensionalGrowthRates, dimensionalPropSpeeds, 
                     eigModesReal, eigModesImag, eigStreamfnsReal, 
                     eigStreamfnsImag, eig_urReal, eig_urImag, eig_uφReal,
                     eig_uφImag, setupString):
    #N.b. if not running from saved data, will need to use r = geom.r[0:(params.halfNr + 1)]
    
    plt.rcParams.update({"text.usetex": True, "font.size": 12})

    makedirs("./Graphs", exist_ok = True) #Make folder if nonexistent
    
    PlotEigvals(params, len(modes), 
                kφs, kzs, dimensionalGrowthRates, dimensionalPropSpeeds, 
                setupString)
     
    PlotEigModeStructures(params, len(modes), kφs, kzs, r, z, eigModesReal, 
                          eigModesImag, setupString)
                          
    PlotStreamfnsAndVelocities(params, geom, len(modes), kφs, kzs, r, φ, z, 
                               eigModesReal, eigModesImag,
                               eigStreamfnsReal, eigStreamfnsImag, eig_urReal,
                               eig_urImag, eig_uφReal, eig_uφImag, setupString)
                               
#def RunVisFromEigSolver(params, geom)

def RunVisFromSavedData(params, geom):

    if params.discretizeVertical:
    
        commonVariables, z, growthDim, propDim = LoadSavedData2D(params, geom)
        kzs                                    = None
        
        setupString = f"Lr{params.Lr:.1E}_Lz{params.Lz:.1E}_Nr{params.Nr}_Nz{params.Nz}_Ro{params.Ro:.1E}_Bu{params.Bu:.1E}_f{params.f0:.1E}"
        
    elif not params.discretizeVertical:
    
        commonVariables, kzs, growthDim, propDim = LoadSavedData1D(params, geom)
        z                                        = None
        
        setupString = f"Lr{params.Lr:.1E}_Nr{params.Nr}_Ro{params.Ro:.1E}_BuInf_f{params.f0:.1E}"
        
    modes, kφs = commonVariables[0], commonVariables[1]
    r, φ       = commonVariables[2], commonVariables[3]
    
    eigModesReal, eigModesImag         = commonVariables[4], commonVariables[5]
    eigStreamfnsReal, eigStreamfnsImag = commonVariables[6], commonVariables[7]
    
    eig_urReal, eig_urImag = commonVariables[8], commonVariables[9]
    eig_uφReal, eig_uφImag = commonVariables[10], commonVariables[11]
    
    RunVisualization(params, geom, modes, kφs, kzs, r, φ, z, growthDim, propDim,
                     eigModesReal, eigModesImag, eigStreamfnsReal, 
                     eigStreamfnsImag, eig_urReal, eig_urImag, eig_uφReal, 
                     eig_uφImag, setupString)