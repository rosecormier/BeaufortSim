import matplotlib.pyplot as plt
import numpy as np

from math import pi
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from netCDF4 import Dataset
from os import makedirs

NormedMappable = ScalarMappable(norm = Normalize(vmin = -1, vmax = 1), 
                                cmap = "RdBu_r")

def Normed(fReal, fImag):
    """
    Normalized components of complex-valued function.
    """
    
    fMaxAmp = np.max(np.sqrt(fReal**2 + fImag**2))
    
    return fReal / fMaxAmp, fImag / fMaxAmp

def vmax(data):
    """
    Get upper limit (i.e., vmax) of range of data.
    """
    return np.max(np.abs(data))

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
    
    ds_var_data = {"modes": ds.variables["mode"][:], 
                   "kφs": ds.variables["kφ"][:], 
                   "r": ds.variables["r"][:],
                   "φ": ds.variables["φ"][:],
                   "eigModesReal": ds.variables["eigMode"][:, :, :, :]["r"],
                   "eigModesImag": ds.variables["eigMode"][:, :, :, :]["i"],
                   "eigStreamfnsReal": None, "eigStreamfnsImag": None,
                   "eig_urReal": None, "eig_urImag": None,
                   "eig_uφReal": None, "eig_uφImag": None}
    
    #Update eigen-streamfunction and -velocity values, if data are present
    
    if "eigStreamfn" in ds.variables:
        ds_var_data["eigStreamfnsReal"] = ds.variables["eigStreamfn"][:, :, :, :, :]["r"]
        ds_var_data["eigStreamfnsImag"] = ds.variables["eigStreamfn"][:, :, :, :, :]["i"]
        
    if ("eig_ur" in ds.variables and "eig_uφ" in ds.variables):
        ds_var_data["eig_urReal"] = ds.variables["eig_ur"][:, :, :, :, :]["r"]
        ds_var_data["eig_urImag"] = ds.variables["eig_ur"][:, :, :, :, :]["i"]
        ds_var_data["eig_uφReal"] = ds.variables["eig_uφ"][:, :, :, :, :]["r"]
        ds_var_data["eig_uφImag"] = ds.variables["eig_uφ"][:, :, :, :, :]["i"]
    
    return ds_var_data

def LoadSavedData1DRadial(params, geom):
    """
    Load results of 1D (r) gen. eig. solver from nc file.
    """

    ds = Dataset(f"./Data/{params.dimString}_Lr{params.Lr:.1E}_Nr{params.Nr}_Ro{params.Ro:.1E}_BuInf_f{params.f0:.1E}.nc")

    commonVariables = LoadCommonVariables(ds)
    kzs             = ds.variables["kz"][:]
    growthDim       = ds.variables["growth_rate"][:, :, :]
    propDim         = ds.variables["prop_speed"][:, :, :]
    
    ds.close()
    
    return commonVariables, kzs, growthDim, propDim

def LoadSavedData1DVertical(params, geom):
    """
    Load results of 1D (z) gen. eig. solver from nc file.
    """

    ds = Dataset(f"./Data/{params.dimString}_Lz{params.Lz:.1E}_Nz{params.Nz}_{params.stratification_kw}Strat_Ro{params.Ro:.1E}_Bu{params.Bu:.1E}_f{params.f0:.1E}.nc")

    commonVariables = LoadCommonVariables(ds)
    z               = ds.variables["z"][:]
    growthDim       = ds.variables["growth_rate"][:, :, :]
    propDim         = ds.variables["prop_speed"][:, :, :]
    
    ds.close()
    
    return commonVariables, z, growthDim, propDim

def LoadSavedData2D(params, geom):
    """
    Load results of 2D gen. eig. solver from nc file.
    """
    
    ds = Dataset(f"./Data/{params.dimString}_Lr{params.Lr:.1E}_Lz{params.Lz:.1E}_Nr{params.Nr}_Nz{params.Nz}_Ro{params.Ro:.1E}_Bu{params.Bu:.1E}_f{params.f0:.1E}.nc")

    commonVariables = LoadCommonVariables(ds)
    z               = ds.variables["z"][:]
    growthDim       = ds.variables["growth_rate"][:, :]
    propDim         = ds.variables["prop_speed"][:, :]
    
    ds.close()
    
    return commonVariables, z, growthDim, propDim
    
def plot_rMaxVelocity_polarGrid(ax, params):
    """
    Plot indication of radial gyre length scale, if it is within computational 
     domain, on polar rφ-grid.
    """
    
    if ((params.discretizeRadial and params.sigmar < params.Lr) 
        or (not params.discretizeRadial and params.sigmar < np.max(params.rs))):
        
        ax.plot(np.linspace(0, (2 * pi), params.Np), 
                (params.sigmar / (2**0.5)) * np.ones(params.Np), color = "k",
                ls = "--", label = "Location of max. U")
                    
def plot_rMaxVelocity_CartesianGrid(ax, params):
    """
    Plot indication of radial gyre length scale, if it is within computational 
     domain and plot domain, on Cartesian grid with horizontal axis r.
    """
    
    if ((params.discretizeRadial and params.sigmar < params.Lr) 
        or (not params.discretizeRadial and np.min(params.rs) < params.sigmar 
            and params.sigmar < np.max(params.rs))):

        ax.axvline(params.sigmar / (2**0.5), color = "k", ls = "--", 
                   label = "Location of max. U")
    
def plot_sigmaz(ax, sigmaz, Lz):
    """
    Plot indication of vertical gyre length scale, if it is within domain.
    """
    if sigmaz < Lz:
        ax.axhline(-sigmaz, color = "k", ls = "--") 
        
def plot_stratification_peaks(ax, params):
    """
    Plot indications of peaks in stratifications, if they are within domain.
    """
    
    if params.stratification_kw in ["doubleTanh", "doubleTanhTWB"]:
    
        if abs(params.z_s) < params.Lz:
            
            ax.axhline(params.z_s, color = "greenyellow", ls = "--")
    
            if abs(params.z_d) < params.Lz:
                ax.axhline(params.z_d, color = "greenyellow", ls = "--")

def PlotEigvals(params, nmodes, kφs, kzs, dimensionalGrowthRates,
                dimensionalPropSpeeds, setupString):
    """
    Visualize growth rates and propagation speeds for different wavenumbers.
    """
        
    if (params.discretizeRadial and params.discretizeVertical):
    
        modeString           = f"first{nmodes}modes"
        dimString, xVariable = params.dimString + "2D", "kphi"

        fig, (axGrowth, axProp) = plt.subplots(1, 2, figsize = (13, 5))
    
        axGrowth.grid(True)
        axProp.grid(True)
        
        if nmodes == 1:
            axGrowth.scatter(kφs, np.ravel(dimensionalGrowthRates[:, 0]),
                             color = "mediumpurple")
            axProp.scatter(kφs, np.ravel(dimensionalPropSpeeds[:, 0]),
                           color = "mediumpurple")
        
        elif nmodes > 1:
            for mode in range(nmodes):
                axGrowth.scatter(kφs, np.ravel(dimensionalGrowthRates[:, mode]),
                                 label = f"Mode {mode}")
                axProp.scatter(kφs, np.ravel(dimensionalPropSpeeds[:, mode]),
                               label = f"Mode {mode}")
            axGrowth.legend()
            axProp.legend()

        axGrowth.set(title = "Growth rate", xlabel = "Azimuthal wavenumber",
                     ylabel = "Growth rate (s$^{{-1}}$)")
        axProp.set(title = "Propagation speed", 
                   xlabel = "Azimuthal wavenumber",
                   ylabel = "Angular velocity (s$^{{-1}}$)")
                   
        fig.savefig(f"./Graphs/omega_vs_{xVariable}_{modeString}_{params.dimString}gyre_{setupString}.png")
        plt.close(fig)
        
    elif (params.discretizeRadial and not params.discretizeVertical):

        for mode in range(nmodes):
    
            modeString           = GetModeString(nmodes, mode)
            dimString, xVariable = params.dimString + "1D", "kz"

            nRows = min(len(kφs), 4)

            fig, axs = plt.subplots(nRows, 2, figsize = (13, (7 * nRows)),
                                    sharex = "col")
            
            for ii in range(nRows):
                
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
            axGrowth.set(xlabel = f"Vertical wavenumber ({params.units[xVariable]})")
            axProp.set(xlabel = f"Vertical wavenumber ({params.units[xVariable]})")
            
            fig.savefig(f"./Graphs/omega_vs_{xVariable}_{modeString}_{params.dimString}gyre_{setupString}.png")
            plt.close(fig)
            
    elif (params.discretizeVertical and not params.discretizeRadial):

        for mode in range(nmodes):
        
            modeString           = GetModeString(nmodes, mode)
            dimString, xVariable = params.dimString + "1D", "r"

            nRows = min(len(kφs), 4)

            fig, axs = plt.subplots(nRows, 2, figsize = (13, (4 * nRows)),
                                    sharex = "col")
            
            for ii in range(nRows):
                
                axGrowth = axs[ii, 0]
                axGrowth.grid(True)
                axGrowth.plot(params.rs, 
                              dimensionalGrowthRates[:, ii, mode],
                              ".-", color = "mediumpurple")
                              
                plot_rMaxVelocity_CartesianGrid(axGrowth, params)
                
                axGrowth.set(title = f"Growth rate; $k_{{\phi}}$ = {kφs[ii]}",
                             ylabel = "Growth rate (s$^{{-1}}$)")

                axProp = axs[ii, 1]
                axProp.grid(True)
                axProp.plot(params.rs,
                            np.ravel(dimensionalPropSpeeds[:, ii, mode]),
                            ".-", color = "mediumpurple")
                            
                plot_rMaxVelocity_CartesianGrid(axProp, params)
                
                axProp.set(title = f"Propagation speed; $k_{{\phi}}$ = {kφs[ii]}",
                           ylabel = "Angular velocity (s$^{{-1}}$)")
            
            #Set x-labels on lowest axes
            axGrowth.set(xlabel = f"$r$ ({params.units[xVariable]})")
            axProp.set(xlabel = f"$r$ ({params.units[xVariable]})")
            
            fig.savefig(f"./Graphs/omega_vs_{xVariable}_{modeString}_{params.dimString}gyre_{setupString}.png")
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
    
            if (params.discretizeRadial and params.discretizeVertical):
                
                #Normalize eigenvector components
                eigModeReal, eigModeImag = Normed(eigModesReal[kφ_idx, :, :,
                                                               mode],
                                                  eigModesImag[kφ_idx, :, :,
                                                               mode])
             
                #Reshape eigenvector to fit rz-grid
                eigModeReal_rz = np.reshape(eigModeReal, (len(r), len(z)))
                eigModeImag_rz = np.reshape(eigModeImag, (len(r), len(z)))
             
                zMesh, rMesh = np.meshgrid(z, r)
             
                fig, axs = plt.subplots(1, 2, figsize = (12, 7), sharey = "row")
                
                for i in range(2):
                    axs[i].grid(False) #Required for pcolormesh
                    
                axs[0].pcolormesh(rMesh, zMesh, eigModeReal_rz, cmap = "RdBu_r",
                                  vmin = -1, vmax = 1)
                axs[0].set(xlabel = f"$r$ ({params.units['r']})", 
                           ylabel = f"$z$ ({params.units['z']})",
                           title = "Re[$\hat{\psi} (r,z)$]")
                axs[1].pcolormesh(rMesh, zMesh, eigModeImag_rz, cmap = "RdBu_r",
                                  vmin = -1, vmax = 1)
                axs[1].set(xlabel = f"$r$ ({params.units['r']})",
                           title = "Im[$\hat{\psi} (r,z)$]")
        
                for i in range(2):
                    
                    #Gyre length scales
                    plot_rMaxVelocity_CartesianGrid(axs[i], params)
                    plot_sigmaz(axs[i], params.sigmaz, params.Lz)
                    
                    plot_stratification_peaks(axs[i], params) #z_s and z_d
                    
                    axs[i].grid(True) #Restore grids for final version
                    
                fig.suptitle(f"Components of fastest-growing eigenmode in $rz$-plane for wavenumber $k_{{\phi}}= {kφ}$\n\n")
                fig.colorbar(NormedMappable, ax = axs.ravel().tolist(),
                             orientation = "horizontal", shrink = 0.8,
                             label = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$")
                fig.savefig(f"./Graphs/eigModeStructure_k{int(kφ)}_{modeString}_{params.dimString}2Dgyre_{setupString}.png")
                plt.close(fig)
                
            elif (params.discretizeRadial and not params.discretizeVertical):
            
                for kz_idx in range(len(kzs)):
            
                    kz = kzs[kz_idx] #Wavenumber to plot for
        
                    #Normalize eigenvector components
                    eigModeReal, eigModeImag = Normed(eigModesReal[kz_idx, 
                                                                   kφ_idx, :,
                                                                   mode],
                                                      eigModesImag[kz_idx, 
                                                                   kφ_idx, :, 
                                                                   mode])
                          
                    fig, ax = plt.subplots(figsize = (10, 8))
        
                    ax.grid(True)
                    ax.plot(r, eigModeReal, "-", color = "mediumpurple", 
                            label = "Re[$\hat{\psi}$]")
                    ax.plot(r, eigModeImag, "--", color = "mediumpurple", 
                            label = "Im[$\hat{\psi}$]")
        
                    plot_rMaxVelocity_CartesianGrid(ax, params)
                
                    ax.set(xlabel = f"$r$ ({params.units['r']})",
                           ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
                           title = f"Components of fastest-growing eigenvector for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz:.4f} {params.units['kz']}")
                    ax.legend()
                    fig.savefig(f"./Graphs/eigModeStructure_k{int(kφ)}_m{kz:.4E}_{modeString}_{params.dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)
                    
            elif (params.discretizeVertical and not params.discretizeRadial):

                for r_plot in params.rs_plot:
                
                    #Get index of point on discrete r-grid closest to r_plot
                    r_idx = np.abs(r - r_plot).argmin()
                    
                    r_int = int(r[r_idx])
                    
                    #Normalize eigenvector components
                    eigModeReal, eigModeImag = Normed(eigModesReal[kφ_idx,
                                                                   r_idx, :,
                                                                   mode],
                                                      eigModesImag[kφ_idx,
                                                                   r_idx, :,
                                                                   mode])

                    fig, ax = plt.subplots(figsize = (10, 8))
        
                    ax.grid(True)
                    ax.plot(eigModeReal, z, "-", color = "mediumpurple", 
                            label = "Re[$\hat{\psi}$]")
                    ax.plot(eigModeImag, z, "--", color = "mediumpurple", 
                            label = "Im[$\hat{\psi}$]")
        
                    plot_sigmaz(ax, params.sigmaz, params.Lz) #Gyre length scale
                    plot_stratification_peaks(ax, params) #z_s and z_d
                
                    ax.set(xlabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
                           ylabel = f"$z$ ({params.units['z']})",
                           title = f"Components of fastest-growing eigenvector for $k_{{\phi}}$ = {kφ} and $r =$ {r_int} {params.units['r']}")
                    ax.legend()
                    fig.savefig(f"./Graphs/eigModeStructure_k{int(kφ)}_r{r_int}_{modeString}_{params.dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)
                    
def PlotStreamfnsAndVelocities(params, geom, nmodes, kφs, kzs, r, φ, z, 
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
            
                φMesh, rMesh = np.meshgrid(φ, r)
            
                z_idx = 1
            
                #Normalize eigen-streamfunction
                StreamReal, StreamImag = Normed(eigStreamfnsReal[kφ_idx, :, :, 
                                                                 :, mode],
                                                eigStreamfnsImag[kφ_idx, :, :, 
                                                                 :, mode])
                
                #Plot eigen-streamfunction on constant-z surface
                
                fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                        subplot_kw = {"projection": "polar"})
    
                for i in range(2):
                    axs[i].grid(False) #Required for pcolormesh
                    
                vmaxStream = max(vmax(StreamReal[:, :, z_idx]), 
                                 vmax(StreamImag[:, :, z_idx]))

                axs[0].pcolormesh(φMesh, rMesh, StreamReal[:, :, z_idx], 
                                  cmap = "RdBu_r", vmin = -vmaxStream, 
                                  vmax = vmaxStream)
                axs[0].set(title = f"Re[$\hat{{\psi}}$ exp($ik\phi$)]")
                axs[1].pcolormesh(φMesh, rMesh, StreamImag[:, :, z_idx], 
                                  cmap = "RdBu_r", vmin = -vmaxStream, 
                                  vmax = vmaxStream)
                axs[1].set(title = f"Im[$\hat{{\psi}}$ exp($ik\phi$)]")
    
                for i in range(2):
                    plot_rMaxVelocity_polarGrid(axs[i], params)
                    axs[i].grid(True) #Restore grids for final version of plot
        
                if params.discretizeRadial:
                    dimString = params.dimString + "2D"
                elif not params.discretizeRadial:
                    dimString = params.dimString + "1D"
        
                fig.subplots_adjust(hspace = 0.75, wspace = 0.5)
                fig.suptitle(f"Components of fastest-growing eigen-streamfunction for $k_{{\phi}}$ = {kφ} in plane $z=$ {z[z_idx]:.1f} {params.units['z']}\n\n\n")
                fig.colorbar(ScalarMappable(norm = Normalize(
                  vmin = -vmaxStream, vmax = vmaxStream), cmap = "RdBu_r"),
                             ax = axs.ravel().tolist(), label = params.units["psi"],
                             orientation = "horizontal", shrink = 0.8)
                fig.savefig(f"./Graphs/streamfn_z{z[z_idx]:.1f}_k{int(kφ)}_{modeString}_{dimString}gyre_{setupString}.png")
                plt.close(fig)
                
                #Normalize eigen-velocity components
                urReal, urImag = Normed(eig_urReal[kφ_idx, :, :, :, mode],
                                        eig_urImag[kφ_idx, :, :, :, mode])
                uφReal, uφImag = Normed(eig_uφReal[kφ_idx, :, :, :, mode],
                                        eig_uφImag[kφ_idx, :, :, :, mode])

                #Plot eigen-velocity on constant-z surface
                
                fig, axs = plt.subplots(2, 2, figsize = (10, 8),
                                        subplot_kw = {"projection": "polar"})

                for i in range(2):
                    for j in range(2):
                        axs[i, j].grid(False) #Required for pcolormesh
                        
                vmax_ur = max(vmax(urReal[:, :, z_idx]), 
                              vmax(urImag[:, :, z_idx]))
                vmax_uφ = max(vmax(uφReal[:, :, z_idx]),
                              vmax(uφImag[:, :, z_idx]))
                            
                axs[0, 0].pcolormesh(φMesh, rMesh, urReal[:, :, z_idx], 
                                     cmap = "RdBu_r", vmin = -vmax_ur, 
                                     vmax = vmax_ur)
                axs[0, 0].set_title(f"Re[$u_r'$]")
                axs[0, 1].pcolormesh(φMesh, rMesh, urImag[:, :, z_idx], 
                                     cmap = "RdBu_r", vmin = -vmax_ur, 
                                     vmax = vmax_ur)
                axs[0, 1].set_title(f"Im[$u_r'$]")
                    
                axs[1, 0].pcolormesh(φMesh, rMesh, uφReal[:, :, z_idx], 
                                     cmap = "RdBu_r", vmin = -vmax_uφ, 
                                     vmax = vmax_uφ)
                axs[1, 0].set_title(f"Re[$u_{{\phi}}'$]")
                axs[1, 1].pcolormesh(φMesh, rMesh, uφImag[:, :, z_idx], 
                                     cmap = "RdBu_r", vmin = -vmax_uφ, 
                                     vmax = vmax_uφ)
                axs[1, 1].set_title(f"Im[$u_{{\phi}}'$]")
                    
                for i in range(2):
                    for j in range(2):
                        plot_rMaxVelocity_polarGrid(axs[i, j], params)
                        axs[i, j].grid(True) #Restore grids for final version

                fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
                fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction \n in plane $z =$ {z[z_idx]:.1f} {params.units['z']} and for $k_{{\phi}} =$ {kφ}")
                fig.colorbar(ScalarMappable(norm = Normalize(
                  vmin = -vmax_ur, vmax = vmax_ur), cmap = "RdBu_r"), 
                             ax = [axs[0, 0], axs[0, 1]], location = "right",
                             shrink = 0.6, label = params.units["u"], pad = 0.1)
                fig.colorbar(ScalarMappable(norm = Normalize(
                  vmin = -vmax_uφ, vmax = vmax_uφ), cmap = "RdBu_r"), 
                             ax = [axs[1, 0], axs[1, 1]], location = "right", 
                             shrink = 0.6, label = params.units["u"], pad = 0.1)
                fig.savefig(f"./Graphs/velocities_z{z[z_idx]:.1f}_k{int(kφ)}_{modeString}_{dimString}gyre_{setupString}.png")
                plt.close(fig)
                
                zMesh, φMesh = np.meshgrid(z, φ)
                
                for r_plot in params.rs_plot:
                
                    #Get index of point on discrete r-grid closest to r_plot
                    r_idx = np.abs(r - r_plot).argmin()
                    
                    r_int = int(r[r_idx])
                
                    #Plot eigen-streamfunction on constant-r surface
                    
                    fig, axs = plt.subplots(1, 2, figsize = (11, 7), 
                                            sharey = "row")
                    
                    for i in range(2):
                        axs[i].grid(False) #Required for pcolormesh
                        
                    vmaxStream = max(vmax(StreamReal[r_idx, :, :]), 
                                     vmax(StreamImag[r_idx, :, :]))
            
                    axs[0].pcolormesh(φMesh, zMesh, StreamReal[r_idx, :, :], 
                                      cmap = "RdBu_r", vmin = -vmaxStream, 
                                      vmax = vmaxStream)
                    axs[0].set(xlabel = "$\phi$",
                               ylabel = f"$z$ ({params.units['z']})",
                               title = "Re[$\hat{{\psi}}$ exp($ik\phi$)]")
                    axs[1].pcolormesh(φMesh, zMesh, StreamImag[r_idx, :, :], 
                                      cmap = "RdBu_r", vmin = -vmaxStream, 
                                      vmax = vmaxStream)
                    axs[1].set(xlabel = "$\phi$",
                               title = "Im[$\hat{{\psi}}$ exp($ik\phi$)]")
            
                    for i in range(2):
                        plot_sigmaz(axs[i], params.sigmaz, params.Lz) #Gyre length scale
                        plot_stratification_peaks(axs[i], params) #z_s and z_d
                        axs[i].grid(True) #Restore grids for final version of plot
                        
                    fig.suptitle(f"Components of fastest-growing eigen-streamfunction for $k_{{\phi}}$ = {kφ} in plane $r=$ {r_int} {params.units['r']}\n\n\n")
                    fig.subplots_adjust(hspace = 0.8)
                    fig.colorbar(ScalarMappable(norm = Normalize(
                      vmin = -vmaxStream, vmax = vmaxStream),
                                                cmap = "RdBu_r"), 
                                 ax = axs.ravel().tolist(),
                                 orientation = "horizontal", shrink = 0.8,
                                 pad = 0.1)
                    fig.savefig(f"./Graphs/streamfn_r{r_int}_k{int(kφ)}_{modeString}_{dimString}gyre_{setupString}.png")
                    plt.close(fig)
    
                    #Plot eigen-velocity on constant-r surface
                    
                    fig, axs = plt.subplots(2, 2, figsize = (11, 8), sharey = "row",
                                            sharex = "col")
    
                    for i in range(2):
                        for j in range(2):
                            axs[i, j].grid(False) #Required for pcolormesh
                            
                    vmax_ur = max(vmax(urReal[r_idx, :, :]), 
                                  vmax(urImag[r_idx, :, :]))
                    vmax_uφ = max(vmax(uφReal[r_idx, :, :]),
                                  vmax(uφImag[r_idx, :, :]))
                                
                    axs[0, 0].pcolormesh(φMesh, zMesh, urReal[r_idx, :, :], 
                                         cmap = "RdBu_r", vmin = -vmax_ur, 
                                         vmax = vmax_ur)
                    axs[0, 0].set(ylabel = f"$z$ ({params.units['z']})", 
                                  title = "Re[$u_r'$]")
                    axs[0, 1].pcolormesh(φMesh, zMesh, urImag[r_idx, :, :], 
                                         cmap = "RdBu_r", vmin = -vmax_ur, 
                                         vmax = vmax_ur)
                    axs[0, 1].set(title = "Im[$u_r'$]")
                        
                    axs[1, 0].pcolormesh(φMesh, zMesh, uφReal[r_idx, :, :], 
                                         cmap = "RdBu_r", vmin = -vmax_uφ, 
                                         vmax = vmax_uφ)
                    axs[1, 0].set(xlabel = "$\phi$",
                                  ylabel = f"$z$ ({params.units['z']})", 
                                  title = "Re[$u_{{\phi}}'$]")
                    axs[1, 1].pcolormesh(φMesh, zMesh, uφImag[r_idx, :, :], 
                                         cmap = "RdBu_r",
                                         vmin = -vmax_uφ, vmax = vmax_uφ)
                    axs[1, 1].set(xlabel = "$\phi$", title = "Im[$u_{{\phi}}'$]")
                        
                    for i in range(2):
                        for j in range(2):
                            plot_sigmaz(axs[i, j], params.sigmaz, params.Lz) #Gyre length scale
                            plot_stratification_peaks(axs[i, j], params) #z_s, z_d
                            axs[i, j].grid(True) #Restore grids for final version
    
                    fig.subplots_adjust(hspace = 0.3)
                    fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction \n on surface $r =$ {int(r[r_idx])} {params.units['r']} for $k_{{\phi}} =$ {kφ}")
                    fig.colorbar(ScalarMappable(norm = Normalize(
                      vmin = -vmax_ur, vmax = vmax_ur), cmap = "RdBu_r"), 
                                 ax = [axs[0, 0], axs[0, 1]], location = "right",
                                 shrink = 0.6, label = params.units["u"], pad = 0.1)
                    fig.colorbar(ScalarMappable(norm = Normalize(
                      vmin = -vmax_uφ, vmax = vmax_uφ), cmap = "RdBu_r"), 
                                 ax = [axs[1, 0], axs[1, 1]], location = "right", 
                                 shrink = 0.6, label = params.units["u"], pad = 0.1)
                    fig.savefig(f"./Graphs/velocities_r{int(r[r_idx])}_k{int(kφ)}_{modeString}_{dimString}gyre_{setupString}.png")
                    plt.close(fig)
            
            elif not params.discretizeVertical:
            
                for kz_idx in range(len(kzs)):
                
                    kz = kzs[kz_idx] #Wavenumber to plot for
                    
                    #Normalize eigen-streamfunction
                    StreamReal, StreamImag = Normed(eigStreamfnsReal[kz_idx,
                                                                      kφ_idx,
                                                                     :, :, 
                                                                     mode],
                                                    eigStreamfnsImag[kz_idx, 
                                                                     kφ_idx, 
                                                                     :, :,
                                                                     mode]
                                                   ) 

                    #Plot eigen-streamfunction in rφ-plane
                    
                    φMesh, rMesh = np.meshgrid(φ, r)

                    fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                            subplot_kw = {"projection": "polar"}
                                           )
                                
                    for i in range(2):
                        axs[i].grid(False) #Required for pcolormesh
                        
                    vmaxStream = max(vmax(StreamReal), vmax(StreamImag))

                    axs[0].pcolormesh(φMesh, rMesh, StreamReal, cmap = "RdBu_r",
                                      vmin = -vmaxStream, vmax = vmaxStream)
                    axs[0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]")
                    axs[1].pcolormesh(φMesh, rMesh, StreamImag, cmap = "RdBu_r",
                                      vmin = -vmaxStream, vmax = vmaxStream)
                    axs[1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]")
                    
                    for i in range(2):
                        plot_rMaxVelocity_polarGrid(axs[i], params)
                        axs[i].grid(True) #Restore grids for final version of plot
            
                    fig.subplots_adjust(hspace = 0.5, wspace = 0.75)
                    fig.suptitle(f"Components of fastest-growing eigen-streamfunction in $r\phi$-plane\n for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz:.4f} {params.units['kz']}\n\n")
                    fig.colorbar(ScalarMappable(norm = Normalize(
                      vmin = -vmaxStream, vmax = vmaxStream), cmap = "RdBu_r"), 
                                 ax = axs.ravel().tolist(), 
                                 orientation = "horizontal", shrink = 0.8)
                    fig.savefig(f"./Graphs/streamfn_k{int(kφ)}_m{kz:.4E}_{modeString}_{params.dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)
                    
                    #Normalize eigen-velocity components
                    urReal, urImag = Normed(eig_urReal[kz_idx, kφ_idx, :, :, 
                                                       mode],
                                            eig_urImag[kz_idx, kφ_idx, :, :, 
                                                       mode]
                                           )
                    uφReal, uφImag = Normed(eig_uφReal[kz_idx, kφ_idx, :, :, 
                                                       mode],
                                            eig_uφImag[kz_idx, kφ_idx, :, :, 
                                                       mode]
                                           )

                    #Plot eigen-velocity in rφ-plane
                    
                    fig, axs = plt.subplots(2, 2, figsize = (8, 8),
                                            subplot_kw = {"projection": 
                                                          "polar"})

                    for i in range(2):
                        for j in range(2):
                            axs[i, j].grid(False) #Required for pcolormesh
                            
                    vmax_ur = max(vmax(urReal), vmax(urImag))
                    vmax_uφ = max(vmax(uφReal), vmax(uφImag))
                            
                    pcm_ur = axs[0, 0].pcolormesh(φMesh, rMesh, urReal,
                                                  cmap = "RdBu_r", 
                                                  vmin = -vmax_ur, 
                                                  vmax = vmax_ur)
                    axs[0, 0].set_title(f"Re[$u_r'(r, \phi)$]")
                    axs[0, 1].pcolormesh(φMesh, rMesh, urImag, cmap = "RdBu_r",
                                         vmin = -vmax_ur, vmax = vmax_ur)
                    axs[0, 1].set_title(f"Im[$u_r'(r, \phi)$]")
                    
                    pcm_uφ = axs[1, 0].pcolormesh(φMesh, rMesh, uφReal,
                                                  cmap = "RdBu_r", 
                                                  vmin = -vmax_uφ, 
                                                  vmax = vmax_uφ)
                    axs[1, 0].set_title(f"Re[$u_{{\phi}}'(r,\phi)$]")
                    axs[1, 1].pcolormesh(φMesh, rMesh, uφImag, cmap = "RdBu_r", 
                                         vmin = -vmax_uφ, vmax = vmax_uφ)
                    axs[1, 1].set_title(f"Im[$u_{{\phi}}'(r,\phi)$]")
                    
                    for i in range(2):
                        for j in range(2):
                            plot_rMaxVelocity_polarGrid(axs[i, j], params)
                            axs[i, j].grid(True) #Restore grids for final version

                    fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
                    fig.suptitle(f"Velocities derived from fastest-growing eigen-streamfunction \n in $r\phi$-plane for wavenumbers $k_{{\phi}} =$ {kφ}, $m =$ {kz:.4f} {params.units['kz']}")
                    fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], 
                                 location = "right", shrink = 0.6,
                                 label = params.units["u"], pad = 0.1)
                    fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], 
                                 location = "right", shrink = 0.6,
                                 label = params.units["u"], pad = 0.1)
                    fig.savefig(f"./Graphs/velocities_k{int(kφ)}_m{kz:.4E}_{modeString}_{params.dimString}1Dgyre_{setupString}.png")
                    plt.close(fig)

def RunVisualization(params, geom, modes, kφs, kzs, r, φ, z, 
                     dimensionalGrowthRates, dimensionalPropSpeeds, 
                     eigModesReal, eigModesImag, eigStreamfnsReal, 
                     eigStreamfnsImag, eig_urReal, eig_urImag, eig_uφReal,
                     eig_uφImag, setupString):
                     
    print("Starting visualization")
    
    plt.rcParams.update({"text.usetex": True, "font.size": 12})

    makedirs("./Graphs", exist_ok = True) #Make folder if nonexistent
    
    PlotEigvals(params, len(modes), 
                kφs, kzs, dimensionalGrowthRates, dimensionalPropSpeeds, 
                setupString)
     
    PlotEigModeStructures(params, len(modes), kφs, kzs, r, z, eigModesReal, 
                          eigModesImag, setupString)
    
    if (eigStreamfnsReal is not None and eig_urReal is not None):      
        PlotStreamfnsAndVelocities(params, geom, 1, kφs, kzs, r, φ, z, 
                                   eigStreamfnsReal, eigStreamfnsImag,
                                   eig_urReal, eig_urImag, eig_uφReal, 
                                   eig_uφImag, setupString)

def RunVisFromSavedData(params, geom):

    if (params.discretizeRadial and not params.discretizeVertical):
    
        commonVariables, kzs, growthDim, propDim = LoadSavedData1DRadial(params,
                                                                         geom)
        z                                        = None
        
        setupString = f"Lr{params.Lr:.1E}_Nr{params.Nr}_{params.stratification_kw}Strat_Ro{params.Ro:.1E}_BuInf_f{params.f0:.1E}"
        
    elif (params.discretizeVertical and not params.discretizeRadial):
    
        commonVariables, z, growthDim, propDim = LoadSavedData1DVertical(params,
                                                                         geom)
        kzs                                    = None
        
        setupString = f"Lz{params.Lz:.1E}_Nz{params.Nz}_{params.stratification_kw}Strat_Ro{params.Ro:.1E}_Bu{params.Bu:.1E}_f{params.f0:.1E}"
        
    elif (params.discretizeRadial and params.discretizeVertical):
    
        commonVariables, z, growthDim, propDim = LoadSavedData2D(params, geom)
        kzs                                    = None
        
        setupString = f"Lr{params.Lr:.1E}_Lz{params.Lz:.1E}_Nr{params.Nr}_Nz{params.Nz}_{params.stratification_kw}Strat_Ro{params.Ro:.1E}_Bu{params.Bu:.1E}_f{params.f0:.1E}"

    RunVisualization(params, geom, commonVariables["modes"], 
                     commonVariables["kφs"], kzs, commonVariables["r"], 
                     commonVariables["φ"], z, growthDim, propDim,
                     commonVariables["eigModesReal"], 
                     commonVariables["eigModesImag"], 
                     commonVariables["eigStreamfnsReal"], 
                     commonVariables["eigStreamfnsImag"], 
                     commonVariables["eig_urReal"], 
                     commonVariables["eig_urImag"], 
                     commonVariables["eig_uφReal"], 
                     commonVariables["eig_uφImag"], setupString)