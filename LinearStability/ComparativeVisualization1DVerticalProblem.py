import argparse
import numpy as np
import matplotlib.pyplot as plt

from VisualizationFunctions import *

parser = argparse.ArgumentParser()
parser.add_argument("-Nz", 
                    help = "Number of grid points", type = int, default = 100)
parser.add_argument("-Lz", 
                    help = "Domain depth", type = float, default = 3.3)
parser.add_argument("--dimString",
                    help = "Use data from which solvers?",
                    choices = ("dimensional", "nondimensional"), type = str,
                    default = "nondimensional")
parser.add_argument("--strat_shape",
                    help = "Stratification type",
                    choices = ("constant", "TWB", "doubleTanh", "doubleTanhTWB"), 
                    default = "doubleTanhTWB")
parser.add_argument("--sigmar",
                    help = "Radial gyre length scale (in m for 'dimensional'; 1 for 'nondimensional')",
                    type = float, default = 1)
parser.add_argument("--sigmaz",
                    help = "Vertical gyre length scale (in m for 'dimensional'; 1 for 'nondimensional')",
                    type = float, default = 1)
parser.add_argument("-r", 
                    help = "r-values at which problem was evaluated (in m for 'dimensional'; dimensionless for 'nondimensional'); enter as (-r start stop step) or as r-values themselves",
                    type = float, default = [1e-5, 10, 1e-2], nargs = "*")
parser.add_argument("-Ro",
                    help = "Rossby number(s) of background profile(s)", 
                    type = float, nargs = "+")
parser.add_argument("-Bu",
                    help = "Burger number(s) of background profile(s)", 
                    type = float, nargs = "+")
parser.add_argument("-f0",
                    help = "Coriolis frequency (Hz)",
                    type = float, default = 1.4e-4)        
args = vars(parser.parse_args())

#Shorthand for linear spacing; parse as [start] [stop] [step]
if len(args["r"]) == 3: 
    rs = np.arange(args["r"][0], args["r"][1], args["r"][2])
                                        
else: #List of r-values is provided
    rs = np.array(args["r"])
                    
for Bu in args["Bu"]:

    fig_eigvecs, ax_eigvecs = plt.subplots(figsize = (10, 10))
    ax_eigvecs.grid(True)
    
    fig_comparative_growths, ax_comparative_growths = plt.subplots(figsize = (10, 10))
    ax_comparative_growths.grid(True)
        
    r_plot = args["sigmar"]
    r_idx  = np.abs(rs - r_plot).argmin()
    r_int  = int(rs[r_idx])
    
    for Ro in args["Ro"]:
        
        #Load and parse saved data
        
        ds_Ro = Dataset(f"./Data/{args['dimString']}_Lz{args['Lz']:.1E}_Nz{args['Nz']}_{args['strat_shape']}Strat_Ro{Ro:.1E}_Bu{Bu:.1E}_f{args['f0']:.1E}.nc")
        
        commonVariables        = LoadCommonVariables(ds_Ro)
        z                      = ds_Ro.variables["z"][:]
        dimensionalGrowthRates = ds_Ro.variables["growth_rate"][:, :, :]
        
        ds_Ro.close()
            
        modes, kfs                 = commonVariables[0], commonVariables[1]
        eigModesReal, eigModesImag = commonVariables[4], commonVariables[5]
        
        #Determine the overall most unstable mode, at sigma_r, for this Ro
        kf_idx           = dimensionalGrowthRates[r_idx, :, 0].argmax()
        most_unstable_kf = kfs[kf_idx]
        
        #Normalize the most unstable eigenmode
        eigModeReal, eigModeImag = Normed(eigModesReal[kf_idx, r_idx, :, 0],
                                          eigModesImag[kf_idx, r_idx, :, 0])
        
        #Plot components of most unstable eigenmode, for this Ro, against z
        ax_eigvecs.plot(eigModeReal, z, "-", label = f"Real component; Ro = {Ro}, k = {most_unstable_kf}")
        ax_eigvecs.plot(eigModeImag, z, "--",label = f"Imaginary component; Ro = {Ro}, k = {most_unstable_kf}")
        
        for mode in range(len(modes)):
            #Plot growth rate for this Ro, at r_plot, against k_phi
            ax_comparative_growths.scatter(kfs, dimensionalGrowthRates[r_idx, :, mode], label = f"Growth rate; mode {mode}")
        
        #Plot max overall growth rate against k_phi and r
        
        kfMesh, rMesh = np.meshgrid(kfs, rs)
        
        fig_growth, ax_growth = plt.subplots(figsize = (10, 5))
        
        ax_growth.grid(False) #Required for pcolormesh
        
        ax_growth.pcolormesh(rMesh, kfMesh, dimensionalGrowthRates[:, :, 0], 
                                      cmap = "Reds")
        
        ax_growth.grid(True) #Restore grids for final version
        ax_growth.set(xlabel = "r", ylabel = "k", title = f"Largest growth rate; Ro = {Ro}")
        fig_growth.colorbar(ScalarMappable(norm = Normalize(
                            vmin = 0, vmax = np.max(dimensionalGrowthRates[:, :, 0])), cmap = "Reds"), 
                             ax = ax_growth, location = "right",
                             shrink = 0.6, label = "Growth rate (s$^{-1})$", pad = 0.1)
        fig_growth.savefig(f"./Graphs/largestGrowthRate_vs_k_and_r_Ro{Ro}_Bu{Bu}_{args['dimString']}1Dgyre.png")
        plt.close(fig_growth)
        
    ax_comparative_growths.set(xlabel = "$k_{\phi}$", ylabel = "Growth rate (s$^{-1}$)",
                               title = f"Growth rate vs. $k_{{\phi}}$ at $r = {r_int}$")
    ax_comparative_growths.legend()
    fig_comparative_growths.savefig(f"./Graphs/growthRate_vs_k_r{r_int}_Nz{args['Nz']}_Bu{Bu}_{args['dimString']}1Dgyre.png")
    plt.close(fig_comparative_growths)
        
    plot_sigmaz(ax_eigvecs, args["sigmaz"], args["Lz"]) #Gyre length scale
    #plot_stratification_peaks(ax, params) #z_s and z_d
    ax_eigvecs.set(xlabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$")#,
                       #ylabel = f"$z$ ({params.units['z']})",
                       #title = f"Components of fastest-growing eigenvector for $k_{{\phi}}$ = {kf} and $r =$ {r_int} {params.units['r']}")
    ax_eigvecs.legend()
    fig_eigvecs.savefig(f"./Graphs/eigModeStructures_variableRo_r{r_int}_fastestgrowing_{args['dimString']}1Dgyre.png")
    plt.close(fig_eigvecs)