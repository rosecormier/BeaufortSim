import numpy as np

from math import pi
from netCDF4 import Dataset
from os import makedirs

from Streamfunctions import *

def SaveToNetCDF(params, geom, dimensionalGrowthRates, dimensionalPropSpeeds,
                 eigVecs):
    
    dφ = 2 * pi / params.Np           #φ-increment
    φ  = dφ * np.arange(0, params.Np) #φ-coords
    
    if params.discretizeRadial:
        Lr, Nr, halfNr = params.Lr, params.Nr, params.halfNr

    if (params.discretizeRadial and not params.discretizeVertical):
        Bu, gridString = "Inf", f"Lr{Lr:.1E}_Nr{Nr}"
        
    elif params.discretizeVertical:
        Lz, Nz, Bu = params.Lz, params.Nz, f"{params.Bu:.1E}"
        
        if params.discretizeRadial:
            gridString = f"Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}"
        elif not params.discretizeRadial:
            gridString = f"Lz{Lz:.1E}_Nz{Nz}"
        
    makedirs("./Data", exist_ok = True) #Create data directory if nonexistent

    ncfile = Dataset(f"./Data/{params.dimString}_{gridString}_Ro{params.Ro:.1E}_Bu{Bu}_f{params.f0:.1E}.nc",
                     mode = "w", auto_complex = True)

    kφs    = params.kps
    nmodes = params.nmodes

    #Create coordinate dimensions
    mode_dim = ncfile.createDimension("mode", nmodes)
    kφ_dim   = ncfile.createDimension("kφ", len(kφs))
    φ_dim    = ncfile.createDimension("φ", len(φ))
    
    if params.discretizeRadial:
        r_dim = ncfile.createDimension("r", (halfNr + 1))
    elif (params.discretizeVertical and not params.discretizeRadial):
        r_dim = ncfile.createDimension("r", len(params.rs))

    #Create variables corresponding to each coordinate dimension
    mode_var = ncfile.createVariable("mode", float, ("mode",))
    kφ_var   = ncfile.createVariable("kφ", float, ("kφ",))
    φ_var    = ncfile.createVariable("φ", float, ("φ",))
    r_var    = ncfile.createVariable("r", float, ("r",))

    #Store data corresponding to each coordinate dimension
    mode_var[:] = range(nmodes)
    kφ_var[:]   = kφs
    φ_var[:]    = φ
    
    if params.discretizeRadial:
        r_var[:] = geom.r[0:(halfNr + 1)] #Physical coordinates only
    elif (params.discretizeVertical and not params.discretizeRadial):
        r_var[:] = params.rs[:]
    
    #Store units
    r_var.units = params.units["r"]
    
    if (params.discretizeRadial and not params.discretizeVertical):
    
        #Also save kz-information
        
        kzs       = params.kzs
        kz_dim    = ncfile.createDimension("kz", len(kzs))
        kz_var    = ncfile.createVariable("kz", float, ("kz",))
        kz_var[:] = kzs
        
        kz_var.units = params.units["kz"]
            
        #Create variables for dimensional growth rates and propagation speeds
        growth_rate = ncfile.createVariable("growth_rate", float,
                                            (kz_dim, kφ_dim, mode_dim))
        prop_speed  = ncfile.createVariable("prop_speed", float,
                                            (kz_dim, kφ_dim, mode_dim))
                                            
        #Save growth-rate and propagation-speed data
        growth_rate[:, :, :] = dimensionalGrowthRates[:, :, :]
        prop_speed[:, :, :]  = dimensionalPropSpeeds[:, :, :]
                                            
        #Create vars for eigenmodes, corresponding streamfunctions/velocities
        eigMode     = ncfile.createVariable("eigMode", complex,
                                            (kz_dim, kφ_dim, r_dim, mode_dim))
        eigStreamfn = ncfile.createVariable("eigStreamfn", complex,
                                            (kz_dim, kφ_dim, r_dim, φ_dim,
                                             mode_dim))
        eig_ur      = ncfile.createVariable("eig_ur", complex,
                                            (kz_dim, kφ_dim, r_dim, φ_dim, 
                                             mode_dim))
        eig_uφ      = ncfile.createVariable("eig_uφ", complex,
                                            (kz_dim, kφ_dim, r_dim, φ_dim,
                                             mode_dim))

    elif params.discretizeVertical:
    
        #Save z-information
        z_dim    = ncfile.createDimension("z", (Nz + 1))
        z_var    = ncfile.createVariable("z", float, ("z",))
        z_var[:] = geom.z[0:(Nz + 1)]
        
        z_var.units = params.units["z"]
        
        if not params.discretizeRadial:
        
            #Create variables for dimensional growth rates and prop. speeds
            growth_rate = ncfile.createVariable("growth_rate", float,
                                                (r_dim, kφ_dim, mode_dim))
            prop_speed  = ncfile.createVariable("prop_speed", float,
                                                (r_dim, kφ_dim, mode_dim))
                                                
            #Save growth-rate and propagation-speed data
            growth_rate[:, :] = dimensionalGrowthRates[:, :, :]
            prop_speed[:, :]  = dimensionalPropSpeeds[:, :, :]
    
        elif params.discretizeRadial:
    
            #Create variables for dimensional growth rates and prop. speeds
            growth_rate = ncfile.createVariable("growth_rate", float,
                                                (kφ_dim, mode_dim))
            prop_speed  = ncfile.createVariable("prop_speed", float,
                                                (kφ_dim, mode_dim))
                                                
            #Save growth-rate and propagation-speed data
            growth_rate[:, :] = dimensionalGrowthRates[:, :]
            prop_speed[:, :]  = dimensionalPropSpeeds[:, :]
        
        #Create variables for eigenmodes and corresponding streamfunctions
        eigMode     = ncfile.createVariable("eigMode", complex,
                                            (kφ_dim, r_dim, z_dim, mode_dim))
        eigStreamfn = ncfile.createVariable("eigStreamfn", complex,
                                            (kφ_dim, r_dim, φ_dim, z_dim, 
                                             mode_dim))
        eig_ur      = ncfile.createVariable("eig_ur", complex,
                                            (kφ_dim, r_dim, φ_dim, z_dim,
                                             mode_dim))
        eig_uφ      = ncfile.createVariable("eig_uφ", complex,
                                            (kφ_dim, r_dim, φ_dim, z_dim,
                                             mode_dim))
        
    #Store units -- always dimensional version
    growth_rate.units = "per second"
    prop_speed.units  = "per second"

    for kφ in kφs:

        kφ_idx = kφs.tolist().index(kφ)
        
        for mode in range(nmodes):
            
            if params.discretizeVertical:
                                                        
                #Evaluate and save eigen-velocities at discrete grid points
                for ell in range(len(φ)):
                
                    if params.discretizeRadial:
                    
                        eigVel = EigvelFrom2DEigmode(params, geom, 
                                                     eigVecs[kφ_idx, :, mode],
                                                     kφ_idx)
                                                     
                        eig_ur[kφ_idx, :, ell, :, mode] = eigVel[0]
                        eig_uφ[kφ_idx, :, ell, :, mode] = eigVel[1]
                        
                        #Reshape and save eigenmode data
                        eigMode[kφ_idx, :, :, 
                                mode] = np.reshape(eigVecs[kφ_idx, :, mode],
                                                   ((halfNr + 1), (Nz + 1))
                                                  )
                                                 
                    elif not params.discretizeRadial:
                    
                        for r_idx in range(len(params.rs)):
                         
                            eigVel = EigvelFrom1DEigvec(params, geom, 
                                                        eigVecs[:, kφ_idx, :, mode],
                                                        kφ_idx)
                                                 
                            eig_ur[kφ_idx, :, ell, :, mode] = eigVel[0]
                            eig_uφ[kφ_idx, :, ell, :, mode] = eigVel[1]
                                                        
                #Evaluate and save eigen-streamfunction at discrete grid points
                for ell in range(len(φ)):
                    eigStreamfn[kφ_idx, :, 
                                ell, :, mode] = Streamfunction(eigMode[kφ_idx,
                                                                       :, :, 
                                                                       mode],
                                                               k = kφ,
                                                               φ = φ[ell])
            
            elif (params.discretizeRadial and not params.discretizeVertical):
            
                for kz in kzs:
                
                    kz_idx = kzs.tolist().index(kz)
                    
                    #Save eigenmode data
                    eigMode[kz_idx, kφ_idx, :, mode] = eigVecs[kz_idx, kφ_idx,
                                                               :, mode]
            
                    #Evaluate streamfunction and velocities at discrete grid points
                    for ell in range(len(φ)):
                    
                        eigStreamfn[kz_idx,
                                    kφ_idx,
                                    :, ell, 
                                    mode]   = Streamfunction(eigMode[kz_idx,
                                                                     kφ_idx, :, 
                                                                     mode],
                                                             k = kφ,
                                                             φ = φ[ell])
                                                             
                        eigVels = EigvelFrom1DEigvec(params, geom, 
                                                     eigVecs[kz_idx, kφ_idx, :,
                                                             mode], 
                                                     kφ_idx)

                        eig_ur[kz_idx, kφ_idx, :, ell, mode] = eigVels[0]
                        eig_uφ[kz_idx, kφ_idx, :, ell, mode] = eigVels[1]

    #Store units (depending on dimensionality of problem)
    eigStreamfn.units = params.units["psi"]
    eig_ur.units      = params.units["u"]
    eig_uφ.units      = params.units["u"]

    ncfile.close()