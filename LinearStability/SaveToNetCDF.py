import numpy as np

from math import pi
from netCDF4 import Dataset
from os import makedirs

from Streamfunctions import Streamfunction

def SaveToNetCDF(params, geom, dimensionalGrowthRates, dimensionalPropSpeeds,
                 eigVecs):
    
    Lr, Nr, halfNr = params.Lr, params.Nr, params.halfNr
    
    dφ = 2 * pi / params.Np           #φ-increment
    φ  = dφ * np.arange(0, params.Np) #φ-coords
    
    f0 = params.f0

    makedirs("./Data", exist_ok = True) #Create data directory if nonexistent
    
    #Whether gen. eig. problem is 2D
    discretizeVertical = params.discretizeVertical
    
    if discretizeVertical:
        Lz, Nz = params.Lz, params.Nz
    
    #Whether variables are nondimensional
    nondimensional = params.nondimensional
    
    if nondimensional:
    
        Ro = params.Ro
    
        if discretizeVertical:
        
            Bu = f"{params.Bu:.1E}"
            
            gridString = f"Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}"
            
        else:
            Bu, gridString = "Inf", f"Lr{Lr:.1E}_Nr{Nr}"
    
        ncfile = Dataset(f"./Data/nondimensional_{gridString}_Ro{Ro:.1E}_Bu{Bu}_f{f0:.1E}.nc", 
                     mode = "w", auto_complex = True)

    else:
    
        Ro = params.Umax / (params.sigmar * f0)
    
        if discretizeVertical:
        
            Bu = (params.Nmax * params.sigmaz / (f0 * params.sigmar))**2
            Bu = f"{Bu:.1E}"
            
            gridString = f"Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}"
            
        else:
            Bu, gridString = "Inf", f"Lr{Lr:.1E}_Nr{Nr}"
        
        ncfile = Dataset(f"./Data/dimensional_{gridString}_Ro{Ro:.1E}_Bu{Bu}_f{f0:.1E}.nc",
                     mode = "w", auto_complex = True)

    kφs    = params.kps
    nmodes = params.nmodes

    #Create coordinate dimensions
    mode_dim = ncfile.createDimension("mode", nmodes)
    kφ_dim   = ncfile.createDimension("kφ", len(kφs))
    r_dim    = ncfile.createDimension("r", (halfNr + 1))
    φ_dim    = ncfile.createDimension("φ", len(φ))

    #Create variables corresponding to each coordinate dimension
    mode_var = ncfile.createVariable("mode", float, ("mode",))
    kφ_var   = ncfile.createVariable("kφ", float, ("kφ",))
    r_var    = ncfile.createVariable("r", float, ("r",))
    φ_var    = ncfile.createVariable("φ", float, ("φ",))

    #Store data corresponding to each coordinate dimension
    mode_var[:] = range(nmodes)
    kφ_var[:]   = kφs
    r_var[:]    = geom.r[0:(halfNr + 1)] #Physical coordinates only
    φ_var[:]    = φ
    
    #Store units
    if nondimensional:
        r_var.units = "times σ_r"
    else:
        r_var.units = "metres"
    
    if not discretizeVertical:
    
        #For 1D gen. eig. problem, also save kz-information
        
        kzs       = params.kzs
        kz_dim    = ncfile.createDimension("kz", len(kzs))
        kz_var    = ncfile.createVariable("kz", float, ("kz",))
        kz_var[:] = kzs
        
        if nondimensional:
            kz_var.units = "per σ_z"
        else:
            kz_var.units = "per metre"
            
        #Create variables for dimensional growth rates and propagation speeds
        growth_rate = ncfile.createVariable("growth_rate", float,
                                            (kz_dim, kφ_dim, mode_dim))
        prop_speed  = ncfile.createVariable("prop_speed", float,
                                            (kz_dim, kφ_dim, mode_dim))
                                            
        #Create variables for eigenmodes and corresponding streamfunctions
        eigMode     = ncfile.createVariable("eigMode", complex,
                                            (kz_dim, kφ_dim, r_dim, mode_dim))
        eigStreamfn = ncfile.createVariable("eigStreamfn", complex,
                                            (kz_dim, kφ_dim, r_dim, φ_dim,
                                             mode_dim))

    elif discretizeVertical:
    
        #For 2D gen. eig. problem, also save z-information
        z_dim    = ncfile.createDimension("z", (Nz + 1))
        z_var    = ncfile.createVariable("z", float, ("z",))
        z_var[:] = geom.z[0:(Nz + 1)]
        
        if nondimensional:
            z_var.units = "times σ_z"
        else:
            z_var.units = "metres"
    
        #Create variables for dimensional growth rates and propagation speeds
        growth_rate = ncfile.createVariable("growth_rate", float,
                                            (kφ_dim, mode_dim))
        prop_speed  = ncfile.createVariable("prop_speed", float,
                                            (kφ_dim, mode_dim))
                                            
        #Create variables for eigenmodes and corresponding streamfunctions
        eigMode     = ncfile.createVariable("eigMode", complex,
                                            (kφ_dim, r_dim, z_dim, mode_dim))
        eigStreamfn = ncfile.createVariable("eigStreamfn", complex,
                                            (kφ_dim, r_dim, φ_dim, z_dim, 
                                             mode_dim))
        
    #Store units -- always dimensional version
    growth_rate.units = "per second"
    prop_speed.units  = "per second"
    
    #Save growth-rate and propagation-speed data
    growth_rate[:, :] = dimensionalGrowthRates[:, :]
    prop_speed[:, :]  = dimensionalPropSpeeds[:, :]

    for kφ in kφs:

        kφ_idx = kφs.tolist().index(kφ)
        
        for mode in range(nmodes):
            
            if discretizeVertical:
                #Save eigenmode data
                eigMode[kφ_idx, :, :, mode] = np.reshape(eigVecs[kφ_idx, :, 
                                                                 mode],
                                                         ((halfNr + 1), 
                                                          (Nz + 1)
                                                         )
                                                        )
                                                        
                #Evaluate streamfunction at discrete grid points and save result
                for ell in range(len(φ)):
                    eigStreamfn[kφ_idx, :, 
                                ell, :, mode] = Streamfunction(eigMode[kφ_idx,
                                                                       :, :, 
                                                                       mode],
                                                               k = kφ,
                                                               φ = φ[ell])
            
            elif not discretizeVertical:
            
                for kz in kzs:
                
                    kz_idx = kzs.tolist().index(kz)
                    
                    #Save eigenmode data
                    eigMode[kz_idx, kφ_idx, :, mode] = eigVecs[kz_idx, kφ_idx,
                                                               :, mode]
            
                    #Evaluate streamfunction at discrete grid points and save result
                    for ell in range(len(φ)):
                        eigStreamfn[kz_idx,
                                    kφ_idx,
                                    :, ell, 
                                    mode]   = Streamfunction(eigMode[kz_idx,
                                                                     kφ_idx, :, 
                                                                     mode],
                                                             k = kφ,
                                                             φ = φ[ell])

    ncfile.close()