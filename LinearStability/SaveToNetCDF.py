import numpy as np
import os 

from netCDF4 import Dataset

from Streamfunctions import Streamfunction

def SaveToNetCDF(params, growth_rates, prop_speeds, eigVecs, φCoords):

    Ro = params.Ro
    Bu = params.Bu
    
    Lr = params.Lr
    Nr = params.Nr

    Nφ = params.Nφ

    Lz = params.Lz
    Nz = params.Nz

    os.makedirs("./Data", exist_ok = True)    
    ncfile = Dataset(f"./Data/data_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.nc", 
                     mode = "w",
                     auto_complex = True)

    kφs    = params.kφs
    nmodes = params.nmodes

    mode_dim = ncfile.createDimension("mode", nmodes)
    kφ_dim   = ncfile.createDimension("kφ", len(kφs))
    r_dim    = ncfile.createDimension("r", Nr // 2)
    φ_dim    = ncfile.createDimension("φ", Nφ)
    z_dim    = ncfile.createDimension("z", Nz)

    growth_rate = ncfile.createVariable("growth_rate", float, (kφ_dim, mode_dim))
    prop_speed  = ncfile.createVariable("prop_speed", float, (kφ_dim, mode_dim))

    eigVec = ncfile.createVariable("eigVec", complex, (kφ_dim, r_dim, z_dim, mode_dim))
    eigStreamfn   = ncfile.createVariable("eigStreamfn", complex, (kφ_dim, r_dim, φ_dim, z_dim, mode_dim))

    growth_rate[:, :] = growth_rates[:, :]
    prop_speed[:, :]  = prop_speeds[:, :]

    for k in range(len(kφs)):
        for mode in range(nmodes):
            eigVec[k, :, :, mode] = np.reshape(eigVecs[k, :, mode], 
                                               (params.halfNr, Nz), 
                                               order = "C")
            for ell in range(Nφ):
                eigStreamfn[k, :, ell, :, mode] = Streamfunction(np.reshape(eigVecs[k, :, mode],
                                                     (params.halfNr, Nz),
                                                     order = "F"),
                                                     k = k,
                                                     φ = φCoords[ell])

    ncfile.close()
