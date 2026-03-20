from netCDF4 import Dataset

def RunDimensionalBaroclinicVis(Lr, Lz, Nr, Nz, Ro, Bu):

    ds = Dataset(f"./Data/data_Lr{Lr:.1E}_Lz{Lz:.1E}_Nr{Nr}_Nz{Nz}_Ro{Ro:.1E}_Bu{Bu:.1E}.nc")

    modes = ds.variables["mode"][:]
    kφs   = ds.variables["kφ"][:]
    r     = ds.variables["r"][:]
    φ     = ds.variables["φ"][:]
    z     = ds.variables["z"][:]
    
    growth_rates = ds.variables["growth_rate"][:, :]
    prop_speeds  = ds.variables["prop_speed"][:, :]

    eigVecs = ds.variables["eigVec"][:, :, :, :]
    eigStreamfns = ds.variables["eigStreamfn"][:, :, :, :]

    ds.close()

RunDimensionalBaroclinicVis(8.0, 3.3, 21, 10, 4.0e-3, 2.5e-3)
