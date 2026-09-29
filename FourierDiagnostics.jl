include("LibraryCoordinates.jl")

using AbstractFFTs
using Interpolations
using Oceananigans.Grids
using OffsetArrays: no_offset_view
using Statistics

function construct_uniform_polar_grid(Lr, Nr, Nφ)
   #=
   Return Nr evenly spaced gridpoints between 'dr' and 'Lr', and Nφ evenly
    spaced gridpoints between 0 and 2π (without double-counting φ % 2π = 0).
   =#

   dr = Lr / (Nr + 1)
   dφ = 2π / Nφ

   r_gridpoints = dr:dr:(Lr - dr)
   φ_gridpoints = 0:dφ:(2π - dφ)

   return r_gridpoints, φ_gridpoints
end

function compute_φFFT_of_sim_data(simField, grid, gridParams, Nr, Nφ;
                                  Hx = nothing, Hy = nothing, Hz = nothing,
                                  rAvg = false, zAvg = false)
   #=
   Compute the Nφ-point azimuthal Fourier transform of 'simField' (which need
    not actually be a prognostic field of a simulation, but does need to be
    defined on an Oceananigans Grid) by first interpolating to a uniformly-
    spaced grid in cylindrical coordinates.
   Default returns a 3D Nr x Nkφ x Nz field, but options 'rAvg' and 'zAvg'
    perform r- and z-averaging, respectively, before returning.
   =#

   #If halos not otherwise specified, set them to 'grid' halo sizes
   if isnothing(Hx)
      Hx = gridParams.Hx
   end
   if isnothing(Hy)
      Hy = gridParams.Hy
   end
   if isnothing(Hz)
      Hz = gridParams.Hz
   end

   #First, construct a uniformly spaced (Nr x Nφ x Nz) grid
   
   #Get the r- and φ-gridpoints, in polar coordinates, of a uniform polar grid
   rCylGrid, φCylGrid = construct_uniform_polar_grid(1, Nr, Nφ) #gridParams.Lr, Nr, Nφ)
   
   #Read in the coordinates of Oceananigans gridpoints
   
   if location(simField)[1] == Center
      xCartVec = no_offset_view(grid.xᶜᵃᵃ)[(Hx + 1):(end - Hx)]
   elseif location(simField)[1] == Face
      xCartVec = no_offset_view(grid.xᶠᵃᵃ)[(Hx + 1):(end - Hx)]
   end

   if location(simField)[2] == Center
      yCartVec = no_offset_view(grid.yᵃᶜᵃ)[(Hy + 1):(end - Hy)]
   elseif location(simField)[2] == Face
      yCartVec = no_offset_view(grid.yᵃᶠᵃ)[(Hy + 1):(end - Hy)]
   end

   if location(simField)[3] == Center
      zCartVec = no_offset_view(grid.z.cᵃᵃᶜ)[(Hz + 1):(end - Hz)]
   elseif location(simField)[3] == Face
      zCartVec = no_offset_view(grid.z.cᵃᵃᶠ)[(Hz + 1):(end - Hz)]
   end
   
   #Tile xCartVec and yCartVec along complementary dimensions, such that each
   # point (x[i], y[j]) on the Oceananigans horizontal grid is represented by
   # (xCartGrid[i, j], yCartGrid[i, j]). Then tile zCartVec in both x- and
   # y-directions to match.
   
   xCartGrid = repeat(xCartVec, inner = [length(yCartVec)])
   yCartGrid = repeat(yCartVec, outer = [length(xCartVec)])
   zCartGrid = repeat(zCartVec, inner = [length(yCartVec)], 
                      outer = [length(xCartVec)])

   #Cartesian coordinates of the points on the cylindrically uniform grid
   xCylGrid, yCylGrid, zCylGrid = compute_Cart_coords(rCylGrid,
                                                      φCylGrid,
                                                      zCartVec)

   #Construct an uninitialized Array to store interpolated data
   interpolated_data = Array{Float64}(undef, size(xCylGrid)[1], 
                                      size(xCylGrid)[2], length(zCylGrid))
   
   #Interpolate data at each z-level. Interpolation only needs to be done in 2D
   # (since zCylGrid == zCartGrid).
   
   for k in 1:1:length(zCylGrid)

      #Create a linear interpolation object without extrapolation
      interp = Interpolations.scale(Interpolations.interpolate(
                              interior(simField)[:, :, k],
                              BSpline(Interpolations.Linear())), 
                        (xCartVec, yCartVec)
                                   )

      #Interpolate 'simField' to CylGrid points
      interpolated_k = interp.(xCylGrid, yCylGrid)
      
      interpolated_data[:, :, k] = interpolated_k
   end

   kφs = (rfftfreq(Nφ, 2π * Nφ)) ./ 2π #The azimuthal wavenumbers of the FFT
   
   #Construct an uninitialized Array to store Fourier-transformed data
   simFieldφFFT = Array{ComplexF64}(undef, Nr, length(kφs), length(zCylGrid))
   
   for k in 1:1:length(zCylGrid) #Loop over z-coords
      for i in 1:1:length(rCylGrid) #Loop over r-coords
         simFieldφFFT[i, :, k] = rfft(interpolated_data[i, :, k]) #[Real] DFT
      end
   end

   if rAvg
   
      if zAvg
         return kφs, mean(simFieldφFFT, dims = (1, 3))
      elseif !zAvg
         return kφs, mean(simFieldφFFT, dims = 1)
      end
      
   elseif (!rAvg && zAvg)
      return kφs, mean(simFieldφFFT, dims = 3)
   
   else
      return kφs, simFieldφFFT
   end
end