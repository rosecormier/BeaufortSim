include("LibraryCoordinates.jl")

using AbstractFFTs, Interpolations
using Oceananigans.BoundaryConditions, Oceananigans.Grids
using Oceananigans
using OffsetArrays: no_offset_view
using Statistics

function ConstructUniformPolarGrid(Lr::Float64, Nr::Int64, Nφ::Int64)
   #=
   Return Nr evenly spaced gridpoints between 'dr' and 'Lr', and Nφ evenly
    spaced gridpoints between 0 and 2π (without double-counting φ % 2π = 0).
   =#

   dr, dφ = Lr / (Nr + 1), 2π / Nφ

   r_gridpoints = dr:dr:(Lr - dr)
   φ_gridpoints = 0:dφ:(2π - dφ)

   return r_gridpoints, φ_gridpoints
end

@inline φWavenumbers(Nφ::Int64) = (rfftfreq(Nφ, 2π * Nφ)) ./ 2π

function compute_φFFT_of_sim_Field(simField::AbstractField, grid, p::GridParams,
                                   Nr::Int64, Nφ::Int64;
                                   rCoordsCylGrid = nothing,
                                   φCoordsCylGrid = nothing, kφs = nothing,
                                   rAvg = false, zAvg = false)
   #=
   Compute the Nφ-point azimuthal Fourier transform of 'simField' (which need
    not actually be a prognostic field of a simulation, but does need to be an
    Oceananigans Field defined on an Oceananigans Grid) by first interpolating
    to a uniformly-spaced grid in cylindrical coordinates.
   Default returns a 3D Nr x Nkφ x Nz field, but options 'rAvg' and 'zAvg'
    perform r- and z-averaging, respectively, before returning.
   =#

   ####################################################
   # Construct a uniformly spaced (Nr x Nφ x Nz) grid #
   ####################################################
   
   #Get the r- and φ-gridpoints, in polar coordinates, of a uniform polar grid
   if (isnothing(rCoordsCylGrid) || isnothing(φCoordsCylGrid))
      rCoordsCylGrid, φCoordsCylGrid = ConstructUniformPolarGrid(p.Lr, Nr, Nφ)
   end
 
   #Read in the coordinates of Oceananigans (Cartesian) gridpoints
   
   if location(simField)[1] == Center
      xCoordsCartGrid = xᶜInt(grid, p)
   elseif location(simField)[1] == Face
      xCoordsCartGrid = xᶠInt(grid, p)
   end

   if location(simField)[2] == Center
      yCoordsCartGrid = yᶜInt(grid, p)
   elseif location(simField)[2] == Face
      yCoordsCartGrid = yᶠInt(grid, p)
   end

   if location(simField)[3] == Center
      zCoordsCartGrid = zᶜInt(grid, p)
   elseif location(simField)[3] == Face
      zCoordsCartGrid = zᶠInt(grid, p)
   end
   
   #Tile xCoordsCartGrid and yCoordsCartGrid along complementary dimensions,
   # such that each point (x[i], y[j]) on the Oceananigans horizontal grid is
   # represented by (xCoordsCartGridTiled[i, j], yCoordsCartGridTiled[i, j]).
   #Then tile zCoordsCartGrid in both x- and y-directions to match.
   
   xCoordsCartGridTiled = repeat(xCoordsCartGrid,
                                 inner = [length(yCoordsCartGrid)])
   yCoordsCartGridTiled = repeat(yCoordsCartGrid, 
                                 outer = [length(xCoordsCartGrid)])
   zCoordsCartGridTiled = repeat(zCoordsCartGrid, 
                                 inner = [length(yCoordsCartGrid)], 
                                 outer = [length(xCoordsCartGrid)])

   #Cartesian coordinates of the points on the cylindrically uniform grid
   xCoordsCylGrid, yCoordsCylGrid, zCoordsCylGrid = ComputeCartesianCoords(rCoordsCylGrid, φCoordsCylGrid, zCoordsCartGrid)

   ###############################################################
   # Perform the interpolation of data to the uniform-polar grid #
   ###############################################################

   #Construct an uninitialized Array to store interpolated data
   interpolated_data = Array{Float64}(undef, size(xCoordsCylGrid)[1], 
                                      size(xCoordsCylGrid)[2], 
                                      length(zCoordsCylGrid))
   
   #Interpolate data at each z-level. Interpolation only needs to be done in 2D
   # (since zCoordsCylGrid == zCoordsCartGrid).
   
   for k in 1:length(zCoordsCylGrid)

      #Create a linear interpolation object without extrapolation
      interp = Interpolations.scale(Interpolations.interpolate(
                              interior(simField)[:, :, k],
                              BSpline(Interpolations.Linear())), 
                        (xCoordsCartGrid, yCoordsCartGrid)
                                   )

      #Interpolate 'simField' to CylGrid points
      interpolated_k = interp.(xCoordsCylGrid, yCoordsCylGrid)
      
      interpolated_data[:, :, k] = interpolated_k
   end
   
   ##########################
   # Take the azimuthal FFT #
   ##########################

   #Get the azimuthal wavenumbers of the FFT
   if isnothing(kφs)
      kφs = φWavenumbers(Nφ)
   end
   
   #Construct an uninitialized Array to store Fourier-transformed data
   simFieldφFFT = Array{ComplexF64}(undef, Nr, length(kφs),
                                    length(zCoordsCylGrid))
   
   for k in 1:length(zCoordsCylGrid)
      for i in 1:length(rCoordsCylGrid)
         simFieldφFFT[i, :, k] = rfft(interpolated_data[i, :, k])
      end
   end

   ############################################
   # Average, if specified, and return result #
   ############################################

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

function compute_φFFT_of_sim_Field(simField::AbstractArray, 
                                   xCoordsCartGrid, yCoordsCartGrid, 
                                   zCoordsCartGrid, p::GridParams, 
                                   Nr::Int64, Nφ::Int64;
                                   rCoordsCylGrid = nothing,
                                   φCoordsCylGrid = nothing, kφs = nothing,
                                   rAvg = false, zAvg = false,
                                   returnModulusSquared = true)
   #=
   Compute the Nφ-point azimuthal Fourier transform of 'simField' (which need
    not actually come from a prognostic field of a simulation, but does need to
    be an Array of data from an Oceananigans Field defined on an Oceananigans 
    Grid) by first interpolating to a uniformly-spaced grid in cylindrical 
    coordinates.
   Default returns a 3D Nr x Nkφ x Nz field, but options 'rAvg' and 'zAvg'
    perform r- and z-averaging, respectively, before returning.
   Default also returns squared modulus of each mode.
   Turning off 'returnModulusSquared' returns the complex-valued modes
    themselves.
   =#

   ####################################################
   # Construct a uniformly spaced (Nr x Nφ x Nz) grid #
   ####################################################
   
   #Get the r- and φ-gridpoints, in polar coordinates, of a uniform polar grid
   if (isnothing(rCoordsCylGrid) || isnothing(φCoordsCylGrid))
      rCoordsCylGrid, φCoordsCylGrid = construct_uniform_polar_grid(p.Lr, Nr, Nφ)
   end
   
   #Tile xCartVec and yCartVec along complementary dimensions, such that each
   # point (x[i], y[j]) on the Oceananigans horizontal grid is represented by
   # (xCartGrid[i, j], yCartGrid[i, j]). Then tile zCartVec in both x- and
   # y-directions to match.
   
   xCoordsCartGridTiled = repeat(xCoordsCartGrid, inner = [length(yCoordsCartGrid)])
   yCoordsCartGridTiled = repeat(yCoordsCartGrid, outer = [length(xCoordsCartGrid)])
   zCoordsCartGridTiled = repeat(zCoordsCartGrid, inner = [length(yCoordsCartGrid)], 
                                 outer = [length(xCoordsCartGrid)])

   #Cartesian coordinates of the points on the cylindrically uniform grid
   xCoordsCylGrid, yCoordsCylGrid, zCoordsCylGrid = ComputeCartesianCoords(rCoordsCylGrid, φCoordsCylGrid, zCoordsCartGrid)

   ###############################################################
   # Perform the interpolation of data to the uniform-polar grid #
   ###############################################################

   #Construct an uninitialized Array to store interpolated data
   interpolated_data = Array{Float64}(undef, size(xCoordsCylGrid)[1], 
                                      size(xCoordsCylGrid)[2], length(zCoordsCylGrid))
   
   #Interpolate data at each z-level. Interpolation only needs to be done in 2D
   # (since zCoordsCylGrid == zCoordsCartGrid).
   
   for k in 1:length(zCoordsCylGrid)

      #Create a linear interpolation object without extrapolation
      interp = Interpolations.scale(Interpolations.interpolate(
                              simField[:, :, k],
                              BSpline(Interpolations.Linear())), 
                        (xCoordsCartGrid, yCoordsCartGrid)
                                   )

      #Interpolate 'simField' to CylGrid points
      interpolated_k = interp.(xCoordsCylGrid, yCoordsCylGrid)
      
      interpolated_data[:, :, k] = interpolated_k
   end
   
   ##########################
   # Take the azimuthal FFT #
   ##########################

   #Get the azimuthal wavenumbers of the FFT
   if isnothing(kφs)
      kφs = φWavenumbers(Nφ)
   end
   
   #Construct an uninitialized Array to store Fourier-transformed data
   simFieldφFFT = Array{ComplexF64}(undef, Nr, length(kφs), length(zCoordsCylGrid))
   
   for k in 1:length(zCoordsCylGrid)
      for i in 1:length(rCoordsCylGrid)
         simFieldφFFT[i, :, k] = rfft(interpolated_data[i, :, k])
      end
   end

   ############################################################################
   # Average if specified and return result with singleton dimensions removed #
   ############################################################################

   if returnModulusSquared
   
      if rAvg
         if zAvg
            return abs2.(mean(simFieldφFFT, dims = (1, 3))[1, :, 1])
         end
      end
   
   elseif !returnModulusSquared #Return the complex-valued modes themselves
   
      if rAvg
   
         if zAvg
            return mean(simFieldφFFT, dims = (1, 3))[1, :, 1]
         elseif !zAvg
            return mean(simFieldφFFT, dims = 1)[1, :, :]
         end
      
      elseif (!rAvg && zAvg)
         return mean(simFieldφFFT, dims = 3)[:, :, 1]
   
      else
         return simFieldφFFT
      end
   end
end