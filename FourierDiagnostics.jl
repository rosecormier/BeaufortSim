include("LibraryCoordinates.jl")

using AbstractFFTs
using Interpolations
using Oceananigans.Grids
using OffsetArrays: no_offset_view

using CairoMakie

function construct_uniform_polar_grid(Lr, Nr, Nφ)
   #=
   Return Nr evenly spaced gridpoints between 'dr' and 'Lr', and Nφ evenly
    spaced gridpoints between 0 and 2π (without double-counting φ = 0).
   =#

   dr = Lr / (Nr + 1)
   dφ = 2π / Nφ

   r_gridpoints = dr:dr:(Lr - dr)
   φ_gridpoints = 0:dφ:(2π - dφ)

   return r_gridpoints, φ_gridpoints
end

function compute_φFFT_of_sim_data(simField, xLoc, yLoc, zLoc, grid, gridParams, Nr, Nφ; Hx = nothing, Hy = nothing, Hz = nothing)
   #=
   Compute the Nφ-point azimuthal Fourier transform of 'simField' (which need
    not actually be a prognostic field of a simulation, but does need to be
    defined on an Oceananigans Grid) by first interpolating to a uniformly-
    spaced grid in cylindrical coordinates.
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
   
   rCylGrid = collect(rCylGrid)
   
   #Read in the coordinates of Oceananigans gridpoints
   
   if (xLoc == "c" || xLoc == "Center")
      xCartVec = no_offset_view(grid.xᶜᵃᵃ)[(Hx + 1):(end - Hx - 1)]
   elseif (xLoc == "f" || xLoc == "Face")
      xCartVec = no_offset_view(grid.xᶠᵃᵃ)[(Hx + 1):(end - Hx - 1)]
   end
   
   if (yLoc == "c" || yLoc == "Center")
      yCartVec = no_offset_view(grid.yᵃᶜᵃ)[(Hy + 1):(end - Hy - 1)]
   elseif (yLoc == "f" || yLoc == "Face")
      yCartVec = no_offset_view(grid.yᵃᶠᵃ)[(Hy + 1):(end - Hy - 1)]
   end

   #=
   if (zLoc == "c" || zLoc == "Center")
      zCartVec = grid.z.cᵃᵃᶜ[Hz:(end - Hz)]
   elseif (zLoc == "f" || zLoc == "Face")
      zCartVec = grid.z.cᵃᵃᶠ[Hz:(end - Hz)]
   end
   =#
   
   #Tile xCartVec and yCartVec along complementary dimensions, such that each
   # point (x[i], y[j]) on the Oceananigans horizontal grid is represented by
   # (xCartGrid[i, j], yCartGrid[i, j]).
   
   xCartGrid = repeat(xCartVec, inner = [length(yCartVec)])
   yCartGrid = repeat(yCartVec, outer = [length(xCartVec)])

   #Cartesian coordinates of the points on the cylindrically uniform grid
   xCylGrid, yCylGrid, zCylGrid = compute_Cart_coords(rCylGrid,
                                                       φCylGrid,
                                                       1) #zCartGrid)
                                                       
   #Create a linear interpolation (CartGrid -> CylGrid) object without
   # extrapolation.
   #Note that the interpolation only needs to be done in 2D (zCart = zCyl).
   
   interp = Interpolations.scale(Interpolations.interpolate(
                        no_offset_view(simField.data)[(Hx + 1):(end - Hx - 1), 
                                                      (Hy + 1):(end - Hy - 1), 
                                                      1], 
                                          BSpline(Interpolations.Linear())), 
                        (xCartVec, yCartVec)
                                )

   #Interpolate 'simField' to CylGrid points
   interpolated = interp.(xCylGrid, yCylGrid)

   fig = Figure(size=(700, 700))
   ax = Axis(fig[1, 1])
   scatter!(ax, xCartGrid, yCartGrid)
   scatter!(ax, vec(xCylGrid), vec(yCylGrid))
   save("testuniformgrid.png", fig)
   
   cylPoints = Point2f.(vec(xCylGrid), vec(yCylGrid))
   
   fig = Figure(size=(700, 700))
   ax = Axis(fig[1, 1])
   sc = scatter!(ax, cylPoints, colormap = :viridis, color = vec(interpolated), colorrange = (0, 1.5), markersize = 50)
   save("testinterpfield.png", fig)
 
   simFieldφFFT = rfft(interpolated)
   print(simFieldφFFT)
   kφs = rfftfreq(Nφ, (Nφ / 2π))
 
   return interpolated
end

testGrid = RectilinearGrid(CPU(),
                          topology = (Bounded, Bounded, Grids.Flat),
                          size = (12, 12), 
                          x = (-1, 1), 
                          y = (-1, 1),
                          halo = (1, 1)
                         )
                         
testField = CenterField(testGrid)

@inline r2function(x, y) = x^2 + y^2

set!(testField, r2function)

fig = Figure(size=(700, 700))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, view(testField, :, :, 1), colorrange = (0, 1.5))
Colorbar(fig[1, 2], hm)
save("testfield.png", fig)

testInterpolated = compute_φFFT_of_sim_data(testField, "c", "c", "c", testGrid, nothing, 1, 5; Hx = 0, Hy = 0, Hz = 0)