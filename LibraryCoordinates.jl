include("Parameters.jl")

using CSV
using Oceananigans
using Oceananigans.AbstractOperations
using Oceananigans.Fields
using OffsetArrays: no_offset_view
using Tables

##############################################
# FUNCTIONS TO BUILD/SAVE OCEANANIGANS GRIDS #
##############################################

function chebyshev_spaced_faces(i, ξ_min, Nξ; ξ_max = 0.0, ξ0 = 0.0)

   Lξ = ξ_max - ξ_min
   
   π_shift = asin(1 + (ξ0 / Lξ))

   N_below_ξ0 = ((Nξ + 1) * π) / (2 * (π - π_shift)) 

   if i == 1
      i_face = ξ_min
   elseif 1 < i <= N_below_ξ0
      i_face = ξ0 + Lξ * (sin((π - π_shift) * (i - 1) / Nξ) - 1)
   elseif i > N_below_ξ0
      i_face = ξ0 - Lξ * (sin((π - π_shift) * (i - 1) / Nξ) - 1)
   end

   return i_face
end

function build_custom_z_grid(z_grid_type, Lz, Nz) 

   custom_z_grids = Dict(:uniform => (-Lz, 0),
                         :chebyshev => k -> chebyshev_spaced_faces(k, -Lz,
                                                                   Nz + 1)
                        )
   
   return custom_z_grids[z_grid_type]
end

function BuildOceananigansRectilinearGrid(p::GridParams)
   #=
   Instantiate an Oceananigans RectilinearGrid that directly uses specified
    GridParams.
   Should be used as, e.g., the model grid for a simulation.
   =#
   
   grid = RectilinearGrid(p.architecture,
                          topology = p.topology,
                          size = (p.Nx, p.Ny, p.Nz), 
                          x = (-p.Lr, p.Lr), 
                          y = (-p.Lr, p.Lr), 
                          z = build_custom_z_grid(p.z_grid_type, p.Lz, p.Nz),
                          halo = (p.Hx, p.Hy, p.Hz)
                         )
end

function BuildTallCopyOceananigansRectilinearGrid(p::GridParams)
   #=
   Instantiate an Oceananigans RectilinearGrid that is in every way identical
    to 'build_Oceananigans_RectilinearGrid(p)' except that one halo point at
    each vertical boundary is converted to an interior point.
   I.e., the tall copy has vertical size 'p.Nz + 2' and vertical halo size
    'p.Hz - 1'.
   =#
   
   regularGrid = BuildOceananigansRectilinearGrid(p)

   #Get zFace-coordinates for tall grid
   zF = @view no_offset_view(adapt(Array, regularGrid.z.cᵃᵃᶠ)
                            )[p.Hz:(end - (p.Hz - 1))]

   tallGrid = RectilinearGrid(p.architecture,
                              topology = p.topology, 
                              size = (p.Nx, p.Ny, p.Nz + 2), 
                              x = (-p.Lr, p.Lr), 
                              y = (-p.Lr, p.Lr), 
                              z = zF, 
                              halo = (p.Hx, p.Hy, p.Hz - 1)
                             )
end

function save_zC_values(z_grid, grid)
   #=
   Save zC values to a csv file, if non-existent for this grid (Chebyshev only).
   =#
   
   if z_grid == "chebyshev"
   
      gridfilepath = joinpath("./Logs", "grid_Nz$(grid.Nz).csv") 

      if !isfile(gridfilepath)
         mkpath(dirname(gridfilepath)) #Make required path
         @views CSV.write(gridfilepath, Tables.table(znodes(grid, Center())),
                          header = false)
      end
   end
end

######################################################################
# FUNCTIONS TO CONVERT BETWEEN CYLINDRICAL AND CARTESIAN COORDINATES #
######################################################################

@inline r_cca(i, j, k, grid) = @inbounds sqrt((grid.xᶜᵃᵃ[i])^2 + (grid.yᵃᶜᵃ[j])^2)
@inline r_cfa(i, j, k, grid) = @inbounds sqrt((grid.xᶜᵃᵃ[i])^2 + (grid.yᵃᶠᵃ[j])^2)
@inline r_fca(i, j, k, grid) = @inbounds sqrt((grid.xᶠᵃᵃ[i])^2 + (grid.yᵃᶜᵃ[j])^2)
@inline r_ffa(i, j, k, grid) = @inbounds sqrt((grid.xᶠᵃᵃ[i])^2 + (grid.yᵃᶠᵃ[j])^2)

@inline φ_cca(i, j, k, grid) = @inbounds atan(grid.yᵃᶜᵃ[j], grid.xᶜᵃᵃ[i])
@inline φ_cfa(i, j, k, grid) = @inbounds atan(grid.yᵃᶠᵃ[j], grid.xᶜᵃᵃ[i])
@inline φ_fca(i, j, k, grid) = @inbounds atan(grid.yᵃᶜᵃ[j], grid.xᶠᵃᵃ[i])
@inline φ_ffa(i, j, k, grid) = @inbounds atan(grid.yᵃᶠᵃ[j], grid.xᶠᵃᵃ[i])

function polar_coords_Fields(grid, xLoc, yLoc, zLoc)
   #=
   Compute polar (r, φ) coordinates for each (x, y)-point at specified 'loc' on
    'grid'; return 'r' and 'φ' as Fields (at appropriate location 'loc') on
    'grid'.
   =#
   
   if (xLoc == "c" || xLoc == "Center")
   
      if (yLoc == "c" || yLoc == "Center")
      
         if (zLoc == "c" || z == "Center")
            r_op = KernelFunctionOperation{Center, Center, Center}(r_cca, grid)
            φ_op = KernelFunctionOperation{Center, Center, Center}(φ_cca, grid)
         elseif (zLoc == "f" || zLoc == "Face")
            r_op = KernelFunctionOperation{Center, Center, Face}(r_cca, grid)
            φ_op = KernelFunctionOperation{Center, Center, Face}(φ_cca, grid)
         end
         
      elseif (yLoc == "f" || yLoc == "Face")
         
         if (zLoc == "c" || zLoc == "Center")
            r_op = KernelFunctionOperation{Center, Face, Center}(r_cfa, grid)
            φ_op = KernelFunctionOperation{Center, Face, Center}(φ_cfa, grid)
         elseif (zLoc == "f" || zLoc == "Face")
            r_op = KernelFunctionOperation{Center, Face, Face}(r_cfa, grid)
            φ_op = KernelFunctionOperation{Center, Face, Face}(φ_cfa, grid)
         end
         
      end
   
   elseif (xLoc == "f" || xLoc == "Face")
   
      if (yLoc == "c" || yLoc == "Center")
      
         if (zLoc == "c" || zLoc == "Center")
            r_op = KernelFunctionOperation{Face, Center, Center}(r_fca, grid)
            φ_op = KernelFunctionOperation{Face, Center, Center}(φ_fca, grid)
         elseif (zLoc == "f" || zLoc == "Face")
            r_op = KernelFunctionOperation{Face, Center, Face}(r_fca, grid)
            φ_op = KernelFunctionOperation{Face, Center, Face}(φ_fca, grid)
         end
         
      elseif (yLoc == "f" || yLoc == "Face")
         
         if (zLoc == "c" || zLoc == "Center")
            r_op = KernelFunctionOperation{Face, Face, Center}(r_ffa, grid)
            φ_op = KernelFunctionOperation{Face, Face, Center}(φ_ffa, grid)
         elseif (zLoc == "f" || zLoc == "Face")
            r_op = KernelFunctionOperation{Face, Face, Face}(r_ffa, grid)
            φ_op = KernelFunctionOperation{Face, Face, Face}(φ_ffa, grid)
         end
      
      end
   end
   
   @compute r = Field(r_op)
   @compute φ = Field(φ_op)

   return r, φ
end

function compute_Cart_coords(r, φ, z)
   #=
   Given cylindrical (r, φ, z) coordinates of point(s), compute representation
   in Cartesian (x, y, z) coordinates.
   =#
   
   x = Array{Float64}(undef, length(r), length(φ))
   y = Array{Float64}(undef, length(r), length(φ))
   
   #If r and φ both 'StepRangeLen', must loop over one of them (can't broadcast)
   for r_idx in 1:1:length(r)
      x[r_idx, :] = r[r_idx] * cos.(φ)
      y[r_idx, :] = r[r_idx] * sin.(φ)
   end

   return x, y, z
end

#r- and φ-components of a vector given in (x, y)-coordinates
@inline vector_r_component_ccc(i, j, k, g, ux, uy) = @inbounds ℑxᶜᵃᵃ(i, j, k, g, ux) * cos(φ_cca(i, j, k, g)) + ℑyᵃᶜᵃ(i, j, k, g, uy) * sin(φ_cca(i, j, k, g))
@inline vector_φ_component_ccc(i, j, k, g, ux, uy) = @inbounds ℑyᵃᶜᵃ(i, j, k, g, uy) * cos(φ_cca(i, j, k, g)) - ℑxᶜᵃᵃ(i, j, k, g, ux) * sin(φ_cca(i, j, k, g))