using CSV
using Oceananigans
using Oceananigans.AbstractOperations
using Oceananigans.Fields
using Tables

##############################################
# FUNCTIONS TO BUILD/SAVE OCEANANIGANS GRIDS #
##############################################

function chebyshev_spaced_faces(i, ξ_min, Nξ; ξ_max = 0.0, ξ0 = 0.0)

   Lξ = ξ_max - ξ_min
   
   pi_shift = asin(1 + (ξ0 / Lξ))

   N_below_ξ0 = ((Nξ + 1) * pi) / (2 * (pi - pi_shift)) 

   if i == 1
      i_face = ξ_min
   elseif 1 < i <= N_below_ξ0
      i_face = ξ0 + Lξ * (sin((pi - pi_shift) * (i - 1) / Nξ) - 1)
   elseif i > N_below_ξ0
      i_face = ξ0 - Lξ * (sin((pi - pi_shift) * (i - 1) / Nξ) - 1)
   end

   return i_face
end

function build_custom_z_grid(z_grid_type, Lz, Nz) 

   custom_z_grids = Dict("uniform" => (-Lz, 0),
                         "chebyshev" => k -> chebyshev_spaced_faces(k, -Lz, 
                                                                    Nz + 1)
                        )
   
   return custom_z_grids[z_grid_type]
end

function build_Oceananigans_RectilinearGrid(gridParams)
   #=
   Instantiate an Oceananigans RectilinearGrid that directly uses 'gridParams'.
   Should be used as, e.g., the model grid for a simulation.
   =#
   
   Nx, Ny, Nz = gridParams.Nx, gridParams.Ny, gridParams.Nz
   Lr, Lz     = gridParams.Lr, gridParams.Lz
   
   grid = RectilinearGrid(gridParams.architecture,
                          topology = gridParams.topology,
                          size = (Nx, Ny, Nz), 
                          x = (-Lr, Lr), 
                          y = (-Lr, Lr), 
                          z = build_custom_z_grid(gridParams.z_grid_type, 
                                                  Lz, Nz),
                          halo = (gridParams.Hx, gridParams.Hy, gridParams.Hz)
                         )
end

function build_tall_copy_Oceananigans_RectilinearGrid(gridParams)
   #=
   Instantiate an Oceananigans RectilinearGrid that is in every way identical
    to 'build_Oceananigans_RectilinearGrid(gridParams)' except that one halo
    point at each vertical boundary is converted to an interior point.
   I.e., the tall copy has vertical size 'gridParams.Nz + 2' and vertical halo
    size 'gridParams.Hz - 1'.
   =#
   
   regularGrid = build_Oceananigans_RectilinearGrid(gridParams)

   Nx, Ny, Nz = gridParams.Nx, gridParams.Ny, gridParams.Nz
   Hx, Hy, Hz = gridParams.Hx, gridParams.Hy, gridParams.Hz
   Lr, Lz     = gridParams.Lr, gridParams.Lz

   #Get zFace-coordinates for tall grid
   zF = @view no_offset_view(adapt(Array, regularGrid.z.cᵃᵃᶠ))[Hz:(end - (Hz - 1))]

   tallGrid = RectilinearGrid(gridParams.architecture, 
                              topology = gridParams.topology, 
                              size = (Nx, Ny, Nz + 2), 
                              x = (-Lr, Lr), 
                              y = (-Lr, Lr), 
                              z = zF, 
                              halo = (Hx, Hy, Hz - 1)
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
   Given cylindrical (r, φ, z) coordinates of a point, compute its Cartesian
    (x, y, z) coordinates.
   =#
   
   x = r * cos(φ)
   y = r * sin(φ)

   return x, y, z
end

function xy_vector_to_rφ(vx, vy, grid, useGPU)
   #=
   Given Cartesian (x, y) components of a vector, stored at 
    {Face, Center, Center} and {Center, Face, Center} respectively, 
    compute the vector's projection to polar (r, φ) coordinates.
   Return both polar components at {Center, Center, Center} on 'grid'.
   =#

   r, φ = polar_coords_Fields(grid, "c", "c", "c")

   if useGPU
   
      @inline interpolate_vx_to_ccc(i, j, k, grid) = @inbounds interpolate((grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.z.cᵃᵃᶜ[k]), vx, (Face(), Center(), Center()), grid)

      @inline interpolate_vy_to_ccc(i, j, k, grid) = @inbounds interpolate((grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.z.cᵃᵃᶜ[k]), vy, (Center(), Face(), Center()), grid)
   
      interpolate_vx_to_ccc_op = KernelFunctionOperation{Center, Center, Center}(interpolate_vx_to_ccc, grid)
      interpolate_vy_to_ccc_op = KernelFunctionOperation{Center, Center, Center}(interpolate_vy_to_ccc, grid)
      
      @compute vx_ccc_vals = Field(interpolate_vx_to_ccc_op)
      @compute vy_ccc_vals = Field(interpolate_vx_to_ccc_op)

      vr_binaryOp = (vx_ccc_vals * cos(φ)) + (vy_ccc_vals * sin(φ))
      vφ_binaryOp = (vy_ccc_vals * cos(φ)) - (vx_ccc_vals * sin(φ))
      
      @compute vr = Field(vr_binaryOp)
      @compute vφ = Field(vφ_binaryOp)
      
   elseif !useGPU

      @compute vr = @at (Center, Center, Center) (vx * cos(φ)) + (vy * sin(φ))
      @compute vφ = @at (Center, Center, Center) (vy * cos(φ)) - (vx * sin(φ))
   end

   return vr, vφ
end