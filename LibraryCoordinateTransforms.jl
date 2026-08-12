using Oceananigans.AbstractOperations, Oceananigans.Fields

@inline r_cca(i, j, k, grid) = @inbounds sqrt((grid.xᶜᵃᵃ[i])^2 + (grid.yᵃᶜᵃ[j])^2)
@inline r_cfa(i, j, k, grid) = @inbounds sqrt((grid.xᶜᵃᵃ[i])^2 + (grid.yᵃᶠᵃ[j])^2)
@inline r_fca(i, j, k, grid) = @inbounds sqrt((grid.xᶠᵃᵃ[i])^2 + (grid.yᵃᶜᵃ[j])^2)
@inline r_ffa(i, j, k, grid) = @inbounds sqrt((grid.xᶠᵃᵃ[i])^2 + (grid.yᵃᶠᵃ[j])^2)

@inline φ_cca(i, j, k, grid) = @inbounds atan(grid.yᵃᶜᵃ[j], grid.xᶜᵃᵃ[i])
@inline φ_cfa(i, j, k, grid) = @inbounds atan(grid.yᵃᶠᵃ[j], grid.xᶜᵃᵃ[i])
@inline φ_fca(i, j, k, grid) = @inbounds atan(grid.yᵃᶜᵃ[j], grid.xᶠᵃᵃ[i])
@inline φ_ffa(i, j, k, grid) = @inbounds atan(grid.yᵃᶠᵃ[j], grid.xᶠᵃᵃ[i])

function polar_coords_Fields(grid, xLoc, yLoc, zLoc; new = false)
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

   return(r, φ)
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