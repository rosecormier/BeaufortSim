using CSV
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using SpecialFunctions, Tables

function chebyshev_spaced_faces(i, ξ_min, Nξ; ξ_max = 0.0, ξ0 = 0.0)

   Lξ = ξ_max - ξ_min
   
   pi_shift = asin(1 + (ξ0/Lξ))

   N_below_ξ0 = (Nξ*pi) / (2*(pi-pi_shift)) 

   if i <= N_below_ξ0
      i_face = ξ0 + Lξ * (sin((pi-pi_shift)*i/Nξ) - 1)
   elseif i > N_below_ξ0
      i_face = ξ0 - Lξ * (sin((pi-pi_shift)*i/Nξ) - 1)
   end

   return i_face
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

function N²DoubleTanh(z, parameters)
   #=
   Evaluate, at z, N² corresponding to double-tanh function.
   =#

   g   = parameters.g
   ρ₀  = parameters.ρ₀
   A_s = parameters.A_s
   C_s = parameters.C_s
   z_s = parameters.z_s
   A_d = parameters.A_d
   C_d = parameters.C_d
   z_d = parameters.z_d

   N² = @. -(g / ρ₀) * ((A_s / C_s) * (sech((z - z_s) / C_s))^2 
                        + (A_d / C_d) * (sech((z - z_d) / C_d))^2)

   return N²   
end

function buoyancyDoubleTanh(z, parameters)
   #=
   Evaluate, at z, buoyancy field from double-tanh function.
   =#

   g   = parameters.g
   ρ₀  = parameters.ρ₀
   A_s = parameters.A_s
   C_s = parameters.C_s
   z_s = parameters.z_s
   A_d = parameters.A_d
   C_d = parameters.C_d
   z_d = parameters.z_d
   
   b = @. -(g / ρ₀) * (A_s * tanh((z - z_s) / C_s)
                       + A_d * tanh((z - z_d) / C_d))
   
   return b
end

function TWB_b_field(grid, f, σr, σz, U; returnAsArray = true)

   @inline b_ccc(i, j, k, g) = @inbounds (-(sqrt(2) * f * U * σr * g.z.cᵃᵃᶜ[k]
                                            / (σz^2)
                                           )
			                                     * exp(0.5 - (g.z.cᵃᵃᶜ[k]/σz)^2)
			                                     * (exp(-(g.xᶜᵃᵃ[i]^2 + g.yᵃᶜᵃ[j]^2) 
                                                   / (σr^2)) 
                                              - 1)
                                         )

   b_op = KernelFunctionOperation{Center, Center, Center}(b_ccc, grid)
   
   @compute b = Field(b_op)

   if returnAsArray
      return adapt(Array, b)
   elseif !returnAsArray
      return b
   end
end

function TWB_b_anon_function(f, σr, σz, U, yFlat)
   #=
   Contribution to background buoyancy from thermal-wind balance with background
    velocity.
   =#
   
   if !yFlat #Return 3D version of function
      TWB_b = (x, y, z) -> (-(sqrt(2) * f * U * σr * z / (σz^2))
                             * exp(0.5 - (z/σz)^2)
                             * (exp(-(x^2 + y^2) / (σr^2)) - 1)
                           )
   elseif yFlat #Return 2D (x, z) version of function, evaluated at y = 0
      TWB_b = (x, z) -> (-(sqrt(2) * f * U * σr * z / (σz^2))
                          * exp((1/2) - (z/σz)^2) * (exp(-(x^2) / (σr^2)) - 1)
                        )
   end
   return TWB_b
end

function TWB_∂b∂z_field(grid, N²_far, f, σr, σz, U; returnAsArray = true)

   @inline ∂b∂z_ccc(i, j, k, g) = @inbounds (N²_far 
                                             - ((sqrt(2) * f * U * σr / (σz^2))
                                                * exp(0.5 
                                                      - (g.z.cᵃᵃᶜ[k] / σz)^2)
                                                * (exp(-(g.xᶜᵃᵃ[i]^2 
                                                         + g.yᵃᶜᵃ[j]^2) 
                                                       / (σr^2))
                                                   - 1) 
                                                * (1 - 2 * (g.z.cᵃᵃᶜ[k]/σz)^2)
                                               )
                                            )

   ∂b∂z_op = KernelFunctionOperation{Center, Center, Center}(∂b∂z_ccc, grid)
                                                    
   @compute ∂b∂z = Field(∂b∂z_op)
   
   if returnAsArray
      return adapt(Array, ∂b∂z)
   elseif !returnAsArray
      return ∂b∂z
   end
end

function buoyancy_BCS(f, σr, σz, U, N²_far, grid, yFlat; doubleTanhParams = nothing)
   #=
   Return appropriate (i.e., preserving gradient-continuity) boundary 
    conditions on buoyancy at top and bottom of simulation domain.
   =#

   if σz == "infinity" #Barotropic case
      b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(constantN²Term),
				                               bottom = GradientBoundaryCondition(constantN²Term))
      
   else #Baroclinic case
   
      #Function to compute contribution from background-state thermal-wind balance
      @inline TWB_∂b∂z_function(x, y, z) = @. (N²_far
                                               - ((sqrt(2) * f * U * σr / (σz^2))
                                                  * exp(0.5 - (z/σz)^2)
                                                  * (exp(-(x^2 + y^2) / (σr^2)) - 1) 
                                                  * (1 - 2 * (z/σz)^2)
                                                 )
                                              )

      z_top = @view grid.z.cᵃᵃᶜ[grid.Nz]
      z_bot = @view grid.z.cᵃᵃᶜ[1]

      if ambientStrat == "constant"
      
         if !yFlat
            b̄z_top = (x, y, t) -> TWB_∂b∂z_function(x, y, z_top)
            b̄z_bot = (x, y, t) -> TWB_∂b∂z_function(x, y, z_bot)
         elseif yFlat #Note: this case needs to be tested
            b̄z_top = (x, t) -> TWB_∂b∂z_function(x, 0, z_top)
            b̄z_bot = (x, t) -> TWB_∂b∂z_function(x, 0, z_bot)
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top),
				                                  bottom = GradientBoundaryCondition(b̄z_bot))

      elseif ambientStrat == "doubleTanh"

         if !yFlat
            b̄z_top = (x, y, t) -> (N²DoubleTanh(z_top, doubleTanhParams)
                                    .+ TWB_∂b∂z_function(x, y, z_top)
                                   )
            b̄z_bot = (x, y, t) -> (N²DoubleTanh(z_bot, doubleTanhParams) 
                                    .+ TWB_∂b∂z_function(x, y, z_bot)
                                   )
         elseif yFlat
            b̄z_top = (x, t) -> (N²DoubleTanh(z_top, doubleTanhParams)
                                 .+ TWB_∂b∂z_function(x, 0, z_top)
                                )
            b̄z_bot = (x, t) -> (N²DoubleTanh(z_bot, doubleTanhParams)
                                 .+ TWB_∂b∂z_function(x, 0, z_bot)
                                )
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top),
				                                  bottom = GradientBoundaryCondition(b̄z_bot))
      end
   end
   return b̄_BCs
end

function bkgd_buoyancy(f, σr, σz, U, N²_far;
                       grid = nothing, 
                       yFlat = false,
                       doubleTanhParams = nothing)

   if σz == "infinity" #Barotropic case
      if !yFlat
         B = (x, y, z) -> N²_far .* z
      elseif yFlat
         B = (x, z) -> N²_far .* z
      end
      
   else #Baroclinic case
   
      TWB_b_function = TWB_b_anon_function(f, σr, σz, U, yFlat)

      if ambientStrat == "constant"
      
         if !yFlat
            B = (x, y, z) -> TWB_b_function(x, y, z) .+ (N²_far .* z)
         elseif yFlat
            B = (x, z) -> TWB_b_function(x, z) .+ (N²_far .* z)
         end
         
      elseif ambientStrat == "doubleTanh"
      
         function B_ccc(i, j, k, g)

            doubleTanh_term = buoyancyDoubleTanh(g.z.cᵃᵃᶜ[k], doubleTanhParams)
            TWB_term        = TWB_b_function(g.xᶜᵃᵃ[i], g.yᵃᶜᵃ[j], g.z.cᵃᵃᶜ[k])
            linear_term     = N²_far * g.z.cᵃᵃᶜ[k]
            
            return doubleTanh_term + TWB_term + linear_term
         end
      
         B_op = KernelFunctionOperation{Center, Center, Center}(B_ccc, grid)
        
         @compute B = Field(B_op)
      end
   end
   return B
end

function bkgd_velocities(σr, σz, U; yFlat = false)

   if yFlat #2D versions, evaluated at y = 0
      
      ū = (x, z) -> 0
   
      if σz == "infinity" #Barotropic case
         v̄ = (x, z) -> -((sqrt(2) * U * x / σr) * exp((1/2) - (x/σr)^2))
      else #Baroclinic case
         v̄ = (x, z) -> -((sqrt(2) * U * x / σr) * exp((1/2) - (x/σr)^2 
                                                       - (z/σz)^2))
      end
   
   elseif !yFlat #3D versions
   
      if σz == "infinity" #Barotropic case
         ū = (x, y, z) -> ((sqrt(2) * U * y / σr)
                            * exp((1/2) - (x^2 + y^2)/(σr^2)))
         v̄ = (x, y, z) -> -((sqrt(2) * U * x / σr)
                             * exp((1/2) - (x^2 + y^2)/(σr^2)))
      else #Baroclinic case
         ū = (x, y, z) -> ((sqrt(2) * U * y / σr)
                            * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))
         v̄ = (x, y, z) -> -((sqrt(2) * U * x / σr)
                             * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))
      end
   end
   return ū, v̄
end

function integrate_N²_upwards(i, j, k, grid, N²)
   #=
   Integrate discrete N²-values from the surface (z[grid.Nz - 1] = 0) to the
    depth z[k], producing a CenterField(grid).
   *Needs to be edited/tested*
   =#

   integral = @views N²[grid.Nz] .* Δzᶜᶜᶜ(i, j, grid.Nz, grid)

   for m in k:grid.Nz:1
      integral += @views N²[m - 1] .* Δzᶜᶜᶜ(i, j, (m - 1), grid)
   end

   return integral
end

function N²FromData(k, parameters)
   #=
   Evaluate N² at z[k] from seasonal data in csv file.
   =#

   file = CSV.File(joinpath("./Data/",
                            "N2_Nz$(parameters.Nz)_MLD$(abs(parameters.d_ML)).csv");
                   header = 2)
   N²   = Tuple(file[parameters.season])
   N²_k = N²[k]
   
   return N²_k
end