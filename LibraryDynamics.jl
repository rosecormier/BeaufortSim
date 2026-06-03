using CSV
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using SpecialFunctions, Tables

function chebyshev_spaced_faces(i, ξ_min, Nξ; ξ_max = 0.0, ξ0 = 0.0)

   i -= 1

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

function save_zC_values(z_grid, d_ML, grid)
   #=
   Save zC values to a csv file, if non-existent for this grid (Chebyshev only).
   =#
   
   if z_grid == "chebyshev"
   
      gridfilepath = joinpath("./Logs", "grid_Nz$(grid.Nz)_MLD$(abs(d_ML)).csv") 

      if !isfile(gridfilepath)
         mkpath(dirname(gridfilepath)) #Make required path
         @views CSV.write(gridfilepath, Tables.table(znodes(grid, Center())),
                          header = false)
      end
   end
end

function TWB_buoyancy_contribution(f, σr, σz, U, yFlat)
   #=
   Contribution to background buoyancy from thermal-wind balance with background
    velocity.
   =#
   
   if !yFlat #Return 3D version of function
      TWB_b̄ = (x, y, z) -> (-(sqrt(2) * f * U * σr * z / (σz^2))
			                        * exp((1/2) - (z/σz)^2)
			                        * (exp(-(x^2 + y^2) / (σr^2)) - 1))
   elseif yFlat #Return 2D (x, z) version of function, evaluated at y = 0
      TWB_b̄ = (x, z) -> (-(sqrt(2) * f * U * σr * z / (σz^2))
                           * exp((1/2) - (z/σz)^2) * (exp(-(x^2) / (σr^2)) - 1))
   end
   return TWB_b̄
end

function N²DoubleTanh(z, parameters)
   #=
   Evaluate double-tanh function for N² at z.
   =#

   g  = parameters.g
   ρ₀ = parameters.ρ₀
   As = parameters.As
   Bs = parameters.Bs
   Cs = parameters.Cs
   Ad = parameters.Ad
   Bd = parameters.Bd
   Cd = parameters.Cd

   N² = -(g / ρ₀) * ((As / Cs) * (sech((z - Bs) / Cs))^2 
                     + (Ad / Cd) * (sech((z - Bd) / Cd))^2)

   return N²   
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

function buoyancy_BCS(σz, constantN²Term, z_top, z_bot, yFlat;
                      parameters = nothing)
   #=
   Compute appropriate (i.e., preserving gradient-continuity) boundary 
    conditions on buoyancy at top and bottom of simulation domain.
   =#

   if σz == "infinity" #Barotropic case
      b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(constantN²Term),
				                               bottom = GradientBoundaryCondition(constantN²Term))
      
   else #Baroclinic case
   
      #Contribution from background-state thermal-wind balance
      @inline TWB_∂b̄∂z(x, y, z, constN²) = constN² .- ((sqrt(2) * f * U * σr / (σz^2))
                                                        * exp((1/2) - (z/σz)^2)
					                                              * (exp(-(x^2 + y^2) / (σr^2)) - 1) 
                                                        * (1 - 2 * (z/σz)^2)
                                                       )

      if isnothing(parameters.additionalN²Top)
      
         if !yFlat
            b̄z_top = (x, y, t) -> TWB_∂b̄∂z(x, y, z_top, constantN²Term)
            b̄z_bot = (x, y, t) -> TWB_∂b̄∂z(x, y, z_bot, constantN²Term)
         elseif yFlat #Note: this case needs to be tested
            b̄z_top = (x, t) -> TWB_∂b̄∂z(x, 0, z_top, constantN²Term)
            b̄z_bot = (x, t) -> TWB_∂b̄∂z(x, 0, z_bot, constantN²Term)
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top),
				                                  bottom = GradientBoundaryCondition(b̄z_bot))

      else

         if !yFlat
            b̄z_top = (x, y, t, p) -> (p.additionalN²Top 
                                       .+ TWB_∂b̄∂z(x, y, z_top, constantN²Term)
                                      )
            b̄z_bot = (x, y, t, p) -> (p.additionalN²Bottom 
                                       .+ TWB_∂b̄∂z(x, y, z_bot, constantN²Term)
                                      )
         elseif yFlat
            b̄z_top = (x, t, p) -> (p.additionalN²Top 
                                    .+ TWB_∂b̄∂z(x, 0, z_top, constantN²Term)
                                   )
            b̄z_bot = (x, t, p) -> (p.additionalN²Bottom 
                                    .+ TWB_∂b̄∂z(x, 0, z_bot, constantN²Term)
                                   )
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top, 
                                                                    parameters = parameters),
				                                  bottom = GradientBoundaryCondition(b̄z_bot, 
                                                                    parameters = parameters))
      end
   end
   return b̄_BCs
end

function integrate_N²_upwards(i, j, k, grid, N²)
   #=
   Integrate discrete N²-values from the surface (z[grid.Nz - 1] = 0) to the 
    depth z[k], producing a CenterField(grid).
   =#

   integral = @views N²[grid.Nz] .* Δzᶜᶜᶜ(i, j, grid.Nz, grid)

   for m in k:grid.Nz:1
      integral += @views N²[m - 1] .* Δzᶜᶜᶜ(i, j, (m - 1), grid)
   end

   return integral
end

function bkgd_buoyancy(f, σr, σz, U;
                       constantN²Term = 0,
                       grid = nothing, 
                       yFlat = false,
                       doubleTanhParams = nothing,
                       N²Data = nothing)

   TWB_b̄ = TWB_buoyancy_contribution(f, σr, σz, U, yFlat)

   if σz == "infinity" #Barotropic case
      if !yFlat
         b̄ = (x, y, z) -> constantN²Term .* z
      elseif yFlat
         b̄ = (x, z) -> constantN²Term .* z
      end
      
   else #Baroclinic case

      if ambientStrat == "constant"
      
         if !yFlat
            b̄ = (x, y, z) -> TWB_b̄(x, y, z) .+ (constantN²Term .* z)      
         elseif yFlat
            b̄ = (x, z) -> TWB_b̄(x, z) .+ (constantN²Term .* z)
         end
         
      elseif ambientStrat == "doubleTanh"

         @inline globalDoubleTanhN²_ccc(i, j, k, grid) = N²DoubleTanh(grid.z.cᵃᵃᶜ[k], doubleTanhParams)

         globalN²_op = KernelFunctionOperation{Center, Center, Center}(globalDoubleTanhN²_ccc, grid)
         
         @compute globalN² = Field(globalN²_op)
      
         function total_b̄_ccc(i, j, k, grid)

            integrated_term = integrate_N²_upwards(i, j, k, grid, globalN²)
            TWB_term        = TWB_b̄(grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.z.cᵃᵃᶜ[k])
            linear_term     = constantN²Term .* grid.z.cᵃᵃᶜ[k]
            
            return integrated_term + TWB_term + linear_term
         end
      
         total_b̄_op = KernelFunctionOperation{Center, Center, Center}(total_b̄_ccc, grid)
        
         @compute b̄ = Field(total_b̄_op)
      end
   end
   return b̄
end

function bkgd_velocities(σr, σz, U, yFlat = false)

   if yFlat #2D versions, evaluated at y = 0
      
      ū = (x, z) -> 0
   
      if σz == "infinity" #Barotropic case
         v̄ = (x, z) -> -((sqrt(2) * U * x / σr) * exp((1/2) - (x^2)/(σr^2)))
      else #Baroclinic case
         v̄ = (x, z) -> -((sqrt(2) * U * x / σr) * exp((1/2) - (x^2)/(σr^2) - (z/σz)^2))
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