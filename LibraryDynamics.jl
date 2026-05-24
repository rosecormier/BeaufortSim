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
   #Save zC values to a csv file, if non-existent for this grid (Chebyshev only)
   
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
   Contribution to background buoyancy from thermal-wind balance with background velocity.
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

function TWB_∂b∂z_contribution(f, σr, σz, U, yFlat)
   #=
   Contribution to z-derivative of background buoyancy from thermal-wind balance with background velocity derivatives.
   =#
   if !yFlat #Return 3D version of function
      TWB_∂b̄∂z = (x, y, z) -> -((sqrt(2) * f * U * σr / (σz^2))
				                         * exp((1/2) - (z/σz)^2)
					                       * (exp(-(x^2 + y^2) / (σr^2)) - 1) 
                                 * (1 - 2 * (z/σz)^2)
                                )
   elseif yFlat #Return 2D (x, z) version of function, evaluated at y = 0
      TWB_∂b̄∂z = (x, z) -> -((sqrt(2) * f * U * σr / (σz^2))
				                      * exp((1/2) - (z/σz)^2)
					                    * (exp(-(x^2) / (σr^2)) - 1) * (1 - 2 * (z/σz)^2)
                             )
   end
   return TWB_∂b̄∂z
end

function N²_from_data(Nz, d_ML, constantN²Term, season)
   file = CSV.File(joinpath("./Data/", "N2_Nz$(Nz)_MLD$(abs(d_ML)).csv"); header = 2)
   N²   = file[season] .+ constantN²Term
end

function buoyancy_BCS(σz,
                      constantN²Term,
                      z_top, 
                      z_bot, 
                      yFlat;
                      parameters=nothing)

   if σz == "infinity" #Barotropic case
      b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(constantN²Term),
				                               bottom = GradientBoundaryCondition(constantN²Term))
      
   else #Baroclinic case
   
      TWB_∂b̄∂z_function = TWB_∂b∂z_contribution(f, σr, σz, U, yFlat)
      
      if isnothing(parameters.N²FromDataTop)
      
         if !yFlat
            
            #Function to compute z-derivative of background buoyancy
            @inline b̄z(x, y, z) = constantN²Term .+ TWB_∂b̄∂z_function(x, y, z)
         
            #Evaluate z-derivatives at top and bottom of domain
            b̄z_top = (x, y, t) -> b̄z(x, y, z_top)
            b̄z_bot = (x, y, t) -> b̄z(x, y, z_bot)
         
         elseif yFlat
         
            #Function to compute z-derivative of background buoyancy
            @inline b̄z(x, z) = constantN²Term .+ TWB_∂b̄∂z_function(x, z)
            
            #Evaluate z-derivatives at top and bottom of domain
            b̄z_top = (x, t) -> b̄z(x, z_top)
            b̄z_bot = (x, t) -> b̄z(x, z_bot)
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top),
				                                  bottom = GradientBoundaryCondition(b̄z_bot))
     
      else
         
         #z-derivatives of buoyancy at top and bottom of domain
         if !yFlat
            b̄z_top = (x, y, t, p) -> p.N²FromDataTop .+ TWB_∂b̄∂z_function(x, y, z_top)
            b̄z_bot = (x, y, t, p) -> p.N²FromDataBottom .+ TWB_∂b̄∂z_function(x, y, z_bot)
         elseif yFlat
            b̄z_top = (x, t, p) -> p.N²FromDataTop .+ TWB_∂b̄∂z_function(x, z_top)
            b̄z_bot = (x, t, p) -> p.N²FromDataBottom .+ TWB_∂b̄∂z_function(x, z_bot)
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top, 
                                                                    parameters = parameters),
				                                  bottom = GradientBoundaryCondition(b̄z_bot, 
                                                                    parameters = parameters))
      end
   end
   return b̄_BCs
end

function integrate_N²_upwards(i, j, k, grid, N², integral)
   #=
   Integrate discrete N² values from the surface (z[grid.Nz - 1] = 0) to the 
    depth z[k], producing a CenterField(grid).
   =#

   integral_m = N²[grid.Nz] .* Δzᶜᶜᶜ(i, j, grid.Nz, grid)

   for m in k:grid.Nz:1
      integral_m += N²[m - 1] .* Δzᶜᶜᶜ(i, j, (m - 1), grid)
   end

   integral[i, j, k] = integral_m
end

function bkgd_buoyancy(f, σr, σz, U;
                       constantN²Term = 0,
                       N²FromData = nothing, 
                       grid = nothing, 
                       yFlat = false)

   TWB_b̄_function = TWB_buoyancy_contribution(f, σr, σz, U, yFlat)

   if σz == "infinity" #Barotropic case
      if !yFlat
         b̄ = (x, y, z) -> constantN²Term .* z
      elseif yFlat
         b̄ = (x, z) -> constantN²Term .* z
      end
      
   else #Baroclinic case

      if isnothing(N²FromData)
      
         if !yFlat
            b̄ = (x, y, z) -> TWB_b̄_function(x, y, z) .+ (constantN²Term .* z)      
         elseif yFlat
            b̄ = (x, z) -> TWB_b̄_function(x, z) .+ (constantN²Term .* z)
         end
         
      else
      
         b̄ = CenterField(grid)
         
         integrate_N²_upwards_op = KernelFunctionOperation{Center, Center, Center}(
                integrate_N²_upwards, grid, N²FromData, b̄)

         compute!(integrate_N²_upwards_op)
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