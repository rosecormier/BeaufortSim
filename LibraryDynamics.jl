include("LibraryCoordinateTransforms.jl")

using CSV
using Oceananigans.AbstractOperations
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using SpecialFunctions, Tables

using Adapt, CairoMakie

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

function TWB_b_anon_function(gyreParams; yFlat = false)
   #=
   Return anonymous function to evaluate contribution to background buoyancy 
    from thermal-wind balance with background velocity.
   =#
   
   f, U, σr, σz = gyreParams.f, gyreParams.U, gyreParams.σr, gyreParams.σz
   
   if !yFlat #Return 3D version of function
      TWB_b = (x, y, z) -> -((sqrt(2) * f * U * σr * z / (σz^2))
                             * exp(0.5 - (z/σz)^2)
                             * (1 - exp(-(x^2 + y^2) / (σr^2)))
                            )
   elseif yFlat #Return 2D (x, z) version of function, evaluated at y = 0
      TWB_b = (x, z) -> -((sqrt(2) * f * U * σr * z / (σz^2))
                          * exp(0.5 - (z/σz)^2) * (1 - exp(-(x^2) / (σr^2)))
                         )
   end
   return TWB_b
end

function TWB_b_field(grid, gyreParams; returnAsArray = true, yFlat = false)
   #=
   Compute thermal-wind-balance contribution to background buoyancy.
   Convert to array before returning, if indicated; otherwise, return as a 
    CenterField on 'grid'.
   =#

   TWB_b_function = TWB_b_anon_function(gyreParams; yFlat = yFlat)

   @inline b_ccc(i, j, k, g) = @inbounds TWB_b_function(g.xᶜᵃᵃ[i], g.yᵃᶜᵃ[j],
                                                        g.z.cᵃᵃᶜ[k])

   b_op = KernelFunctionOperation{Center, Center, Center}(b_ccc, grid)
   
   @compute b = Field(b_op)

   if returnAsArray
      return adapt(Array, b)
   elseif !returnAsArray
      return b
   end
end

function TWB_∂b∂z_anon_function(gyreParams; yFlat = false)
   #=
   Return anonymous function to evaluate contribution to background z-derivative
    of buoyancy from thermal-wind balance with background velocity.
   =#
   
   f, U, σr, σz = gyreParams.f, gyreParams.U, gyreParams.σr, gyreParams.σz

   if !yFlat #Return 3D version of function
      TWB_∂b∂z = (x, y, z) -> @. (-(sqrt(2) * f * U * σr / (σz^2))
                                  * exp(0.5 - (z/σz)^2)
                                  * (1 - exp(-(x^2 + y^2) / (σr^2))) 
                                  * (1 - 2 * (z/σz)^2)
                                 )
   elseif yFlat #Return 2D (x, z) version of function, evaluated at y = 0
      TWB_∂b∂z = (x, z) -> @. (-(sqrt(2) * f * U * σr / (σz^2))
                               * exp(0.5 - (z/σz)^2)
                               * (1 - exp(-(x/σr)^2)) 
                               * (1 - 2 * (z/σz)^2)
                              )
   end
   return TWB_∂b∂z
end

function TWB_∂b∂z_field(grid, gyreParams; returnAsArray = true, yFlat = false)
   #=
   Compute thermal-wind-balance contribution to background z-derivative of 
    buoyancy.
   Convert to array before returning, if indicated; otherwise, return as a 
    ZFaceField on 'grid'.
   =#                        

   TWB_∂b∂z_function = TWB_∂b∂z_anon_function(gyreParams; yFlat = yFlat)

   @inline ∂b∂z_ccf(i, j, k, g) = @inbounds TWB_∂b∂z_function(g.xᶜᵃᵃ[i],
                                                              g.yᵃᶜᵃ[j],
                                                              g.z.cᵃᵃᶠ[k])
   
   ∂b∂z_op = KernelFunctionOperation{Center, Center, Face}(∂b∂z_ccf, grid)
                                                    
   @compute ∂b∂z = Field(∂b∂z_op)
   
   if returnAsArray
      return adapt(Array, ∂b∂z)
   elseif !returnAsArray
      return ∂b∂z
   end
end

function buoyancy_BCs(gyreParams, grid, yFlat, ambientStrat; 
                      Hz = 3, doubleTanhParams = nothing, includeDefaultBCs = false)
   #=
   Return appropriate (i.e., preserving gradient-continuity) boundary 
    conditions on buoyancy at top and bottom of simulation domain.
   =#
   
   N²_far = gyreParams.N²_far

   if σz == "infinity" #Barotropic case
      b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(N²_far),
				                               bottom = GradientBoundaryCondition(N²_far))
      
   else #Baroclinic case
   
      TWB_∂b∂z_function = TWB_∂b∂z_anon_function(gyreParams; yFlat = yFlat)

      @inline TWB_plus_const_∂b∂z(x, y, z) = @inbounds N²_far .+ TWB_∂b∂z_function(x, y, z)

      z_top = @view no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ))[end - Hz]
      z_bot = @view no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ))[Hz + 1]

      if ambientStrat == "constant"
      
         if !yFlat
            b̄z_top = (x, y, t) -> TWB_plus_const_∂b∂z(x, y, z_top)
            b̄z_bot = (x, y, t) -> TWB_plus_const_∂b∂z(x, y, z_bot)
         elseif yFlat #Note: this case needs to be tested
            b̄z_top = (x, t) -> TWB_plus_const_∂b∂z(x, 0, z_top)
            b̄z_bot = (x, t) -> TWB_plus_const_∂b∂z(x, 0, z_bot)
         end
         
         if includeDefaultBCs #Return conditions on all 6 boundaries
            b̄_BCs = FieldBoundaryConditions(grid, 
                                             (Center(), Center(), Center()), 
                                             east = PeriodicBoundaryCondition(), 
                                             west = PeriodicBoundaryCondition(), 
                                             north = PeriodicBoundaryCondition(), 
                                             south = PeriodicBoundaryCondition(), 
                                             top = GradientBoundaryCondition(b̄z_top), 
                                             bottom = GradientBoundaryCondition(b̄z_bot))
         else #Return only non-default conditions
            b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top),
                                          bottom = GradientBoundaryCondition(b̄z_bot))
         end

      elseif ambientStrat == "doubleTanh"

         if !yFlat
            b̄z_top = (x, y, t) -> (N²DoubleTanh(z_top, doubleTanhParams)
                                    .+ TWB_plus_const_∂b∂z(x, y, z_top)
                                   )
            b̄z_bot = (x, y, t) -> (N²DoubleTanh(z_bot, doubleTanhParams) 
                                    .+ TWB_plus_const_∂b∂z(x, y, z_bot)
                                   )
         elseif yFlat
            b̄z_top = (x, t) -> (N²DoubleTanh(z_top, doubleTanhParams)
                                 .+ TWB_plus_const_∂b∂z(x, 0, z_top)
                                )
            b̄z_bot = (x, t) -> (N²DoubleTanh(z_bot, doubleTanhParams)
                                 .+ TWB_plus_const_∂b∂z(x, 0, z_bot)
                                )
         end
         
         b̄_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_top),
				                                  bottom = GradientBoundaryCondition(b̄z_bot))
      end
   end
   return b̄_BCs
end

function bkgd_buoyancy(gyreParams, ambientStrat;
                       grid = nothing, yFlat = false, 
                       doubleTanhParams = nothing)

   if σz == "infinity" #Barotropic case
      if !yFlat
         B = (x, y, z) -> gyreParams.N²_far * z
      elseif yFlat
         B = (x, z) -> gyreParams.N²_far * z
      end
      
   else #Baroclinic case
      
      TWB_b_function = TWB_b_anon_function(gyreParams; yFlat = yFlat)

      if ambientStrat == "constant"
      
         if !yFlat
            B = (x, y, z) -> TWB_b_function(x, y, z) + (gyreParams.N²_far * z)
         elseif yFlat
            B = (x, z) -> TWB_b_function(x, z) + (gyreParams.N²_far * z)
         end
         
      elseif ambientStrat == "doubleTanh"
      
         function B_ccc(i, j, k, g)

            doubleTanh_term = buoyancyDoubleTanh(g.z.cᵃᵃᶜ[k], doubleTanhParams)
            TWB_term        = TWB_b_function(g.xᶜᵃᵃ[i], g.yᵃᶜᵃ[j], g.z.cᵃᵃᶜ[k])
            linear_term     = gyreParams.N²_far * g.z.cᵃᵃᶜ[k]
            
            return doubleTanh_term + TWB_term + linear_term
         end
      
         B_op = KernelFunctionOperation{Center, Center, Center}(B_ccc, grid)
        
         @compute B = Field(B_op)
      end
   end
   return B
end

function discrete_Cartesian_TWB_ICs(simGrid, tallGrid, gyreParams, 
                                    cylindrical_Ψ_anon_function, ambientStrat;
                                    Hz = 3, includeDefaultBCs = false)

   cylindrical_Ψ_function = cylindrical_Ψ_anon_function(gyreParams)

   @inline Ψ_ccc(i, j, k, g) = @inbounds cylindrical_Ψ_function(r_cca(i, j, k, g), g.z.cᵃᵃᶜ[k])
   @inline Ψ_ccf(i, j, k, g) = @inbounds cylindrical_Ψ_function(r_cca(i, j, k, g), g.z.cᵃᵃᶠ[k])
   @inline Ψ_ffc(i, j, k, g) = @inbounds cylindrical_Ψ_function(r_ffa(i, j, k, g), g.z.cᵃᵃᶜ[k])
   
   Ψ_ccc_op = KernelFunctionOperation{Center, Center, Center}(Ψ_ccc, tallGrid)
   Ψ_ccf_op = KernelFunctionOperation{Center, Center, Face}(Ψ_ccf, tallGrid)
   Ψ_ffc_op = KernelFunctionOperation{Face, Face, Center}(Ψ_ffc, tallGrid)
   
   #Compute Ψ Fields at necessary locations
   @compute tempTall_Ψ_ccc_Field = Field(Ψ_ccc_op)
   @compute tempTall_Ψ_ccf_Field = Field(Ψ_ccf_op)
   @compute tempTall_Ψ_ffc_Field = Field(Ψ_ffc_op)
   
   #Compute Ψ-derivative Fields at necessary locations
   @compute tempTall_∂Ψ∂x_cfc_Field    = Field(∂x(tempTall_Ψ_ffc_Field))
   @compute tempTall_∂Ψ∂y_fcc_Field    = Field(∂y(tempTall_Ψ_ffc_Field))
   @compute tempTall_∂Ψ∂z_ccc_Field    = Field(∂z(tempTall_Ψ_ccf_Field))
   #@compute tempTall_∂2Ψ∂z∂x_cff_Field = Field(∂z(tempTall_∂Ψ∂x_cfc_Field))
   #@compute tempTall_∂2Ψ∂z∂y_fcf_Field = Field(∂z(tempTall_∂Ψ∂y_fcc_Field))
   @compute tempTall_∂2Ψ∂z2_ccf_Field  = Field(∂z(tempTall_∂Ψ∂z_ccc_Field))
   
   #We now need to impose constancy of the first derivatives of Ψ across the top
   # and bottom boundaries (of 'simGrid').
   #It's not actually important that the Fields involving second derivatives of
   # Ψ satisfy these boundary conditions, because we never use the second-
   # derivative Fields directly, but they do need to be applied to the first-
   # derivative Fields before we can diagnose the thermal-wind-balanced buoyancy
   # and velocity.
   
   ###
   @compute tempTall_∂Ψ∂x_fcf_Field    = Field(∂x(tempTall_Ψ_ccf_Field))
   @compute tempTall_∂Ψ∂y_cff_Field    = Field(∂y(tempTall_Ψ_ccf_Field))
   
   @compute tempTall_∂2Ψ∂z∂x_fcc_Field = Field(∂z(tempTall_∂Ψ∂x_fcf_Field))
   @compute tempTall_∂2Ψ∂z∂y_cfc_Field = Field(∂z(tempTall_∂Ψ∂y_cff_Field))
   ###
   
   ∂2Ψ∂z∂x_east = @view tempTall_∂2Ψ∂z∂x_fcc_Field[(end - Hz), :, :]
   ∂2Ψ∂z∂x_west = @view tempTall_∂2Ψ∂z∂x_fcc_Field[2, :, :]

   ∂2Ψ∂z∂y_north = @view tempTall_∂2Ψ∂z∂y_cfc_Field[:, (end - Hz), :]
   ∂2Ψ∂z∂y_south = @view tempTall_∂2Ψ∂z∂y_cfc_Field[:, 2, :]
   
   ∂2Ψ∂z2_top = @view tempTall_∂2Ψ∂z2_ccf_Field[:, :, (end - Hz)]
   ∂2Ψ∂z2_bot = @view tempTall_∂2Ψ∂z2_ccf_Field[:, :, 2]

   fig = Figure(size = (1400, 700))
   ax0 = Axis(fig[1, 1], xlabel = "Psi", ylabel = L"$z$ [m]")
   ax1 = Axis(fig[1, 2], xlabel = "dPsi/dz", ylabel = L"$z$ [m]")
   ax2 = Axis(fig[1, 3], xlabel = "d2Psi/dz2", ylabel = L"$z$ [m]")

   lines!(ax0, no_offset_view(adapt(Array, tempTall_Ψ_ccf_Field[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), label = "Ψ_ccf_Field")
   scatter!(ax0, no_offset_view(adapt(Array, tempTall_Ψ_ccf_Field[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), label = "Ψ_ccf_Field (at gridpoints)")
   
   lines!(ax1, no_offset_view(adapt(Array, tempTall_∂Ψ∂z_ccc_Field[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶜ)), label = "∂Ψ∂z_ccc_Field")
   scatter!(ax1, no_offset_view(adapt(Array, tempTall_∂Ψ∂z_ccc_Field[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶜ)), label = "∂Ψ∂z_ccc_Field (at gridpoints)")

   lines!(ax2, no_offset_view(adapt(Array, tempTall_∂2Ψ∂z2_ccf_Field[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), label = "∂2Ψ∂z2_ccf_Field")
   scatter!(ax2, no_offset_view(adapt(Array, tempTall_∂2Ψ∂z2_ccf_Field[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), label = "∂2Ψ∂z2_ccf_Field (at gridpoints)")
   lines!(ax2, no_offset_view(adapt(Array, Field(∂z(tempTall_∂Ψ∂z_ccc_Field))[28, 28, :])), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)))
   
   save(joinpath("./Plots", "testPsi.png"), fig)
   
   @inline u_TWB_fcc(i, j, k, g) = @inbounds tempTall_∂Ψ∂y_fcc_Field[i, j, (k + 1)]
   @inline v_TWB_cfc(i, j, k, g) = @inbounds -tempTall_∂Ψ∂x_cfc_Field[i, j, (k + 1)]
   
   u_TWB_op = KernelFunctionOperation{Face, Center, Center}(u_TWB_fcc, simGrid)
   v_TWB_op = KernelFunctionOperation{Center, Face, Center}(v_TWB_cfc, simGrid)

   @compute u_TWB = Field(u_TWB_op)
   @compute v_TWB = Field(v_TWB_op)
   
   @compute ∂u∂x_TWB = Field(∂x(u_TWB))
   @compute ∂v∂y_TWB = Field(∂y(v_TWB))

   @inline ∂w∂z_TWB_ccc(i, j, k, g) = @inbounds -∂u∂x_TWB[i, j, k] - ∂v∂y_TWB[i, j, k]
   @inline w_TWB_ccf(i, j, k, g) = @inbounds ∂w∂z_TWB_ccc(i, j, k, g) * Δzᶜᶜᶜ(i, j, k, g)
   
   w_TWB_op = KernelFunctionOperation{Center, Center, Face}(w_TWB_ccf, simGrid)
   
   @compute w_TWB = Field(w_TWB_op)
   
   print(adapt(Array, w_TWB)[25, :, 1])
   
   @inline b_linear_ccc(i, j, k, g) = @inbounds gyreParams.N²_far * g.z.cᵃᵃᶜ[k]
   @inline b_TWB_ccc(i, j, k, g) = @inbounds gyreParams.f * tempTall_∂Ψ∂z_ccc_Field[i, j, (k + 1)]

   if ambientStrat == "constant"
   
      b_linear_op = KernelFunctionOperation{Center, Center, Center}(b_linear_ccc, simGrid)
      b_TWB_op    = KernelFunctionOperation{Center, Center, Center}(b_TWB_ccc, simGrid)
      
      @compute b_total = Field(b_linear_op + b_TWB_op)

      ∂b∂z_TWB_top   = @. gyreParams.f * ∂2Ψ∂z2_top + gyreParams.N²_far
      ∂b∂z_TWB_bot   = @. gyreParams.f * ∂2Ψ∂z2_bot + gyreParams.N²_far
      ∂b∂x_TWB_east  = @. gyreParams.f * ∂2Ψ∂z∂x_east
      ∂b∂x_TWB_west  = @. gyreParams.f * ∂2Ψ∂z∂x_west
      ∂b∂y_TWB_north = @. gyreParams.f * ∂2Ψ∂z∂y_north
      ∂b∂y_TWB_south = @. gyreParams.f * ∂2Ψ∂z∂y_south
   end
   
   if includeDefaultBCs #Return conditions, even defaults, on all boundaries
   
      b_TWB_BCs = FieldBoundaryConditions(grid, 
                                          (Center(), Center(), Center()), 
                                          east = PeriodicBoundaryCondition(), 
                                          west = PeriodicBoundaryCondition(), 
                                          north = PeriodicBoundaryCondition(), 
                                          south = PeriodicBoundaryCondition(), 
                                          top = GradientBoundaryCondition(∂b∂z_TWB_top), 
                                          bottom = GradientBoundaryCondition(∂b∂z_TWB_bot))
                                          
      return b_TWB_BCs
      
   else #Return the prognostic variables and only the non-default BCs
   
      b_TWB_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(∂b∂z_TWB_top),
                                          bottom = GradientBoundaryCondition(∂b∂z_TWB_bot))
      #u_TWB_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(∂2Ψ∂z∂y_top), 
      #                                    bottom = GradientBoundaryCondition(∂2Ψ∂z∂y_bot))
      #v_TWB_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(-∂2Ψ∂z∂x_top),
      #                                    bottom = GradientBoundaryCondition(-∂2Ψ∂z∂x_bot))
      
      return b_total, u_TWB, v_TWB, w_TWB, b_TWB_BCs #, u_TWB_BCs, v_TWB_BCs
   end
end

function bkgd_velocities(gyreParams; yFlat = false)

   U, σr, σz = gyreParams.U, gyreParams.σr, gyreParams.σz

   if yFlat #2D versions, evaluated at y = 0
      
      ū= (x, z) -> 0
   
      if σz == "infinity" #Barotropic case
         v̄ = (x, z) -> -((sqrt(2) * U * x / σr) * exp(0.5 - (x/σr)^2))
      else #Baroclinic case
         v̄ = (x, z) -> -((sqrt(2) * U * x / σr) 
                          * exp(0.5 - (x/σr)^2 - (z/σz)^2))
      end
   
   elseif !yFlat #3D versions
   
      if σz == "infinity" #Barotropic case
         ū = (x, y, z) -> ((sqrt(2) * U * y / σr)
                            * exp(0.5 - (x^2 + y^2)/(σr^2)))
         v̄ = (x, y, z) -> -((sqrt(2) * U * x / σr)
                             * exp(0.5 - (x^2 + y^2)/(σr^2)))
      else #Baroclinic case
         ū = (x, y, z) -> ((sqrt(2) * U * y / σr)
                            * exp(0.5 - (x^2 + y^2)/(σr^2) - (z/σz)^2))
         v̄ = (x, y, z) -> -((sqrt(2) * U * x / σr)
                            * exp(0.5 - (x^2 + y^2)/(σr^2) - (z/σz)^2))
      end
   end
   
   return ū, v̄
end

function bkgd_B_cylindrical_coords(gyreParams, doubleTanhParams, 
                                   ambientStrat)
   #=
   Return anonymous function to evaluate B at cylindrical coords (r, z).
   =#

   TWB_b_function = TWB_b_anon_function(gyreParams; yFlat = true)

   if ambientStrat == "constant"
      B = (r, z) -> TWB_b_function(r, z) + (gyreParams.N²_far .* z)
   elseif ambientStrat == "doubleTanh"
      B = (r, z) -> (TWB_b_function(r, z) 
                     + buoyancyDoubleTanh(z, doubleTanhParams) 
                     .+ (gyreParams.N²_far .* z)
                    )
   end
   
   return B
end

function bkgd_Q_cylindrical_coords(gyreParams, doubleTanhParams, 
                                   ambientStrat)
   #=
   Return anonymous function to evaluate Q at cylindrical coords (r, z).
   THIS FUNCTION IS WRONG AND NEEDS TO BE UPDATED.
   =#

   TWB_N²_cyl_coords_function = TWB_∂b∂z_anon_function(gyreParams; yFlat = true)
   
   function ∂N²∂z_doubleTanh(z)
   
      g   = doubleTanhParams.g
      ρ₀  = doubleTanhParams.ρ₀
      A_s = doubleTanhParams.A_s
      C_s = doubleTanhParams.C_s
      z_s = doubleTanhParams.z_s
      A_d = doubleTanhParams.A_d
      C_d = doubleTanhParams.C_d
      z_d = doubleTanhParams.z_d
   
      ∂N²∂z = @. ((2 * g / ρ₀) 
                  * ((A_s / C_s^2) * (sech((z - z_s) / C_s))^2 
                     * tanh((z - z_s) / C_s)
                     + (A_d / C_d^2) * (sech((z - z_d) / C_d))^2 
                     * tanh((z - z_d) / C_d)
                    )
                 )
   end
   
   f, U   = gyreParams.f, gyreParams.U
   σr, σz = gyreParams.σr, gyreParams.σz
   
   ∂N²∂z_TWB = (r, z) -> (-((sqrt(2) * f * U * σr / (σz^2))
                              .* (1 .- exp.(-(r./σr).^2)))
                              .* exp(0.5 .- (z./σz).^2)
                              .* ((2 / σz^2) .* z .* (2 .* (z./σz).^2 .- 3))
                         )
   
   
   if ambientStrat == "constant"
   
      N2    = (r, z) -> gyreParams.N²_far ##TWB_N²_cyl_coords_function(r, z) + gyreParams.N²_far
      ∂N2∂z = (r, z) -> 0##∂N²∂z_TWB(r, z)
      
   elseif ambientStrat == "doubleTanh"
   
      N2    = (r, z) -> (TWB_N²_cyl_coords_function(r, z) 
                         + N²DoubleTanh(z, doubleTanhParams) 
                         .+ gyreParams.N²_far
                        )
      ∂N2∂z = (r, z) -> ∂N²∂z_TWB(r, z) + ∂N²∂z_doubleTanh(z)
   end
   
   Q = (r, z) -> (1/f) * ((sqrt(8) * U / σr) 
                             .* exp.(0.5 .- (z./σz).^2)
                             .* ((1 .- (r./σr).^2) .* exp.(-(r./σr).^2)
                                 + (f^2 / (2 * N2(r, z) * σz^2)) 
                                    .* (2 .* (z./σz).^2 .- 1) 
                                    .* (1 .- exp.(-(r./σr).^2))
                                 )
                             )
                             
   return Q
end

function bkgd_Uφ_cylindrical_coords(gyreParams)
   #=
   Return anonymous function to evaluate Uφ at cylindrical coords (r, z). If r 
    is provided as a negative value, it will be treated as abs(r).
   =#
   
   U, σr, σz = gyreParams.U, gyreParams.σr, gyreParams.σz

   Uφ = (r, z) -> (-(sqrt(2) * U / σr) .* abs.(r) .* exp.(0.5 .- (r./σr).^2) 
                   * exp.(-(z./σz).^2))

   return Uφ
end

function bkgd_Ψ_cylindrical_coords(gyreParams)
   #=
   Return anonymous function to evaluate Ψ at cylindrical coords (r, z).
   =#
   
   U, σr, σz = gyreParams.U, gyreParams.σr, gyreParams.σz
   
   Ψ = (r, z) -> ((U * σr / sqrt(2)) .* (1 .- exp.(-(r./σr).^2)) 
                   .* exp.(0.5 .- (z./σz).^2))

   return Ψ
end

function compute_Q_QG_Cartesian(grid, gyreParams, Ux, Uy; Hz = 3)
   #=
   Return background-state QG potential vorticity, computed at cell centres from
    Cartesian components of velocity.
   Equivalent to cylindrical formulation.
   =#

   rcoords = CenterField(grid)
   Ψ       = CenterField(grid)
   
   set!(rcoords, polar_coords_Fields(grid, "c", "c", "c")[1])

   Ψ_function = bkgd_Ψ_cylindrical_coords(gyreParams)
   
   zcoords = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ)
                           )[(Hz + 1):(length(grid.z.cᵃᵃᶜ) - Hz)]
   
   set!(Ψ, Ψ_function(rcoords, 
                      reshape(zcoords, 1, 1, length(zcoords))
                     )
       )
       
   f, N²_far = gyreParams.f, gyreParams.N²_far

   ∂z2Ψ = ∂z(∂z(Ψ))
   ζa   = ∂x(Uy) - ∂y(Ux) + f

   return (ζa + f^2 * ∂z2Ψ / N²_far) / f
end

function compute_Q_QG_cylindrical(grid, gyreParams, Uφ, φcoords; Hz = 3)
   #=
   Return background-state QG potential vorticity, computed at cell centres from
    cylindrical components of velocity.
   Equivalent to Cartesian formulation.
   =#

   rcoords = CenterField(grid)
   Ψ       = CenterField(grid)
   
   set!(rcoords, polar_coords_Fields(grid, "c", "c", "c")[1])

   Ψ_function = bkgd_Ψ_cylindrical_coords(gyreParams)
   
   zcoords = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ))
   
   set!(Ψ, Ψ_function(rcoords,
                      reshape(zcoords[Hz:(length(zcoords) - Hz - 1)],
                              1, 1, (length(zcoords) - 2 * Hz)
                             )
                     )
       )
       
   f, N²_far = gyreParams.f, gyreParams.N²_far

   ∂z2Ψ = ∂z(∂z(Ψ))
   ζa   = (-(1/rcoords) * (cos(φcoords) * ∂x(rcoords * Uφ) 
           + sin(φcoords) * ∂y(rcoords * Uφ)) + f)

   return (ζa + f^2 * ∂z2Ψ / N²_far) / f
end

function compute_Q_Ertel_Cartesian(grid, gyreParams, Ux, Uy, Uz, B)
   #=
   Return background-state Ertel potential vorticity, computed at cell centres
    from Cartesian components of velocity.
   =#
   
   ωx_initial, ωy_initial, ωz_initial = CenterFields_ω(grid, Ux, Uy, Uz)

   f, N²_far = gyreParams.f, gyreParams.N²_far

   return (ωx_initial * ∂x(B) + ωy_initial * ∂y(B) 
           + (f + ωz_initial) * ∂z(B)) / (N²_far * f)
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