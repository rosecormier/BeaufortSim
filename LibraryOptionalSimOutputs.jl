include("LibraryCoordinates.jl")
include("LibraryVisualization.jl")

################################################################################
# FROM PKE BUDGET: THE INDIVIDUAL TERMS ui′uj′∂jUi (CARTESIAN INDICES) THAT    #
#  CONTRIBUTE TO THE TOTAL BACKGROUND-ENERGY CONVERSION.                       # 
# FOR ANY GENERIC BACKGROUND FLOW, THE SUM OF THESE NINE TERMS SHOULD EQUAL    #
# THE SUM OF THE BTI- AND BCI-CONVERSION TERMS.                                #
################################################################################

@inline ux′_fcc(i, j, k, g, ux, Ux) = @inbounds ux[i, j, k] - Ux[i, j, k]
@inline ux′_ccc(i, j, k, g, ux, Ux) = @inbounds ℑxᶜᵃᵃ(i, j, k, g, ux′_fcc, ux, Ux)
@inline uy′_cfc(i, j, k, g, uy, Uy) = @inbounds uy[i, j, k] - Uy[i, j, k]
@inline uy′_ccc(i, j, k, g, uy, Uy) = @inbounds ℑyᵃᶜᵃ(i, j, k, g, uy′_cfc, uy, Uy)
@inline uz′_ccf(i, j, k, g, uz, Uz) = @inbounds uz[i, j, k] - Uz[i, j, k]
@inline uz′_ccc(i, j, k, g, uz, Uz) = @inbounds ℑzᵃᵃᶜ(i, j, k, g, uz′_ccf, uz, Uz)

function bkgd_transfer_ux′ux′∂xUx(simulation, Ux_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities

   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(ux′_ccc(i, j, k, g, u, Ux_bkgd) * ux′_ccc(i, j, k, g, u, Ux_bkgd) * ∂xᶜᶜᶜ(i, j, k, g, Ux_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_ux′uy′∂yUx(simulation, Ux_bkgd, Uy_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities

   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(ux′_ccc(i, j, k, g, u, Ux_bkgd) * uy′_ccc(i, j, k, g, v, Uy_bkgd) * ∂yᶜᶜᶜ(i, j, k, g, Ux_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_ux′uz′∂zUx(simulation, Ux_bkgd, Uz_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(ux′_ccc(i, j, k, g, u, Ux_bkgd) * uz′_ccc(i, j, k, g, w, Uz_bkgd) * ∂zᶜᶜᶜ(i, j, k, g, Ux_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_uy′ux′∂xUy(simulation, Ux_bkgd, Uy_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(uy′_ccc(i, j, k, g, v, Uy_bkgd) * ux′_ccc(i, j, k, g, u, Ux_bkgd) * ∂xᶜᶜᶜ(i, j, k, g, Uy_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_uy′uy′∂yUy(simulation, Uy_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities

   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(uy′_ccc(i, j, k, g, v, Uy_bkgd) * uy′_ccc(i, j, k, g, v, Uy_bkgd) * ∂yᶜᶜᶜ(i, j, k, g, Uy_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_uy′uz′∂zUy(simulation, Uy_bkgd, Uz_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities

   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(uy′_ccc(i, j, k, g, v, Uy_bkgd) * uz′_ccc(i, j, k, g, w, Uz_bkgd) * ∂zᶜᶜᶜ(i, j, k, g, Uy_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_uz′ux′∂xUz(simulation, Ux_bkgd, Uz_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities

   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(uz′_ccc(i, j, k, g, w, Uz_bkgd) * ux′_ccc(i, j, k, g, u, Ux_bkgd) * ∂xᶜᶜᶜ(i, j, k, g, Uz_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_uz′uy′∂yUz(simulation, Uy_bkgd, Uz_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(uz′_ccc(i, j, k, g, w, Uz_bkgd) * uy′_ccc(i, j, k, g, v, Uy_bkgd) * ∂yᶜᶜᶜ(i, j, k, g, Uz_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

function bkgd_transfer_uz′uz′∂zUz(simulation, Uz_bkgd)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities

   bkgd_transfer_ccc(i, j, k, g) = @inbounds -(uz′_ccc(i, j, k, g, w, Uz_bkgd) * uz′_ccc(i, j, k, g, w, Uz_bkgd) * ∂zᶜᶜᶜ(i, j, k, g, Uz_bkgd))
   bkgd_transfer_op = KernelFunctionOperation{Center, Center, Center}(bkgd_transfer_ccc, grid)
   
   compute!(Integral(Field(bkgd_transfer_op)))
end

######################
# VERTICAL INTEGRALS #
######################

function integrate_field_upwards(i, j, k, grid, ψ, integral) #This guy might be broken
   #=
   Integrate an arbitrary field ψ, defined at (centre, centre, face)-points,
    from the surface (z[grid.Nz - 1] = 0) to the depth z[k].
   =#

   integral_m = ψ[i, j, (grid.Nz + 1)] .* Δzᶜᶜᶠ(i, j, (grid.Nz + 1), grid)

   for m in (grid.Nz - 1):(-1):k
      integral_m += ψ[i, j, (m + 1)] .* Δzᶜᶜᶠ(i, j, (m + 1), grid)
   end

   integral[i, j, k] = integral_m
end

function compute_∫bdz_0_z(simulation)
   #=
   Evalute the integral of b w.r.t. z over the interval [0, z].
   =#

   int_bdz = CenterField(simulation.model.grid)

   integrate_b_upwards_op = KernelFunctionOperation{Center, Center, Center}(
                integrate_field_upwards, simulation.model.grid, simulation.model.tracers.b, int_bdz)

   compute!(integrate_b_upwards_op)
end

########################################
# DERIVATIVES W.R.T. POLAR COORDINATES #
########################################

@inline φ_ccc_vals(grid) = @inbounds polar_coords_Fields(grid, "c", "c", "c")[2]

@inline ∂p∂φ_over_r(grid, p) = -sin(φ_ccc_vals(grid)) * ∂x(p) + cos(φ_ccc_vals(grid)) * ∂y(p)

function compute_∂p∂φ_over_r(simulation)
   #=
   Evaluate (1/r) times the derivative of kinematic pressure deviation (p) 
    w.r.t. φ.
   =#

   compute!(∂p∂φ_over_r(simulation.model.grid, 
			                  (simulation.model.pressures.pNHS + simulation.model.pressures.pHY′)
	       	             )
	         )
end

@inline ∂p∂r(grid, p) = cos(φ_ccc_vals(grid)) * ∂x(p) + sin(φ_ccc_vals(grid)) * ∂y(p)

function compute_∂p∂r(simulation)
   #=
   Evaluate derivative of kinematic pressure deviation (p) w.r.t. r.
   =#

   compute!(∂p∂r(simulation.model.grid, 
                 (simulation.model.pressures.pNHS + simulation.model.pressures.pHY′)
	            	)
	         )
end

@inline ∂φ_∫bdz_over_r(grid, int_bdz) = (-sin(φ_ccc_vals(grid)) * ∂x(int_bdz) 
				         + cos(φ_φ_ccc_vals(grid)) * ∂y(int_bdz))

function compute_∂φ_∫bdz_over_r(simulation)
   #=
   Evaluate (1/r) times the derivative, w.r.t. φ, of the integral of b w.r.t.
    z, where the integral is over the interval [0, z].
   =#

   int_bdz = compute_∫bdz_0_z(simulation)

   compute!(∂φ_∫bdz_over_r(simulation.model.grid, int_bdz))
end

@inline ∂r_∫bdz(grid, int_bdz) = (cos(φ_ccc_vals(grid)) * ∂x(int_bdz) 
				                          + sin(φ_ccc_vals(grid)) * ∂y(int_bdz))

function compute_∂r_∫bdz(simulation)
   #=
   Evaluate the derivative, w.r.t. r, of the integral of b w.r.t. z, where the
    integral is over the interval [0, z].
   =#

   int_bdz = compute_∫bdz_0_z(simulation)

   compute!(∂r_∫bdz(simulation.model.grid, int_bdz))
end