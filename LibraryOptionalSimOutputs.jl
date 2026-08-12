include("LibraryCoordinateTransforms.jl")
include("LibraryVisualization.jl")

function integrate_field_upwards(i, j, k, grid, ψ, integral)
   #=
   Integrate an arbitrary field ψ, defined at (centre, centre, face)-points,
    from the surface (z[grid.Nz - 1] = 0) to the depth z[k]
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