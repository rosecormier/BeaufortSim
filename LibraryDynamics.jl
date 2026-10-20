include("LibraryCoordinates.jl")
include("Parameters.jl")

using Adapt, CSV, CUDA
using Oceananigans.AbstractOperations
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Random, SpecialFunctions

#Background-state streamfunction
@inline Ψ(r, z, p::GyreParams) = (p.U * p.σr / sqrt(2)) .* (1 .- exp.(-(r./p.σr).^2)) .* exp.(0.5 .- (z./p.σz).^2)

#Double-tanh buoyancy function
@inline buoyancyDoubleTanh(z, p::DoubleTanhStratParams) = @. -(p.gravity / p.ρ₀) * (p.A_s * tanh((z - p.z_s) / p.C_s) + p.A_d * tanh((z - p.z_d) / p.C_d))

#N² corresponding to double-tanh buoyancy function
@inline N²DoubleTanh(z, p::DoubleTanhStratParams) = @. -(p.gravity / p.ρ₀) * ((p.A_s / p.C_s) * (sech((z - p.z_s) / p.C_s))^2 + (p.A_d / p.C_d) * (sech((z - p.z_d) / p.C_d))^2)

function default_xy_tracerBCs(gridParams::GridParams)
   #=
   Default boundary conditions in x- and y- directions for "tracer-located"
    fields (i.e., at cell centers in x and y respectively).
   =#
   
   tracerBCs = Dict{DataType, BoundaryCondition}(Periodic => PeriodicBoundaryCondition(),
                                                 Bounded => GradientBoundaryCondition(0))

   x_BC = tracerBCs[gridParams.topology[1]]
   y_BC = gridParams.topology[2] == Flat ? nothing : tracerBCs[gridParams.topology[2]]
   
   return x_BC, y_BC
end

#Given appropriate derivatives of a thermal-wind-balanced Ψ, these functions
# compute the thermal-wind-balanced zonal and meridional velocity components.
@inline u_TWB_fcc(i, j, k, g, ∂Ψ∂y_fcc_tallField) = @inbounds ∂Ψ∂y_fcc_tallField[i, j, (k + 1)]
@inline v_TWB_cfc(i, j, k, g, ∂Ψ∂x_cfc_tallField) = @inbounds -∂Ψ∂x_cfc_tallField[i, j, (k + 1)]

#Given appropriate derivatives of a thermal-wind-balanced Ψ, these functions
# compute the horizontal divergence of flow velocity and diagnose the balanced
# vertical velocity, respectively.
@inline divH_u_ccf(i, j, k, g, ∂Ψ∂y_fcc_Field, ∂Ψ∂x_cfc_Field) = @inbounds ∂xᶜᶜᶠ(i, j, k, g, u_TWB_fcc, ∂Ψ∂y_fcc_Field) + ∂yᶜᶜᶠ(i, j, k, g, v_TWB_cfc, ∂Ψ∂x_cfc_Field)
@inline w_TWB_ccf(i, j, k, g, ∂Ψ∂y_fcc_Field, ∂Ψ∂x_cfc_Field) = @inbounds -divH_u_ccf(i, j, k, g, ∂Ψ∂y_fcc_Field, ∂Ψ∂x_cfc_Field) * Δzᶜᶜᶠ(i, j, k, g)

#These functions compute the linear, thermal-wind-balanced (i.e., dynamic) and
# double-tanh-stratification contributions to buoyancy, respectively.
@inline b_linear_ccc(i, j, k, g, p::GyreParams) = @inbounds p.N²_far * g.z.cᵃᵃᶜ[k]
@inline b_TWB_ccc(i, j, k, g, p::GyreParams, ∂Ψ∂z_ccc_tallField) = @inbounds p.f * ∂Ψ∂z_ccc_tallField[i, j, (k + 1)]
@inline b_double_tanh_ccc(i, j, k, g, p::DoubleTanhStratParams) = @inbounds buoyancyDoubleTanh(g.z.cᵃᵃᶜ[k], p)

function GridAndDiscretelyBalancedBackgroundFields(gridParams::GridParams,
                                                   gyreParams::GyreParams,
                                                   doubleTanhParams::DoubleTanhStratParams,
                                                   perturbParams::PerturbationParams;
                                                   visualizeΨ::Bool = true)
   #=
   Build a model grid using 'gridParams'.
   Construct discretely thermal-wind-balanced background Fields (buoyancy and
    components of velocity) on this grid by first defining a streamfunction on
    a tall copy of the grid (which has one extra interior gridpoint at each
    vertical boundary, necessary to properly define buoyancy boundary
    conditions).
   Return grid, background Fields, and appropriate (i.e., preserving continuity
    of background-field gradient) boundary conditions on buoyancy at top and
    bottom of simulation domain.
   =#
   
   #########
   # Grids #
   #########
   
   simGrid  = BuildOceananigansRectilinearGrid(gridParams)
   tallGrid = BuildTallCopyOceananigansRectilinearGrid(gridParams)
   
   ################################
   # Compute the background state #
   ################################
   
   @inline Ψ_ccf(i, j, k, g) = @inbounds Ψ(r_cca(i, j, k, g), g.z.cᵃᵃᶠ[k], gyreParams)
   @inline Ψ_ffc(i, j, k, g) = @inbounds Ψ(r_ffa(i, j, k, g), g.z.cᵃᵃᶜ[k], gyreParams)
   
   Ψ_ccf_op = KernelFunctionOperation{Center, Center, Face}(Ψ_ccf, tallGrid)
   Ψ_ffc_op = KernelFunctionOperation{Face, Face, Center}(Ψ_ffc, tallGrid)
   
   #Compute Ψ Fields at necessary locations
   @compute tall_Ψ_ccf_Field = Field(Ψ_ccf_op)
   @compute tall_Ψ_ffc_Field = Field(Ψ_ffc_op)

   #Compute Ψ-derivative Fields at necessary locations to diagnose buoyancy
   @compute tall_∂Ψ∂z_ccc_Field   = Field(∂z(tall_Ψ_ccf_Field))
   @compute tall_∂2Ψ∂z2_ccf_Field = Field(∂z(tall_∂Ψ∂z_ccc_Field))
   
   if visualizeΨ
      visualize_Ψ_and_z_derivs(simGrid, tall_Ψ_ccf_Field, tall_∂Ψ∂z_ccc_Field,
                               tall_∂2Ψ∂z2_ccf_Field)
   end
   
   if gridParams.architecture == GPU()

      #On GPU, must compute boundary values using CUDA.jl functionality

      ∂2Ψ∂z2_top = CuArray{Float64}(undef, (gridParams.Nx, gridParams.Ny))
      ∂2Ψ∂z2_bot = CuArray{Float64}(undef, (gridParams.Nx, gridParams.Ny))
      
      @views ∂2Ψ∂z2_top[:, :] .= tall_∂2Ψ∂z2_ccf_Field[:, :, 
                                                       (end - gridParams.Hz)]
      @views ∂2Ψ∂z2_bot[:, :] .= tall_∂2Ψ∂z2_ccf_Field[:, :, 2]
      
   elseif gridParams.architecture == CPU()
   
      #On CPU, can compute boundary values by direct scalar indexing
   
      ∂2Ψ∂z2_top = tall_∂2Ψ∂z2_ccf_Field[:, :, (end - gridParams.Hz)]
      ∂2Ψ∂z2_bot = tall_∂2Ψ∂z2_ccf_Field[:, :, 2]
   end
   
   #Compute Ψ-derivative Fields at necessary locations to diagnose velocities
   @compute tall_∂Ψ∂x_cfc_Field = Field(∂x(tall_Ψ_ffc_Field))
   @compute tall_∂Ψ∂y_fcc_Field = Field(∂y(tall_Ψ_ffc_Field))
   
   u_TWB_op = KernelFunctionOperation{Face, Center, Center}(u_TWB_fcc, simGrid, tall_∂Ψ∂y_fcc_Field)
   v_TWB_op = KernelFunctionOperation{Center, Face, Center}(v_TWB_cfc, simGrid, tall_∂Ψ∂x_cfc_Field)
   w_TWB_op = KernelFunctionOperation{Center, Center, Face}(w_TWB_ccf, simGrid, tall_∂Ψ∂y_fcc_Field, tall_∂Ψ∂x_cfc_Field)
   
   b_x_BC, b_y_BC = default_xy_tracerBCs(gridParams)

   if gyreParams.doubleTanhAmbientStrat

      @inline b_total_doubleTanhAmbientStrat_ccc(i, j, k, g) = @inbounds b_linear_ccc(i, j, k, g, gyreParams) + b_TWB_ccc(i, j, k, g, gyreParams, tall_∂Ψ∂z_ccc_Field) + b_double_tanh_ccc(i, j, k, g, doubleTanhParams)
      
      @compute b_total = Field(KernelFunctionOperation{Center, Center, Center}(b_total_doubleTanhAmbientStrat_ccc, simGrid))
      
      N²DoubleTanh_top = @views N²DoubleTanh(simGrid.z.cᵃᵃᶠ[gridParams.Nz + 1], doubleTanhParams)
      N²DoubleTanh_bot = @views N²DoubleTanh(simGrid.z.cᵃᵃᶠ[1], doubleTanhParams)
      
      ∂b∂z_top = @views @. gyreParams.N²_far + gyreParams.f * ∂2Ψ∂z2_top + N²DoubleTanh_top
      ∂b∂z_bot = @views @. gyreParams.N²_far + gyreParams.f * ∂2Ψ∂z2_bot + N²DoubleTanh_bot
      
   else #Defaults to constant ambient stratification
      
      @inline b_total_TWBambientStrat_ccc(i, j, k, g) = @inbounds b_linear_ccc(i, j, k, g, gyreParams) + b_TWB_ccc(i, j, k, g, gyreParams, tall_∂Ψ∂z_ccc_Field)
      
      @compute b_total = Field(KernelFunctionOperation{Center, Center, Center}(b_total_TWBambientStrat_ccc, simGrid))
      
      ∂b∂z_top = @views @. gyreParams.N²_far + gyreParams.f * ∂2Ψ∂z2_top
      ∂b∂z_bot = @views @. gyreParams.N²_far + gyreParams.f * ∂2Ψ∂z2_bot
   end
   
   b_total_BCs = FieldBoundaryConditions(top = GradientBoundaryCondition(∂b∂z_top),
                                         bottom = GradientBoundaryCondition(∂b∂z_bot))

   @compute b_TWB    = Field(KernelFunctionOperation{Center, Center, Center}(b_TWB_ccc, tallGrid, gyreParams, tall_∂Ψ∂z_ccc_Field))
   @compute b_linear = Field(KernelFunctionOperation{Center, Center, Center}(b_linear_ccc, tallGrid, gyreParams))
   
   ###################################
   # Now compute the perturbed state #
   ###################################
   
   randArray = adapt(CuArray, 2 * (rand(Float64, (gridParams.Nx, 
                                                  gridParams.Ny, 
                                                  gridParams.Nz + 3)) 
                                   .- 0.5)
                    )
   
   Random.seed!(perturbParams.seed) #For reproducibility

   @inline Ψ_perturbed_ccf(i, j, k, g) = @inbounds ((perturbParams.max_ψ′ * randArray[i, j, k]  + 1) * Ψ(r_cca(i, j, k, g), g.z.cᵃᵃᶠ[k], gyreParams))
   
   @compute tall_Ψ_perturbed_ccf_Field = Field(KernelFunctionOperation{Center, Center, Face}(Ψ_perturbed_ccf, tallGrid))
   
   #Interpolate...done in steps for GPU compatibility...
   @inline Ψ_perturbed_ccc(i, j, k, g) = @inbounds ℑzᵃᵃᶜ(i, j, k, g, Ψ_perturbed_ccf)
   @inline Ψ_perturbed_cfc(i, j, k, g) = @inbounds ℑyᵃᶠᵃ(i, j, k, g, Ψ_perturbed_ccc)
   @inline Ψ_perturbed_ffc(i, j, k, g) = @inbounds ℑxᶠᵃᵃ(i, j, k, g, Ψ_perturbed_cfc) 
   
   @compute tall_Ψ_perturbed_ffc_Field = Field(KernelFunctionOperation{Face, Face, Center}(Ψ_perturbed_ffc, tallGrid))
   
   #Compute Ψ-derivative Fields at necessary locations to diagnose buoyancy
   @compute tall_∂Ψ∂z_perturbed_ccc_Field   = Field(∂z(tall_Ψ_perturbed_ccf_Field))
   @compute tall_∂2Ψ∂z2_perturbed_ccf_Field = Field(∂z(tall_∂Ψ∂z_perturbed_ccc_Field))
   
   #Compute Ψ-derivative Fields at necessary locations to diagnose velocities
   @compute tall_∂Ψ∂x_perturbed_cfc_Field = Field(∂x(tall_Ψ_perturbed_ffc_Field))
   @compute tall_∂Ψ∂y_perturbed_fcc_Field = Field(∂y(tall_Ψ_perturbed_ffc_Field))
   
   #Operations representing the thermal-wind-balanced, perturbed velocity
   # components.
   u_perturbed_TWB_op = KernelFunctionOperation{Face, Center, Center}(u_TWB_fcc, simGrid, tall_∂Ψ∂y_perturbed_fcc_Field)
   v_perturbed_TWB_op = KernelFunctionOperation{Center, Face, Center}(v_TWB_cfc, simGrid, tall_∂Ψ∂x_perturbed_cfc_Field)
   w_perturbed_TWB_op = KernelFunctionOperation{Center, Center, Face}(w_TWB_ccf, simGrid, tall_∂Ψ∂y_perturbed_fcc_Field, tall_∂Ψ∂x_perturbed_cfc_Field)
   
   if gyreParams.doubleTanhAmbientStrat

      @inline b_perturbed_doubleTanhAmbientStrat_ccc(i, j, k, g) = @inbounds b_linear_ccc(i, j, k, g, gyreParams) + b_TWB_ccc(i, j, k, g, gyreParams, tall_∂Ψ∂z_perturbed_ccc_Field) + b_double_tanh_ccc(i, j, k, g, doubleTanhParams)
      
      @compute b_perturbed = Field(KernelFunctionOperation{Center, Center, Center}(b_perturbed_doubleTanhAmbientStrat_ccc, simGrid))
      
   else #Defaults to constant ambient stratification
      
      @inline b_perturbed_TWBambientStrat_ccc(i, j, k, g) = @inbounds b_linear_ccc(i, j, k, g, gyreParams) + b_TWB_ccc(i, j, k, g, gyreParams, tall_∂Ψ∂z_perturbed_ccc_Field)
      
      @compute b_perturbed = Field(KernelFunctionOperation{Center, Center, Center}(b_perturbed_TWBambientStrat_ccc, simGrid))
   end
   
   @compute return simGrid, (b_total, Field(u_TWB_op), Field(v_TWB_op), Field(w_TWB_op)), b_total_BCs, (Field(b_perturbed), Field(u_perturbed_TWB_op), Field(v_perturbed_TWB_op), Field(w_perturbed_TWB_op))
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