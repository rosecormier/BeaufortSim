include("LibraryCoordinates.jl")
include("Parameters.jl")

using Random

#Background-state streamfunction
@inline Ψ(r::Float64, z::Float64, p::GyreParams) = (p.U * p.σr / sqrt(2)) .* (exp.(-(r./p.σr).^2)) .* exp.(0.5 .- (z./p.σz).^2)

#Function to evaluate Ψ on a 2D slice at z = constant_z
@inline Ψ_2D(i, j, g, constant_z::Float64, p::GyreParams) = @inbounds Ψ(r_cc_2D(i, j, g), constant_z, p)

#Functions to compute, directly by centered-differencing, approximations for
# ∂F/∂z at z-Center and z-Face, respectively.
#In future, I will make them more concise using Oceananigans' Δz operators (but
# first, I need to clarify whether, e.g., Δzᵃᵃᶠ(i, j, k, g) returns the spacing
# between k^th and (k+1)^th z-Faces or the spacing between z-Centers adjacent to
# k^th z-Face).
#For a uniformly spaced grid, though, it shouldn't matter.

@inline ΔFΔzᵃᵃᶜ_2D(i, j, g, k::Float64, Fᵃᵃᶠ_minus, Fᵃᵃᶠ_plus) = @inbounds (Fᵃᵃᶠ_plus[i, j, convert(Int64, k + 0.5)] - Fᵃᵃᶠ_minus[i, j, convert(Int64, k - 0.5)]) / (g.z.cᵃᵃᶠ[convert(Int64, k + 0.5)] - g.z.cᵃᵃᶠ[convert(Int64, k - 0.5)])
@inline ΔFΔzᵃᵃᶠ_2D(i, j, g, k::Float64, Fᵃᵃᶜ_minus, Fᵃᵃᶜ_plus) = @inbounds (Fᵃᵃᶜ_plus[i, j, convert(Int64, k + 0.5)] - Fᵃᵃᶜ_minus[i, j, convert(Int64, k - 0.5)]) / (g.z.cᵃᵃᶜ[convert(Int64, k + 0.5)] - g.z.cᵃᵃᶜ[convert(Int64, k - 0.5)])

function GridAndDiscretelyBalancedBackgroundFields(gridParams::GridParams,
                                                   gyreParams::GyreParams,
                                                   perturbParams::PerturbationParams,
                                                   visParams::VisParams;
                                                   x_idx = nothing, 
                                                   y_idx = nothing)
                                                   
   grid = BuildOceananigansRectilinearGrid(gridParams)
   
   #The background streamfunction as a ZFaceField
   @inline Ψᶜᶜᶠ(i, j, k, g) = @inbounds Ψ(r_cca(i, j, k, g), g.z.cᵃᵃᶠ[k], gyreParams)
   @compute Ψ_ccf_Field     = Field(KernelFunctionOperation{Center, Center, Face}(Ψᶜᶜᶠ, grid))
   
   #The same streamfunction field, interpolated to various staggered locations
   @inline Ψᶜᶠᶠ(i, j, k, g) = @inbounds ℑyᵃᶠᵃ(i, j, k, g, Ψ_ccf_Field)
   @compute Ψ_cff_Field     = Field(KernelFunctionOperation{Center, Face, Face}(Ψᶜᶠᶠ, grid))
   @inline Ψᶠᶠᶠ(i, j, k, g) = @inbounds ℑxᶠᵃᵃ(i, j, k, g, Ψ_cff_Field)
   @compute Ψ_fff_Field     = Field(KernelFunctionOperation{Face, Face, Face}(Ψᶠᶠᶠ, grid))
   @inline Ψᶠᶠᶜ(i, j, k, g) = @inbounds ℑzᵃᵃᶜ(i, j, k, g, Ψ_fff_Field)
   @compute Ψ_ffc_Field     = Field(KernelFunctionOperation{Face, Face, Center}(Ψᶠᶠᶜ, grid))

   #First z-derivative of Ψ at necessary location to diagnose buoyancy
   @compute ∂Ψ∂z_ccc_Field = Field(∂z(Ψ_ccf_Field))
   
   #The halo values of ∂Ψ∂z_ccc_Field are 0 by default.
   #However, the value of ∂Ψ/∂z at the innermost halo cell *Center* at each vertical boundary needs to be computed before ∂2Ψ∂z2_ccf_Field can be computed correctly.
   
   ##########################################################
   # Explicitly compute the innermost halo-face values of Ψ #
   ##########################################################
   
   #The z-coordinates
   bot_inner_halo_zF = no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ))[gridParams.Hz]
   top_inner_halo_zF = no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ)
                                     )[gridParams.Nz + gridParams.Hz + 2]

   #These slices are 2D grids
   bot_inner_halo_zF_slice = slice(grid, :, :, 0, z = bot_inner_halo_zF)
   top_inner_halo_zF_slice = slice(grid, :, :, gridParams.Nz + 2, 
                                   z = top_inner_halo_zF)
   
   #Compute these as Fields on the 2D grids
   @compute bot_inner_halo_Ψ_ccf = Field(KernelFunctionOperation{Center, Center, Nothing}(Ψ_2D, bot_inner_halo_zF_slice, bot_inner_halo_zF, gyreParams))
   @compute top_inner_halo_Ψ_ccf = Field(KernelFunctionOperation{Center, Center, Nothing}(Ψ_2D, top_inner_halo_zF_slice, top_inner_halo_zF, gyreParams))

   ##############################################################
   # Explicitly compute the innermost halo-cell values of ∂Ψ/∂z #
   ##############################################################

   #Explicitly compute ∂Ψ/∂z at the innermost halo cell *Center* point at each
   # vertical boundary, by finite differencing.
   #Note that we do not impose these values as BoundaryConditions on
   # ∂Ψ∂z_ccc_Field because:
   # (a) we can't: Oceananigans can only apply a ValueBoundaryCondition at a
   #     cell Face; and 
   # (b) we don't need to dynamically evolve ∂Ψ∂z_ccc_Field anyway.
   
   #The z-coordinates
   bot_inner_halo_zC = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ))[gridParams.Hz]
   top_inner_halo_zC = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ)
                                     )[gridParams.Nz + gridParams.Hz + 1]
   
   #More 2D grids
   bot_inner_halo_zC_slice = slice(grid, :, :, 0, z = bot_inner_halo_zC)
   top_inner_halo_zC_slice = slice(grid, :, :, gridParams.Nz + 1, 
                                   z = top_inner_halo_zC)

   #Compute these as Fields on the 2D grids
   @compute bot_inner_halo_∂Ψ∂z_ccc = Field(KernelFunctionOperation{Center, Center, Nothing}(ΔFΔzᵃᵃᶜ_2D, grid, 0.5, bot_inner_halo_Ψ_ccf, Ψ_ccf_Field))
   @compute top_inner_halo_∂Ψ∂z_ccc = Field(KernelFunctionOperation{Center, Center, Nothing}(ΔFΔzᵃᵃᶜ_2D, grid, (gridParams.Nz + 1.5), Ψ_ccf_Field, top_inner_halo_Ψ_ccf))

   ################################################################
   # Explicitly compute the innermost halo-face values of ∂2Ψ/∂z2 #
   ################################################################
   
   @compute bot_∂2Ψ∂z2_ccf = Field(KernelFunctionOperation{Center, Center, Nothing}(ΔFΔzᵃᵃᶠ_2D, grid, 0.5, bot_inner_halo_∂Ψ∂z_ccc, ∂Ψ∂z_ccc_Field))
   @compute top_∂2Ψ∂z2_ccf = Field(KernelFunctionOperation{Center, Center, Nothing}(ΔFΔzᵃᵃᶠ_2D, grid, (gridParams.Nz + 1.5), ∂Ψ∂z_ccc_Field, top_inner_halo_∂Ψ∂z_ccc))

   #############################################
   # Diagnose thermal-wind velocity components #
   #############################################
   
   @compute u_TWB = Field(-∂y(Ψ_ffc_Field))
   @compute v_TWB = Field(∂x(Ψ_ffc_Field))
   @compute w_TWB = Field((-∂x(∂y(Ψ_fff_Field)) + ∂y(∂x(Ψ_fff_Field))))
   
   ##############################################################
   # Diagnose thermal-wind buoyancy and its boundary conditions #
   ##############################################################

   @inline LinearBᵃᵃᶜ(i, j, k, g)                 = @inbounds gyreParams.N²_far * g.z.cᵃᵃᶜ[k]
   @inline ThermalWindBᶜᶜᶜ(i, j, k, g)            = @inbounds gyreParams.f * ∂Ψ∂z_ccc_Field[i, j, k]
   @inline TotalConstAmbientStratBᶜᶜᶜ(i, j, k, g) = @inbounds LinearBᵃᵃᶜ(i, j, k, g) + ThermalWindBᶜᶜᶜ(i, j, k, g)
      
   @compute b_TWB = Field(KernelFunctionOperation{Center, Center, Center}(TotalConstAmbientStratBᶜᶜᶜ, grid))
      
   @inline ∂b∂z(i, j, g, ∂2Ψ∂z2_ccf_Field) = @inbounds gyreParams.N²_far + gyreParams.f * ∂2Ψ∂z2_ccf_Field[i, j]
   
   @compute ∂b∂z_bot = Field(KernelFunctionOperation{Center, Center, Nothing}(∂b∂z, grid, bot_∂2Ψ∂z2_ccf))
   @compute ∂b∂z_top = Field(KernelFunctionOperation{Center, Center, Nothing}(∂b∂z, grid, top_∂2Ψ∂z2_ccf))
   
   ##########################
   # Optional visualization #
   ##########################
   
   if visParams.vis_B_and_N²
      visualize_B_and_N²_vs_z(b_TWB, 
                              Field(KernelFunctionOperation{Center, Center, Center}(ThermalWindBᶜᶜᶜ, grid)), 
                              Field(KernelFunctionOperation{Center, Center, Center}(LinearBᵃᵃᶜ, grid)), 
                              grid, x_idx, y_idx, gyreParams, nothing, 
                              gridParams.Hz)
   end

   #######################################
   # Compute the perturbed initial state #
   #######################################
   
   Random.seed!(perturbParams.seed) #For reproducibility
   
   randArray = adapt(CuArray, 2 * (rand(Float64, (gridParams.Nx, 
                                                  gridParams.Ny, 
                                                  gridParams.Nz + 1)) 
                                   .- 0.5)
                    )

   @inline Ψperturbedᶜᶜᶠ(i, j, k, g) = @inbounds (perturbParams.max_ψ′ * randArray[i, j, k] + 1) * Ψᶜᶜᶠ(i, j, k, g)
   @compute Ψperturbed_ccf_Field = Field(KernelFunctionOperation{Center, Center, Face}(Ψperturbedᶜᶜᶠ, grid))
   
   @inline Ψperturbedᶜᶠᶠ(i, j, k, g) = @inbounds ℑyᵃᶠᵃ(i, j, k, g, Ψperturbed_ccf_Field)
   @compute Ψperturbed_cff_Field = Field(KernelFunctionOperation{Center, Face, Face}(Ψperturbedᶜᶠᶠ, grid))
   @inline Ψperturbedᶠᶠᶠ(i, j, k, g) = @inbounds ℑxᶠᵃᵃ(i, j, k, g, Ψperturbed_cff_Field)
   @compute Ψperturbed_fff_Field = Field(KernelFunctionOperation{Face, Face, Face}(Ψperturbedᶠᶠᶠ, grid))
   @inline Ψperturbedᶠᶠᶜ(i, j, k, g) = @inbounds ℑzᵃᵃᶜ(i, j, k, g, Ψperturbed_fff_Field)
   @compute Ψperturbed_ffc_Field = Field(KernelFunctionOperation{Face, Face, Center}(Ψperturbedᶠᶠᶜ, grid))
   
   @compute ∂Ψ∂z_perturbed_ccc_Field = Field(∂z(Ψperturbed_ccf_Field))
   
   @inline ThermalWind_b_perturbedᶜᶜᶜ(i, j, k, g) = @inbounds gyreParams.f * ∂Ψ∂z_perturbed_ccc_Field[i, j, k]
   @inline perturbed_bᶜᶜᶜ(i, j, k, g) = @inbounds LinearBᵃᵃᶜ(i, j, k, g) + ThermalWind_b_perturbedᶜᶜᶜ(i, j, k, g)
   
   @compute b_perturbed = Field(KernelFunctionOperation{Center, Center, Center}(perturbed_bᶜᶜᶜ, grid))
   @compute u_perturbed = Field(-∂y(Ψperturbed_ffc_Field))
   @compute v_perturbed = Field(∂x(Ψperturbed_ffc_Field))
   @compute w_perturbed = Field((-∂x(∂y(Ψperturbed_fff_Field)) + ∂y(∂x(Ψperturbed_fff_Field))))
   
   return (grid, 
           (b_TWB, u_TWB, v_TWB, w_TWB), 
           FieldBoundaryConditions(top = GradientBoundaryCondition(∂b∂z_top),
                                   bottom = GradientBoundaryCondition(∂b∂z_bot)),
           (b_perturbed, u_perturbed, v_perturbed, w_perturbed)
          )
end

@inline b_double_tanh_ccc(i, j, k, g, p::DoubleTanhStratParams) = @inbounds buoyancyDoubleTanh(g.z.cᵃᵃᶜ[k], p)

#Double-tanh buoyancy function
@inline buoyancyDoubleTanh(z, p::DoubleTanhStratParams) = @. -(p.gravity / p.ρ₀) * (p.A_s * tanh((z - p.z_s) / p.C_s) + p.A_d * tanh((z - p.z_d) / p.C_d))

#N² corresponding to double-tanh buoyancy function
@inline N²DoubleTanh(z, p::DoubleTanhStratParams) = @. -(p.gravity / p.ρ₀) * ((p.A_s / p.C_s) * (sech((z - p.z_s) / p.C_s))^2 + (p.A_d / p.C_d) * (sech((z - p.z_d) / p.C_d))^2)