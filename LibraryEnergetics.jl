using Adapt
using LinearAlgebra: norm
using Oceananigans
using Oceananigans.Architectures
#using Oceananigans.Fields
using Oceananigans.Grids

#Functions to compute cylindrical components of perturbation velocity
@inline ur′(i, j, k, grid, ur, Ur) = @inbounds ur[i, j, k] - Ur[i, j, k]
@inline uφ′(i, j, k, grid, uφ, Uφ) = @inbounds uφ[i, j, k] - Uφ[i, j, k]
@inline uz′(i, j, k, grid, uz, Uz) = @inbounds uz[i, j, k] - Uz[i, j, k]

#Function to compute 2-norm of a perturbation field
@inline perturbation_norm(field, bkgd_field) = norm(field - bkgd_field)

#Function to compute square of a perturbation field
@inline ψ′²(i, j, k, grid, ψ, ψ̄) = @inbounds (ψ[i, j, k] - ψ̄[i, j, k])^2

#Function to compute PKE in single control volume
@inline PKE_ccc(i, j, k, grid, u, v, w, Ux, Uy, Uz) = @inbounds (
     		      		ℑxᶜᵃᵃ(i, j, k, grid, ψ′², u, Ux) + 
     		      		ℑyᵃᶜᵃ(i, j, k, grid, ψ′², v, Uy) +
                      		ℑzᵃᵃᶜ(i, j, k, grid, ψ′², w, Uz)) / 2

#Function to compute product of b′ and uz′ in single control volume
@inline b′uz′_ccc(i, j, k, grid, b, uz, B, Uz) = @inbounds (
		  (b[i, j, k] - B[i, j, k]) * ℑzᵃᵃᶜ(i, j, k, grid, uz′, uz, Uz))

function ∂rUφ_ur′_uφ′_ccc(i, j, k, grid, ux, uy, Ur, Uφ, ∂rUφ)
   #=
   Function to compute (ur′ times uφ′ times r-derivative of Uφ) in single 
    control volume.
   =#

   φ = @inbounds atan(ℑyᵃᶜᵃ(i, j, k, grid, ynodes(grid, Center())),
                      ℑxᶜᵃᵃ(i, j, k, grid, xnodes(grid, Center()))
                     )

   ux_ccc = @inbounds ℑxᶜᵃᵃ(i, j, k, grid, ux)
   uy_ccc = @inbounds ℑyᵃᶜᵃ(i, j, k, grid, uy)
   ur_ccc = @inbounds (ux_ccc * cos(φ)) + (uy_ccc * sin(φ))
   uφ_ccc = @inbounds (uy_ccc * cos(φ)) - (ux_ccc * sin(φ))

   ur′_ccc  = @inbounds ur′(i, j, k, grid, ur_ccc, Ur)
   uφ′_ccc  = @inbounds uφ′(i, j, k, grid, uφ_ccc, Uφ)
   ∂rUφ_ccc = @inbounds ∂rUφ[i, j, k]

   return @inbounds -(∂rUφ_ccc * ur′_ccc * uφ′_ccc)
end

function ∂zUφ_uφ′_uz′_ccc(i, j, k, grid, ux, uy, uz, Uφ, Uz, ∂zUφ)
   #=
   Function to compute (negative uφ′ times uz′ times z-derivative of Uφ) in 
    single control volume.
   =#

   φ = @inbounds atan(ℑyᵃᶜᵃ(i, j, k, grid, ynodes(grid, Center())),
                      ℑxᶜᵃᵃ(i, j, k, grid, xnodes(grid, Center()))
                     )

   ux_ccc = @inbounds ℑxᶜᵃᵃ(i, j, k, grid, ux)
   uy_ccc = @inbounds ℑyᵃᶜᵃ(i, j, k, grid, uy)
   uφ_ccc = @inbounds (uy_ccc * cos(φ)) - (ux_ccc * sin(φ)) 

   uφ′_ccc  = @inbounds uφ′(i, j, k, grid, uφ_ccc, Uφ)
   uz′_ccc  = @inbounds ℑzᵃᵃᶜ(i, j, k, grid, uz′, uz, Uz) 
   ∂zUφ_ccc = @inbounds ∂zUφ[i, j, k]

   return @inbounds -(∂zUφ_ccc * uφ′_ccc * uz′_ccc)
end
                         
function IsWithinGyreRegion(i, j, k, grid, integrand; parameters) 
   #=
   Return true if gridpoint coords are r ≤ (2 * σr) and (only if baroclinic 
    background state) -(2 * σz) ≤ z; return false if gridpoint is outside this
    region.
   =#
   
   σr, σz       = parameters.σr, parameters.σz
   architecture = parameters.architecture

   if architecture == GPU()
      xnodes = adapt(CuArray, xnodes(grid, Center()))
      ynodes = adapt(CuArray, ynodes(grid, Center()))
      znodes = adapt(CuArray, znodes(grid, Center()))
   elseif architecture == CPU()
      xnodes = grid.xᶜᵃᵃ
      ynodes = grid.yᵃᶜᵃ
      znodes = grid.z.cᵃᵃᶜ
   end

   isWithinGyreRegion = @views (-(2 * σr) ≤ ynodes[j] ≤ (2 * σr)
                                && -sqrt((2 * σr)^2 - ynodes[j]^2) ≤ 
                                    xnodes[i] ≤ sqrt((2 * σr)^2 - ynodes[j]^2)
                                && (σz == "infinity"
                                    || -(2 * σz) ≤ znodes[k] ≤ 0
                                   )
                               )
end
   
function total_PKE(simulation; Ux, Uy, Uz)
   #=
   Return computed integral of PKE over entire domain.
   =#
   
   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   pointwise_PKE_op = KernelFunctionOperation{Center, Center, Center}(PKE_ccc, 
                                       grid, u, v, w, Ux, Uy, Uz)
   
   compute!(Integral(Field(pointwise_PKE_op)))
end

function gyre_PKE(simulation; gyreParameters)
   #=
   Return computed integral of PKE, taken over gyre region only.
   =#

   @inline IsWithinModelGyreRegion(i, j, k, grid, integrand) = @inbounds IsWithinGyreRegion(i, j, k, grid, integrand; parameters = gyreParameters)
   
   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   pointwise_PKE_op = KernelFunctionOperation{Center, Center, Center}(PKE_ccc, 
                                       grid, u, v, w, Ux, Uy, Uz)

   compute!(Integral(Field(pointwise_PKE_op), condition = IsWithinModelGyreRegion))
end

function total_PAPE_to_PKE(simulation; B, Uz)
   #=
   Return computed integral, over entire domain, of conversion from perturbation 
    APE to PKE.
   =#
   
   grid = simulation.model.grid
   b, w = simulation.model.tracers.b, simulation.model.velocities.w

   PAPE_to_PKE_op = KernelFunctionOperation{Center, Center, Center}(b′uz′_ccc,
                                           grid, b, w, B, Uz)

   compute!(Integral(Field(PAPE_to_PKE_op)))
end

function gyre_PAPE_to_PKE(simulation; gyreParameters)
   #=
   Return computed integral, taken over gyre region only, of conversion from 
    perturbation APE to PKE.
   =#
   
   @inline IsWithinModelGyreRegion(i, j, k, grid, integrand) = @inbounds IsWithinGyreRegion(i, j, k, grid, integrand; parameters = gyreParameters)
   
   grid = simulation.model.grid
   b, w = simulation.model.tracers.b, simulation.model.velocities.w

   PAPE_to_PKE_op = KernelFunctionOperation{Center, Center, Center}(b′uz′_ccc,
                                           grid, b, w, B, Uz)

   compute!(Integral(Field(PAPE_to_PKE_op), condition = IsWithinModelGyreRegion))
end

function BTI_transfer(simulation; bkgdParameters)
   #=
   Return computed integral, over entire domain, of BTI-transfer term in PKE
    budget.
   =#
   
   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   Ur   = bkgdParameters.Ur
   Uφ   = bkgdParameters.Uφ
   ∂rUφ = bkgdParameters.∂rUφ
   
   BTI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂rUφ_ur′_uφ′_ccc, grid, u, v, Ur, Uφ, ∂rUφ)
   
   compute!(Integral(Field(BTI_transfer_op)))
end

function gyre_BTI_transfer(simulation; bkgdParameters, gyreParameters)
   #=
   Return computed integral, taken over gyre region only, of BTI-transfer term
    in PKE budget.
   =#

   @inline IsWithinModelGyreRegion(i, j, k, grid, integrand) = @inbounds IsWithinGyreRegion(i, j, k, grid, integrand; parameters = gyreParameters)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   Ur   = bkgdParameters.Ur
   Uφ   = bkgdParameters.Uφ
   ∂rUφ = bkgdParameters.∂rUφ

   BTI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂rUφ_ur′_uφ′_ccc, grid, u, v, Ur, Uφ, ∂rUφ)

   compute!(Integral(Field(BTI_transfer_op), condition = IsWithinModelGyreRegion))
end

function BCI_transfer(simulation; bkgdParameters)
   #=
   Return computed integral, over entire domain, of BCI-transfer term in PKE
    budget.
   =#
   
   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   Uφ   = bkgdParameters.Uφ
   Uz   = bkgdParameters.Uz
   ∂zUφ = bkgdParameters.∂zUφ
   
   BCI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂zUφ_uφ′_uz′_ccc, grid, u, v, w, Uφ, Uz, ∂zUφ)

   compute!(Integral(Field(BCI_transfer_op)))
end

function gyre_BCI_transfer(simulation; bkgdParameters, gyreParameters)
   #=
   Return computed integral, taken over gyre region only, of BCI-transfer term
    in PKE budget.
   =#
   
   @inline IsWithinModelGyreRegion(i, j, k, grid, integrand) = @inbounds IsWithinGyreRegion(i, j, k, grid, integrand; parameters = gyreParameters)

   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   Uφ   = bkgdParameters.Uφ
   Uz   = bkgdParameters.Uz
   ∂zUφ = bkgdParameters.∂zUφ
   
   BCI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂zUφ_uφ′_uz′_ccc, grid, u, v, w, Uφ, Uz, ∂zUφ)

   compute!(Integral(Field(BCI_transfer_op), condition = IsWithinModelGyreRegion))
end