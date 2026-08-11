using Adapt
using LinearAlgebra: norm
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Fields
using Oceananigans.Grids
using Oceanostics.KineticEnergyEquation

#Function to compute total potential energy in single control volume
@inline totalPE_ccc(i, j, k, grid, b, g) = @inbounds ((g - b[i, j, k]) * grid.z.cᵃᵃᶜ[k])

function totalKE(simulation)
   #=
   Return computed integral, over entire domain, of total kinetic energy.
   =#

   u, v, w = simulation.model.velocities

   totalKE_op = KineticEnergy(simulation.model, u, v, w)
   
   compute!(Integral(Field(totalKE_op)))
end

function totalKEadvFlux(simulation; useOceanostics = false)
   #=
   Return computed integral, over entire domain, of total advective KE-flux.
   =#
   
   if useOceanostics #Note this requires a NonhydrostaticModel
      
      totalKEadvFlux_op = KineticEnergyAdvection(simulation.model)
      
   elseif !useOceanostics #Can run with Nonhydrostatic or HydrostaticFreeSurface
   
      u, v, w = simulation.model.velocities
   
      totalKE_op = KineticEnergy(simulation.model, u, v, w)
      
      @compute KE = Field(totalKE_op)
      
      @compute ∂xKE = ∂x((Center, Center, Center), KE)
      @compute ∂yKE = ∂y((Center, Center, Center), KE)
      @compute ∂zKE = ∂z((Center, Center, Center), KE)
      
      @inline KEadvFlux_ccc(i, j, k, grid) = @inbounds (u[i, j, k] * ∂xKE[i, j, k] + v[i, j, k] * ∂yKE[i, j, k] + w[i, j, k] * ∂zKE[i, j, k]) / 2
      #Note the factor of 1/2 is required for consistency with Oceanostics
      # functions.
   
      totalKEadvFlux_op = KernelFunctionOperation{Center, Center, Center}(KEadvFlux_ccc, simulation.model.grid)
   end
   
   compute!(Integral(Field(-totalKEadvFlux_op)))
end

function totalPressureWork(simulation; useNHS = nothing, useOceanostics = false)
   #=
   Return computed integral, over entire domain, of total pressure work.
   =#
   
   if useOceanostics #Note this requires a NonhydrostaticModel
      
      totalPressureWork_op = KineticEnergyPressureRedistribution(simulation.model)
   
   elseif !useOceanostics #Can run with Nonhydrostatic or HydrostaticFreeSurface
   
      if useNHS
         p = simulation.model.pressures.pNHS #Note this is kinematic pressure
      elseif !useNHS
         p = simulation.model.pressure.pHY′ #Note this is kinematic pressure
      end

      u, v, w = simulation.model.velocities
      
      @compute ∂xp = ∂x(p)
      @compute ∂yp = ∂y(p)
      @compute ∂zp = ∂z(p)
      
      @inline pressureWork_ccc(i, j, k, grid) = @inbounds (-(u[i, j, k] * ∂xp[i, j, k] + v[i, j, k] * ∂yp[i, j, k] + w[i, j, k] * ∂zp[i, j, k]))
      
      totalPressureWork_op = KernelFunctionOperation{Center, Center, Center}(pressureWork_ccc, simulation.model.grid)
   end
   
   compute!(Integral(Field(totalPressureWork_op)))
end

function totalPE(simulation, g)
   #=
   Return computed integral, over entire domain, of total kinetic energy.
   =#
   
   grid = simulation.model.grid
   b    = simulation.model.tracers.b

   totalPE_op = KernelFunctionOperation{Center, Center, Center}(totalPE_ccc, grid, b, g)

   compute!(Integral(Field(totalPE_op)))
end

function totalGravityWork(simulation, g)
   #=
   Return computed integral, over entire domain, of total gravity work.
   =#

   uz = simulation.model.velocities.w
   
   @inline gravityWork_ccc(i, j, k, grid) = @inbounds uz[i, j, k]
      
   totalGravityWork_op = KernelFunctionOperation{Center, Center, Center}(gravityWork_ccc, simulation.model.grid)

   compute!(Integral(Field(g * totalGravityWork_op)))
end

function totalBuoyancyAdvFlux(simulation)
   #=
   Return computed integral, over entire domain, of total advective KE-flux.
   =#
   
   b       = simulation.model.tracers.b
   u, v, w = simulation.model.velocities
      
   @compute ∂xb = ∂x((Center, Center, Center), b)
   @compute ∂yb = ∂y((Center, Center, Center), b)
   @compute ∂zb = ∂z((Center, Center, Center), b)
      
   @inline buoyancyAdvFlux_ccc(i, j, k, grid) = @inbounds ((u[i, j, k] * ∂xb[i, j, k] + v[i, j, k] * ∂yb[i, j, k] + w[i, j, k] * ∂zb[i, j, k]) * grid.z.cᵃᵃᶜ[k])
   
   totalBuoyancyAdvFlux_op = KernelFunctionOperation{Center, Center, Center}(buoyancyAdvFlux_ccc, simulation.model.grid)

   compute!(Integral(Field(totalBuoyancyAdvFlux_op)))
end

#Functions to compute cylindrical components of perturbation velocity
@inline ur′(i, j, k, grid, ur, Ur) = @inbounds ur[i, j, k] - Ur[i, j, k]
@inline uφ′(i, j, k, grid, uφ, Uφ) = @inbounds uφ[i, j, k] - Uφ[i, j, k]
@inline uz′(i, j, k, grid, uz, Uz) = @inbounds uz[i, j, k] - Uz[i, j, k]

#Function to compute 2-norm of a perturbation field
@inline perturbation_norm(field, bkgd_field) = norm(field - bkgd_field)

#Function to compute square of a perturbation field
@inline ψ′²(i, j, k, grid, ψ, ψ̄) = @inbounds (ψ[i, j, k] - ψ̄[i, j, k])^2

#Function to compute PKE in single control volume
@inline PKE_ccc(i, j, k, grid, u, v, w, Ux, Uy, Uz) = @inbounds (ℑxᶜᵃᵃ(i, j, k, grid, ψ′², u, Ux) +	ℑyᵃᶜᵃ(i, j, k, grid, ψ′², v, Uy) + ℑzᵃᵃᶜ(i, j, k, grid, ψ′², w, Uz)) / 2

#Function to compute product of b′ and uz′ in single control volume
@inline b′uz′_ccc(i, j, k, grid, b, uz, B, Uz) = @inbounds ((b[i, j, k] - B[i, j, k]) * ℑzᵃᵃᶜ(i, j, k, grid, uz′, uz, Uz))

function ur′_uφ′_Uφ_over_r_ccc(i, j, k, grid, ux, uy, Uφ) ##Ur, Uφ, ∂rUφ)
   #=
   Function to compute (1/r times ur′ times uφ′ times Uφ) in single control
    volume.
   =#

   r = @inbounds sqrt(xnodes(grid, Center())[i, j, k]^2 + ynodes(grid, Center())[i, j, k]^2)
   φ = @inbounds atan(ℑyᵃᶜᵃ(i, j, k, grid, ynodes(grid, Center())),
                      ℑxᶜᵃᵃ(i, j, k, grid, xnodes(grid, Center()))
                     )

   ux_ccc = @inbounds ℑxᶜᵃᵃ(i, j, k, grid, ux)
   uy_ccc = @inbounds ℑyᵃᶜᵃ(i, j, k, grid, uy)
   ur_ccc = @inbounds (ux_ccc * cos(φ)) + (uy_ccc * sin(φ))
   uφ_ccc = @inbounds (uy_ccc * cos(φ)) - (ux_ccc * sin(φ))

   ur′_ccc       = @inbounds ur′(i, j, k, grid, ur_ccc, 0) #Ur)
   uφ′_ccc       = @inbounds uφ′(i, j, k, grid, uφ_ccc, Uφ)
   Uφ_over_r_ccc = @inbounds Uφ[i, j, k] / r[i, j, k]

   return @inbounds -(Uφ_over_r_ccc * ur′_ccc * uφ′_ccc)
end

function neg∂rUφ_ur′_uφ′_ccc(i, j, k, grid, ux, uy, Uφ, ∂rUφ) ##Ur, Uφ, ∂rUφ)
   #=
   Function to compute (negative ur′ times uφ′ times r-derivative of Uφ) in
    single control volume.
   =#

   φ = @inbounds atan(ℑyᵃᶜᵃ(i, j, k, grid, ynodes(grid, Center())),
                      ℑxᶜᵃᵃ(i, j, k, grid, xnodes(grid, Center()))
                     )

   ux_ccc = @inbounds ℑxᶜᵃᵃ(i, j, k, grid, ux)
   uy_ccc = @inbounds ℑyᵃᶜᵃ(i, j, k, grid, uy)
   ur_ccc = @inbounds (ux_ccc * cos(φ)) + (uy_ccc * sin(φ))
   uφ_ccc = @inbounds (uy_ccc * cos(φ)) - (ux_ccc * sin(φ))

   ur′_ccc  = @inbounds ur′(i, j, k, grid, ur_ccc, 0) #Ur)
   uφ′_ccc  = @inbounds uφ′(i, j, k, grid, uφ_ccc, Uφ)
   ∂rUφ_ccc = @inbounds ∂rUφ[i, j, k]

   return @inbounds -(∂rUφ_ccc * ur′_ccc * uφ′_ccc)
end

function neg∂zUφ_uφ′_uz′_ccc(i, j, k, grid, ux, uy, uz, Uφ, Uz, ∂zUφ)
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
   
   σr, σz = parameters.σr, parameters.σz
   
   isWithinGyreRegion = @inbounds (-(2 * σr) ≤ grid.yᵃᶜᵃ[j] ≤ (2 * σr)
                                   && -sqrt((2 * σr)^2 - grid.yᵃᶜᵃ[j]^2) ≤
                                      grid.xᶜᵃᵃ[i] ≤ 
                                      sqrt((2 * σr)^2 - grid.yᵃᶜᵃ[j]^2)
                                   && (σz == "infinity"
                                       || -(2 * σz) ≤ grid.z.cᵃᵃᶜ[k] ≤ 0
                                      )
                                  )
   return isWithinGyreRegion
end
   
function PKE(simulation, Ux, Uy, Uz)
   #=
   Return computed integral of PKE over entire domain.
   =#
   
   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   PKE_op = KernelFunctionOperation{Center, Center, Center}(PKE_ccc, 
                                       grid, u, v, w, Ux, Uy, Uz)
   
   compute!(Integral(Field(PKE_op)))
end

function gyre_PKE(simulation; gyreParameters)
   #=
   Return computed integral of PKE, taken over gyre region only.
   =#

   @inline IsWithinModelGyreRegion(i, j, k, grid, integrand) = @inbounds IsWithinGyreRegion(i, j, k, grid, integrand; parameters = gyreParameters)
   
   grid    = simulation.model.grid
   u, v, w = simulation.model.velocities
   
   PKE_op = KernelFunctionOperation{Center, Center, Center}(PKE_ccc, grid, u, v, w, Ux, Uy, Uz)
   PKE    = Field(PKE_op)

   #Mask areas far from gyre
   @compute mask = Field(KernelFunctionOperation{Center, Center, Center}(
                         IsWithinModelGyreRegion, grid, PKE))

   compute!(Integral(PKE, mask = mask))
end

function PAPE_to_PKE(simulation, B, Uz)
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
   PAPE_to_PKE    = Field(PAPE_to_PKE_op)
   
   #Mask areas far from gyre
   @compute mask = Field(KernelFunctionOperation{Center, Center, Center}(
                         IsWithinModelGyreRegion, grid, PAPE_to_PKE))

   compute!(Integral(PAPE_to_PKE, mask = mask))
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
   
   BTI_transfer_ccc(i, j, k, grid) = ur′_uφ′_Uφ_over_r_ccc(i, j, k, grid, u, v, Uφ) + neg∂rUφ_ur′_uφ′_ccc(i, j, k, grid, u, v, Uφ, ∂rUφ)
   
   BTI_transfer_op = KernelFunctionOperation{Center, Center, Center}(BTI_transfer_ccc, grid)#, u, v, Ur, Uφ, ∂rUφ)
   
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

   #BTI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂rUφ_ur′_uφ′_ccc, 
   #                                                      grid, u, v, Uφ, ∂rUφ) ## Ur, Uφ, ∂rUφ)

   BTI_transfer_ccc(i, j, k, grid) = ur′_uφ′_Uφ_over_r_ccc(i, j, k, grid, u, v, Uφ) + neg∂rUφ_ur′_uφ′_ccc(i, j, k, grid, u, v, Uφ, ∂rUφ)
   
   BTI_transfer_op = KernelFunctionOperation{Center, Center, Center}(BTI_transfer_ccc, grid)#, u, v, Ur, Uφ, ∂rUφ)
   BTI_transfer    = Field(BTI_transfer_op)

   #Mask areas far from gyre
   @compute mask = Field(KernelFunctionOperation{Center, Center, Center}(
                         IsWithinModelGyreRegion, grid, BTI_transfer))

   compute!(Integral(BTI_transfer_op, mask = mask))
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
   
   BCI_transfer_op = KernelFunctionOperation{Center, Center, Center}(neg∂zUφ_uφ′_uz′_ccc, grid, u, v, w, Uφ, Uz, ∂zUφ)

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
   
   BCI_transfer_op = KernelFunctionOperation{Center, Center, Center}(neg∂zUφ_uφ′_uz′_ccc, 
                                                        grid, u, v, w, Uφ, Uz, ∂zUφ)
   BCI_transfer    = Field(BCI_transfer_op)
   
   #Mask areas far from gyre
   @compute mask = Field(KernelFunctionOperation{Center, Center, Center}(
                         IsWithinModelGyreRegion, grid, BCI_transfer))

   compute!(Integral(BCI_transfer, mask = mask))
end