using Adapt, CUDA
using Dates: canonicalize, format, now
using KernelAbstractions: @index, @kernel
using NCDatasets
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.Units
using Oceananigans.Utils
using Oceanostics
using OffsetArrays: no_offset_view
using Printf, Random

#Numbers of gridpoints
const Nx = 12
const Ny = 12
const Nz = 12

#Lengths of axes
const Lx = 2.5e3 * kilometer
const Ly = 2.5e3 * kilometer
const Lz = 1 * kilometer

#Latitude (deg. N)
const lat = 74.0

#f-plane and Coriolis frequency
fPlane  = FPlane(latitude = lat)
const f = fPlane.f

#Characteristic scales for basic state
const σr  = 250 * kilometer
const U   = 1.5e-1 * (meter/second)
const N²₀ = 3e-3 * (second^(-2))

#Time-stepping parameters
const Δt      = 10 * second
const tf      = 3 * hour
const Δt_save = 1 * hour

#Max. relative magnitude of initial u-perturbations
const max_u′ = 1e-8


const use_GPU = false

use_GPU ? architecture = GPU() : architecture = CPU()

grid = RectilinearGrid(architecture,
		       topology = (Bounded, Bounded, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lx/2, Lx/2), 
                       y = (-Ly/2, Ly/2), 
                       z = (-Lz, 0.0),
		       halo = (3, 3, 3))

b̄ = (x, y, z) -> N²₀ * z
ū = (x, y, z) -> (sqrt(2)*U*y/σr) * exp((1/2) - (x^2 + y^2)/(σr^2))
v̄ = (x, y, z) -> -(sqrt(2)*U*x/σr) * exp((1/2) - (x^2 + y^2)/(σr^2))

b̄z_top = (x, y, t) -> N²₀
b̄z_bot = (x, y, t) -> N²₀
b̄_BCs  = FieldBoundaryConditions(top    = GradientBoundaryCondition(b̄z_top),
         			 bottom = GradientBoundaryCondition(b̄z_bot))

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :RungeKutta3,
                            advection = UpwindBiasedFifthOrder(), 
                            coriolis = fPlane,
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
			    boundary_conditions = (; b = b̄_BCs,))

set!(model, u = ū, v = v̄)
fill_halo_regions!(model.velocities)
set!(model, b = b̄)

datetimenow = format(now(), "yymmdd-HHMMSS")

Ux = XFaceField(model.grid)
Uy = YFaceField(model.grid)
Uz = ZFaceField(model.grid)

set!(Ux, model.velocities.u)
set!(Uy, model.velocities.v)
set!(Uz, model.velocities.w)
fill_halo_regions!(Ux)
fill_halo_regions!(Uy)
fill_halo_regions!(Uz)

@kernel function pKE_ccc!(pKE, grid, u, v, w, Ux, Uy, Uz)
   i, j, k = @index(Global, NTuple)
   @inbounds pKE[i, j, k] = (
                              ℑxᶜᵃᵃ(i, j, k, grid, ψ′², u, Ux) +
                              ℑyᵃᶜᵃ(i, j, k, grid, ψ′², v, Uy) +
                              ℑzᵃᵃᶜ(i, j, k, grid, ψ′², w, Uz)
                             ) / 2
end

pKE = CenterField(model.grid)

#pKE_ccc!(pKE, model.grid, model.velocities.u, model.velocities.v, model.velocities.w, Ux, Uy, Uz)

pKE_op = KernelFunctionOperation{Center, Center, Center}(pKE_ccc!, 
							 grid, 
							 model.velocities.u, 
							 model.velocities.v, 
							 model.velocities.w,
							 Ux, Uy, Uz)

#pKE_field = Field(pKE_op)
compute!(pKE_op)#field)

#Perturb velocity components to trigger instability
@inline u_perturbed(x, y, z) = (ū(x, y, z) + (2 * (rand() - 0.5)) * (max_u′ / sqrt(2)))
@inline v_perturbed(x, y, z) = (v̄(x, y, z) + (2 * (rand() - 0.5)) * (max_u′ / sqrt(2)))
set!(model, u = u_perturbed, v = v_perturbed)
fill_halo_regions!(model.velocities, model.tracers.b)

simulation = Simulation(model, Δt = Δt, stop_time = tf)

energy_diagnostics = (; pKE = pKE_op)
energyfilepath = joinpath("./Output", "energetics_$(datetimenow).nc")
mkpath(dirname(energyfilepath)) #Make path if nonexistent

energy_writer = NetCDFOutputWriter(model, 
				   energy_diagnostics,
                                   with_halos = true,
				   filename = energyfilepath,
				   schedule = TimeInterval(Δt_save),
				   file_splitting = FileSizeLimit(30GiB))

simulation.output_writers[:energy_writer] = energy_writer

run!(simulation)
