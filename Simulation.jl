include("LibraryDynamics.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Visualization.jl")

using Dates: canonicalize, format, now
using LinearAlgebra: norm
using NCDatasets
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using Oceananigans.Utils
using Oceanostics
using Printf, Random

######################
# SPECIFY PARAMETERS #
######################

#Numbers of gridpoints
const Nx = 1200
const Ny = 1200
const Nz = 12

#Lengths of axes
const Lx = 2.5e3 * kilometer
const Ly = 2.5e3 * kilometer
const Lz = 1 * kilometer

#Eddy viscosities and diffusivities
const νh = 0 * (meter^2/second)
const νv = 0 * (meter^2/second)
const κh = 0 * (meter^2/second)
const κv = 0 * (meter^2/second)

#Latitude (deg. N)
const lat = 74.0

#f-plane and Coriolis frequency
fPlane  = FPlane(latitude = lat)
const f = fPlane.f

#Gyre scales
const σr = 250 * kilometer
const σz = "infinity" ##300 * meter

#Speed and buoyancy frequency at surface of gyre
const U   = 1.5e-1 * (meter/second)
const N²₀ = 3e-3 * (second^(-2))

#Max buoyancy frequency (equal to N²₀ for uniform stratification)
const N²_max = 3e-3 * (second^(-2))

#Mixed-layer depth
const d_ML = -50 * meter

#Time-stepping parameters
const Δti     = 10 * second
const Δt_max  = 3 * hour
const CFL     = 0.2
const tf      = 20 * day #30 * day
const Δt_save = 6 * hour

#Architecture
const use_GPU = true

#Max. relative magnitude of initial u-perturbations
const max_u′ = 1e-8

#Whether to run visualization functions
const vis_const_x = false
const vis_const_y = false
const vis_const_z = false
const vis_norms   = true
const vis_z_grid  = false #Can only be done on CPU

#Indices at which to plot fields
const x_idx      = 259
const y_idx      = 259
const z_idx      = 11
const t_idx_skip = 1

#Seeds for 2 random-number generators
const seed1 = 12345
const seed2 = 67890

if !isnothing(seed1)
   Random.seed!(seed1)
end

#########################
# SET UP GRID AND MODEL #
#########################

use_GPU ? architecture = GPU() : architecture = CPU()

#z_grid_spacing(k) = chebyshev_spaced_faces(k, -Lz, Nz; ξ_centre = d_ML)

grid = RectilinearGrid(architecture,
		       topology = (Bounded, Bounded, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lx/2, Lx/2), 
                       y = (-Ly/2, Ly/2), 
                       z = (-Lz, 0.0),
		       halo = (3, 3, 3))
#                       z = z_grid_spacing)

closure = (HorizontalScalarDiffusivity(ν = νh, κ = κh), 
	   VerticalScalarDiffusivity(ν = νv, κ = κv))

const bkgd_N²_top = N²₀ #lognormal_strat(N²₀, N²_max, d_ML, 0)[1]
const bkgd_N²_bot = N²₀ #lognormal_strat(N²₀, N²_max, d_ML, -Lz)[1]

b̄, ū, v̄, ūφ_abs, b̄_BCs = bkgd_fields(f, σr, σz, U, 
				     bkgd_N²_top, bkgd_N²_bot)
model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :RungeKutta3,
                            advection = UpwindBiasedFifthOrder(), 
                            coriolis = fPlane,
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
			    boundary_conditions = (; b = b̄_BCs,))#,
#			    closure = closure)

set!(model, u = ū, v = v̄)
set!(model, b = b̄)

#Prints warnings if the respective instabilities are present
check_inert_stability(model.grid, f, model.velocities.u, model.velocities.v;
		      plot_ζz_abs = false, z_idx = z_idx)
check_grav_stability(model.tracers.b; plot_∂b∂z = false, grid = model.grid,
                     x_idx = x_idx)

#########################
# SAVE BACKGROUND STATE #
#########################

datetimestart = now()
datetimenow   = "250816-182636" #format(datetimestart, "yymmdd-HHMMSS")
print("Date-time label: $(datetimenow)", "\n")

#Create fields to store background state
Ux = Field{Face, Center, Center}(model.grid)
Uy = Field{Center, Face, Center}(model.grid)
Ur = Field{Center, Center, Center}(model.grid)
Uφ = Field{Center, Center, Center}(model.grid)
Uz = Field{Center, Center, Face}(model.grid)
B  = Field{Center, Center, Center}(model.grid)

Ur_vals, Uφ_vals = xy_vector_to_rφ(model.velocities.u, 
				   model.velocities.v, model.grid)

#Save background-state data
Ux .= model.velocities.u
Uy .= model.velocities.v
Ur .= Ur_vals
Uφ .= Uφ_vals
Uz .= model.velocities.w
B  .= model.tracers.b

@inline perturbation_norm(field, bkgd_field) = norm(field - bkgd_field)

#############################
# SET UP AND RUN SIMULATION #
#############################

#Perturb velocity components to trigger BCI
#=
@inline u_perturbed(x, y, z) = (ū(x, y, z) 
				+ (2 * (rand() - 0.5)) * (max_u′ / sqrt(2)))

if !isnothing(seed2)
   Random.seed!(seed2) #Update seed so next random number is independent
end

@inline v_perturbed(x, y, z) = (v̄(x, y, z) 
				+ (2 * (rand() - 0.5)) * (max_u′ / sqrt(2)))

set!(model, u = u_perturbed, v = v_perturbed) #Update initial condition to trigger BCI

simulation = Simulation(model, Δt = Δti, stop_time = tf)

#wizard = TimeStepWizard(cfl = CFL, max_Δt = Δt_max)
#simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

function progress(sim)
   umax = maximum(abs, sim.model.velocities.u)
   wmax = maximum(abs, sim.model.velocities.w)
   bmax = maximum(abs, sim.model.tracers.b)
   @info @sprintf("Iter: %d; time: %.2e days; Δt: %s",
		  iteration(sim), (time(sim)/day),  prettytime(sim.Δt))
   @info @sprintf("max|u|: %.2e; max|w|: %.2e; max|b|: %.2e",
		  umax, wmax, bmax)
   @info @sprintf("norm u' = %.10e", norm(sim.model.velocities.u - Ux))
   return nothing
end

add_callback!(simulation, progress, TimeInterval(Δt_save))

ur, uφ = xy_vector_to_rφ(model.velocities.u, model.velocities.v, model.grid)

outputs = (ur = ur,
	   uφ = uφ,
	   ux = model.velocities.u,
	   uy = model.velocities.v,
	   uz = model.velocities.w,
	   b  = model.tracers.b)

outfilepath = joinpath("./Output", "output_$(datetimenow).nc")
mkpath(dirname(outfilepath)) #Make path if nonexistent

field_writer = NetCDFOutputWriter(model, 
				  outputs, 
                                  with_halos = true,
		                  filename = outfilepath, 
                                  schedule = TimeInterval(Δt_save),
				  file_splitting = FileSizeLimit(30GiB))

ux_perturbation_norm(model) = perturbation_norm(model.velocities.u, Ux)
uy_perturbation_norm(model) = perturbation_norm(model.velocities.v, Uy)
ur_perturbation_norm(model) = perturbation_norm(ur, Ur)
uφ_perturbation_norm(model) = perturbation_norm(uφ, Uφ)
uz_perturbation_norm(model) = perturbation_norm(model.velocities.w, Uz)
b_perturbation_norm(model)  = perturbation_norm(model.tracers.b, B)

scalar_diagnostics = (ux′_norm = ux_perturbation_norm,
		      uy′_norm = uy_perturbation_norm,
		      ur′_norm = ur_perturbation_norm,
		      uφ′_norm = uφ_perturbation_norm,
		      uz′_norm = uz_perturbation_norm,
		      b′_norm = b_perturbation_norm)

scalarfilepath = joinpath("./Output", "scalars_$(datetimenow).nc")
mkpath(dirname(scalarfilepath)) #Make path if nonexistent

scalar_writer = NetCDFOutputWriter(model, 
				   scalar_diagnostics, 
				   with_halos = true, 
				   filename = scalarfilepath, 
				   schedule = TimeInterval(Δt_save),
                                   file_splitting = FileSizeLimit(30GiB),
				   dimensions = (ux′_norm = (),
						 uy′_norm = (),
						 ur′_norm = (),
						 uφ′_norm = (),
						 uz′_norm = (),
						 b′_norm = ()))

simulation.output_writers[:field_writer] = field_writer
simulation.output_writers[:scalar_writer] = scalar_writer

run!(simulation)

duration = canonicalize(now() - datetimestart)

pad_filenames(datetimenow)

###############################
# SAVE PARAMETERS TO LOG FILE #
###############################

logfilepath = joinpath("./Logs", "log_$(datetimenow).txt")
mkpath(dirname(logfilepath)) #Make path if nonexistent

open(logfilepath, "w") do file
   write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
   write(file, "Lx, Ly, Lz = $(Lx), $(Ly), $(Lz) \n\n")
   write(file, "νh, νv, κh, κv = $(νh), $(νv), $(κh), $(κv) \n\n")
   write(file, "lat = $(lat) \n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "U, N²₀ = $(U), $(N²₀) \n")
   write(file, "Max. u' = $(max_u′) \n")
   write(file, "Random-number seeds = $(seed1), $(seed2) \n\n")
   #write(file, "Δti, Δt_max, Δt_save = $(Δti), $(Δt_max), $(Δt_save) \n")
   write(file, "Δt = $(Δti) \n")
   write(file, "CFL = $(CFL) \n")
   write(file, "tf = $(tf) \n\n")
   write(file, "Total number of iterations = $(iteration(simulation)) \n")
   write(file, "Δtf = $(prettytime(simulation.Δt)) \n\n")
   write(file, "Simulation runtime = $(duration) \n")
   write(file, "Output filesize = $(pretty_filesize(filesize(outfilepath)))")
end
=#
###################################
# RUN VISUALIZATION, IF INDICATED #
###################################

if vis_const_x
   visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny; 
		      bkgd_datetime = bkgd_datetime,
		      x_idx = x_idx,
		      plot_animation = true,
                      t_idx_skip = t_idx_skip)
   visualize_fields_const_x(datetimenow, x_idx; 
			    bkgd_datetime = bkgd_datetime,
                            plot_animation = true, 
			    t_idx_skip = t_idx_skip)
end

if vis_const_y
   visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny;
		      bkgd_datetime = bkgd_datetime,
                      y_idx = y_idx, 
		      plot_animation = true,
                      t_idx_skip = t_idx_skip)
end

if vis_const_z
   #visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny; 
   #		      bkgd_datetime = bkgd_datetime,
   #                   z_idx = z_idx,
   #		      plot_animation = true, 
   #		      t_idx_skip = t_idx_skip)
   visualize_fields_const_z(datetimenow, z_idx,  B, Uφ; 
			    plot_animation = true, t_idx_skip = t_idx_skip)
end

if vis_norms
   visualize_norms(datetimenow)
end

if vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end
