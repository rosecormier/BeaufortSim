include("LibraryDynamics.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Visualization.jl")

using Dates: canonicalize, format, now
using NCDatasets
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using Oceananigans.Utils
using Printf, Random

######################
# SPECIFY PARAMETERS #
######################

#Numbers of gridpoints
const Nx = 512 #1024
const Ny = 512 #1024
const Nz = 16 #256

#Lengths of axes
const Lx = 4e3 * kilometer #2e3 * kilometer
const Ly = 4e3 * kilometer #2e3 * kilometer
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
const σz = "infinity" #300 * meter

#Speed and buoyancy frequency at surface of gyre
const U   = 1.5e-1 * (meter/second)
const N²₀ = 3e-3 * (second^(-2))

#Max buoyancy frequency (equal to N²₀ for uniform stratification)
const N²_max = 3e-3 * (second^(-2))

#Mixed-layer depth
const d_ML = -50 * meter

#Time-stepping parameters
const Δti     = 5 * second
const Δt_max  = 1 * hour
const CFL     = 0.2
const tf      = 60 * day
const Δt_save = 6 * hour

#Architecture
const use_GPU = true

#Max. relative magnitude of initial u-perturbations
const max_u′ = 1e-5

const save_bkgd = true #Whether to save background state to a NetCDF file
const bkgd_datetime = nothing #If save_bkgd == true, must == nothing

#Whether to run visualization functions
const vis_const_x = false
const vis_const_y = false
const vis_const_z = true
const vis_norms   = true
const vis_z_grid  = false #Can only be done on CPU

#Indices at which to plot fields
const x_idx      = 259
const y_idx      = 259
const z_idx      = 9 #253
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
#                       z = z_grid_spacing,
#                       halo = (3, 3, 3))

closure = (HorizontalScalarDiffusivity(ν = νh, κ = κh), 
	   VerticalScalarDiffusivity(ν = νv, κ = κv))

#const bkgd_N²_top    = lognormal_strat(N²₀, N²_max, d_ML, 0)[1]
#const bkgd_N²_bottom = lognormal_strat(N²₀, N²_max, d_ML, -Lz)[1]

const bkgd_N²_top = N²₀
const bkgd_N²_bot = N²₀

b̄, ū, v̄, ūφ_abs, b̄_BCs = bkgd_fields(f, σr, σz, U, 
				     bkgd_N²_top, bkgd_N²_bot)

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :RungeKutta3,
                            advection = UpwindBiasedFifthOrder(),
                            closure = closure, 
                            coriolis = fPlane,
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
			    boundary_conditions = (; b = b̄_BCs,))

b       = model.tracers.b
u, v, w = model.velocities

set!(model, u = ū, v = v̄, b = b̄)

#Prints warnings if the respective instabilities are present
check_inert_stability(model.grid, f, model.velocities.u, model.velocities.v;
		      plot_ζz_abs = false, z_idx = z_idx)
check_grav_stability(model.tracers.b; plot_∂b∂z = false, grid = model.grid,
                     x_idx = x_idx)

#######################################
# SAVE BACKGROUND STATE, IF INDICATED #
#######################################

datetimestart = now()
datetimenow   = format(datetimestart, "yymmdd-HHMMSS")
print("Date-time label: $(datetimenow)", "\n")

if save_bkgd
   
   bkgd_simulation = Simulation(model, Δt = Δti, stop_iteration = 1)

   ur, uφ = xy_vector_to_rφ(model.velocities.u, 
			    model.velocities.v, 
			    model.grid)

   bkgd_outputs = (Ur = ur,
                   Uφ = uφ,
                   Ux = model.velocities.u,
                   Uy = model.velocities.v,
                   Uz = model.velocities.w,
                   B  = model.tracers.b)

   bkgd_filepath = joinpath("./Output", "bkgd_$(datetimenow).nc")
   mkpath(dirname(bkgd_filepath)) #Make path if nonexistent

   bkgd_writer = NetCDFOutputWriter(model,
                                    bkgd_outputs,
                                    with_halos = true,
                                    filename = bkgd_filepath,
                                    schedule = SpecifiedTimes([0]))

   bkgd_simulation.output_writers[:field_writer] = bkgd_writer
   run!(bkgd_simulation)
end

#############################
# SET UP AND RUN SIMULATION #
#############################

@inline speed_perturb(x, y, z) = (2 * (rand() - 0.5)) * max_u′ * ūφ_abs(x, y, z)

if !isnothing(seed2)
   Random.seed!(seed2)
end

@inline direction_perturb(x, y, z) = 2pi * rand()

@inline u_perturbed(x, y, z) = ū(x, y, z) + (speed_perturb(x, y, z) * cos(direction_perturb(x, y, z)))
@inline v_perturbed(x, y, z) = v̄(x, y, z) + (speed_perturb(x, y, z) * sin(direction_perturb(x, y, z)))

set!(model, u = u_perturbed, v = v_perturbed) #Update initial condition to trigger BCI

simulation = Simulation(model, Δt = Δti, stop_time = tf)

wizard = TimeStepWizard(cfl = CFL, max_Δt = Δt_max)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

function progress(sim)
   umax = maximum(abs, sim.model.velocities.u)
   wmax = maximum(abs, sim.model.velocities.w)
   bmax = maximum(abs, sim.model.tracers.b)
   @info @sprintf("Iter: %d; time: %.2e days; Δt: %s",
		  iteration(sim), (time(sim)/day),  prettytime(sim.Δt))
   @info @sprintf("max|u|: %.2e; max|w|: %.2e; max|b|: %.2e",
		  umax, wmax, bmax)
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

simulation.output_writers[:field_writer] = field_writer

run!(simulation)
duration = canonicalize(now() - datetimestart)

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
   write(file, "Δti, Δt_max, Δt_save = $(Δti), $(Δt_max), $(Δt_save) \n")
   write(file, "CFL = $(CFL) \n")
   write(file, "tf = $(tf) \n\n")
   write(file, "Total number of iterations = $(iteration(simulation)) \n")
   write(file, "Δtf = $(prettytime(simulation.Δt)) \n\n")
   write(file, "Simulation runtime = $(duration) \n")
   write(file, "Output filesize = $(pretty_filesize(filesize(outfilepath)))")
end

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
   visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny; 
		      bkgd_datetime = bkgd_datetime,
                      z_idx = z_idx,
		      plot_animation = true, 
		      t_idx_skip = t_idx_skip)
   visualize_fields_const_z(datetimenow, z_idx; 
			    bkgd_datetime = bkgd_datetime,
                            plot_animation = true, 
			    t_idx_skip = t_idx_skip)
end

if vis_norms
   visualize_norms(datetimenow, model.grid; bkgd_datetime = bkgd_datetime)
end

if vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end
