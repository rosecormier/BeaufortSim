include("LibraryDynamics.jl")
include("LibraryStability.jl")
include("Visualization.jl")

using Dates: canonicalize, format, now
using Oceananigans
#using Oceananigans.AbstractOperations, Oceananigans.Fields 
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Coriolis
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using Printf, Random

######################
# SPECIFY PARAMETERS #
######################

#Numbers of gridpoints
const Nx = 512
const Ny = 512
const Nz = 256

#Lengths of axes
const Lx = 2000 * kilometer
const Ly = 2000 * kilometer
const Lz = 1000 * meter

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
@printf("f = %.2e Hz \n", f)

#Gyre scales
const σr = 250 * kilometer
const σz = 300 * meter

#Speed and buoyancy frequency at surface of gyre
const U   = 1.5e-1 * (meter/second)
const N²₀ = 3e-3 * (second^(-2))
@printf("Bu = %.2e \n", compute_Bu(σr, σz, f, N²₀))

#Max buoyancy frequency (equal to N²₀ for uniform stratification)
const N²_max = 3e-3 * (second^(-2))

#Mixed-layer depth
const d_ML = -50 * meter

#Time-stepping parameters
const Δti     = 5 * second
const Δt_max  = 1200 * second 
const CFL     = 0.2
const tf      = 3 * day
const Δt_save = 12 * hour

#Architecture
const use_GPU = true

#Max. magnitude of initial b-perturbations (0 for no perturbation)
const max_b′ = 4e-3 * (meter/(second^2))

#Whether to run visualization functions
const do_vis_const_x = false
const do_vis_const_y = false
const do_vis_const_z = false
const do_vis_norms   = false
const do_vis_z_grid  = false #Can only be done on CPU

#Indices at which to plot fields
const x_idx      = 259
const y_idx      = 259
const z_idx      = 253
const t_idx_skip = 1

#Seed for random-number generator
const seed = 12345

if !isnothing(seed)
   Random.seed!(seed)
end

##############################
# INSTANTIATE GRID AND MODEL #
##############################

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
const bkgd_N²_bottom = N²₀

@inline dbdz_top(x, y, t)    = (bkgd_N²_top
				.+ (sqrt(2) * f * U * σr / (σz^2)
				   * exp(1/2) 
				   * (1 - exp(-(x^2 + y^2)/(σr^2)))))
@inline dbdz_bottom(x, y, t) = (bkgd_N²_bottom
				.+ (sqrt(2) * f * U * σr / (σz^2)
				   * exp((1/2) - (Lz/σz)^2) 
			      	   * (1 - exp(-(x^2 + y^2)/(σr^2))) 
			           * (1 - 2 * (Lz/σz)^2)))

b_top_BC    = GradientBoundaryCondition(dbdz_top)
b_bottom_BC = GradientBoundaryCondition(dbdz_bottom)

b_BCs = FieldBoundaryConditions(top = b_top_BC, bottom = b_bottom_BC)

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :RungeKutta3,
                            advection = UpwindBiasedFifthOrder(),
                            closure = closure, 
                            coriolis = fPlane,
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
			    boundary_conditions = (; b = b_BCs,))

#phi_KernOp = KernelFunctionOperation{Center, Center, Center}(φ, model.grid)
#phi = Field(phi_KernOp)
#compute!(phi)
compute_polar_coords(model.grid)
##########################
# SET INITIAL CONDITIONS #
##########################

b       = model.tracers.b
u, v, w = model.velocities

ū(x, y, z) = ((sqrt(2)*U*y/σr) 
	       * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))
v̄(x, y, z) = -((sqrt(2)*U*x/σr) 
	        * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))
b̄(x, y, z) = (lognormal_strat(N²₀, N²_max, d_ML, z)[2]
	      + ((sqrt(2)*f*U*σr*z/(σz^2)) 
	          * exp((1/2) - (z/σz)^2) 
		  * (1 - exp(-(x^2 + y^2)/(σr^2)))
		  * (1 - ((sqrt(2)*U/(f*σr)) * exp((1/2) - (z/σz)^2)
			  * (1 + exp(-(x^2 + y^2)/(σr^2)))
			  )
		     )
		  )
	      )

set!(model, u = ū, v = v̄, b = b̄)

#Prints warnings if the respective instabilities are present
check_inert_stability(model.grid, f, model.velocities.u, model.velocities.v;
		      plot_ζz_abs = false, z_idx = z_idx)
check_grav_stability(model.tracers.b; plot_∂b∂z = false, grid = model.grid,
                     x_idx = x_idx)

b_perturbed(x, y, z) = (max_b′ * rand()) + b̄(x, y, z) 

set!(model, b = b_perturbed)

#############################
# SET UP AND RUN SIMULATION #
#############################
#=
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

outputs = (u = model.velocities.u,
	   v = model.velocities.v,
	   w = model.velocities.w,
	   b = model.tracers.b)

datetimestart = now()
=#
datetimenow   = "250425-143740" #format(datetimestart, "yymmdd-HHMMSS")
outfilename   = "output_$(datetimenow).nc"
outfilepath   = joinpath("./Output", outfilename)
#=mkpath(dirname(outfilepath)) #Make path if nonexistent

outputwriter = NetCDFOutputWriter(model, 
				  outputs, 
                                  with_halos = true,
		                  filename = outfilepath, 
                                  schedule = TimeInterval(Δt_save))

simulation.output_writers[:field_writer] = outputwriter

run!(simulation)
print("Date-time label: $(datetimenow)", "\n")

duration = canonicalize(now() - datetimestart)

###############################
# SAVE PARAMETERS TO LOG FILE #
###############################

logfilename = "log_$(datetimenow).txt"
logfilepath = joinpath("./Logs", logfilename)
mkpath(dirname(logfilepath)) #Make path if nonexistent

open(logfilepath, "w") do file
   write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
   write(file, "Lx, Ly, Lz = $(Lx), $(Ly), $(Lz) \n\n")
   write(file, "νh, νv, κh, κv = $(νh), $(νv), $(κh), $(κv) \n\n")
   write(file, "lat = $(lat) \n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "U, N²₀ = $(U), $(N²₀) \n")
   write(file, "Computed Bu = $(compute_Bu(σr, σz, f, N²₀)) \n\n")
   write(file, "Max. b', random-number seed = $(max_b′), $(seed) \n\n")
   write(file, "Δti, Δt_max, Δt_save = $(Δti), $(Δt_max), $(Δt_save) \n")
   write(file, "CFL = $(CFL) \n")
   write(file, "tf = $(tf) \n\n")
   write(file, "Total number of iterations = $(iteration(simulation)) \n")
   write(file, "Δtf = $(prettytime(simulation.Δt)) \n\n")
   write(file, "Simulation runtime = $(duration) \n")
   write(file, "Output filesize = $(filesize(outfilepath)) bytes")
end
=#
###################################
# RUN VISUALIZATION, IF INDICATED #
###################################

if do_vis_const_x
   visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny;
                      x_idx = x_idx, plot_animation = true,
                      t_idx_skip = t_idx_skip)
   #plot_background_ζa(datetimenow, U, f, σr, σz; x_idx = x_idx)
end

if do_vis_const_y
   visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny;
                      y_idx = y_idx, plot_animation = true,
                      t_idx_skip = t_idx_skip)
end

if do_vis_const_z
   visualize_b_and_ωz(datetimenow, Lx/Nx, Ly/Ny;
                      z_idx = z_idx, plot_animation = true, 
		      t_idx_skip = t_idx_skip)
end

if do_vis_norms
   visualize_norms(datetimenow; f = f)
end

if do_vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end
