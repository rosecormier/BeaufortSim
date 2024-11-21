include("LibraryStability.jl")
include("Visualization.jl")

using Dates: canonicalize, format, now
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Coriolis
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using Printf
using .Stability

######################
# SPECIFY PARAMETERS #
######################

#Numbers of gridpoints
const Nx = 256
const Ny = 256
const Nz = 256

#Lengths of axes
const Lx = 2000 * kilometer
const Ly = 2000 * kilometer
const Lz = 1000 * meter

#Eddy viscosities
const νh = 5e-2 * (meter^2/second)
const νv = 5e-5 * (meter^2/second)
const κh = 5e-2 * (meter^2/second)
const κv = 5e-5 * (meter^2/second)

#Latitude (deg. N)
const lat = 74.0

#f-plane and Coriolis frequency
fPlane  = FPlane(latitude = lat)
const f = fPlane.f
@printf("f = %.2e Hz \n", f)

#Gyre scales
const σr = 250 * kilometer
const σz = 300 * meter

#Gyre speed and buoyancy frequency
const U  = 0.1 * (meter/second)
const N2 = 1e-3 * (second^(-2))
@printf("Bu = %.2e \n", compute_Bu(σr, σz, f, N2))

#Time-stepping parameters
const Δti     = 1 * second
const Δt_max  = 120 * second 
const CFL     = 0.1
const tf      = 5 * day
const Δt_save = 1 * hour

#Architecture
const use_GPU = true

#Max. magnitude of initial b-perturbations (0 for no perturbation)
const max_b′ = 0 # 1e-3

#Whether to run visualization functions
const do_vis_const_x = true
const do_vis_const_y = false
const do_vis_const_z = true

#Indices at which to plot fields
const x_idx      = 131
const y_idx      = 259
const z_idx      = 252
const t_idx_skip = 1

##############################
# INSTANTIATE GRID AND MODEL #
##############################

use_GPU ? architecture = GPU() : architecture = CPU()

grid = RectilinearGrid(architecture,
		       topology = (Bounded, Bounded, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lx/2, Lx/2), 
                       y = (-Ly/2, Ly/2), 
                       z = (-Lz, 0),
                       halo = (3, 3, 3))

closure = (HorizontalScalarDiffusivity(ν = νh, κ = κh), 
	   VerticalScalarDiffusivity(ν = νv, κ = κv))

@inline dbdz_top(x, y, t)    = (N2 
				+ (sqrt(2) * f * U * σr / (σz^2)
				   * exp(1/2) 
				   * (1 - exp(-(x^2 + y^2)/(σr^2)))))
@inline dbdz_bottom(x, y, t) = (N2 
				+ (sqrt(2) * f * U * σr / (σz^2)
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

#Prints warnings if the respective instabilities are present
check_inert_stability(σr, σz, f, U,
                      xnodes(model.grid, Face(), Face(), Face()),
                      ynodes(model.grid, Face(), Face(), Face()),
                      znodes(model.grid, Face(), Face(), Face()))
check_grav_stability(σr, σz, f, U, N2,
		     xnodes(model.grid, Face(), Face(), Face()),
                     ynodes(model.grid, Face(), Face(), Face()),
                     znodes(model.grid, Face(), Face(), Face()))

##########################
# SET INITIAL CONDITIONS #
##########################

b       = model.tracers.b
u, v, w = model.velocities

ū(x,y,z) = ((sqrt(2) * U * y / σr) 
	    * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))
v̄(x,y,z) = -((sqrt(2) * U * x / σr) 
	     * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))

b′(x,y,z) = max_b′ * rand() * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2)
b̄(x,y,z)  = (N2 * z 
	     + (sqrt(2) * f * U * σr * z / (σz^2) 
	        * exp((1/2) - (z/σz)^2) 
		* (1 - exp(-(x^2 + y^2)/(σr^2))))
	     + b′(x,y,z))

set!(model, u = ū, v = v̄, b = b̄)

#############################
# SET UP AND RUN SIMULATION #
#############################

simulation = Simulation(model, Δt = Δti, stop_time = tf)

wizard = TimeStepWizard(cfl = CFL, max_Δt = Δt_max)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1)) #(8))

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
datetimenow   = format(datetimestart, "yymmdd-HHMMSS")
outfilename   = "output_$(datetimenow).nc"
outfilepath   = joinpath("./Output", outfilename)
mkpath(dirname(outfilepath)) #Make path if nonexistent

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
   write(file, "U, N2 = $(U), $(N2) \n")
   write(file, "Computed Bu = $(compute_Bu(σr, σz, f, N2)) \n\n")
   write(file, "Max. b' = $(max_b′) \n\n")
   write(file, "Δti, Δt_max, Δt_save = $(Δti), $(Δt_max), $(Δt_save) \n")
   write(file, "CFL = $(CFL) \n")
   write(file, "tf = $(tf) \n\n")
   write(file, "Total number of iterations = $(iteration(simulation)) \n")
   write(file, "Δtf = $(prettytime(simulation.Δt)) \n\n")
   write(file, "Simulation runtime = $(duration) \n")
   write(file, "Output filesize = $(filesize(outfilepath)) bytes")
end

###################################
# RUN VISUALIZATION, IF INDICATED #
###################################

if do_vis_const_x
   visualize_fields_const_x(datetimenow, x_idx; t_idx_skip = t_idx_skip)
   visualize_fields_const_x(datetimenow, x_idx; plot_animation = false)
   #visualize_q_const_x(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, x_idx)
   #plot_background_ζa(datetimenow, U, f, σr, σz; x_idx = x_idx)
end

if do_vis_const_y
   visualize_fields_const_y(datetimenow, y_idx; t_idx_skip = t_idx_skip)
   visualize_fields_const_y(datetimenow, y_idx; plot_animation = false)
   #visualize_q_const_y(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, y_idx)
   #plot_background_ζa(datetimenow, U, f, σr, σz; y_idx = y_idx)
end

if do_vis_const_z
   visualize_fields_const_z(datetimenow, z_idx; t_idx_skip = t_idx_skip)
   visualize_fields_const_z(datetimenow, z_idx; plot_animation = false)
   #visualize_q_const_z(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, z_idx)
end
