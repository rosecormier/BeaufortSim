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
const Nx = 1024
const Ny = 1024
const Nz = 512

#Lengths of axes
const Lx = 4000 * kilometer
const Ly = 4000 * kilometer
const Lz = 2000 * meter

#Eddy viscosities
const νh = (5e-2) * (meter^2/second)
const νv = (5e-5) * (meter^2/second)

#Latitude (deg. N)
const lat = 74.0

#f-plane and Coriolis frequency
fPlane  = FPlane(latitude = lat)
const f = fPlane.f

#Gyre scales
const σr  = 250 * kilometer
const σz  = 300 * meter

#Gyre speed and buoyancy frequency
const U  = 0.75 * U_upper_bound(σr, f) * (meter/second)
const N2 = 1.1 * N2_lower_bound(σr, σz, f, U) * (second^(-2))

#Time increments
const Δti     = 0.005 * second
const Δt_max  = 200 * second 
const CFL     = 0.1
const tf      = 5 * day
const Δt_save = 30 * minute

#Architecture
const use_GPU = true

#Whether to run visualization functions
const do_vis_const_x = true
const do_vis_const_y = false
const do_vis_const_z = true

#Indices at which to plot fields
const x_idx = 512
const y_idx = 256
const z_idx = 500

##############################
# INSTANTIATE GRID AND MODEL #
##############################

use_GPU ? architecture = GPU() : architecture = CPU()

grid = RectilinearGrid(architecture,
		       topology = (Periodic, Periodic, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lx/2, Lx/2), 
                       y = (-Ly/2, Ly/2), 
                       z = (-Lz, 0),
                       halo = (3,3,3))

closure = (HorizontalScalarDiffusivity(ν = νh), 
	   VerticalScalarDiffusivity(ν = νv))

b_BCs = FieldBoundaryConditions(bottom = GradientBoundaryCondition(N2))

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :QuasiAdamsBashforth2, 
                            advection = UpwindBiasedFifthOrder(),
                            closure = closure, 
                            coriolis = fPlane,
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
			    boundary_conditions = (b = b_BCs,))

#Prints warnings if the respective instabilities are present
check_inert_stability(σr, σz, f, U,
                      xnodes(model.grid, Face(), Face(), Face()),
                      ynodes(model.grid, Face(), Face(), Face()),
                      znodes(model.grid, Face(), Face(), Face()))
check_grav_stability(σr, σz, f, U, N2)

##########################
# SET INITIAL CONDITIONS #
##########################

#Prints warnings if the respective instabilities are present
check_inert_stability(σr, σz, f, U, 
		      xnodes(model.grid, Face(), Face(), Face()), 
		      ynodes(model.grid, Face(), Face(), Face()), 
		      znodes(model.grid, Face(), Face(), Face()))
check_grav_stability(σr, σz, f, U, N2)

b       = model.tracers.b
u, v, w = model.velocities

ū(x,y,z)  = (U*y/σr) * exp(1 - (x^2 + y^2)/(σr^2) - (z/σz)^2)
v̄(x,y,z)  = -(U*x/σr) * exp(1 - (x^2 + y^2)/(σr^2) - (z/σz)^2)
bʹ(x,y,z) = (1e-6) * rand()
b̄(x,y,z)  = (N2*z - (U*f*σr/(σz^2)) * z * exp(1 - (x^2 + y^2)/(σr^2) 
					      -(z/σz)^2)) #+ bʹ(x,y,z))

set!(model, u = ū, v = v̄, b = b̄)

#############################
# SET UP AND RUN SIMULATION #
#############################

simulation = Simulation(model, Δt = Δti, stop_time = tf)

wizard = TimeStepWizard(cfl = CFL, max_Δt = Δt_max)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

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
   write(file, "νh, νv = $(νh), $(νv) \n\n")
   write(file, "lat = $(lat) \n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "U, N2 = $(U), $(N2) \n\n")
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
   visualize_fields_const_x(datetimenow, x_idx)
   #visualize_q_const_x(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, x_idx)
   #plot_background_ζa(datetimenow, U, f, σr, σz; x_idx = x_idx)
end

if do_vis_const_y
   visualize_fields_const_y(datetimenow, y_idx)
   #visualize_q_const_y(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, y_idx)
   #plot_background_ζa(datetimenow, U, f, σr, σz; y_idx = y_idx)
end

if do_vis_const_z
   visualize_fields_const_z(datetimenow, z_idx)
   #visualize_q_const_z(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, z_idx)
end
