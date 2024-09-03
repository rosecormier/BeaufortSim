include("Visualization.jl")

using Dates: format, now
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Coriolis
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using Printf

######################
# SPECIFY PARAMETERS #
######################

#Numbers of gridpoints
const Nx = 512
const Ny = 512
const Nz = 128

#Lengths of axes
const Lx = 2000 * kilometer
const Ly = 2000 * kilometer
const Lz = 1000 * meter

#Eddy viscosities
const νh = (5e-2) * (meter^2/second)
const νv = (5e-5) * (meter^2/second)

#Gyre scales
const lat = 74.0  #Degrees N
const U   = 1.0 * meter/second
const σr  = 250 * kilometer
const σz  = 300 * meter
const N2  = (3e-4) * (1/second^2)

#Time increments
const Δti     = 1 * second
const Δt_max  = 100 * second 
const CFL     = 0.2
const tf      = 1 * day
const Δt_save = 1 * hour

#Architecture
const use_GPU = true

#Whether to run visualization functions
const do_vis_const_x = true
const do_vis_const_y = true
const do_vis_const_z = true

#Indices at which to plot fields
const x_idx = 256
const y_idx = 256
const z_idx = 20

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

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :QuasiAdamsBashforth2, 
                            advection = UpwindBiasedFifthOrder(),
                            closure = closure, 
                            coriolis = FPlane(latitude = lat),
                            tracers = (:b),
                            buoyancy = BuoyancyTracer())

##########################
# SET INITIAL CONDITIONS #
##########################

f       = model.coriolis.f
b       = model.tracers.b
u, v, w = model.velocities

ū(x,y,z)  = (U*y/σr) * exp(-(x^2 + y^2)/(σr^2) - (z/σz)^2)
v̄(x,y,z)  = (-U*x/σr) * exp(-(x^2 + y^2)/(σr^2) - (z/σz)^2)
bʹ(x,y,z) = (1e-4) * rand()
b̄(x,y,z)  = (N2*z 
	     - (σr*f*U*z/(σr^2)) * exp(-(x^2 + y^2)/(σr^2) - (z/σz)^2)
	     + bʹ(x,y,z))

set!(model, u = ū, v = v̄, b = b̄)

#############################
# SET UP AND RUN SIMULATION #
#############################

simulation = Simulation(model, Δt = Δti, stop_time = tf)

wizard = TimeStepWizard(cfl = CFL, max_Δt = Δt_max)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

function progress(sim)
   umax = maximum(abs, sim.model.velocities.u)
   wmax = maximum(abs, sim.model.velocities.w)
   bmax = maximum(abs, sim.model.tracers.b)
   @info @sprintf("Iter: %d; time: %.2e; Δt: %s",
		  iteration(sim), time(sim), prettytime(sim.Δt))
   @info @sprintf("max|u|: %.2e; max|w|: %.2e; max|b|: %.2e",
		  umax, wmax, bmax)
   return nothing
end

add_callback!(simulation, progress, TimeInterval(Δt_save))

outputs = (u = model.velocities.u,
	   v = model.velocities.v,
	   w = model.velocities.w,
	   b = model.tracers.b)

datetimenow = format(now(), "yymmdd-HHMMSS")
outfilename = "output_$(datetimenow).nc"
outfilepath = joinpath("./Output", outfilename)
mkpath(dirname(outfilepath)) #Make path if nonexistent

outputwriter = NetCDFOutputWriter(model, 
				  outputs, 
                                  with_halos = true,
		                  filename = outfilepath, 
                                  schedule = TimeInterval(Δt_save))

simulation.output_writers[:field_writer] = outputwriter

run!(simulation)
print("Date-time label: $(datetimenow)", "\n")

###############################
# SAVE PARAMETERS TO LOG FILE #
###############################

logfilename = "log_$(datetimenow).txt"
logfilepath = joinpath("./Logs", logfilename)
mkpath(dirname(logfilepath)) #Make path if nonexistent

open(logfilepath, "w") do file
   write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
   write(file, "Lx, Ly, Lz = $(Lx), $(Ly), $(Lz) \n")
   write(file, "νh, νv = $(νh), $(νv) \n")
   write(file, "lat = $(lat) \n")
   write(file, "U = $(U) \n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "N2 = $(N2) \n")
   write(file, "Δti, Δt_max, Δt_save = $(Δti), $(Δt_max), $(Δt_save) \n")
   write(file, "CFL = $(CFL) \n")
   write(file, "tf = $(tf)")
end

###################################
# RUN VISUALIZATION, IF INDICATED #
###################################

if do_vis_const_x
   visualize_perturbs_const_x(datetimenow, x_idx)
   visualize_q_const_x(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, x_idx)
end

if do_vis_const_y
   visualize_perturbs_const_y(datetimenow, y_idx)
   visualize_q_const_y(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, y_idx)
end

if do_vis_const_z
   visualize_perturbs_const_z(datetimenow, z_idx)
   visualize_q_const_z(datetimenow, Lx/Nx, Ly/Ny, Lz/Nz, f, z_idx)
end
