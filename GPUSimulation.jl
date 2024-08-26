include("Library.jl")

using Dates: format, now
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Coriolis
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.TurbulenceClosures
using Oceananigans.Units

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
const tf      = 20 * day
const Δt_save = 1 * hour

#Architecture
const use_GPU = true

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

#####################
# SET UP SIMULATION #
#####################

simulation = Simulation(model, Δt = Δti, stop_time = tf)

wizard = TimeStepWizard(cfl = CFL, max_Δt = Δt_max)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

progress(sim) = @info string(
    "Iteration: $(iteration(sim)), time: $(time(sim)), Δt: $(sim.Δt)")
add_callback!(simulation, progress, IterationInterval(100))

outputs = Dict("u" => model.velocities.u,
	       "v" => model.velocities.v,
	       "w" => model.velocities.w,
	       "b" => model.tracers.b)
timenow = format(now(), "yymmdd-HHMMSS")

output_filename = "output_$(timenow).nc"
output_filepath = joinpath("./Output", output_filename)
mkpath(dirname(output_filepath)) #Make path if nonexistent

outputwriter = NetCDFOutputWriter(model, 
				  outputs, 
                                  with_halos = true,
		                  filename = output_filepath, 
                                  schedule = TimeInterval(Δt_save))

simulation.output_writers[:field_writer] = outputwriter

run!(simulation)
print(timenow, "\n")
