include("Library.jl")

using Dates
using Oceananigans
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Coriolis
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.TurbulenceClosures
#using Oceananigans.Units

# GET INPUTS

params = ReadInputFile("InputSimulation.txt")

const Nx = params["Nx"]
const Ny = params["Ny"]
const Nz = params["Nz"]

const Lx = params["Lx"] * 1e3
const Ly = params["Ly"] * 1e3
const Lz = params["Lz"] #Units m

const νh = params["h_diffusivity"]
const νv = params["v_diffusivity"]

const latitude = params["latitude"]

const p̃0 = params["p̃0"]
const σr = params["σr"]
const σz = params["σz"]
const N2 = params["N"]^2

const Δti         = params["Δt_initial"]
const stop_time   = params["stop_time"]
const save_interv = params["save_interval"] 

# SIMULATION GEOMETRY

grid = RectilinearGrid(GPU(), 
                  size = (Nx, Ny, Nz), 
                  x = (-Lx/2, Lx/2), 
                  y = (-Ly/2, Ly/2), 
                  z = (-Lz, 0),
                  topology = (Periodic, Periodic, Bounded),
                  halo = (3,3,3))

# INSTANTIATE MODEL WITH ROTATION

model = NonhydrostaticModel(; 
           grid = grid, 
           timestepper = :QuasiAdamsBashforth2, 
           advection = UpwindBiasedFifthOrder(),
           closure = (HorizontalScalarDiffusivity(ν = νh), 
              VerticalScalarDiffusivity(ν = νv)), 
           coriolis = FPlane(latitude = latitude),
           tracers = (:b),
           buoyancy = BuoyancyTracer())

# SET INITIAL CONDITIONS

f = model.coriolis.f

b = model.tracers.b
u, v, w = model.velocities

p̃(x,y,z) = p̃0 * exp(-(x^2+y^2)/σr^2 - (z/σz)^2)
ū(x,y,z) = 2 * y * p̃(x,y,z) / (f*σr^2)
v̄(x,y,z) = -2 * x * p̃(x,y,z) / (f*σr^2)
b̄(x,y,z) = N2*z - (2 * z * p̃(x,y,z) / σz^2)

bi(x,y,z) = b̄(x,y,z) + (1e-5)*rand()

set!(model, u = ū, v = v̄, b = bi)

simulation = Simulation(model, Δt = Δti, stop_time = stop_time)

wizard = TimeStepWizard(cfl = 0.2, max_Δt = 100)
simulation.callbacks[:wizard] = Callback(wizard, 
				         IterationInterval(100))

progress(sim) = @info string(
    "Iteration: $(iteration(sim)), time: $(time(sim)), Δt: $(sim.Δt)")
add_callback!(simulation, progress, IterationInterval(100))

ζz_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
ζz = Field(ζz_op)
#b_perturb = b - b̄

outputs = Dict("u" => model.velocities.u,
	       "v" => model.velocities.v,
	       "w" => model.velocities.w,
	       "b" => model.tracers.b,
	       "ζz" => ζz)
timenow = Dates.format(now(), "yymmdd-HHMMSS")

output_filename = "output_$(timenow).nc"
output_filepath = joinpath("./Output", output_filename)
mkpath(dirname(output_filepath)) #Make directory if nonexistent

simulation.output_writers[:field_writer] = NetCDFOutputWriter(
                model, 
		outputs, 
                with_halos = true,
		filename = output_filepath, 
                schedule = TimeInterval(save_interv))

run!(simulation)
print(timenow, "\n")
