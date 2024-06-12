using Dates
using Oceananigans, Oceananigans.Coriolis, Oceananigans.TurbulenceClosures
using Oceananigans.Units

const Nx = 256
const Ny = 256
const Nz = 64

const Lx = 2000kilometer
const Ly = 2000kilometer
const Lz = 1000meter

const νₕ = 5e-1                  # m^2/s
const νᵥ = 5e-5                  # m^2/s

const latitude = 74.0            # Degrees N

const p̃0 = 100                   # Reference [reduced] pressure (m^2/s^2)
const σₕ = 250                   # Radial gyre lengthscale (km)
const σᵥ = 300                   # Vertical gyre lengthscale (m)
const N  = 5e-3                  # BV frequency (1/s); cf. Jackson et al, 2012

const Δtᵢ           = 60         # Initial timestep (s)
const stop_time     = 432000     # Wall time (s; simulation time)
const save_interval = 3600       # Interval at which to save model fields (s; simulation time)

grid = RectilinearGrid(GPU(), 
                      size = (Nx, Ny, Nz), 
                         x = (-Lx/2, Lx/2), 
                         y = (-Ly/2, Ly/2), 
                         z = (-Lz, 0),
                  topology = (Periodic, Periodic, Bounded),
                      halo = (3,3,3))

model = NonhydrostaticModel(; grid = grid, 
    timestepper = :QuasiAdamsBashforth2, 
      advection = UpwindBiasedFifthOrder(),
        closure = (
            HorizontalScalarDiffusivity(ν = νₕ), 
            VerticalScalarDiffusivity(ν = νᵥ)), 
      coriolis = FPlane(latitude = latitude),
       tracers = (:b),
      buoyancy = BuoyancyTracer())

N² = N^2
f  = model.coriolis.f
b  = model.tracers.b
u, v, w = model.velocities

p̃(x,y,z) =  p̃0 * exp(-(x^2+y^2)/σₕ^2 - (z/σᵥ)^2)
ū(x,y,z) =   2 * y / (f * σₕ^2) * p̃(x,y,z)
v̄(x,y,z) = - 2 * x / (f * σₕ^2) * p̃(x,y,z)
b̄(x,y,z) = - 2 * z / σᵥ^2 * p̃(x,y,z) + N² * z

bᵢ(x,y,z) = b̄(x,y,z) + 1e-5*rand()

set!(model, u = ū, v = v̄, b = bᵢ)

simulation = Simulation(model, Δt =Δtᵢ, stop_time = stop_time)

wizard = TimeStepWizard(cfl = 0.2, max_Δt = 60)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

progress(sim) = @info string(
    "Iteration: $(iteration(sim)), time: $(time(sim)), Δt: $(sim.Δt)")
add_callback!(simulation, progress, IterationInterval(100))

#ωz = ∂x(v) - ∂y(u)
#b′ = b - b̄

outputs = (; u, v, w, b) #ωz, b′)
timenow = Dates.format(now(), "yymmdd-HHMMSS")

output_filename = "output_$(timenow).nc"
output_filepath = joinpath("./Output", output_filename)
mkpath(dirname(output_filepath)) #Make directory if nonexistent

simulation.output_writers[:field_writer] = NetCDFOutputWriter(
                model, outputs, 
                filename = output_filepath, 
                schedule = TimeInterval(save_interval))

run!(simulation)
print(timenow, "\n")
