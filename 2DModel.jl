using Oceananigans, Distributions, Random

grid = RectilinearGrid(size=(64, 64), extent=(2π, 2π), topology=(Periodic, Periodic, Flat))

#Instantiate the model
model = NonhydrostaticModel(; 
                            grid, 
                            timestepper=:QuasiAdamsBashforth2, 
                            advection=UpwindBiasedFifthOrder(), 
                            closure=ScalarDiffusivity(ν=2e-3)
)

u, v, w = model.velocities

#Set the target RMS speed
u_RMS_target = 2.23488

#Set up Gaussian distribution
dist = Normal(1, u_RMS_target^2)

uᵢ= rand(dist, size(u))
vᵢ= rand(dist, size(v))

set!(model, u=uᵢ, v=vᵢ)

simulation = Simulation(model, Δt=0.005, stop_time=5)
#Reduced timestep from what was in the Oceananigans tutorial

progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", time(sim))
add_callback!(simulation, progress, IterationInterval(1))

#Set up an output writer for the simulation that saves vorticity, speed, u, v regularly

u, v, w = model.velocities
ω = ∂x(v) - ∂y(u)
s = sqrt(u^2 +v^2)

output_fields = Dict("u" => u, "v" => v, "ω" => ω, "s" => s)

simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, output_fields, filename="two_dimensional_turbulence.nc", schedule=TimeInterval(0.1))

run!(simulation)
