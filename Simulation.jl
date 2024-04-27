include("Library.jl")

using Dates
using Oceananigans, Oceananigans.Coriolis, Oceananigans.TurbulenceClosures

# GET INPUTS

input_params = ReadInputFile("InputSimulation.txt")

# SIMULATION PARAMETERS AND GEOMETRY

Nx, Ny, Nz = input_params["Nx"], input_params["Ny"], input_params["Nz"]

#Convert horizontal extents to m
Lx, Ly = input_params["Lx"] * 1e3, input_params["Ly"] * 1e3
Lz = input_params["Lz"] #Units m

grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                        x = (-Lx/2, Lx/2), 
                        y = (-Ly/2, Ly/2), 
                        z = (-Lz, 0),
                        topology = (Periodic, Periodic, Bounded),
                        halo = (3,3,3))

# INSTANTIATE THE MODEL WITH ROTATION AND CONSTANT ANISOTROPIC DIFFUSIVITIES

model = NonhydrostaticModel(; grid = grid, 
    timestepper = :QuasiAdamsBashforth2, 
    advection = UpwindBiasedFifthOrder(),
    closure = (
        HorizontalScalarDiffusivity(
            ν = input_params["h_diffusivity"]), 
        VerticalScalarDiffusivity(
            ν = input_params["v_diffusivity"])), 
    coriolis = FPlane(latitude = input_params["latitude"]),
    tracers = (:b),
    buoyancy = BuoyancyTracer())

# SET INITIAL CONDITIONS

#Physical parameters
σr = input_params["σr"] * 1e3 #Radial lengthscale - convert to m
σz = input_params["σz"] #Vertical lengthscale [m]
N2 = input_params["N"]^2
p̃0 = input_params["p̃0"] #Reference value for reduced pressure 

f = 2 * 7.2921 * 1e-5 * sin(π * input_params["latitude"] / 180)

#Function to compute initial buoyancy profile
b_initial(x,y,z) = (N2*z 
    - (2*p̃0/σz^2) * z * exp(-(x^2+y^2)/σr^2 - (z/σz)^2))

#Functions to compute initial velocities
u_initial(x,y,z) = (2*p̃0/(f*σr^2)) * y * exp(-(x^2+y^2)/σr^2 - (z/σz)^2)
v_initial(x,y,z) = (-2*p̃0/(f*σr^2)) * x * exp(-(x^2+y^2)/σr^2 - (z/σz)^2)
w_initial(x,y,z) = 0

u, v, w = model.velocities

set!(model, 
    u = u_initial, 
    v = v_initial, 
    w = w_initial, 
    b = b_initial)

# DEFINE SIMULATION; SET UP CALLBACK AND OUTPUT WRITER

simulation = Simulation(model, 
                Δt = input_params["timestep"], 
                stop_time = input_params["stop_time"])

progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", time(sim))
add_callback!(simulation, progress, IterationInterval(1))

b = model.tracers.b
u, v, w = model.velocities
ωx = ∂y(w) - ∂z(v)
ωy = ∂z(u) - ∂x(w)
ωz = ∂x(v) - ∂y(u)
s = sqrt(u^2 + v^2 + w^2)

output_fields = Dict("u" => u, "v" => v, "w" => w, 
                    "ωx" => ωx, "ωy" => ωy, "ωz" => ωz,
                    "s" => s, "b" => b)

timenow = Dates.format(now(), "yymmdd-HHMMSS")

output_filename = "output_$(timenow).nc"
output_filepath = joinpath("./Output", output_filename)
mkpath(dirname(output_filepath)) #Make directory if nonexistent

simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, 
                output_fields, 
                filename=output_filepath, 
                schedule=IterationInterval(1))

# RUN SIMULATION

run!(simulation)
print(timenow, "\n")

# SAVE PARAMETERS TO LOG FILE

log_filename = "log_$(timenow).txt"
log_filepath = joinpath("./Logs", log_filename)
mkpath(dirname(log_filepath)) #Make directory if nonexistent

open(log_filepath, "w") do file
    write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
    write(file, "Lx, Ly, Lz = $(Lx), $(Ly), $(Lz) \n")
    write(file, "Horizontal diffusivity = 
        $(input_params["h_diffusivity"]) \n")
    write(file, "Vertical diffusivity = 
        $(input_params["v_diffusivity"]) \n")
    write(file, "f-Plane latitude = 
        $(input_params["latitude"]) degrees N \n")
    write(file, "p̃0, σr, σz = $(p̃0), $(σr), $(σz) \n")
    write(file, "N^2 = $(N2) \n")
    write(file, "Timestep, stop time = 
        $(input_params["timestep"]), 
        $(input_params["stop_time"])")
end