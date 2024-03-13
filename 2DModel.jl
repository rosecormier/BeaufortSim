using Dates
using Oceananigans, Oceananigans.Coriolis, Oceananigans.TurbulenceClosures

#=
SIMULATION PARAMETERS
=#

Nx, Ny, Nz = 16, 16, 16 #Numbers of gridpoints
Lx, Ly, Lz = 2π * 100, 2π * 100, 500 #Domain extents

#=
SIMULATION SETUP
=#

grid = RectilinearGrid(size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz), topology=(Bounded, Bounded, Bounded))

#Diffusivities (constant anisotropic)

h_diffusivity, v_diffusivity = 1e-3, 1e-6
h_closure = HorizontalScalarDiffusivity(ν=h_diffusivity)
v_closure = VerticalScalarDiffusivity(ν=v_diffusivity)

#Effect of rotation

latitude = 73 #Degrees
BG_f = 2 * 7.2921 * 1e-5 * sin(π * latitude / 180) #Units s^-1
BG_f_plane = FPlane(f=BG_f)

#Instantiate the model
model = NonhydrostaticModel(; 
    grid=grid, 
    timestepper=:QuasiAdamsBashforth2, 
    closure=(h_closure, v_closure), 
    coriolis=BG_f_plane,
    tracers=(:b),
    buoyancy=BuoyancyTracer())

u, v, w = model.velocities

uᵢ = 10 * (1 .- rand(size(u)...))
vᵢ = 10 * (1 .- rand(size(v)...))
wᵢ = (1 .- rand(size(w)...))

set!(model, u=uᵢ, v=vᵢ, w=wᵢ)

simulation = Simulation(model, Δt=10, stop_time=500)
#Reduced timestep from what was in the Oceananigans tutorial

progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", time(sim))
add_callback!(simulation, progress, IterationInterval(10))

#Set up an output writer for the simulation that saves vorticity, speed, u, v regularly

u, v, w = model.velocities
ωx = ∂y(w) - ∂z(v)
ωy = ∂z(u) - ∂x(w)
ωz = ∂x(v) - ∂y(u)
s = sqrt(u^2 + v^2 + w^2)

output_fields = Dict("u" => u, "v" => v, "w" => w, 
                    "ωx" => ωx, "ωy" => ωy, "ωz" => ωz,
                    "s" => s)
output_filename = "output.nc"
rm("output.nc", force=true) #Remove file if already existing

simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, output_fields, filename=output_filename, schedule=IterationInterval(1))

#=
RUN SIMULATION
=#

run!(simulation)

#=
SAVE PARAMETERS TO LOG FILE
=#

timenow=Dates.format(now(), "yymmdd-HHMMSS")
log_filename = "log$(timenow).txt"
log_filepath = joinpath("./Logs", log_filename)
mkpath(dirname(log_filepath)) #Make log directory if it doesn't exist

open(log_filepath, "w") do file
    write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
    write(file, "Lx, Ly, Lz = $(Lx), $(Ly), $(Lz) \n")
    write(file, "Horizontal diffusivity = $(h_diffusivity) \n")
    write(file, "Vertical diffusivity = $(v_diffusivity) \n")
    write(file, "f-Plane latitude = $(latitude) degrees N \n")
end