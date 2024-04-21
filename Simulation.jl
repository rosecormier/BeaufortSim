using Dates, DelimitedFiles
using Oceananigans, Oceananigans.Coriolis, Oceananigans.TurbulenceClosures

#=
GET INPUTS
=#

#print(readdlm("InputSimulation.txt", Any))

#=
SIMULATION PARAMETERS AND GEOMETRY
=#

Nx, Ny, Nz = 128, 128, 2 #Numbers of gridpoints
Lx, Ly, Lz = 2000 * 1e3, 2000 * 1e3, 1000 #Domain extents

grid = RectilinearGrid(size=(Nx, Ny, Nz), 
    x = (-Lx/2, Lx/2),
    y = (-Ly/2, Ly/2),
    z = (-Lz, 0),
    topology=(Periodic, Periodic, Bounded))

#=
DIFFUSIVITIES
=#

#Constant anisotropic diffusivities
h_diffusivity, v_diffusivity = 5e-1, 5e-5 

h_closure = HorizontalScalarDiffusivity(ν=h_diffusivity)
v_closure = VerticalScalarDiffusivity(ν=v_diffusivity)

#=
EFFECT OF ROTATION
=#

latitude = 73 #Degrees N
f = 2 * 7.2921 * 1e-5 * sin(π * latitude / 180) #Units s^-1
f_plane = FPlane(f=f)

#=
INSTANTIATE THE MODEL
=#

model = NonhydrostaticModel(; 
    grid=grid, 
    timestepper=:QuasiAdamsBashforth2, 
    advection=UpwindBiasedFifthOrder(),
    closure=(h_closure, v_closure), 
    coriolis=f_plane,
    tracers=(:b),
    buoyancy=BuoyancyTracer())

#=
SET INITIAL CONDITIONS
=#

σr, σz = 2.5e5, 3e2 #Radial and vertical gyre lengthscales

N = 5e-3 #Approx. BV frequency (1/s); cf. Jackson et al, 2012
N2 = N^2

p0 = 100 #Reference value for reduced pressure (m^2/s^2)

#Function to compute initial buoyancy profile
function b_initial(x,y,z)
    mean_b = N2*z - (2*p0/σz^2) * z * exp(-(x^2+y^2)/σr^2 - (z/σz)^2)
    b_perturbation = 0
    buoyancy = mean_b + b_perturbation
end

#Functions to compute initial velocities
u_initial(x,y,z) = (2*p0/(f*σr^2)) * y * exp(-(x^2+y^2)/σr^2 - (z/σz)^2)
v_initial(x,y,z) = (-2*p0/(f*σr^2)) * x * exp(-(x^2+y^2)/σr^2 - (z/σz)^2)
w_initial(x,y,z) = 0

u, v, w = model.velocities
b = model.tracers.b

set!(model, u=u_initial, v=v_initial, w=w_initial, b=b_initial)

#=
DEFINE SIMULATION
=#

timestep = 10
simulation = Simulation(model, Δt=timestep, stop_time=5*60)

#=
SET UP CALLBACK AND OUTPUT WRITER FOR SIMULATION
=#

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
mkpath(dirname(output_filepath)) #Make output directory if nonexistent

simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, 
    output_fields, 
    filename=output_filepath, 
    schedule=IterationInterval(1))

#=
RUN SIMULATION
=#

run!(simulation)
print(timenow, "\n")

#=
SAVE PARAMETERS TO LOG FILE
=#

log_filename = "log_$(timenow).txt"
log_filepath = joinpath("./Logs", log_filename)
mkpath(dirname(log_filepath)) #Make log directory if nonexistent

open(log_filepath, "w") do file
    write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
    write(file, "Lx, Ly, Lz = $(Lx), $(Ly), $(Lz) \n")
    write(file, "Horizontal diffusivity = $(h_diffusivity) \n")
    write(file, "Vertical diffusivity = $(v_diffusivity) \n")
    write(file, "f-Plane latitude = $(latitude) degrees N \n")
    write(file, "p0, σr, σz = $(p0), $(σr), $(σz) \n")
    write(file, "N^2 = $(N2) \n")
    write(file, "Timestep = $(timestep)")
end