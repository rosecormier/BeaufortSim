using Dates
using Oceananigans, Oceananigans.Coriolis, Oceananigans.TurbulenceClosures

#=
SIMULATION PARAMETERS AND GEOMETRY
=#

Nx, Ny, Nz = 32, 32, 32 #Numbers of gridpoints
Lx, Ly, Lz = 2000 * 1e3, 2000 * 1e3, 1400 #Domain extents

grid = RectilinearGrid(size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz), topology=(Periodic, Periodic, Bounded))

#=
DIFFUSIVITIES
=#

#Constant anisotropic diffusivities
h_diffusivity, v_diffusivity = 5e-1, 5e-5 #5e-3, 5e-7
h_closure = HorizontalScalarDiffusivity(ν=h_diffusivity)
v_closure = VerticalScalarDiffusivity(ν=v_diffusivity)

#=
EFFECT OF ROTATION
=#

latitude = 73 #Degrees
BG_f = 2 * 7.2921 * 1e-5 * sin(π * latitude / 180) #Units s^-1
BG_f_plane = FPlane(f=BG_f)

#=
INSTANTIATE THE MODEL
=#

model = NonhydrostaticModel(; 
    grid=grid, 
    timestepper=:QuasiAdamsBashforth2, 
    closure=(h_closure, v_closure), 
    coriolis=BG_f_plane,
    tracers=(:b),
    buoyancy=BuoyancyTracer())

#=
SET INITIAL CONDITIONS
=#

σr, σz = 3e5, 4e2 #Standard deviations of pressure in r and z
p_surf_max = 1.3e5 #Maximum surface pressure
p0 = p_surf_max / exp(1) # 2π * σr * σz * p_surf_max / exp(1)
x0, y0 = Lx / 2, Ly / 2 #Centres of horizontal axes

#Function to compute (mean) pressure
p_bar(x,y,z) = p0 * exp(-((x-x0)/σr)^2) * exp(-((y-y0)/σr)^2) * exp(-(z/σz)^2)
    
grav = 9.80665 #Gravitational acceleration; units m/s^2
rho_sw = 1020 #Surface seawater density

#Function to compute mean density from mean pressure
rho_bar(x,y,z) = rho_sw - (2 * p_bar(x,y,z) * z / (grav * (σz^2)))

N2 = 1e-3 #Squared BV frequency

#Function to compute initial buoyancy profile with perturbation
function initial_b(x,y,z)
    mean_rho = rho_bar(x,y,z)
    mean_b = mean_rho .* N2
    #add perturbation
    buoyancy = mean_b
end

#Function to compute initial zonal velocity with perturbation
function initial_u(x,y,z)
    pbar, rhobar = p_bar(x,y,z), rho_bar(x,y,z)
    mean_u = 2 * pbar * (y-y0) / (rhobar * BG_f * (σr^2))
    #add perturbation
    u = mean_u
end

#Function to compute initial zonal velocity with perturbation
function initial_v(x,y,z)
    pbar, rhobar = p_bar(x,y,z), rho_bar(x,y,z)
    mean_v = -2 * pbar * (x-x0) / (rhobar * BG_f * (σr^2))
    #add perturbation
    v = mean_v
end

#Function to compute initial vertical velocity with perturbation
function initial_w(x,y,z)
    mean_w = 0
    #add perturbation
    w = mean_w
end

u, v, w = model.velocities
b = model.tracers.b

set!(model, u=initial_u, v=initial_v, w=initial_w, b=initial_b)

#=
DEFINE SIMULATION
=#

timestep = 10
simulation = Simulation(model, Δt=timestep, stop_time=5*60*60)

#=
SET UP CALLBACK AND OUTPUT WRITER FOR SIMULATION
=#

progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", time(sim))
add_callback!(simulation, progress, IterationInterval(10))

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

simulation.output_writers[:field_writer] = NetCDFOutputWriter(model, 
    output_fields, 
    filename=output_filename, 
    schedule=IterationInterval(1))

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
    write(file, "Max. surf. p, σr, σz = $(p_surf_max), $(σr), $(σz) \n")
    write(file, "x0, y0 = $(x0), $(y0) \n")
    write(file, "g = $(grav) \n")
    write(file, "N^2 = $(N2) \n")
    write(file, "Timestep = $(timestep)")
end