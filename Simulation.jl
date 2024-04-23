using Dates
using Oceananigans, Oceananigans.Coriolis, Oceananigans.TurbulenceClosures

#=
GET INPUTS
=#

input_params = Dict()

for (i, line) in enumerate(eachline("InputSimulation.txt"))
    
    var_len, line_len = 0, length(line)
    if (line_len > 1 && line[1] != '#')
 
        for char in line
            if char == '='
                break
            else
                var_len += 1
            end
        end
        
        #Variable name (whitespace stripped)
        var = strip(line[1:var_len])
        #Strip leading/trailing whitespace from rest of line
        line = strip(line[var_len:end])
        
        val_len = 0
        for char in line
            if char == '#' 
                break #Isolate parameter value from trailing comment
            else
                val_len += 1
            end
        end
        
        #Strip leading/trailing whitespace from parameter value
        val_str = strip(line[2:val_len]) #Index 1 is "="
        
        val = tryparse(Int, val_str)
        if isnothing(val)
            val = parse(Float64, val_str)
        end
        input_params[var] = val
        
    end    
end

#=
SIMULATION PARAMETERS AND GEOMETRY
=#

#Numbers of gridpoints
Nx, Ny, Nz = input_params["Nx"], input_params["Ny"], input_params["Nz"]

#Numbers of gridpoints
Lx, Ly, Lz = input_params["Lx"], input_params["Ly"], input_params["Lz"]

grid = RectilinearGrid(size=(Nx, Ny, Nz), 
    x = (-Lx/2, Lx/2), y = (-Ly/2, Ly/2), z = (-Lz, 0),
    topology=(Periodic, Periodic, Bounded))

#=
DIFFUSIVITIES
=#

#Constant anisotropic diffusivities
h_diffusivity = input_params["h_diffusivity"]
v_diffusivity = input_params["v_diffusivity"]

h_closure = HorizontalScalarDiffusivity(ν=h_diffusivity)
v_closure = VerticalScalarDiffusivity(ν=v_diffusivity)

#=
EFFECT OF ROTATION
=#

latitude = input_params["latitude"]
f = 2 * 7.2921 * 1e-5 * sin(π * latitude / 180) #Units s^-1
f_plane = FPlane(f=f)

#=
INSTANTIATE THE MODEL
=#

model = NonhydrostaticModel(; grid=grid, 
    timestepper=:QuasiAdamsBashforth2, 
    advection=UpwindBiasedFifthOrder(),
    closure=(h_closure, v_closure), 
    coriolis=f_plane,
    tracers=(:b),
    buoyancy=BuoyancyTracer())

#=
SET INITIAL CONDITIONS
=#

#Radial and vertical gyre lengthscales
σr, σz = input_params["σr"], input_params["σz"]

#BV frequency
N = input_params["N"]
N2 = N^2

#Reference value for reduced pressure 
p0 = input_params["p0"]

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

timestep, stop_time = input_params["timestep"], input_params["stop_time"]
simulation = Simulation(model, Δt=timestep, stop_time=stop_time)

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
    write(file, "Timestep, stop time = $(timestep), $(stop_time)")
end