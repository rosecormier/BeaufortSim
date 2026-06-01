include("LibraryDynamics.jl")
include("LibraryEnergetics.jl")
include("LibraryOptionalSimOutputs.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Visualization.jl")

using Adapt, CSV, CUDA
using Dates: canonicalize, format, now
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Units
using Oceananigans.Utils 
using Printf, Random

######################
# SPECIFY PARAMETERS #
######################

const Nx = 252 #x-grid size
const Ny = 252 #y-grid size
const Nz = 20  #z-grid size

const Hx = 3 #Number of x halo cells per boundary
const Hy = 3 #Number of y halo cells per boundary
const Hz = 3 #Number of z halo cells per boundary

const Lr = 2.5e3 * kilometer #[Minimum] domain radius
const Lz = 1 * kilometer     #z-axis length

const lat = 74.0     #Latitude (deg. N)
fPlane    = FPlane(latitude = lat)
const f   = fPlane.f #Coriolis frequency

const U  = 3.5 * (meter/second) #Maximum gyre velocity scale (at surface)
const σr = 250 * kilometer      #Radial gyre length scale
const σz = 300 * meter 	        #Vertical gyre length scale

#Ambient (i.e., excluding gyre's thermal-wind contribution) N²-value at z = 0
const constantN²Term = 4e-3 * second^(-2)

#Type of ambient stratification to construct ('fromData', 'doubleTanh', or 'constant')
ambientStrat = "doubleTanh"

#Season to use data (Timmermans & Toole, 2023) from when constructing stratification.
# Note this stratification can only be constructed on CPU; fails on GPU.
const season = "Ma06"

#Parameters for double-tanh stratification (defined as in Kosty et al., 2026)
const g  = 9.81 * meter * (second^2)
const ρ₀ = 1000 * meter^(-3)
const As = 1.7e-1 * meter^(-3)
const Bs = 40 * meter
const Cs = 75 * meter
const Ad = 2e-2 * meter^(-3)
const Bd = 200 * meter
const Cd = 150 * meter

const z_grid = "uniform" #Either 'uniform' or 'chebyshev'
const d_ML   = -30 * meter #Mixed-layer depth (<= 0); only necessary for Chebyshev grid

const Δt         = 600 * second #parse(Float64, ARGS[1]) #Simulation timestep (s)
const tf         = 2400 * second #parse(Float64, ARGS[2]) #Simulation stop time (s)
const Δt_checkpt = 250 * day   		         #Checkpoint interval
#=
#Set save interval
if parse(Float64, ARGS[3]) < tf / 250
   print("Save interval too small for given duration. Using tf/250 instead.")
   const Δt_save = tf / 250
else
   const Δt_save = parse(Float64, ARGS[3])
end
=#
const Δt_save = 600 * second

const useGPU = true #Whether to use GPU
const useNHS = true #Whether to use NonhydrostaticModel

const max_u′ = 1e-10 #Max. relative magnitude of initial velocity perturbation

#Whether to run visualization functions
const vis_const_x    = true
const vis_const_y    = false
const vis_const_z    = false
const vis_norms      = false
const vis_energetics = false
const vis_z_grid     = false #Note: currently can only be done on CPU

const x_idx      = Nx ÷ 2 #Visualize yz-slice at this x-index
const y_idx      = Ny ÷ 2 #Visualize xz-slice at this y-index
const z_idx      = Nz - 1 #Visualize xy-slice at this z-index
const t_idx_skip = 1      #Step size for animations and timeseries

#Seeds for 2 random-number generators
const seed1 = 12345
const seed2 = 67890

if !isnothing(seed1)
   Random.seed!(seed1)
end

#########################
# SET UP GRID AND MODEL #
#########################

useGPU ? architecture = GPU() : architecture = CPU()

z_grids = Dict("uniform"   => (-Lz, 0),
               "chebyshev" => k -> chebyshev_spaced_faces(k, -Lz, Nz; ξ0 = d_ML))

grid = RectilinearGrid(architecture,
		                   topology = (Bounded, Bounded, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lr, Lr),
		                   y = (-Lr, Lr),
                       z = z_grids[z_grid],
		                   halo = (Hx, Hy, Hz))

save_zC_values(z_grid, d_ML, grid) #If Chebyshev z-grid, save values to csv file

if ambientStrat == "doubleTanh"

   doubleTanhParams = (g = g, ρ₀ = ρ₀,
                       As = As, Bs = Bs, Cs = Cs, 
                       Ad = Ad, Bd = Bd, Cd = Cd)

   #Construct ambient stratification using double-tanh function
   N²DoubleTanhFunction = N²DoubleTanh(doubleTanhParams)
   discreteN²           = @. N²DoubleTanhFunction(grid.z.cᵃᵃᶜ)
   
   const additionalN²Top    = @view discreteN²[end-1]
   const additionalN²Bottom = @view discreteN²[1]
   
elseif ambientStrat == "fromData"
   #Construct ambient stratification from seasonal data
   const discreteN²         = N²FromData(Nz, d_ML, constantN²Term, season)
   const additionalN²Top    = @view N²_from_data(Nz, d_ML, constantN²Term, season)[end-1]
   const additionalN²Bottom = @view N²_from_data(Nz, d_ML, constantN²Term, season)[1]
elseif ambientStrat == "constant"
   const discreteN²         = nothing
   const additionalN²Top    = nothing
   const additionalN²Bottom = nothing
end

b̄_BCs = buoyancy_BCS(σz, constantN²Term, 0, -Lz, false;
                      parameters = (additionalN²Top = additionalN²Top, 
                                    additionalN²Bottom = additionalN²Bottom)
                     )

if useNHS
   model = NonhydrostaticModel(; 
                               grid = grid, 
                               timestepper = :RungeKutta3,
                               advection = WENO(),
                               coriolis = fPlane,
                               tracers = (:b),
                               buoyancy = BuoyancyTracer(),
                               boundary_conditions = (; b = b̄_BCs))
elseif !useNHS
   model = HydrostaticFreeSurfaceModel(;
                                       grid = grid,
                                       momentum_advection = WENO(),
                                       tracer_advection = WENO(),
                                       coriolis = fPlane,
                                       tracers = (:b),
                                       buoyancy = BuoyancyTracer(),
                                       boundary_conditions = (; b = b̄_BCs))
end

b̄     = bkgd_buoyancy(f, σr, σz, U;
                       constantN²Term = constantN²Term,
                       discreteN² = discreteN², 
                       grid = model.grid)
ū, v̄ = bkgd_velocities(σr, σz, U)

set!(model, u = ū, v = v̄, b = b̄)

#Print warnings if the respective instabilities are present
check_inert_stability(model.grid, f, model.velocities.u, model.velocities.v; 
                      z_idx = z_idx)
check_grav_stability(model.tracers.b)

#########################################################
# SAVE BACKGROUND STATE AND DEFINE DIAGNOSTIC FUNCTIONS #
#########################################################

datetimestart = now()
datetimenow   = format(datetimestart, "yymmdd-HHMMSS")

print("Date-time label: $(datetimenow)", "\n")

Ur_vals, Uφ_vals = xy_vector_to_rφ(model.velocities.u, model.velocities.v, 
                                   model.grid, useGPU)

#Create fields to store background state
Ux = XFaceField(model.grid)
Uy = YFaceField(model.grid)
Ur = CenterField(model.grid)
Uφ = CenterField(model.grid)
Uz = ZFaceField(model.grid)
B  = CenterField(model.grid)

#Prescribe background values to those fields
set!(Ux, ū)
set!(Uy, v̄)
set!(Ur, Ur_vals)
set!(Uφ, Uφ_vals)
set!(Uz, model.velocities.w)
set!(B, b̄ )
fill_halo_regions!(Ux)
fill_halo_regions!(Uy)
fill_halo_regions!(Uz)
fill_halo_regions!(B)

#Create fields that are used in computing PKE budget terms
φcoords = CenterField(model.grid)
∂rUφ    = CenterField(model.grid)
∂zUφ    = CenterField(model.grid)

#Prescribe initial values to those fields
set!(φcoords, compute_polar_coords(grid)[2])
set!(∂rUφ, cos(φcoords) * ∂x(Uφ) + sin(φcoords) * ∂y(Uφ))
set!(∂zUφ, ∂z(Uφ))

#############################
# SET UP AND RUN SIMULATION #
#############################

#Add random perturbations to horizontal velocity components

@inline u_perturbed(x, y, z) = @inbounds ū(x, y, z)*(1 + 2*(rand()-0.5) * max_u′/(U*sqrt(2)))

if !isnothing(seed2)
   Random.seed!(seed2) #Update seed so next random number is independent
end

@inline v_perturbed(x, y, z) = @inbounds v̄(x, y, z)*(1 + 2*(rand()-0.5) * max_u′/(U*sqrt(2)))

set!(model, u = u_perturbed, v = v_perturbed) #Set perturbed ICs

initial_KE = CenterField(model.grid)
set!(initial_KE, (model.velocities.u^2 + model.velocities.v^2
                  + model.velocities.w^2) / 2)
total_initial_KE = Field(Integral(initial_KE))
compute!(total_initial_KE)

simulation = Simulation(model;
                        Δt = Δt,
                        stop_time = tf, 
                        align_time_step = false, 
                        minimum_relative_step = 1e-9)

function progress(sim)
   umax = maximum(abs, sim.model.velocities.u)
   wmax = maximum(abs, sim.model.velocities.w)
   bmax = maximum(abs, sim.model.tracers.b)
   @info @sprintf("Iter: %d; time: %.2e days; Δt: %s",
		  iteration(sim), (time(sim)/day),  prettytime(sim.Δt))
   @info @sprintf("max|u|: %.2e; max|w|: %.2e; max|b|: %.2e", umax, wmax, bmax)
   @info @sprintf("Norm of u' = %.10e", norm(sim.model.velocities.u - Ux))
   return nothing
end

add_callback!(simulation, progress, TimeInterval(Δt_save))

ur, uφ = xy_vector_to_rφ(model.velocities.u, model.velocities.v, model.grid, useGPU)

outputs = (ur = ur,
	         uφ = uφ,
	         ux = model.velocities.u,
	         uy = model.velocities.v,
	         uz = model.velocities.w,
	         b  = model.tracers.b)

#Define output filepaths
outfilepath    = joinpath("./Output", "output_$(datetimenow).nc")
scalarfilepath = joinpath("./Output", "scalars_$(datetimenow).nc")
energyfilepath = joinpath("./Output", "energetics_$(datetimenow).nc")
logfilepath    = joinpath("./Logs", "log_$(datetimenow).txt")

#Make required paths if nonexistent
mkpath(dirname(outfilepath))
mkpath(dirname(scalarfilepath))
mkpath(dirname(energyfilepath))
mkpath(dirname(logfilepath))
mkpath("./Checkpoints")

field_writer = NetCDFWriter(model, outputs, with_halos = true,
                            filename = outfilepath, 
                            schedule = TimeInterval(Δt_save),
                            file_splitting = FileSizeLimit(30GiB))

ux_perturbation_norm(model) = perturbation_norm(model.velocities.u, Ux)
uy_perturbation_norm(model) = perturbation_norm(model.velocities.v, Uy)
ur_perturbation_norm(model) = perturbation_norm(ur, Ur)
uφ_perturbation_norm(model) = perturbation_norm(uφ, Uφ)
uz_perturbation_norm(model) = perturbation_norm(model.velocities.w, Uz)
b_perturbation_norm(model)  = perturbation_norm(model.tracers.b, B)

scalar_diagnostics = (ux′_norm = ux_perturbation_norm,
		                  uy′_norm = uy_perturbation_norm,
		                  ur′_norm = ur_perturbation_norm,
		                  uφ′_norm = uφ_perturbation_norm,
		                  uz′_norm = uz_perturbation_norm,
		                  b′_norm  = b_perturbation_norm)

scalar_writer = NetCDFWriter(model, scalar_diagnostics,
		                         filename = scalarfilepath, 
				                     schedule = TimeInterval(1*hour),
                             file_splitting = FileSizeLimit(30GiB),
		                         dimensions = (ux′_norm = (),
					                                 uy′_norm = (),
					                                 ur′_norm = (),
					                                 uφ′_norm = (),
					                                 uz′_norm = (),
					                                 b′_norm  = ()))

bkgdParameters = (Ur = Ur, Uφ = Uφ, Uz = Uz, ∂rUφ = ∂rUφ, ∂zUφ = ∂zUφ)
gyreParameters = (σr = σr, σz = σz)

energy_diagnostics = (; 
   integrated_pKE = total_PKE(simulation; Ux, Uy, Uz),
   integrated_pAPE_to_pKE = total_PAPE_to_PKE(simulation; B, Uz),
   integrated_BTI_transfer = BTI_transfer(simulation; 
                                          bkgdParameters = bkgdParameters),
   integrated_BCI_transfer = BCI_transfer(simulation; 
                                          bkgdParameters = bkgdParameters),
   gyre_integrated_pKE = gyre_PKE(simulation; 
                                  gyreParameters = gyreParameters),
   gyre_integrated_pAPE_to_pKE = gyre_PAPE_to_PKE(simulation; 
                                                  gyreParameters = gyreParameters),
   gyre_BTI_transfer = gyre_BTI_transfer(simulation; 
                                         bkgdParameters = bkgdParameters, 
                                         gyreParameters = gyreParameters),
   gyre_BCI_transfer = gyre_BCI_transfer(simulation; 
                                         bkgdParameters = bkgdParameters, 
                                         gyreParameters = gyreParameters)
                     )

energy_writer = NetCDFWriter(model, energy_diagnostics,
				   filename = energyfilepath,
			 	   schedule = TimeInterval(Δt_save),
			     file_splitting = FileSizeLimit(30GiB))

checkpointer = Checkpointer(model; 
                            schedule = TimeInterval(Δt_checkpt),
                            dir = "Checkpoints", 
			    	                prefix = "checkpoint_$(datetimenow)", 
		    	                  properties = [:grid, :clock, :timestepper,
					                                :velocities, :tracers])

simulation.output_writers[:field_writer]  = field_writer
simulation.output_writers[:scalar_writer] = scalar_writer
simulation.output_writers[:energy_writer] = energy_writer
simulation.output_writers[:checkpointer]  = checkpointer

run!(simulation; pickup = false) 
     #pickup = joinpath("../Checkpoints", "checkpoint_260310-084141_iteration720000.jld2"))

duration = canonicalize(now() - datetimestart)

#Append zeros to filenames so they can be accessed in chronological order
pad_filenames(datetimenow)
pad_filenames(datetimenow; prefix = "energetics")
pad_filenames(datetimenow; prefix = "scalars")

#Save parameters to logfile
open(logfilepath, "w") do file
   write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
   write(file, "Lr, Lz = $(Lr), $(Lz) \n\n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "U = $(U) \n")
   write(file, "Constant N² term = $(constantN²Term) \n")
   write(file, "Season = $(season) \n")
   write(file, "Max. u' = $(max_u′) \n")
   write(file, "Random-number seeds = $(seed1), $(seed2) \n\n")
   write(file, "Δt, tf = $(Δt), $(tf) \n\n")
   write(file, "Total number of iterations = $(iteration(simulation)) \n")
   write(file, "Simulation runtime = $(duration) \n")
   write(file, "Output filesize = $(pretty_filesize(filesize(outfilepath)))")
end

#####################
# RUN VISUALIZATION #
#####################

if vis_const_x
   visualize_fields_2D_slice(datetimenow, "x", x_idx, B, Uφ, Hx, Hy, Hz; 
			     t_idx_skip = t_idx_skip)
end

if vis_const_y
   visualize_fields_2D_slice(datetimenow, "y", y_idx, B, Uφ, Hx, Hy, Hz; 
			     t_idx_skip = t_idx_skip) 
end

if vis_const_z
   visualize_fields_2D_slice(datetimenow, "z", z_idx, B, Uφ, Hx, Hy, Hz; 
			     t_idx_skip = t_idx_skip)
end

if vis_norms
   visualize_norms(datetimenow, idxStartLinGrowth_b = 24, idxEndLinGrowth_b = 37,
                idxStartLinGrowth_ur = 100, idxEndLinGrowth_ur = 103,
                idxStartLinGrowth_uφ = 100, idxEndLinGrowth_uφ = 103,
                idxStartLinGrowth_ux = 1, idxEndLinGrowth_ux = 5,
                idxStartLinGrowth_uy = 1, idxEndLinGrowth_uy = 5,
                idxStartLinGrowth_uz = 24, idxEndLinGrowth_uz = 37,
                idxStartPlot = 2, idxEndPlot = -1)
end

if vis_energetics
   visualize_energetics(datetimenow, model.grid, total_initial_KE.data)
end

if vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end