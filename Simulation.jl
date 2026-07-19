include("LibraryDynamics.jl")
include("LibraryEnergetics.jl")
include("LibraryOptionalSimOutputs.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Visualization.jl")

using Adapt, CUDA
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

const Nx = 400 #x-grid size
const Ny = 400 #y-grid size
const Nz = 100 #z-grid size

const Hx = 3 #Number of x halo cells per boundary
const Hy = 3 #Number of y halo cells per boundary
const Hz = 3 #Number of z halo cells per boundary

const Lr = 2.5e3 * kilometer #[Minimum] domain radius
const Lz = 1 * kilometer     #z-axis length

const lat = 74.0     #Latitude (deg. N)
fPlane    = FPlane(latitude = lat)
const f   = fPlane.f #Coriolis frequency

const U  = 5e-2 * (meter/second) #Maximum gyre velocity scale (at surface)
const σr = 250 * kilometer       #Radial gyre length scale
const σz = 300 * meter 	         #Vertical gyre length scale

#Ambient (i.e., excluding gyre's TWB contribution) N²-value at z --> -infty
const N²_far = 2e-5 * second^(-2)

const z_grid = "chebyshev" #Either 'uniform' or 'chebyshev'

#Type of ambient stratification to construct ('doubleTanh' or 'constant')
# (note that a TWB contribution will always be included)
const ambientStrat = "doubleTanh"

#Parameters for double-tanh stratification (defined as in Kosty et al., 2026)
const g   = -9.81 * meter * (second^2)
const ρ₀  = 1025.5 * meter^(-3)
const A_s = 2.5 * meter^(-3)
const z_s = -40 * meter
const C_s = 15 * meter
const A_d = 1.05 * meter^(-3)
const z_d = -200 * meter
const C_d = 60 * meter

doubleTanhParams = (g = g, ρ₀ = ρ₀, A_s = A_s, C_s = C_s, z_s = z_s,
                    A_d = A_d, C_d = C_d, z_d = z_d)

const Δt         = 600 #parse(Float64, ARGS[1]) #Simulation timestep (s)
const tf         = 600 #parse(Float64, ARGS[2]) #Simulation stop time (s)
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
const Δt_save = 600

const useGPU = false  #Whether to use GPU
const useNHS = true #Whether to use NonhydrostaticModel

const max_u′ = 1e-10 #Max. relative magnitude of initial velocity perturbation

#Whether to run visualization functions
const vis_const_x       = false
const vis_const_y       = false
const vis_const_z       = false
const vis_norms         = true
const vis_energetics    = false
const vis_z_grid        = false #Note: currently can only be done on CPU
const vis_bkgd_profiles = false

const x_idx      = Nx ÷ 2 #Visualize yz-slice at this x-index
const y_idx      = Ny ÷ 2 #Visualize xz-slice at this y-index
const z_idx      = Nz - 1 #Visualize xy-slice at this z-index
const t_idx_skip = 10      #Step size for animations and timeseries

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

custom_z_grids = Dict("uniform"   => (-Lz, 0),
                      "chebyshev" => k -> chebyshev_spaced_faces(k, -Lz, Nz + 1)
                     )

grid = RectilinearGrid(architecture,
		                   topology = (Periodic, Periodic, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lr, Lr),
		                   y = (-Lr, Lr),
                       z = custom_z_grids[z_grid],
		                   halo = (Hx, Hy, Hz)
                      )

b̄_BCs = buoyancy_BCS(f, σr, σz, U, N²_far, grid, false;
                      doubleTanhParams = doubleTanhParams)

#box_sponge = Relaxation(rate = 1, mask = PiecewiseLinearMask{:x}(center = 9 * σr, width = σr))

if useNHS
   model = NonhydrostaticModel(;
                               grid = grid, 
                               timestepper = :RungeKutta3,
                               advection = WENO(),
                               coriolis = fPlane,
                               tracers = (:b),
                               buoyancy = BuoyancyTracer(),
                               boundary_conditions = (; b = b̄_BCs)
                              )
elseif !useNHS
   model = HydrostaticFreeSurfaceModel(;
                                       grid = grid,
                                       momentum_advection = WENO(),
                                       tracer_advection = WENO(),
                                       coriolis = fPlane,
                                       tracers = (:b),
                                       buoyancy = BuoyancyTracer(),
                                       boundary_conditions = (; b = b̄_BCs)
                                      )
end

b̄     = bkgd_buoyancy(f, σr, σz, U, N²_far;
                       grid = model.grid,
                       doubleTanhParams = doubleTanhParams)
ū, v̄ = bkgd_velocities(σr, σz, U)

set!(model, u = ū, v = v̄, b = b̄)

#Print warnings if the respective instabilities are present
check_inert_stability(model.grid, f, model.velocities.u, model.velocities.v; 
                      z_idx = z_idx)
check_grav_stability(model.tracers.b, model.grid)

######################################################
# SAVE BACKGROUND STATE AND DEFINE DIAGNOSTIC FIELDS #
######################################################

datetimestart = now()
datetimenow   = "260706-163110" #format(datetimestart, "yymmdd-HHMMSS")

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

@inline u_perturbed(x, y, z) = @inbounds (ū(x, y, z)
                                          * (1 + 2 * (rand() - 0.5) * max_u′
                                                 / (U * sqrt(2))
                                            )
                                         )

if !isnothing(seed2)
   Random.seed!(seed2) #Update seed so next random number is independent
end

@inline v_perturbed(x, y, z) = @inbounds (v̄(x, y, z)
                                          * (1 + 2 * (rand() - 0.5) * max_u′ 
                                                 / (U * sqrt(2))
                                            )
                                         )

set!(model, u = u_perturbed, v = v_perturbed) #Set perturbed ICs
#=
simulation = Simulation(model;
                        Δt = Δt,
                        stop_time = tf, 
                        align_time_step = false, 
                        minimum_relative_step = 1e-9
                       )

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
	         b = model.tracers.b)
=#
#Define output filepaths
outfilepath    = joinpath("./Output", "output_$(datetimenow).nc")
scalarfilepath = joinpath("./Output", "scalars_$(datetimenow).nc")
energyfilepath = joinpath("./Output", "energetics_$(datetimenow).nc")
logfilepath    = joinpath("./Logs", "log_$(datetimenow).txt")
#=
#Make required paths if nonexistent
mkpath(dirname(outfilepath))
mkpath(dirname(scalarfilepath))
mkpath(dirname(energyfilepath))
mkpath(dirname(logfilepath))
mkpath("./Checkpoints")

field_writer = NetCDFWriter(model, outputs,
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
				                     schedule = TimeInterval(1 * hour),
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
   total_KE            = totalKE(simulation),
   total_KE_adv_flux   = totalKEadvFlux(simulation),
   total_KE_production = totalProduction(simulation; useNHS = useNHS),
   total_pressure_work = totalPressureWork(simulation; useNHS = useNHS),
   total_PE            = totalPE(simulation, g),
   total_b_adv_flux    = totalBuoyancyAdvFlux(simulation),
   total_gravity_work  = totalGravityWork(simulation, g),
   total_PKE           = PKE(simulation, Ux, Uy, Uz),
   total_PAPE_to_PKE   = PAPE_to_PKE(simulation, B, Uz),
   total_BTI_transfer  = BTI_transfer(simulation; 
                                      bkgdParameters = bkgdParameters),
   total_BCI_transfer  = BCI_transfer(simulation; 
                                      bkgdParameters = bkgdParameters),
   gyre_PKE            = gyre_PKE(simulation; 
                                  gyreParameters = gyreParameters),
   gyre_PAPE_to_PKE    = gyre_PAPE_to_PKE(simulation; 
                                          gyreParameters = gyreParameters),
   gyre_BTI_transfer   = gyre_BTI_transfer(simulation; 
                                           bkgdParameters = bkgdParameters, 
                                           gyreParameters = gyreParameters),
   gyre_BCI_transfer   = gyre_BCI_transfer(simulation; 
                                           bkgdParameters = bkgdParameters, 
                                           gyreParameters = gyreParameters)
                     )

energy_writer = NetCDFWriter(model, energy_diagnostics,
                             filename       = energyfilepath,
                             schedule       = TimeInterval(Δt_save),
                             file_splitting = FileSizeLimit(30GiB)
                            )

checkpointer = Checkpointer(model; 
                            schedule   = TimeInterval(Δt_checkpt),
                            dir        = "Checkpoints", 
			    	                prefix     = "checkpoint_$(datetimenow)", 
		    	                  properties = [:grid, :clock, :timestepper,
					                                :velocities, :tracers]
                           )

simulation.output_writers[:field_writer]  = field_writer
simulation.output_writers[:scalar_writer] = scalar_writer
simulation.output_writers[:energy_writer] = energy_writer
simulation.output_writers[:checkpointer]  = checkpointer

run!(simulation; pickup = false)
    #pickup = joinpath("./Checkpoints", "checkpoint_260605-075732_iteration6.jld2"))

duration = canonicalize(now() - datetimestart)

#Append zeros to filenames so they can be accessed in chronological order
pad_filenames(datetimenow)
pad_filenames(datetimenow; prefix = "energetics")
pad_filenames(datetimenow; prefix = "scalars")
=#
#Save parameters to logfile
open(logfilepath, "w") do file
   write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
   write(file, "Lr, Lz = $(Lr), $(Lz) \n\n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "U = $(U) \n")
   write(file, "Far-field N² term = $(N²_far) \n")
   write(file, "Stratification type = $(ambientStrat) \n")
   write(file, "Stratification parameters: g = $(g), ρ₀ = $(ρ₀), A_s = $(A_s),
                C_s = $(C_s), z_s = $(z_s), A_d = $(A_d), C_d = $(C_d), 
                z_d = $(z_d) \n")
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
   visualize_fields_2D_slice(datetimenow, "x", x_idx, B, Uφ; 
                             t_idx_skip = t_idx_skip, plot_speed_animation = true)
end

if vis_const_y
   visualize_fields_2D_slice(datetimenow, "y", y_idx, B, Uφ; 
                             t_idx_skip = t_idx_skip) 
end

if vis_const_z
   visualize_fields_2D_slice(datetimenow, "z", z_idx, B, Uφ; 
                             t_idx_skip = t_idx_skip, plot_speed_animation = true)
end

if vis_norms
   visualize_norms(datetimenow; 
                   idxStartLinGrowth_b = 24, idxEndLinGrowth_b = 37,
                   idxStartLinGrowth_ur = 100, idxEndLinGrowth_ur = 103,
                   idxStartLinGrowth_uφ = 100, idxEndLinGrowth_uφ = 103,
                   idxStartLinGrowth_ux = 1, idxEndLinGrowth_ux = 5,
                   idxStartLinGrowth_uy = 1, idxEndLinGrowth_uy = 5,
                   idxStartLinGrowth_uz = 24, idxEndLinGrowth_uz = 37,
                   idxStartPlot = 540, idxEndPlot = -1, 
                   growth_rate = "timeseries")
end

if vis_energetics
   visualize_total_energy_budgets(datetimenow, model.grid)
   visualize_PKE(datetimenow, model.grid)
end

if vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end

if vis_bkgd_profiles
   visualize_B_U_Q_Ψ_vs_r_and_z(U, model.grid, f, σr, σz, N²_far, 
                                doubleTanhParams, ambientStrat, Nx ÷ 2, Nz, 
                                1e6, Lz)
   visualize_B_and_N²_vs_z(B, model.grid, x_idx, y_idx, doubleTanhParams, f, 
                           σr, σz, U, N²_far; Hz = Hz)
   
   ωx_initial = CenterField(model.grid)
   ωy_initial = CenterField(model.grid)
   ωz_initial = CenterField(model.grid)
   Q_Ertel    = CenterField(model.grid)
   Q_QG       = CenterField(model.grid)
   ∂rQ_Ertel  = CenterField(model.grid)
   ∂rQ_QG     = CenterField(model.grid)
   
   rcoords = CenterField(model.grid)
   set!(rcoords, compute_polar_coords(model.grid)[1])
   
   Ψ = CenterField(model.grid)
   Ψ_function = bkgd_Ψ_cylindrical_coords(σr, σz, U)
   
   set!(Ψ, Ψ_function(rcoords, reshape(no_offset_view(adapt(Array, model.grid.z.cᵃᵃᶜ))[4:length(model.grid.z.cᵃᵃᶜ)-3], 1, 1, Nz)))
   
   set!(ωx_initial, ∂y(Uz) - ∂z(Uy))
   set!(ωy_initial, ∂z(Ux) - ∂x(Uz))
   set!(ωz_initial, ∂x(Uy) - ∂y(Ux))
   
   set!(Q_Ertel, (ωx_initial * ∂x(B) + ωy_initial * ∂y(B) + (f + ωz_initial) * ∂z(B))/f)
   set!(Q_QG, -(1/rcoords) * (cos(φcoords) * ∂x(rcoords * Uφ) + sin(φcoords) * ∂y(rcoords * Uφ)) + f^2 * ∂z(∂z(Ψ) / ∂z(B)) + f)
   #set!(Q_QG, ∂x(Uy) - ∂y(Ux) + (f^2 * ∂z(∂z(Ψ) / ∂z(B))) + f)
   set!(∂rQ_Ertel, cos(φcoords) * ∂x(Q_Ertel) + sin(φcoords) * ∂y(Q_Ertel))
   set!(∂rQ_QG, cos(φcoords) * ∂x(Q_QG) + sin(φcoords) * ∂y(Q_QG))
   
   visualize_Q_and_∂Q∂r_from_ICs(datetimenow, Q_Ertel, Q_QG, ∂rQ_Ertel, ∂rQ_QG, model.grid.xᶜᵃᵃ, model.grid.z.cᵃᵃᶜ, y_idx)
end