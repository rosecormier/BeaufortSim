include("LibraryCoordinates.jl")
include("LibraryDynamics.jl")
include("LibraryEnergetics.jl")
include("LibraryOptionalSimOutputs.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Parameters.jl")
include("SimulationSetup.jl")
include("Visualization.jl")
include("WriteLogFile.jl")

using CUDA
using Dates: canonicalize, format, now
#using Oceananigans
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Solvers
using Oceananigans.Units
using Oceananigans.Utils 
using Printf

######################
# SPECIFY PARAMETERS #
######################

gridParams            = GridParams()
modelParams           = ModelParams()
gyreParams            = cudaconvert(GyreParams())
doubleTanhStratParams = cudaconvert(DoubleTanhStratParams())
perturbParams         = cudaconvert(PerturbationParams())
simParams             = SimulationParams()
visParams             = VisParams()

const Δt = parse(Float64, ARGS[1]) #Simulation timestep (s)
const tf = parse(Float64, ARGS[2]) #Simulation stop time (s)

#Set save interval
if parse(Float64, ARGS[3]) < tf / 250
   print("Save interval too small for given duration. Using tf/250 instead.")
   const Δt_save = tf / 250
else
   const Δt_save = parse(Float64, ARGS[3])
end

#Plot 2D slices at these indices
const x_idx = gridParams.Hx + (gridParams.Nx ÷ 2)
const y_idx = gridParams.Hy + (gridParams.Ny ÷ 2)
const z_idx = gridParams.Hz + gridParams.Nz

#Step size for animations and timeseries
const t_idx_skip = 1

#########################
# SET UP GRID AND MODEL #
#########################

grid, backgroundFieldInteriors, B_BCs, initialFields = GridAndDiscretelyBalancedBackgroundFields(gridParams, gyreParams, perturbParams, visParams; x_idx = x_idx, y_idx = y_idx)

if modelParams.useNonhydrostaticModel
   model = NonhydrostaticModel(grid;
                        timestepper = modelParams.timestepper,
                        advection = WENO(),
                        coriolis = FPlane(latitude = gyreParams.lat),
                        pressure_solver = FourierTridiagonalPoissonSolver(grid),
                        hydrostatic_pressure_anomaly = CenterField(grid),
                        tracers = (:b),
                        buoyancy = BuoyancyTracer(),
                        boundary_conditions = (; b = B_BCs)
                              )
else
   model = HydrostaticFreeSurfaceModel(;
                                       grid = grid,
                                       momentum_advection = WENO(),
                                       tracer_advection = WENO(),
                                       coriolis = FPlane(latitude = gyreParams.lat),
                                       tracers = (:b),
                                       buoyancy = BuoyancyTracer(),
                                       boundary_conditions = (; b = B_BCs)
                                      )
end

BInterior, UxInterior, UyInterior, UzInterior = backgroundFieldInteriors

#Fill prognostic fields.
#Note: Must 'set!' in separate lines because we are setting with Fields, not
# functions.
set!(model.velocities.u, UxInterior)
set!(model.velocities.v, UyInterior)
set!(model.velocities.w, UzInterior)
set!(model.tracers.b, BInterior)

fill_halo_regions!(model.tracers.b)
fill_halo_regions!(model.velocities.u)
fill_halo_regions!(model.velocities.v)
fill_halo_regions!(model.velocities.w)

#Print warnings if the respective instabilities are present
check_inertial_stability(model, gyreParams.f)
check_gravitational_stability(model)

######################################################
# SAVE BACKGROUND STATE AND DEFINE DIAGNOSTIC FIELDS #
######################################################

datetimestart = now()
datetimenow   = format(datetimestart, "yymmdd-HHMMSS")

print("Date-time label: $(datetimenow)", "\n")

#Azimuthal coordinates of model.grid-points
φ = Field(KernelFunctionOperation{Center, Center, Center}(φ_cca, model.grid))

#Create fields to store background-state quantities
Ux = XFaceField(model.grid; boundary_conditions = model.velocities.u.boundary_conditions)
Uy = YFaceField(model.grid; boundary_conditions = model.velocities.v.boundary_conditions)
Uz = ZFaceField(model.grid; boundary_conditions = model.velocities.w.boundary_conditions)
B  = CenterField(model.grid; boundary_conditions = model.tracers.b.boundary_conditions)

#Note: these have technically incorrect halo values (0) but should be
# approximately accurate for a wide enough domain that accommodates the whole
# gyre.
Ur = Field(KernelFunctionOperation{Center, Center, Center}(vector_r_component_ccc, model.grid, UxInterior, UyInterior))
Uφ = Field(KernelFunctionOperation{Center, Center, Center}(vector_φ_component_ccc, model.grid, UxInterior, UyInterior))

#Prescribe background values to those fields
set!(Ux, UxInterior)
set!(Uy, UyInterior)
set!(Uz, UzInterior)
set!(B, BInterior)

fill_halo_regions!(Ux)
fill_halo_regions!(Uy)
fill_halo_regions!(Uz)
fill_halo_regions!(B)

ζ_bkgd = CenterField(model.grid)
set!(ζ_bkgd, ∂x(Uy) - ∂y(Ux))

#Create fields that are used in computing PKE budget terms
∂rUφ = CenterField(model.grid)
∂zUφ = CenterField(model.grid)

#Prescribe initial values to those fields
set!(∂rUφ, cos(φ) * ∂x(Uφ) + sin(φ) * ∂y(Uφ))
set!(∂zUφ, ∂z(Uφ))

#Store everything in one tuple to pass to diagnostic functions
backgroundFields = (Ux = Ux, Uy = Uy, Ur = Ur, Uφ = Uφ, Uz = Uz, B = B,
                    ∂rUφ = ∂rUφ, ∂zUφ = ∂zUφ)

#############################
# SET UP AND RUN SIMULATION #
#############################

b_initial, u_initial, v_initial, w_initial = initialFields

#To do: check that this is equiv to setting in separate lines
set!(model, b = b_initial, u = u_initial, v = v_initial, w = w_initial) #Set perturbed ICs

fill_halo_regions!(model.tracers.b)
fill_halo_regions!(model.velocities.u)
fill_halo_regions!(model.velocities.v)
fill_halo_regions!(model.velocities.w)
fill_halo_regions!(model.pressures)

simulation = CreateSimulation(model, Δt, tf, Δt_save, simParams, modelParams, gyreParams, backgroundFields, datetimenow)

#Make path if nonexistent (all other paths will be generated ad hoc)
mkpath("./Checkpoints")

simulation.output_writers[:checkpointer] = Checkpointer(model; 
                                                        schedule = TimeInterval(simParams.Δt_checkpt),
                                                        dir = "Checkpoints", 
			    	                                            prefix = "checkpoint_$(datetimenow)"
                                                       )

run!(simulation; pickup = false)
#     pickup = joinpath("./Checkpoints", "checkpoint_261026-145018_iteration36000.jld2"))

duration = canonicalize(now() - datetimestart)

#Append zeros to filenames so they can be accessed in chronological order
pad_filenames(datetimenow)
pad_filenames(datetimenow; prefix = "energetics")
pad_filenames(datetimenow; prefix = "scalars")

#Save parameters to logfile - need to fix path args
WriteLogFile(simulation, duration, datetimenow, gridParams, 
             gyreParams, doubleTanhStratParams, perturbParams, simParams)

#####################
# RUN VISUALIZATION #
#####################

#Fouriertesting(datetimenow, model, gridParams)

if visParams.vis_const_x
   visualize_fields_2D_slice(datetimenow, "x", x_idx, B, Ur, Uφ, Uz, gridParams.Hx, gridParams.Hy, gridParams.Hz;
                             t_idx_skip = t_idx_skip, 
                             plot_speed_animation = false, 
                             plot_animation = true)
   visualize_b_and_ζ(datetimenow, "x", x_idx, model.grid, B, ζ_bkgd,
                     gridParams.Hx, gridParams.Hy, gridParams.Hz; 
                     plot_animation = false)
end

if visParams.vis_const_y
   visualize_fields_2D_slice(datetimenow, "y", y_idx, B, Ur, Uφ, Uz, gridParams.Hx, gridParams.Hy, gridParams.Hz; 
                             t_idx_skip = t_idx_skip) 
end

if visParams.vis_const_z
   visualize_fields_2D_slice(datetimenow, "z", z_idx, B, Ur, Uφ, Uz, gridParams.Hx, gridParams.Hy, gridParams.Hz; 
                             t_idx_skip = t_idx_skip, 
                             plot_speed_animation = false, 
                             plot_animation = true)
   visualize_b_and_ζ(datetimenow, "z", z_idx, model.grid, B, ζ_bkgd,
                     gridParams.Hx, gridParams.Hy, gridParams.Hz;
                     plot_animation = false)
end

if visParams.vis_norms
   visualize_norms(datetimenow; 
                   idxStartLinGrowth_b = 2, idxEndLinGrowth_b = -1,
                   idxStartLinGrowth_ur = 2, idxEndLinGrowth_ur = 5,
                   idxStartLinGrowth_uφ = 2, idxEndLinGrowth_uφ = 5,
                   idxStartLinGrowth_ux = 2, idxEndLinGrowth_ux = 5,
                   idxStartLinGrowth_uy = 2, idxEndLinGrowth_uy = 5,
                   idxStartLinGrowth_uz = 2, idxEndLinGrowth_uz = 5,
                   idxStartPlot = 2, idxEndPlot = -1, 
                   growth_rate = "timeseries", t_skip_idx = t_idx_skip)
end

if visParams.vis_energetics
   visualize_total_QG_energy_budgets(datetimenow, model.grid)
   visualize_PKE(datetimenow, model.grid; visualize_individual_bkgd_terms = true)
end

if visParams.vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -gridParams.Lz)
end