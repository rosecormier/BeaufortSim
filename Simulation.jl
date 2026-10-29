include("LibraryCoordinates.jl")
include("LibraryDynamics.jl")
include("LibraryEnergetics.jl")
include("LibraryOptionalSimOutputs.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Parameters.jl")
include("Visualization.jl")
include("WriteLogFile.jl")

using CUDA
using Dates: canonicalize, format, now
using Oceananigans
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

simulation = Simulation(model;
                        Δt = Δt,
                        stop_time = tf, 
                        align_time_step = simParams.align_time_step, 
                        minimum_relative_step = simParams.minimum_relative_step
                       )

function progress(sim)

   u_absmax             = maximum(abs, sim.model.velocities.u)
   w_absmax             = maximum(abs, sim.model.velocities.w)
   b_absmax             = maximum(abs, sim.model.tracers.b)
   pHY′_absmax          = maximum(abs, sim.model.pressures.pHY′)
   pNHS_absmax          = maximum(abs, sim.model.pressures.pNHS)
   interior_pNHS_absmax = maximum(abs, interior(sim.model.pressures.pNHS))
   interior_pNHS_absmin = minimum(abs, interior(sim.model.pressures.pNHS))
   interior_hdiv_absmax = maximum(abs, interior(∂x(sim.model.velocities.u) 
                                                + ∂y(sim.model.velocities.v)
                                               )
                                 )
   interior_vdiv_absmax = maximum(abs, interior(∂z(sim.model.velocities.w)))
   
   @info @sprintf("Iter: %d; time: %.2e days; Δt: %s",
		  iteration(sim), (time(sim)/day),  prettytime(sim.Δt))
   @info @sprintf("max|u|: %.2e; max|w|: %.2e; max|b|: %.2e", 
                  u_absmax, w_absmax, b_absmax)
   @info @sprintf("max|pHY′|: %.2e", pHY′_absmax)
   @info @sprintf("max|pNHS|: %.2e", pNHS_absmax)
   @info @sprintf("max|pNHS| in interior: %.2e", interior_pNHS_absmax)
   @info @sprintf("min|pNHS| in interior: %.2e", interior_pNHS_absmin)
   @info @sprintf("max|horizontal divergence of velocity| in interior: %.2e",
                  interior_hdiv_absmax)
   @info @sprintf("max|vertical divergence of velocity| in interior: %.2e", 
                  interior_vdiv_absmax)
   
   return nothing
end

add_callback!(simulation, progress, TimeInterval(Δt_save))

ur = KernelFunctionOperation{Center, Center, Center}(vector_r_component_ccc, model.grid, model.velocities.u, model.velocities.v)
uφ = KernelFunctionOperation{Center, Center, Center}(vector_φ_component_ccc, model.grid, model.velocities.u, model.velocities.v)

outputs = (ur = ur,
           uφ = uφ,
	         ux = model.velocities.u,
	         uy = model.velocities.v,
	         uz = model.velocities.w,
	         b = model.tracers.b,
           ζ = CenterField_ζ′(model.grid, model.velocities.u, model.velocities.v, Ux, Uy)
          )

#Define output filepaths
outfilepath    = joinpath("./Output", "output_$(datetimenow).nc")
scalarfilepath = joinpath("./Output", "scalars_$(datetimenow).nc")
energyfilepath = joinpath("./Output", "energetics_$(datetimenow).nc")
logfilepath    = joinpath("./Logs", "log_$(datetimenow).txt")

#Make paths if nonexistent (all other paths will be generated ad hoc)
mkpath(dirname(logfilepath))
mkpath("./Checkpoints")

field_writer = NetCDFWriter(model, 
                            outputs,
                            filename = outfilepath, 
                            schedule = TimeInterval(Δt_save),
                            array_type = Array{Float64},
                            with_halos = true,
                            file_splitting = FileSizeLimit(30GiB)
                           )
                           
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
		                  b′_norm  = b_perturbation_norm
                     )

scalar_writer = NetCDFWriter(model, 
                             scalar_diagnostics,
		                         filename = scalarfilepath, 
				                     schedule = TimeInterval(Δt_save),
                             array_type = Array{Float64},
                             file_splitting = FileSizeLimit(30GiB),
		                         dimensions = (ux′_norm = (),
					                                 uy′_norm = (),
					                                 ur′_norm = (),
					                                 uφ′_norm = (),
					                                 uz′_norm = (),
					                                 b′_norm  = ()
                                          )
                            )

bkgdParameters = (Ur = Ur, Uφ = Uφ, Uz = Uz, ∂rUφ = ∂rUφ, ∂zUφ = ∂zUφ)

energy_diagnostics = (; 
   total_KE            = total_KE(simulation),
   total_KE_adv_flux   = total_KE_advective_flux(simulation),
   total_pressure_work = total_pressure_work(simulation, modelParams),
   total_KE_production = total_KE_production(simulation, modelParams),
   total_PE            = total_PE(simulation),
   total_int_PE_source = total_internal_PE_sources(simulation),
   total_gravity_work  = total_gravity_work(simulation, GravitationalAcceleration),
   total_PKE           = PKE(simulation, Ux, Uy, Uz),
   total_PAPE_to_PKE   = PAPE_to_PKE(simulation, B, Uz),
   total_BTI_transfer  = BTI_transfer(simulation; 
                                      bkgdParameters = bkgdParameters),
   total_BCI_transfer  = BCI_transfer(simulation; 
                                      bkgdParameters = bkgdParameters),
   gyre_PKE            = gyre_PKE(simulation, gyreParams),
   gyre_PAPE_to_PKE    = gyre_PAPE_to_PKE(simulation, gyreParams),
   gyre_BTI_transfer   = gyre_BTI_transfer(simulation, gyreParams; 
                                           bkgdParameters = bkgdParameters),
   gyre_BCI_transfer   = gyre_BCI_transfer(simulation, gyreParams; 
                                           bkgdParameters = bkgdParameters),
   total_bkgd_transfer = bkgd_to_PKE_transfer(simulation, Ux, Uy, Uz),
   total_bkgd_xxx = bkgd_transfer_ux′ux′∂xUx(simulation, Ux),
   total_bkgd_xyx = bkgd_transfer_ux′uy′∂yUx(simulation, Ux, Uy),
   total_bkgd_xzx = bkgd_transfer_ux′uz′∂zUx(simulation, Ux, Uz),
   total_bkgd_yxy = bkgd_transfer_uy′ux′∂xUy(simulation, Ux, Uy),
   total_bkgd_yyy = bkgd_transfer_uy′uy′∂yUy(simulation, Uy),
   total_bkgd_yzy = bkgd_transfer_uy′uz′∂zUy(simulation, Uy, Uz),
   total_bkgd_zxz = bkgd_transfer_uz′ux′∂xUz(simulation, Ux, Uz),
   total_bkgd_zyz = bkgd_transfer_uz′uy′∂yUz(simulation, Uy, Uz),
   total_bkgd_zzz = bkgd_transfer_uz′uz′∂zUz(simulation, Uz),
   total_TKE_prod = TKEProduction(simulation, Ux, Uy, Uz)
                     )

energy_writer = NetCDFWriter(model, energy_diagnostics,
                             filename = energyfilepath,
                             schedule = TimeInterval(Δt_save),
                             array_type = Array{Float64},
                             file_splitting = FileSizeLimit(30GiB)
                            )

checkpointer = Checkpointer(model; 
                            schedule = TimeInterval(simParams.Δt_checkpt),
                            dir = "Checkpoints", 
			    	                prefix = "checkpoint_$(datetimenow)"
                           )

simulation.output_writers[:field_writer]  = field_writer
simulation.output_writers[:scalar_writer] = scalar_writer
simulation.output_writers[:energy_writer] = energy_writer
simulation.output_writers[:checkpointer]  = checkpointer

run!(simulation; pickup = false)
#     pickup = joinpath("./Checkpoints", "checkpoint_261026-145018_iteration36000.jld2"))

duration = canonicalize(now() - datetimestart)

#Append zeros to filenames so they can be accessed in chronological order
pad_filenames(datetimenow)
pad_filenames(datetimenow; prefix = "energetics")
pad_filenames(datetimenow; prefix = "scalars")

#Save parameters to logfile
WriteLogFile(simulation, duration, logfilepath, outfilepath, gridParams, 
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