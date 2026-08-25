include("LibraryCoordinateTransforms.jl")
include("LibraryDynamics.jl")
include("LibraryEnergetics.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Visualization.jl")

using CUDA
using Dates: canonicalize, format, now
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Units
using Oceananigans.Utils 
using Printf, Random

using CairoMakie
using Oceananigans.Solvers

######################
# SPECIFY PARAMETERS #
######################

const Nx = 50 #x-grid size
const Ny = 50 #y-grid size
const Nz = 10 #z-grid size

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

#Ambient (i.e., excluding gyre's TWB contribution) N²-value at z -> -infty
const N²_far = 5e-5 * second^(-2)

gyreScaleParams = (f = f, U = U, σr = σr, σz = σz, N²_far = N²_far)

const z_grid = "uniform" #Either 'uniform' or 'chebyshev' 

#Type of ambient stratification to construct ('doubleTanh' or 'constant')
# (note that a TWB contribution will always be included)
const ambientStrat = "constant"

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

const Δt         = 3 #parse(Float64, ARGS[1]) #Simulation timestep (s)
const tf         = 3 #parse(Float64, ARGS[2]) #Simulation stop time (s)
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
const Δt_save = 3

const useGPU = false #Whether to use GPU
const useNHS = true  #Whether to use NonhydrostaticModel

const max_u′ = 1e-10 #Max. relative magnitude of initial velocity perturbation

#Whether to run visualization functions
const vis_const_x       = true
const vis_const_y       = false
const vis_const_z       = true
const vis_norms         = false
const vis_energetics    = false
const vis_z_grid        = false #Note: currently can only be done on CPU
const vis_bkgd_profiles = true
const vis_q_timeseries  = false

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
                      
#Note -- will need to fix zgrid (put this all into a function) and update for Chebyshev case
tallGrid = RectilinearGrid(architecture, 
                           topology = (Periodic, Periodic, Bounded), 
                           size = (Nx, Ny, Nz + 2), 
                           x = (-Lr, Lr), 
                           y = (-Lr, Lr), 
                           z = (-Lz - (Lz/Nz), Lz/Nz), 
                           halo = (Hx, Hy, Hz - 1)
                          )
                      
size(grid.yᵃᶜᵃ)[1] > 1 ? yFlat = false : yFlat = true

B_vals, Ux_vals, Uy_vals, Uz_vals, B_BCs = discrete_Cartesian_TWB_ICs(grid, tallGrid, gyreScaleParams, bkgd_Ψ_cylindrical_coords, ambientStrat; Hz = Hz)

#box_sponge = Relaxation(rate = 1, mask = PiecewiseLinearMask{:x}(center = 9 * σr, width = σr))

if useNHS
   model = NonhydrostaticModel(;
                               grid = grid, 
                               timestepper = :RungeKutta3,
                               advection = WENO(),
                               coriolis = fPlane,
                               pressure_solver = FourierTridiagonalPoissonSolver(grid),
                               tracers = (:b),
                               buoyancy = BuoyancyTracer(),
                               boundary_conditions = (; b = B_BCs)
                              )
elseif !useNHS
   model = HydrostaticFreeSurfaceModel(;
                                       grid = grid,
                                       momentum_advection = WENO(),
                                       tracer_advection = WENO(),
                                       coriolis = fPlane,
                                       tracers = (:b),
                                       buoyancy = BuoyancyTracer(),
                                       boundary_conditions = (; b = B_BCs)
                                      )
end

#Note: we must 'set!' in separate lines because we are setting with Fields, not
# functions.
set!(model.velocities.u, Ux_vals)
set!(model.velocities.v, Uy_vals)
set!(model.velocities.w, Uz_vals)
set!(model.tracers.b, B_vals)

pHY_initial = no_offset_view(adapt(Array, model.pressures.pHY′))[:, 25, :]
pNHS_initial = no_offset_view(adapt(Array, model.pressures.pNHS))[:, 25, :]

fig  = Figure(size = (1800, 600))
ax_HY = Axis(fig[1, 1], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = "Initial hydrostatic pressure anomaly")
ax_NHS = Axis(fig[1, 2], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = "Initial NHS pressure")
ax_b = Axis(fig[1, 3], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", title = "Initial buoyancy")

hm_HY = heatmap!(ax_HY, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), pHY_initial, colormap = :balance)
hm_NHS = heatmap!(ax_NHS, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), pNHS_initial, colormap = :balance)
hm_b = heatmap!(ax_b, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), no_offset_view(adapt(Array, model.tracers.b))[:, 25, :], colormap = :balance)

Colorbar(fig[2, 1], hm_HY, tickformat = "{:.1e}", label = "", 
            vertical = false, width = Relative(3/4))
Colorbar(fig[2, 2], hm_NHS, tickformat = "{:.1e}", label = "", 
            vertical = false, width = Relative(3/4))
Colorbar(fig[2, 3], hm_b, tickformat = "{:.1e}", label = "", vertical = false, width = Relative(3/4))
            
save(joinpath("./Plots", "ptest_initial.png"), fig)

@compute hdiv_initial = Field(∂x(model.velocities.u) + ∂y(model.velocities.v))
@compute vdiv_initial = Field(∂z(model.velocities.w))

fig_div = Figure(size = (1200, 600))
ax_hdiv = Axis(fig_div[1, 1], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = L"Initial $\nabla_h \cdot \vec{u}$ [1/s]")
ax_vdiv = Axis(fig_div[1, 2], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = L"Initial $\partial_z w$ [1/s]")
               
hm_hdiv = heatmap!(ax_hdiv, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), no_offset_view(adapt(Array, hdiv_initial))[:, 25, :], colormap = :balance)
hm_vdiv = heatmap!(ax_vdiv, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), no_offset_view(adapt(Array, vdiv_initial))[:, 25, :], colormap = :balance)

Colorbar(fig_div[2, 1], hm_hdiv, vertical = false)
Colorbar(fig_div[2, 2], hm_vdiv, vertical = false)

save(joinpath("./Plots", "div_initial.png"), fig_div)

#Print warnings if the respective instabilities are present
check_inertial_stability(model.grid, f, model.velocities.u, model.velocities.v)
check_gravitational_stability(model.tracers.b, model.grid)

######################################################
# SAVE BACKGROUND STATE AND DEFINE DIAGNOSTIC FIELDS #
######################################################

datetimestart = now()
datetimenow   = format(datetimestart, "yymmdd-HHMMSS")

print("Date-time label: $(datetimenow)", "\n")

Ur_vals, Uφ_vals = xy_vector_to_rφ(model.velocities.u, model.velocities.v, 
                                   model.grid, useGPU)

#Create fields to store background-state primitive variables and PV quantities
Ux        = XFaceField(model.grid)
Uy        = YFaceField(model.grid)
Ur        = CenterField(model.grid)
Uφ        = CenterField(model.grid)
Uz        = ZFaceField(model.grid)
B         = CenterField(model.grid;
                        boundary_conditions = discrete_Cartesian_TWB_ICs(grid,
                           tallGrid, gyreScaleParams, bkgd_Ψ_cylindrical_coords,
                           ambientStrat;
                           Hz = Hz, includeDefaultBCs = true)
                       )

Q_Ertel   = CenterField(model.grid)
Q_QG      = CenterField(model.grid)
∂rQ_Ertel = CenterField(model.grid)
∂rQ_QG    = CenterField(model.grid)

#Prescribe background values to those fields
set!(Ux, Ux_vals)
set!(Uy, Uy_vals)
set!(Ur, Ur_vals)
set!(Uφ, Uφ_vals)
set!(Uz, Uz_vals)
set!(B, B_vals)

fill_halo_regions!(B, model.clock, B)
#Note: this syntax is necessary because 'B' isn't a prognostic field of 'model'

#Create fields that are used in computing PKE budget terms
φ_ccc_vals = CenterField(model.grid)
∂rUφ       = CenterField(model.grid)
∂zUφ       = CenterField(model.grid)

#Prescribe initial values to those fields
set!(φ_ccc_vals, polar_coords_Fields(model.grid, "c", "c", "c")[2])
set!(∂rUφ, cos(φ_ccc_vals) * ∂x(Uφ) + sin(φ_ccc_vals) * ∂y(Uφ))
set!(∂zUφ, ∂z(Uφ))

### testing ###
Q_QG_fn = bkgd_Q_cylindrical_coords(gyreScaleParams, doubleTanhParams, 
                                   ambientStrat)
Q_QG_xyz(x, y, z) = Q_QG_fn(sqrt(x^2 + y^2), z)

#Compute background PVs (and their r-derivatives) and prescribe to fields
set!(Q_Ertel, 
     compute_Q_Ertel_Cartesian(model.grid, gyreScaleParams, Ux, Uy, Uz, B)
    )
set!(Q_QG, Q_QG_xyz) ##compute_Q_QG_Cartesian(model.grid, gyreScaleParams, Ux, Uy; Hz = Hz))
set!(∂rQ_Ertel, cos(φ_ccc_vals) * ∂x(Q_Ertel) + sin(φ_ccc_vals) * ∂y(Q_Ertel))
set!(∂rQ_QG, cos(φ_ccc_vals) * ∂x(Q_QG) + sin(φ_ccc_vals) * ∂y(Q_QG))

#############################
# SET UP AND RUN SIMULATION #
#############################
#=
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
=#
simulation = Simulation(model;
                        Δt = Δt,
                        stop_time = tf, 
                        align_time_step = false, 
                        minimum_relative_step = 1e-9
                       )

function progress(sim)

   u_absmax             = maximum(abs, sim.model.velocities.u)
   w_absmax             = maximum(abs, sim.model.velocities.w)
   b_absmax             = maximum(abs, sim.model.tracers.b)
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
   @info @sprintf("max|pNHS|: %.2e", pNHS_absmax)
   @info @sprintf("max|pNHS| in interior: %.2e", interior_pNHS_absmax)
   @info @sprintf("min|pNHS| in interior: %.2e", interior_pNHS_absmin)
   @info @sprintf("max|horizontal divergence of velocity| in interior: %.2e",
                  interior_hdiv_absmax)
   @info @sprintf("max|vertical divergence of velocity| in interior: %.2e", 
                  interior_vdiv_absmax)
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
                   
#if vis_q_timeseries
#   merge(outputs, (q_Ertel = CenterFields_q_Ertel(model.tracers.b, model.velocities.u, model.velocities.v, model.velocities.w, model.grid, f), ))
#end

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

field_writer = NetCDFWriter(model, 
                            outputs,
                            filename = outfilepath, 
                            schedule = TimeInterval(Δt_save),
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
   total_KE            = totalKE(simulation),
   total_KE_adv_flux   = totalKEadvFlux(simulation),
   total_pressure_work = totalPressureWork(simulation; useNHS = useNHS),
   total_KE_production = totalProduction(simulation; useNHS = useNHS),
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
                                  gyreParameters = gyreScaleParams),
   gyre_PAPE_to_PKE    = gyre_PAPE_to_PKE(simulation; 
                                          gyreParameters = gyreScaleParams),
   gyre_BTI_transfer   = gyre_BTI_transfer(simulation; 
                                           bkgdParameters = bkgdParameters, 
                                           gyreParameters = gyreScaleParams),
   gyre_BCI_transfer   = gyre_BCI_transfer(simulation; 
                                           bkgdParameters = bkgdParameters, 
                                           gyreParameters = gyreScaleParams)
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

pHY_final = no_offset_view(adapt(Array, interior(model.pressures.pHY′)))[:, 25, :]
pNHS_final = no_offset_view(adapt(Array, interior(model.pressures.pNHS)))[:, 25, :]

fig  = Figure(size = (1200, 600))
ax_HY = Axis(fig[1, 1], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = "Final hydrostatic pressure anomaly")
ax_NHS = Axis(fig[1, 2], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = "Final NHS pressure")

hm_HY = heatmap!(ax_HY, no_offset_view(model.grid.xᶜᵃᵃ)[4:end-3], no_offset_view(model.grid.z.cᵃᵃᶜ)[4:end-3], pHY_final, colormap = :balance)
hm_NHS = heatmap!(ax_NHS, no_offset_view(model.grid.xᶜᵃᵃ)[4:end-3], no_offset_view(model.grid.z.cᵃᵃᶜ)[4:end-3], pNHS_final, colormap = :balance)

Colorbar(fig[2, 1], hm_HY, tickformat = "{:.1e}", label = "", 
            vertical = false, width = Relative(3/4))
Colorbar(fig[2, 2], hm_NHS, tickformat = "{:.1e}", label = "", 
            vertical = false, width = Relative(3/4))
            
save(joinpath("./Plots", "ptest_final.png"), fig)

@compute hdiv_final = Field(∂x(model.velocities.u) + ∂y(model.velocities.v))
@compute vdiv_final = Field(∂z(model.velocities.w))

fig_div = Figure(size = (1200, 600))
ax_hdiv = Axis(fig_div[1, 1], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = L"Final $\nabla_h \cdot \vec{u}$ [1/s]")
ax_vdiv = Axis(fig_div[1, 2], xlabel = L"$x$ [m]", ylabel = L"$z$ [m]", 
               title = L"Final $\partial_z w$ [1/s]")
               
hm_hdiv = heatmap!(ax_hdiv, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), no_offset_view(adapt(Array, hdiv_final))[:, 25, :], colormap = :balance)
hm_vdiv = heatmap!(ax_vdiv, no_offset_view(model.grid.xᶜᵃᵃ), no_offset_view(model.grid.z.cᵃᵃᶜ), no_offset_view(adapt(Array, vdiv_final))[:, 25, :], colormap = :balance)

Colorbar(fig_div[2, 1], hm_hdiv, vertical = false)
Colorbar(fig_div[2, 2], hm_vdiv, vertical = false)

save(joinpath("./Plots", "div_final.png"), fig_div)

if vis_const_x
   visualize_fields_2D_slice(datetimenow, "x", x_idx, B, Uφ; 
                             t_idx_skip = t_idx_skip, 
                             plot_speed_animation = false, 
                             plot_animation = false)
end

if vis_const_y
   visualize_fields_2D_slice(datetimenow, "y", y_idx, B, Uφ; 
                             t_idx_skip = t_idx_skip) 
end

if vis_const_z
   visualize_fields_2D_slice(datetimenow, "z", z_idx, B, Uφ; 
                             t_idx_skip = t_idx_skip, 
                             plot_speed_animation = false, 
                             plot_animation = false)
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
                   growth_rate = "timeseries", t_skip_idx = t_idx_skip)
end

if vis_energetics
   visualize_total_QG_energy_budgets(datetimenow, model.grid)
   visualize_PKE(datetimenow, model.grid)
end

if vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end

if vis_bkgd_profiles
   visualize_B_U_Q_Ψ_vs_r_and_z(model.grid, gyreScaleParams, 
                                doubleTanhParams, ambientStrat, Nx ÷ 2, Nz, 
                                1e6, Lz) 
   visualize_B_and_N²_vs_z(B, model.grid, x_idx, y_idx, gyreScaleParams, 
                           doubleTanhParams; yFlat = yFlat, Hz = Hz)
   visualize_Q_and_∂Q∂r(Q_Ertel, Q_QG, ∂rQ_Ertel, ∂rQ_QG,
                        model.grid.xᶜᵃᵃ, model.grid.z.cᵃᵃᶜ, y_idx)
end

if vis_q_timeseries
   visualize_q_2D_slice(datetimenow, "z", z_idx, f; plot_animation = true)
end