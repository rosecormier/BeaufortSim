include("LibraryDynamics.jl")
include("LibraryOptionalSimOutputs.jl")
include("LibraryStability.jl")
include("LibraryVisualization.jl")
include("Visualization.jl")

using Adapt, CSV, CUDA
using Dates: canonicalize, format, now
using LinearAlgebra: norm
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Units
using Oceananigans.Utils 
using Printf, Random, Tables

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

const U   = 3.5 * (meter/second) #Gyre velocity scale (at surface)
const N²₀ = 4e-3 * (second^(-2)) #1e-4 * (second^(-2)) #Buoyancy frequency squared (at surface)

const σr = 250 * kilometer #Radial gyre length scale
const σz = 300 * meter 	   #Vertical gyre length scale

#Max buoyancy frequency (equal to N²₀ for uniform stratification)
const N²_max = 4e-3 * (second^(-2)) #1e-4 * (second^(-2))

const d_ML = -30 * meter #Mixed-layer depth

const Δt         = parse(Float64, ARGS[1]) #Simulation timestep (s)
const tf         = parse(Float64, ARGS[2]) #Simulation stop time (s)
const Δt_checkpt = 250 * day   		         #Checkpoint interval

#Set save interval
if parse(Float64, ARGS[3]) < tf / 200
   print("Save interval too small for given duration. Using tf/200 instead.")
   const Δt_save = tf / 200
else
   const Δt_save = parse(Float64, ARGS[3])
end

const useGPU = false #Whether to use GPU

const max_u′ = 1e-10 #Max. relative magnitude of initial velocity perturbation

#Whether to run visualization functions
const vis_const_x    = false
const vis_const_y    = false
const vis_const_z    = false
const vis_norms      = false
const vis_energetics = false
const vis_z_grid     = true #Note: currently can only be done on CPU

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

z_grid_spacing(k) = chebyshev_spaced_faces(k, -Lz, Nz; ξ0 = d_ML)

grid = RectilinearGrid(architecture,
		                   topology = (Bounded, Bounded, Bounded),
                       size = (Nx, Ny, Nz), 
                       x = (-Lr, Lr),
		                   y = (-Lr, Lr),
                       z = z_grid_spacing,
		                   halo = (Hx, Hy, Hz))
                       #z = (-Lz, 0.0),

gridfilepath = joinpath("./Logs", "grid_Nz$(Nz)_MLD$(abs(d_ML)).csv") 

mkpath(dirname(gridfilepath)) #Make required path if nonexistent
CSV.write(gridfilepath, 
          Tables.table(znodes(grid, Center())), header = false) #Save zC values

const bkgd_N²_top = N²₀ #lognormal_strat(N²₀, N²_max, d_ML, 0)[1]
const bkgd_N²_bot = N²₀ #lognormal_strat(N²₀, N²_max, d_ML, -Lz)[1]

b̄, ū, v̄, b̄_BCs = bkgd_fields_3D(f, σr, σz, U, bkgd_N²_top, bkgd_N²_bot, 0, -Lz)

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            coriolis = fPlane,
                            tracers = (:b),
               		    buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; b = b̄_BCs))

set!(model, u = ū, v = v̄, b = b̄)

#Print warnings if the respective instabilities are present
check_inert_stability(model.grid, f, model.velocities.u, model.velocities.v;
		      z_idx = z_idx)
check_grav_stability(model.tracers.b; grid = model.grid, x_idx = x_idx)

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
set!(Ux, model.velocities.u)
set!(Uy, model.velocities.v)
set!(Ur, Ur_vals)
set!(Uφ, Uφ_vals)
set!(Uz, model.velocities.w)
set!(B, model.tracers.b)
fill_halo_regions!(Ux)
fill_halo_regions!(Uy)
fill_halo_regions!(Uz)
fill_halo_regions!(B)

∂z_Uφ = CenterField(model.grid)
set!(∂z_Uφ, ∂z(Uφ))

φcoords = CenterField(model.grid)
set!(φcoords, compute_polar_coords(grid)[2])

∂r_Uφ = CenterField(model.grid)
set!(∂r_Uφ, cos(φcoords) * ∂x(Uφ) + sin(φcoords) * ∂y(Uφ))

@inline ur′(i, j, k, grid, ur, Ur) = @inbounds ur[i, j, k] - Ur[i, j, k]

@inline uφ′(i, j, k, grid, uφ, Uφ) = @inbounds uφ[i, j, k] - Uφ[i, j, k]

@inline uz′(i, j, k, grid, uz, Uz) = @inbounds uz[i, j, k] - Uz[i, j, k]

@inline perturbation_norm(field, bkgd_field) = norm(field - bkgd_field)

@inline ψ′²(i, j, k, grid, ψ, ψ̄) = @inbounds (ψ[i, j, k] - ψ̄[i, j, k])^2

@inline pKE_ccc(i, j, k, grid, u, v, w, Ux, Uy, Uz) = @inbounds (
     		      		ℑxᶜᵃᵃ(i, j, k, grid, ψ′², u, Ux) + 
     		      		ℑyᵃᶜᵃ(i, j, k, grid, ψ′², v, Uy) +
                      		ℑzᵃᵃᶜ(i, j, k, grid, ψ′², w, Uz)) / 2

pointwise_pKE_op = KernelFunctionOperation{Center, Center, Center}(pKE_ccc,
   grid, model.velocities.u, model.velocities.v, model.velocities.w, Ux, Uy, Uz)

function compute_integrated_pKE(sim)
   compute!(Integral(Field(pointwise_pKE_op)))
end

function ∂rUφ_ur′_uφ′_ccc(i, j, k, grid, ux, uy, Ur, Uφ, ∂r_Uφ)

   φ = atan(ℑyᵃᶜᵃ(i, j, k, grid, ynodes(grid, Center())), ℑxᶜᵃᵃ(i, j, k, grid, xnodes(grid, Center())))

   ux_ccc = @inbounds ℑxᶜᵃᵃ(i, j, k, grid, ux)
   uy_ccc = @inbounds ℑyᵃᶜᵃ(i, j, k, grid, uy)
   ur_ccc = @inbounds (ux_ccc * cos(φ)) + (uy_ccc * sin(φ))
   uφ_ccc = @inbounds (uy_ccc * cos(φ)) - (ux_ccc * sin(φ))

   ur′_ccc = ur′(i, j, k, grid, ur_ccc, Ur)
   uφ′_ccc = uφ′(i, j, k, grid, uφ_ccc, Uφ)
   ∂r_Uφ_ccc = @inbounds ∂r_Uφ[i, j, k]

   return @inbounds -(∂r_Uφ_ccc * ur′_ccc * uφ′_ccc)
end

BTI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂rUφ_ur′_uφ′_ccc, grid, model.velocities.u, model.velocities.v, Ur, Uφ, ∂r_Uφ)

function compute_BTI_transfer(sim)
   compute!(Integral(Field(BTI_transfer_op)))
end

@inline b′uz′_ccc(i, j, k, grid, b, uz, B, Uz) = @inbounds (
		  (b[i, j, k] - B[i, j, k]) * ℑzᵃᵃᶜ(i, j, k, grid, uz′, uz, Uz))

pAPE_to_pKE_op = KernelFunctionOperation{Center, Center, Center}(b′uz′_ccc,
   grid, model.tracers.b, model.velocities.w, B, Uz)

function compute_integrated_pAPE_to_pKE(sim)
   compute!(Integral(Field(pAPE_to_pKE_op)))
end

function ∂zUφ_uφ′_uz′_ccc(i, j, k, grid, ux, uy, uz, Uφ, Uz, ∂z_Uφ)

   φ = atan(ℑyᵃᶜᵃ(i, j, k, grid, ynodes(grid, Center())), ℑxᶜᵃᵃ(i, j, k, grid, xnodes(grid, Center())))

   ux_ccc = @inbounds ℑxᶜᵃᵃ(i, j, k, grid, ux)
   uy_ccc = @inbounds ℑyᵃᶜᵃ(i, j, k, grid, uy)
   uφ_ccc = @inbounds (uy_ccc * cos(φ)) - (ux_ccc * sin(φ)) 

   uφ′_ccc = uφ′(i, j, k, grid, uφ_ccc, Uφ)
   uz′_ccc = @inbounds ℑzᵃᵃᶜ(i, j, k, grid, uz′, uz, Uz) 
   ∂z_Uφ_ccc = @inbounds ∂z_Uφ[i, j, k]

   return @inbounds -(∂z_Uφ_ccc * uφ′_ccc * uz′_ccc)
end

BCI_transfer_op = KernelFunctionOperation{Center, Center, Center}(∂zUφ_uφ′_uz′_ccc, grid, model.velocities.u, model.velocities.v, model.velocities.w, Uφ, Uz, ∂z_Uφ)

function compute_BCI_transfer(sim)
   compute!(Integral(Field(BCI_transfer_op)))
end

#############################
# SET UP AND RUN SIMULATION #
#############################

#Add random perturbations to horizontal velocity components

@inline u_perturbed(x, y, z) = ū(x, y, z)*(1 + 2*(rand()-0.5) * max_u′/(U*sqrt(2)))

if !isnothing(seed2)
   Random.seed!(seed2) #Update seed so next random number is independent
end

@inline v_perturbed(x, y, z) = v̄(x, y, z)*(1 + 2*(rand()-0.5) * max_u′/(U*sqrt(2)))

set!(model, u = u_perturbed, v = v_perturbed) #Set perturbed ICs

initial_KE = CenterField(model.grid)
set!(initial_KE, (model.velocities.u^2 + model.velocities.v^2
                  + model.velocities.w^2) / 2)
total_initial_KE = Field(Integral(initial_KE))
compute!(total_initial_KE)

simulation = Simulation(model; Δt = Δt,
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

energy_diagnostics = (; integrated_pKE = compute_integrated_pKE(simulation),
	integrated_pAPE_to_pKE = compute_integrated_pAPE_to_pKE(simulation),
	integrated_BTI_transfer = compute_BTI_transfer(simulation),
	integrated_BCI_transfer = compute_BCI_transfer(simulation))

energy_writer = NetCDFWriter(model, energy_diagnostics,
				   filename = energyfilepath,
			 	   schedule = TimeInterval(Δt_save),
			     file_splitting = FileSizeLimit(30GiB))

checkpointer = Checkpointer(model; schedule = TimeInterval(Δt_checkpt),
			    		dir = "Checkpoints", 
			    	     prefix = "checkpoint_$(datetimenow)", 
			    	 properties = [:grid, :clock, :timestepper,
					       :velocities, :tracers])

simulation.output_writers[:field_writer]  = field_writer
simulation.output_writers[:scalar_writer] = scalar_writer
simulation.output_writers[:energy_writer] = energy_writer
simulation.output_writers[:checkpointer]  = checkpointer

run!(simulation)#; pickup = joinpath("./Checkpoints", "checkpoint_260306-200453_iteration360000.jld2"))

duration = canonicalize(now() - datetimestart)

#Append zeros to filenames so they can be accessed in chronological order
pad_filenames(datetimenow)
pad_filenames(datetimenow; prefix = "energetics")

#Save parameters to logfile
open(logfilepath, "w") do file
   write(file, "Nx, Ny, Nz = $(Nx), $(Ny), $(Nz) \n")
   write(file, "Lr, Lz = $(Lr), $(Lz) \n\n")
   write(file, "σr, σz = $(σr), $(σz) \n")
   write(file, "U, N²₀ = $(U), $(N²₀) \n")
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
                idxStartLinGrowth_uz = 24, idxEndLinGrowth_uz = 37)
end

if vis_energetics
   visualize_energetics(datetimenow, model.grid, total_initial_KE.data)
end

if vis_z_grid
   visualize_z_grid(datetimenow, model.grid, -Lz)
end