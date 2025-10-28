include("LibraryVisualization.jl")
include("LibraryStability.jl")
include("Visualization.jl")

using Dates: canonicalize, format, now
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.Units
using Printf, Random

Nx, Nz = 40, 12

Lx, Lz = 2.5e6, 1e3

N²₀ = 3e-3

Δt = 10*minutes
Δt_save = 1*hour
tf = 100*day

grid = RectilinearGrid(topology = (Bounded, Flat, Bounded),
                           size = (Nx, Nz), 
                              x = (-Lx/2, Lx/2), 
                              z = (-Lz, 0.0),
		                   halo = (3, 3))

b̄_BCs = FieldBoundaryConditions(top    = GradientBoundaryCondition(N²₀),
				                bottom = GradientBoundaryCondition(N²₀))

model = NonhydrostaticModel(; 
                            grid = grid, 
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            coriolis = FPlane(latitude = 74.0),
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; b = b̄_BCs))

Ur_vals, Uφ_vals = xy_vector_to_rφ(model.velocities.u,
                                   model.velocities.v, model.grid)

u_perturbed(x, z)  = (2 * (rand() - 0.5)) * 1e-8
v_perturbed(x, z)  = (2 * (rand() - 0.5)) * 1e-8
b_background(x, z) = N²₀*z

set!(model, u = u_perturbed, v = v_perturbed, b = b_background)

#Prints warnings if the respective instabilities are present
check_inert_stability(model.grid, model.coriolis.f, model.velocities.u, model.velocities.v;
		      plot_ζz_abs = false, z_idx = 6)
check_grav_stability(model.tracers.b; plot_∂b∂z = false, grid = model.grid,
                     x_idx = 20)

datetimestart = now()
datetimenow   = format(datetimestart, "yymmdd-HHMMSS")
print("Date-time label: $(datetimenow)", "\n")

#Create fields to store background state

Ux = XFaceField(model.grid)
Uy = YFaceField(model.grid)
Ur = CenterField(model.grid)
Uφ = CenterField(model.grid)
Uz = ZFaceField(model.grid)
B  = CenterField(model.grid)

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

simulation = Simulation(model, Δt = Δt, stop_time = tf)

function progress(sim)
   umax = maximum(abs, sim.model.velocities.u)
   wmax = maximum(abs, sim.model.velocities.w)
   bmax = maximum(abs, sim.model.tracers.b)
   @info @sprintf("Iter: %d; max|u|: %.2e; max|w|: %.2e; max|b|: %.2e", 
                    iteration(sim), umax, wmax, bmax)
   return nothing
end

add_callback!(simulation, progress, TimeInterval(Δt_save))

ur, uφ = xy_vector_to_rφ(model.velocities.u, model.velocities.v, model.grid)

outputs = (ur = ur,
           uφ = uφ,
           ux = model.velocities.u,
           uy = model.velocities.v,
           uz = model.velocities.w,
           b  = model.tracers.b)

outfilepath = joinpath("./Output", "output_$(datetimenow).nc")
mkpath(dirname(outfilepath)) #Make path if nonexistent

field_writer = NetCDFOutputWriter(model,
                                  outputs,
                                  with_halos = true,
                                  filename = outfilepath,
                                  schedule = TimeInterval(Δt_save),
                                  file_splitting = FileSizeLimit(30GiB))

simulation.output_writers[:field_writer] = field_writer

run!(simulation)

visualize_fields_const_y(datetimenow, nothing, B, Uφ;
                            plot_animation = true, t_idx_skip = 1)
