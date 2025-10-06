using Adapt
using Oceananigans
using Oceananigans.Advection
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Grids
using Oceananigans.Units
using Printf

Nz = 8
Hz = 3

grid = RectilinearGrid(CPU(), topology = (Flat, Flat, Bounded), size = (Nz), z = (-1000 * meter, 0))

N²       = 3e-3 * (second^(-2))
b̄z_BC(t) = N²
b̄(z)     = N² * z
b̄_BCs    = FieldBoundaryConditions(top = GradientBoundaryCondition(b̄z_BC),
                                bottom = GradientBoundaryCondition(b̄z_BC))


model = NonhydrostaticModel(;
                            grid = grid,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            tracers = (:b),
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; b = b̄_BCs,))

set!(model, b = b̄)

b_prt = parent(model.tracers.b)
z_prt = parent(grid.zᵃᵃᶜ)
view(b_prt, 1, 1, 1:Hz) .= b̄(view(z_prt, 1:Hz))
view(b_prt, 1, 1, Nz+Hz+1:Nz+2Hz) .= b̄(view(z_prt, Nz+Hz+1:Nz+2Hz))
#fill_halo_regions!(model.tracers.b)

simulation = Simulation(model, Δt = 2.5 * second, stop_time = 2.5 * second)

function progress(sim)
   @info @sprintf("Iteration %d; time = %.2e days",
                  iteration(sim), (time(sim)/day))
   print("Buoyancy slice: ", sim.model.tracers.b[1, 1, :], "\n")
   return nothing
end

add_callback!(simulation, progress, IterationInterval(1))

run!(simulation)
