using Adapt
using Oceananigans
using Oceananigans.Advection
using Oceananigans.Architectures
using Oceananigans.BoundaryConditions
using Oceananigans.Grids
using Oceananigans.Units
using Printf, Random

Nz = 8
Hz = 1 #3

grid = RectilinearGrid(CPU(), 
		       topology = (Flat, Bounded, Bounded), 
		       size = (Nz, Nz),
		       y = (-1000 * meter, 1000 * meter),
		       z = (-1000 * meter, 0),
		       halo = (Hz, Hz))

N²       = 3e-3 * (second^(-2))
b̄z_BC(y,t) = N²
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

Random.seed!(12345)

b_perturbed(y, z) = b̄(z) + rand()
v_initial(y, z)   = 0.0001 * y^2

set!(model, b = b_perturbed, v = v_initial)

b_prt = parent(model.tracers.b)
z_prt = parent(grid.zᵃᵃᶜ)

for jj = 1:Nz+2Hz
   view(b_prt, 1, jj, 1:Hz) .= b̄(view(z_prt, 1:Hz))
   view(b_prt, 1, jj, Nz+Hz+1:Nz+2Hz) .= b̄(view(z_prt, Nz+Hz+1:Nz+2Hz))
end

simulation = Simulation(model, Δt = 2.5 * second, stop_time = 1*hour) #10*second)

function progress(sim)
   @info @sprintf("Iteration %d; time = %.2e days",
                  iteration(sim), (time(sim)/day))
   print("Buoyancy slice: ", sim.model.tracers.b[1, 1, :], "\n")
   return nothing
end

add_callback!(simulation, progress, IterationInterval(200))

run!(simulation)
