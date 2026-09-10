include("../LibraryCoordinates.jl")
include("../LibraryDynamics.jl")
include("../LibraryStability.jl")
include("../LibraryVisualization.jl")

using CairoMakie, CUDA, LinearAlgebra
using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Coriolis
using Oceananigans.Fields
using Oceananigans.Solvers
using Oceananigans.Units
using Oceananigans.Utils 
using Printf

######################
# SPECIFY PARAMETERS #
######################

Nxs = Vector{Int64}([12, 25, 50, 100, 200, 400, 800])

const Nz = 10 #z-grid size

const Hx = 3 #Number of x halo cells per boundary
const Hy = 3 #Number of y halo cells per boundary
const Hz = 3 #Number of z halo cells per boundary

const Lr = 2.5e3 * kilometer #[Minimum] domain radius
const Lz = 1 * kilometer     #z-axis length

const lat = 74.0     #Latitude (deg. N)
fPlane    = FPlane(latitude = lat)
const f   = fPlane.f #Coriolis frequency

const U  = 5e-2 * (meter/second) #Maximum gyre speed (at surface)
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

const useGPU = false #Whether to use GPU
const useNHS = true #Whether to use NonhydrostaticModel

#########################
# SET UP GRID AND MODEL #
#########################

Ur_norms  = Vector{Float64}()

for Nx in Nxs

   print("computing with Nx = $(Nx)", "\n")
   
   Ny    = Nx
   yFlat = false

   gridParams = (architecture = CPU(),
              Tx = Periodic, Ty = Periodic, Tz = Bounded,
              Nx = Nx, Ny = Ny, Nz = Nz,
              Hx = Hx, Hy = Hy, Hz = Hz,
              Lr = Lr, Lz = Lz,
              z_grid_type = z_grid)

   grid = build_Oceananigans_RectilinearGrid(gridParams)

   B_vals, Ux_vals, Uy_vals, Uz_vals, B_BCs = discrete_Cartesian_TWB_ICs(
       grid, gridParams, gyreScaleParams, bkgd_Ψ_cylindrical_coords, ambientStrat, useGPU;
       Hz = Hz, visualizePsi = false)

   if useNHS
      model = NonhydrostaticModel(;
                        grid = grid, 
                        timestepper = :RungeKutta3,
                        advection = WENO(),
                        coriolis = fPlane,
                        pressure_solver = FourierTridiagonalPoissonSolver(grid),
                        hydrostatic_pressure_anomaly = CenterField(grid),
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

   #Must 'set!' in separate lines because we are setting with Fields, not functions
   set!(model.velocities.u, Ux_vals)
   set!(model.velocities.v, Uy_vals)
   set!(model.velocities.w, Uz_vals)
   set!(model.tracers.b, B_vals)
   fill_halo_regions!(model.tracers.b)
   fill_halo_regions!(model.velocities.u)
   fill_halo_regions!(model.velocities.v)

   #Print warnings if the respective instabilities are present
   check_inertial_stability(model.grid, f, model.velocities.u, model.velocities.v)
   check_gravitational_stability(model.tracers.b, model.grid)

   Ur_vals = xy_vector_to_rφ(Ux_vals, Uy_vals, model.grid, useGPU)[1]

   Ur = CenterField(model.grid)
   
   set!(Ur, Ur_vals)
   
   #Note: this syntax is necessary because Ur isn't a prognostic model field
   fill_halo_regions!(Ur, model.clock, Ur_vals)
   
   Ur_norm = norm(Ur)
   
   push!(Ur_norms, Ur_norm)
end

fig = Figure(size = (600, 600))
ax  = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = L"$N_x=N_y$", ylabel = L"$||U_r||_2$")

scatter!(ax, Nxs, Ur_norms)
save("convergencerates.png", fig)