using Oceananigans.Architectures
using Oceananigans.Grids
using Oceananigans.Units

#Geophysical constants
const GravitationalAcceleration = -9.80665 * meter * (second^2)
const EarthRotationRate         = 7.2921e-5 * (second^(-1))

@kwdef struct GridParams
   
   architecture::AbstractArchitecture = GPU()
   
   topology::Tuple = (Bounded, Bounded, Bounded)
   
   #Numbers of gridpoints
   Nx::Int64 = 512
   Ny::Int64 = 512
   Nz::Int64 = 100
   
   #Halo sizes
   Hx::Int64 = 3
   Hy::Int64 = 3
   Hz::Int64 = 3
   
   #Axis sizes
   Lr::Float64 = 1e3 * kilometer
   Lz::Float64 = 1 * kilometer
   
   z_grid_type::Symbol = :uniform
end

@kwdef struct ModelParams
   useNonhydrostaticModel::Bool = true
   timestepper::Symbol          = :RungeKutta3
end

@kwdef struct GyreParams

   lat::Float64 = 74.0                              #Latitude (deg. N)
   f::Float64   = 2 * EarthRotationRate * sind(lat) #Coriolis frequency (s^{-1})
   U::Float64   = 5e-2 * (meter/second^2)           #Maximum gyre speed (at surface)
   σr::Float64  = 250 * kilometer                   #Radial gyre length scale
   σz::Float64  = 300.0 * meter                     #Vertical gyre length scale
   
   #Ambient (i.e., excluding gyre's TWB contribution) N²-value at z -> -infty
   N²_far::Float64 = 1e-4 * second^(-2)
   
   doubleTanhAmbientStrat::Bool = false
   #If false, defaults to constant ambient stratification (note that a TWB
   # contribution will always be included).
end

@kwdef struct PerturbationParams

   #Maximum relative amplitude of initial streamfunction perturbation
   max_ψ′::Float64 = 1e-2
   
   seed::Int64 = 12345 #Random number seed for reproducibility
end

@kwdef struct DoubleTanhStratParams
   #=
   Defined as in Kosty et al. (2026).
   =#
   
   gravity::Float64 = GravitationalAcceleration 
   ρ₀::Float64      = 1025.5 * meter^(-3)
   A_s::Float64     = 2.5 * meter^(-3)
   z_s::Float64     = -40 * meter
   C_s::Float64     = 15 * meter
   A_d::Float64     = 1.05 * meter^(-3)
   z_d::Float64     = -200 * meter
   C_d::Float64     = 60 * meter
end

@kwdef struct SimulationParams

   Δt_checkpt::Float64            = 250 * day #Checkpoint interval
   align_time_step::Bool          = false
   minimum_relative_step::Float64 = 1e-9

end

@kwdef struct VisParams
   
   #Whether to run visualization functions
   vis_const_x    = true
   vis_const_y    = false
   vis_const_z    = true
   vis_norms      = true
   vis_energetics = true
   vis_z_grid     = false #Note: currently can only be done on CPU
end