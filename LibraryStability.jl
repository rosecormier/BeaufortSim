include("LibraryVisualization.jl")
using .VisFunctions

using Adapt, CairoMakie, OffsetArrays
using Oceananigans.AbstractOperations, Oceananigans.Fields, Oceananigans.Operators 

####################

module Stability
   export check_inert_stability, check_grav_stability, compute_Bu
end

####################

function ζz_abs_ffc(i, j, k, grid, f, u, v)
   return f + ζ₃ᶠᶠᶜ(i, j, k, grid, u, v)
end

function check_inert_stability(grid, f, u, v; 
		               z_idx = nothing, plot_ζz_abs = false)
  
   ζz_abs_KernOp = KernelFunctionOperation{Face, Face, Center}(ζz_abs_ffc,
							       grid, f, u, v)
   ζz_abs        = Field(ζz_abs_KernOp)
   
   compute!(ζz_abs)

   if any(z -> z <= 0, ζz_abs)
      print("Warning: system is inertially unstable.")
   end

   if plot_ζz_abs

      mkpath("./Plots") #Make visualization directory if nonexistent
      
      x = xnodes(grid, Face()) ./ 1000 #Convert to km for readability
      y = ynodes(grid, Face()) ./ 1000 #Convert to km for readability
      z = znodes(grid, Center())

      if !isnothing(z_idx)

	 ζz_abs_z     = @views OffsetArrays.no_offset_view(
						      ζz_abs[:, :, z_idx].data)
	 ζz_abs_slice = @views adapt(Array, ζz_abs_z)[:, :, 1]

         idx_kwargs           = (z_idx = z_idx,)
         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
							   z_idx = z_idx)

         h_dim, v_dim, const_dim, units = x, y, "z", "m"
      end

      fig = Figure(size = (500, 400))
      ax  = Axis(fig[2, 1];
                 title = "Absolute vorticity; t = 0", axis_kwargs...)
      hm  = heatmap!(ax, h_dim, v_dim, ζz_abs_slice,
	             colorrange = get_range_lims(ζz_abs_slice), 
		     colormap = :balance)

      Colorbar(fig[2, 2], hm, tickformat = "{:.1e}", label = "1/s")

      title = @sprintf("Absolute vorticity at %s = %i %s; t = 0.0 days",
		       const_dim, nearest, units)

      fig[1, 1:2] = Label(fig, title, fontsize = 18, tellwidth = false)

      save(joinpath("./Plots", "zeta_abs_$(const_dim)$(nearest)_t0.png"), fig)
   end
end

function check_grav_stability(b; plot_∂b∂z = false)
   
   if any(n -> n <= 0, ∂z(b))
      print("Warning: system is gravitationally unstable.")
   end

   #if plot_∂b∂z

   #end
end

function compute_Bu(σr, σz, f, N²)
   Bu = N² * (σz / (f * σr))^2
end
