include("LibraryVisualization.jl")

using Adapt, CairoMakie
using Oceananigans.AbstractOperations
using Oceananigans.Fields
using Oceananigans.Operators 
using OffsetArrays, Printf

####################

function ζz_abs_ffc(i, j, k, g, f, u, v)
   return f + ζ₃ᶠᶠᶜ(i, j, k, g, u, v)
end

function check_inertial_stability(grid, f, u, v; 
		                              plot_ζz_abs = false, 
                                  x_idx = nothing, 
			                            y_idx = nothing, 
                                  z_idx = nothing)
  
   ζz_abs_KernOp = KernelFunctionOperation{Face, Face, Center}(ζz_abs_ffc, grid, f, u, v)
   ζz_abs        = Field(ζz_abs_KernOp)
   
   compute!(ζz_abs)

   if any(z -> z < 0, ζz_abs)
      print("Warning: system is inertially unstable.\n")
   end

   if plot_ζz_abs

      mkpath("./Plots") #Make visualization directory if nonexistent
      
      x = xnodes(grid, Face()) ./ 1000 #Convert to km for readability
      y = ynodes(grid, Face()) ./ 1000 #Convert to km for readability
      z = znodes(grid, Center())

      if !isnothing(x_idx)

         ζz_abs_x     = @views OffsetArrays.no_offset_view(ζz_abs[x_idx, :, :].data)
         ζz_abs_slice = @views adapt(Array, ζz_abs_x)[1, :, :]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "x";
                                                           x_idx = x_idx)

         h_dim, v_dim, const_dim, units = y, z, "x", "km"

      elseif !isnothing(y_idx)

         ζz_abs_y     = @views OffsetArrays.no_offset_view(ζz_abs[:, y_idx, :].data)
         ζz_abs_slice = @views adapt(Array, ζz_abs_y)[:, y, :]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "y";
                                                           y_idx = y_idx)

         h_dim, v_dim, const_dim, units = x, z, "y", "km"

      elseif !isnothing(z_idx)

         ζz_abs_z     = @views OffsetArrays.no_offset_view(ζz_abs[:, :, z_idx].data)
	       ζz_abs_slice = @views adapt(Array, ζz_abs_z)[:, :, 1]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "z";
							                                             z_idx = z_idx)

         h_dim, v_dim, const_dim, units = x, y, "z", "m"
      end

      fig = Figure(size = (500, 400))
      ax  = Axis(fig[2, 1]; axis_kwargs...)
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

function check_gravitational_stability(b, grid; 
                                       plot_∂b∂z = false,
		                                   x_idx = nothing, 
                                       y_idx = nothing, 
                                       z_idx = nothing,
                                       Hx = 3, Hy = 3, Hz = 3)
   
   @compute @at (Center, Center, Center) ∂b∂z = Field(∂z(b))
   ∂b∂z_array = no_offset_view(adapt(Array, ∂b∂z))
   
   ∂b∂z_interior = ∂b∂z_array[Hx + 1:length(no_offset_view(grid.xᶜᵃᵃ)) - Hx,
                              Hy + 1:length(no_offset_view(grid.yᵃᶜᵃ)) - Hy, 
                              Hz + 1:length(grid.z.cᵃᵃᶜ) - Hz]
   
   if any(n -> n < 0, ∂b∂z_interior)
      print("Warning: system is gravitationally unstable.\n")
   end
   
   if plot_∂b∂z

      mkpath("./Plots") #Make visualization directory if nonexistent

      x = adapt(Array, grid.xᶜᵃᵃ)[Hx+1:length(grid.xᶜᵃᵃ)-Hx] ./ 1000 #Convert to km for readability
      y = adapt(Array, grid.yᵃᶜᵃ)[Hy+1:length(grid.yᵃᶜᵃ)-Hy] ./ 1000 #Convert to km for readability
      z = adapt(Array, grid.z.cᵃᵃᶜ)[1:length(grid.z.cᵃᵃᶜ)-2*Hz]

      if !isnothing(x_idx)

         ∂b∂z_slice = ∂b∂z_interior[x_idx, :, :] #@views adapt(Array, ∂z(b))[x_idx, :, :]
         
         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "x";
                                                           x_idx = x_idx)

         h_dim, v_dim, const_dim, units = y, z, "x", "km"

      elseif !isnothing(y_idx)

         ∂b∂z_slice = @views adapt(Array, ∂z(b))[:, y_idx, :]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "y";
                                                           y_idx = y_idx)

         h_dim, v_dim, const_dim, units = x, z, "y", "km"

      elseif !isnothing(z_idx)

         ∂b∂z_slice = @views adapt(Array, ∂z(b))[:, :, z_idx]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "z";
                                                           z_idx = z_idx)

         h_dim, v_dim, const_dim, units = x, y, "z", "m"
      end

      fig = Figure(size = (500, 400))
      ax  = Axis(fig[2, 1]; axis_kwargs...)
      hm  = heatmap!(ax, h_dim, v_dim, ∂b∂z_slice,
                     colorrange = get_range_lims(∂b∂z_slice),
                     colormap = :balance)

      Colorbar(fig[2, 2], hm, tickformat = "{:.1e}", label = "1/s²")

      title = @sprintf("∂b/∂z at %s = %i %s; t = 0.0 days",
                       const_dim, nearest, units)

      fig[1, 1:2] = Label(fig, title, fontsize = 18, tellwidth = false)

      save(joinpath("./Plots", "dbdz_$(const_dim)$(nearest)_t0.png"), fig)
   end
end

@inline compute_Bu(σr, σz, f, N²) = N² * (σz / (f * σr))^2