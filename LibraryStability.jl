include("LibraryVisualization.jl")

using Adapt, CairoMakie
using Oceananigans.AbstractOperations
using Oceananigans.Fields
using Oceananigans.Operators 
using OffsetArrays, Printf

@inline ζz_abs_ffc(i, j, k, g, f, u, v) = @inbounds f + ζ₃ᶠᶠᶜ(i, j, k, g, u, v)

function check_inertial_stability(grid, f, u, v; 
		                              plot_ζz_abs = false, 
                                  x_idx = nothing, 
			                            y_idx = nothing, 
                                  z_idx = nothing)
  
   ζz_abs_KernOp = KernelFunctionOperation{Face, Face, Center}(ζz_abs_ffc, grid, f, u, v)
   
   @compute ζz_abs = Field(ζz_abs_KernOp)

   if any(z -> z < 0, ζz_abs)
      print("Warning: system is inertially unstable.\n")
   end

   if plot_ζz_abs

      mkpath("./Plots") #Make visualization directory if nonexistent
      
      #Convert spatial coords to km for readability
      x = no_offset_view(adapt(Array, grid.xᶠᵃᵃ)
                        )[(Hx + 1):(length(grid.xᶠᵃᵃ) - Hx)] ./ 1000
      y = no_offset_view(adapt(Array, grid.yᵃᶠᵃ)
                        )[(Hy + 1):(length(grid.yᵃᶠᵃ) - Hy)] ./ 1000
      z = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ)
                        )[(Hz + 1):(length(grid.z.cᵃᵃᶜ) - Hz)] ./ 1000

      if !isnothing(x_idx)

         ζz_abs_slice = interior(ζz_abs)[x_idx, :, :]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "x";
                                                           x_idx = x_idx)

         h_dim, v_dim, const_dim, units = y, z, "x", "km"

      elseif !isnothing(y_idx)

         ζz_abs_slice = interior(ζz_abs)[:, y_idx, :]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "y";
                                                           y_idx = y_idx)

         h_dim, v_dim, const_dim, units = x, z, "y", "km"

      elseif !isnothing(z_idx)

         ζz_abs_slice = interior(ζz_abs)[:, :, z_idx]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "z";
							                                             z_idx = z_idx)

         h_dim, v_dim, const_dim, units = x, y, "z", "m"
      end

      fig = Figure(size = (500, 400))
      ax  = Axis(fig[2, 1]; axis_kwargs...)
      hm  = heatmap!(ax, h_dim, v_dim, ζz_abs_slice,
	                   colorrange = get_symm_range_lims(ζz_abs_slice), colormap = :balance)

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
   
   ∂b∂z = ZFaceField(grid)
   set!(∂b∂z, ∂z(b))

   if any(n -> n < 0, interior(∂b∂z))
      print("Warning: system is gravitationally unstable.\n")
   end
   
   if plot_∂b∂z

      mkpath("./Plots") #Make visualization directory if nonexistent

      #Load interior spatial coords; convert them to km for readability
      x = no_offset_view(adapt(Array, grid.xᶜᵃᵃ)
                        )[(Hx + 1):(length(grid.xᶜᵃᵃ) - Hx)] ./ 1000
      y = no_offset_view(adapt(Array, grid.yᵃᶜᵃ)
                        )[(Hy + 1):(length(grid.yᵃᶜᵃ) - Hy)] ./ 1000
      z = no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ)
                        )[(Hz + 1):(length(grid.z.cᵃᵃᶠ) - Hz)] ./ 1000
      
      if !isnothing(x_idx)

         ∂b∂z_slice = interior(∂b∂z)[x_idx, :, :]
         
         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "x";
                                                           x_idx = x_idx)

         h_dim, v_dim, const_dim, units = y, z, "x", "km"

      elseif !isnothing(y_idx)

         ∂b∂z_slice = interior(∂b∂z)[:, y_idx, :]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "y";
                                                           y_idx = y_idx)

         h_dim, v_dim, const_dim, units = x, z, "y", "km"

      elseif !isnothing(z_idx)

         ∂b∂z_slice = interior(∂b∂z)[:, :, z_idx]

         nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z, "z";
                                                           z_idx = z_idx)

         h_dim, v_dim, const_dim, units = x, y, "z", "m"
      end

      fig = Figure(size = (500, 400))
      ax  = Axis(fig[2, 1]; axis_kwargs...)
      hm  = heatmap!(ax, h_dim, v_dim, ∂b∂z_slice,
                     colorrange = get_symm_range_lims(∂b∂z_slice), 
                     colormap = :balance)

      Colorbar(fig[2, 2], hm, tickformat = "{:.1e}", label = "1/s²")

      title = @sprintf("∂b/∂z at %s = %i %s; t = 0.0 days",
                       const_dim, nearest, units)

      fig[1, 1:2] = Label(fig, title, fontsize = 18, tellwidth = false)

      save(joinpath("./Plots", "dbdz_$(const_dim)$(nearest)_t0.png"), fig)
   end
end

@inline compute_Bu(σr, σz, f, N²) = N² * (σz / (f * σr))^2