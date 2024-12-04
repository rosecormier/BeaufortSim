include("LibraryVisualization.jl")

using CairoMakie, CommonDataModel, DataStructures, LaTeXStrings, NCDatasets, Oceananigans, Printf
using .ComputeSecondaries

function open_dataset(datetime)

   outfilepath = joinpath("./Output", "output_$(datetime).nc")
   
   ds = NCDataset(outfilepath, "r")
   x  = ds["xC"][:] ./ 1000 #Convert to km for readability
   y  = ds["yC"][:] ./ 1000 #Convert to km for readability
   z  = ds["zC"][:]
   t  = ds["time"][1:end-1] #Drop the last index, in case it contains NaN
   Nt = length(t)

   return ds, x, y, z, t, Nt
end

function get_background_fields(ds)
   bb = ds["b"][:, :, :, 1]
   ub = ds["u"][:, :, :, 1]
   vb = ds["v"][:, :, :, 1]
   wb = ds["w"][:, :, 1:end-1, 1]
   return bb, ub, vb, wb
end

function get_range_lims(final_field; prescribed_max = 0)
   field_max  = max(maximum(abs.(final_field)), prescribed_max)
   field_lims = [-field_max, field_max]
end

function get_2D_spatial_axis_kwargs(x, y, z; 
		                    x_idx = nothing, y_idx = nothing, 
				    z_idx = nothing)
   if !isnothing(x_idx)
      nearest     = round(Int, x[x_idx])
      axis_kwargs = (xlabel = "y [km]", ylabel = "z [m]")
   elseif !isnothing(y_idx)
      nearest     = round(Int, y[y_idx])
      axis_kwargs = (xlabel = "x [km]", ylabel = "z [m]")
   elseif !isnothing(z_idx)
      nearest     = round(Int, -z[z_idx])
      axis_kwargs = (xlabel = "x [km]", ylabel = "y [km]")
   end
   return nearest, axis_kwargs
end

function visualize_growth_rate(datetime)

   ds, x, y, z, times, Nt = open_dataset(datetime)

   fig_gr    = Figure(size = (1200, 700))
   fig_norms = Figure(size = (1200, 700))

   ax_b_gr = Axis(fig_gr[2, 1]; title = "Growth rate of b'", 
		  xlabel = L"$t$ [days]", 
		  ylabel = L"Growth rate of $||b'||$ [1/s]")
   ax_w_gr = Axis(fig_gr[2, 2]; title = "Growth rate of w'", 
	          xlabel = L"$t$ [days]", 
		  ylabel = L"Growth rate of $||w'||$ [1/s]")
   ax_u_gr = Axis(fig_gr[3, 1]; title = "Growth rate of u'",
	          xlabel = L"$t$ [days]", 
		  ylabel = L"Growth rate of $||u'||$ [1/s]")
   ax_v_gr = Axis(fig_gr[3, 2]; title = "Growth rate of v'",
	          xlabel = L"$t$ [days]", 
		  ylabel = L"Growth rate of $||v'||$ [1/s]")

   ax_b_norm = Axis(fig_norms[2, 1]; title = "Norm of b'",
                    xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]")
   ax_w_norm = Axis(fig_norms[2, 2]; title = "Norm of w'",
                    xlabel = L"$t$ [days]", ylabel = L"$||w'||$ [m/s]")
   ax_u_norm = Axis(fig_norms[3, 1]; title = "Norm of u'",
                    xlabel = L"$t$ [days]", ylabel = L"$||u'||$ [m/s]")
   ax_v_norm = Axis(fig_norms[3, 2]; title = "Norm of v'",
                    xlabel = L"$t$ [days]", ylabel = L"$||v'||$ [m/s]")

   b_initial = growth_rate(ds["b"], 1, times)[2]
   u_initial = growth_rate(ds["u"], 1, times)[2]
   v_initial = growth_rate(ds["v"], 1, times)[2]
   w_initial = growth_rate(ds["w"], 1, times)[2]
	
   n      = Observable(2)
   b_gr   = @lift growth_rate(ds["b"], $n, times)[1] / b_initial
   b_norm = @lift growth_rate(ds["b"], $n, times)[2]
   w_gr   = @lift growth_rate(ds["w"], $n, times)[1] / w_initial
   w_norm = @lift growth_rate(ds["w"], $n, times)[2]
   u_gr   = @lift growth_rate(ds["u"], $n, times)[1] / u_initial
   u_norm = @lift growth_rate(ds["u"], $n, times)[2]
   v_gr   = @lift growth_rate(ds["v"], $n, times)[1] / v_initial
   v_norm = @lift growth_rate(ds["v"], $n, times)[2]
   
   @lift scatter!(ax_b_gr, times[$n]/86400, $b_gr, color = :black)
   @lift scatter!(ax_w_gr, times[$n]/86400, $w_gr, color = :black)
   @lift scatter!(ax_u_gr, times[$n]/86400, $u_gr, color = :black)
   @lift scatter!(ax_v_gr, times[$n]/86400, $v_gr, color = :black)

   @lift scatter!(ax_b_norm, times[$n]/86400, $b_norm, color = :black)
   @lift scatter!(ax_w_norm, times[$n]/86400, $w_norm, color = :black)
   @lift scatter!(ax_u_norm, times[$n]/86400, $u_norm, color = :black)
   @lift scatter!(ax_v_norm, times[$n]/86400, $v_norm, color = :black)

   for i = 2:Nt-1
      yield()
      n[] = i
   end

   fig_gr[1, 1:2]    = Label(fig_gr, "Growth rates of perturbation fields", 
		          fontsize = 24, tellwidth = false)
   fig_norms[1, 1:2] = Label(fig_norms, "Norms of perturbation fields",
			  fontsize = 24, tellwidth = false)

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "growth_rates_$(datetime).png"), fig_gr)
   save(joinpath("./Plots", "norm_fields_$(datetime).png"), fig_norms)
   close(ds)
end

function visualize_b_and_ωz(datetime, z_idx, Δx, Δy; 
		            plot_animation = false, t_idx_skip = 1)

   ds, x, y, z, times, Nt = open_dataset(datetime)
   bb, ub, vb, wb         = get_background_fields(ds)

   ωb_xy = ωz(ds["u"][:, :, z_idx, 1], ds["v"][:, :, z_idx, 1], Δx, Δy) 

   b_total_f_xy = ds["b"][:, :, z_idx, Nt]
   ω_total_f_xy = ωz(ds["u"][:, :, z_idx, Nt], ds["v"][:, :, z_idx, Nt], 
		     Δx, Δy)

   Δb_f_xy = b_total_f_xy .- bb[:, :, z_idx]
   Δω_f_xy = ω_total_f_xy .- ωb_xy

   lims_b_total = get_range_lims(b_total_f_xy)
   lims_ω_total = get_range_lims(ω_total_f_xy)

   lims_Δb = get_range_lims(Δb_f_xy; prescribed_max = 1e-16)
   lims_Δω = get_range_lims(Δω_f_xy; prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   depth_nearest, axis_kwargs_xy = get_2D_spatial_axis_kwargs(x, y, z; z_idx = z_idx)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      n = Observable(1)

      b_total_xy = @lift ds["b"][:, :, z_idx, $n]
      ω_total_xy = @lift ωz(ds["u"][:, :, z_idx, $n], ds["v"][:, :, z_idx, $n], 
			    Δx, Δy)

      Δb_xy = @lift $b_total_xy .- bb[:, :, z_idx]
      Δω_xy = @lift $ω_total_xy .- ωb_xy 
   
      fig_total   = Figure(size = (1200, 500))
      fig_perturb = Figure(size = (1200, 500))

      ax_b_total = Axis(fig_total[2, 1];
                        title = "Total buoyancy (b)", axis_kwargs_xy...)
      ax_ω_total = Axis(fig_total[2, 3];
                        title = "Total vertical vorticity (ζ)", 
			axis_kwargs_xy...)

      ax_b_perturb = Axis(fig_perturb[2, 1];
                          title = "Buoyancy perturbation (b')", 
			  axis_kwargs_xy...)
      ax_ω_perturb = Axis(fig_perturb[2, 3];
                          title = "Vertical vorticity perturbation (ζ')",
                          axis_kwargs_xy...)
      
      hm_b_total = heatmap!(ax_b_total, x, y, b_total_xy,
                            colorrange = lims_b_total, colormap = :balance)
      hm_ω_total = heatmap!(ax_ω_total, x, y, ω_total_xy,
                            colorrange = lims_ω_total, colormap = :balance)

      hm_b_perturb = heatmap!(ax_b_perturb, x, y, Δb_xy,
                              colorrange = lims_Δb, colormap = :balance)
      hm_ω_perturb = heatmap!(ax_ω_perturb, x, y, Δω_xy,
                              colorrange = lims_Δω, colormap = :balance)

      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", label = "1/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
               label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}",
               label = "1/s")

      title_total = @lift @sprintf("Fields at %i-m depth; t = %.2f days",
                                   depth_nearest, times[$n]/(3600*24))
      title_perturb = @lift @sprintf(
                            "Perturbation fields at %i-m depth; t = %.2f days",
                                     depth_nearest, times[$n]/(3600*24))

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                                  tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                                  tellwidth = false)
   
      frames = 1:Nt
   
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
         recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end
   
      save(joinpath("./Plots", "bzeta_total_z-$(depth_nearest)_$(datetime).mp4"), 
	            video_total)
      save(joinpath("./Plots", "bzeta_perturbs_z-$(depth_nearest)_$(datetime).mp4"),
	            video_perturb)
   end

   #Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 500))
   fig_perturb = Figure(size = (1200, 500))

   ax_b_total = Axis(fig_total[2, 1];
                     title = "Total buoyancy (b)", axis_kwargs_xy...)
   ax_ω_total = Axis(fig_total[2, 3];
                     title = "Total vertical vorticity (ζ)", 
		     axis_kwargs_xy...)

   ax_b_perturb = Axis(fig_perturb[2, 1];
                       title = "Buoyancy perturbation (b')",
                       axis_kwargs_xy...)
   ax_ω_perturb = Axis(fig_perturb[2, 3];
                       title = "Vertical vorticity perturbation (ζ')",
                       axis_kwargs_xy...)

   hm_b_total = heatmap!(ax_b_total, x, y, b_total_f_xy,
                         colorrange = lims_b_total, colormap = :balance)
   hm_ω_total = heatmap!(ax_ω_total, x, y, ω_total_f_xy,
                         colorrange = lims_ω_total, colormap = :balance)

   hm_b_perturb = heatmap!(ax_b_perturb, x, y, Δb_f_xy,
                           colorrange = lims_Δb, colormap = :balance)
   hm_ω_perturb = heatmap!(ax_ω_perturb, x, y, Δω_f_xy,
                           colorrange = lims_Δω, colormap = :balance)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ω_total, tickformat = "{:.1e}", label = "1/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ω_perturb, tickformat = "{:.1e}",
            label = "1/s")

   title_total   = @sprintf("Fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt]/(3600*24))
   title_perturb = @sprintf(
                          "Perturbation fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt]/(3600*24))

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "bzeta_total_z-$(depth_nearest)_tf_$(datetime).png"),
        fig_total)
   save(joinpath("./Plots", "bzeta_perturbs_z-$(depth_nearest)_tf_$(datetime).png"),
        fig_perturb)
   close(ds)
end

function visualize_fields_const_x(datetime, x_idx; 
		                  plot_animation = false, t_idx_skip = 1)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)
   bb, ub, vb, wb         = get_background_fields(ds)

   z_plt = 1 #div(length(z[:]), 2) #z-index to start plotting at

   b_total_f_yz = ds["b"][x_idx, :, z_plt:end, Nt]
   u_total_f_yz = ds["u"][x_idx, :, z_plt:end, Nt]
   v_total_f_yz = ds["v"][x_idx, :, z_plt:end, Nt]
   w_total_f_yz = ds["w"][x_idx, :, z_plt:end-1, Nt]

   Δb_f_yz = b_total_f_yz .- bb[x_idx, :, z_plt:end]
   Δu_f_yz = u_total_f_yz .- ub[x_idx, :, z_plt:end]
   Δv_f_yz = v_total_f_yz .- vb[x_idx, :, z_plt:end]
   Δw_f_yz = w_total_f_yz .- wb[x_idx, :, z_plt:end]

   lims_b_total = get_range_lims(b_total_f_yz)
   lims_u_total = get_range_lims(u_total_f_yz; prescribed_max = 1e-16)
   lims_v_total = get_range_lims(v_total_f_yz; prescribed_max = 1e-16)
   lims_w_total = get_range_lims(w_total_f_yz; prescribed_max = 1e-16)

   lims_Δb = get_range_lims(Δb_f_yz; prescribed_max = 1e-16)
   lims_Δu = get_range_lims(Δu_f_yz; prescribed_max = 1e-16)
   lims_Δv = get_range_lims(Δv_f_yz; prescribed_max = 1e-16)
   lims_Δw = get_range_lims(Δw_f_yz; prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   x_nearest, axis_kwargs_yz = get_2D_spatial_axis_kwargs(x, y, z; x_idx = x_idx)
   
   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      n = Observable(1)

      b_total_yz = @lift ds["b"][x_idx, :, z_plt:end, $n]
      u_total_yz = @lift ds["u"][x_idx, :, z_plt:end, $n]
      v_total_yz = @lift ds["v"][x_idx, :, z_plt:end, $n]
      w_total_yz = @lift ds["w"][x_idx, :, z_plt:end-1, $n]

      Δb_yz = @lift $b_total_yz .- bb[x_idx, :, z_plt:end]
      Δu_yz = @lift $u_total_yz .- ub[x_idx, :, z_plt:end]
      Δv_yz = @lift $v_total_yz .- vb[x_idx, :, z_plt:end]
      Δw_yz = @lift $w_total_yz .- wb[x_idx, :, z_plt:end]
      
      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total = Axis(fig_total[2, 1];
                        title = "Total buoyancy (b)", axis_kwargs_yz...)
      ax_w_total = Axis(fig_total[2, 3];
                        title = "Total vertical velocity (w)", 
			axis_kwargs_yz...)
      ax_u_total = Axis(fig_total[3, 1];
                        title = "Total zonal velocity (u)", axis_kwargs_yz...)
      ax_v_total = Axis(fig_total[3, 3];
                        title = "Total meridional velocity (v)", 
			axis_kwargs_yz...)

      ax_b_perturb = Axis(fig_perturb[2, 1];
                          title = "Buoyancy perturbation (b')", 
			  axis_kwargs_yz...)
      ax_w_perturb = Axis(fig_perturb[2, 3];
                          title = "Vertical velocity perturbation (w')", 
		          axis_kwargs_yz...)
      ax_u_perturb = Axis(fig_perturb[3, 1];
                          title = "Zonal velocity perturbation (u')", 
		          axis_kwargs_yz...)
      ax_v_perturb = Axis(fig_perturb[3, 3];
                          title = "Meridional velocity perturbation (v')", 
	                  axis_kwargs_yz...)

      hm_b_total = heatmap!(ax_b_total, y, z[z_plt:end], b_total_yz, 
			    colorrange = lims_b_total, colormap = :balance)
      hm_w_total = heatmap!(ax_w_total, y, z[z_plt:end], w_total_yz,
                            colorrange = lims_w_total, colormap = :balance)
      hm_u_total = heatmap!(ax_u_total, y, z[z_plt:end], u_total_yz,
                            colorrange = lims_u_total, colormap = :balance)
      hm_v_total = heatmap!(ax_v_total, y, z[z_plt:end], v_total_yz,
                            colorrange = lims_v_total, colormap = :balance)

      hm_b_perturb = heatmap!(ax_b_perturb, y, z[z_plt:end], Δb_yz,
                              colorrange = lims_Δb, colormap = :balance)
      hm_w_perturb = heatmap!(ax_w_perturb, y, z[z_plt:end], Δw_yz,
                              colorrange = lims_Δw, colormap = :balance)
      hm_u_perturb = heatmap!(ax_u_perturb, y, z[z_plt:end], Δu_yz,
                              colorrange = lims_Δu, colormap = :balance)
      hm_v_perturb = heatmap!(ax_v_perturb, y, z[z_plt:end], Δv_yz,
                              colorrange = lims_Δv, colormap = :balance)

      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", 
	       label = "m/s²")
      Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", 
	       label = "m/s")
      Colorbar(fig_total[3, 2], hm_u_total, tickformat = "{:.1e}", 
	       label = "m/s")
      Colorbar(fig_total[3, 4], hm_v_total, tickformat = "{:.1e}", 
	       label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}", 
               label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}", 
               label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_u_perturb, tickformat = "{:.1e}", 
               label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_v_perturb, tickformat = "{:.1e}", 
	       label = "m/s")

      title_total   = @lift @sprintf("Fields at x = %i km; t = %.2f days",
                                     x_nearest, times[$n]/(3600*24))
      title_perturb = @lift @sprintf(
			    "Perturbation fields at x = %i km; t = %.2f days",
                            x_nearest, times[$n]/(3600*24))

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24, 
			          tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24, 
			          tellwidth = false)

      frames = 1:Nt
   
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
         recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end

      save(joinpath("./Plots", "fields_x$(x_nearest)_$(datetime).mp4"),
                    video_total)
      save(joinpath("./Plots", "perturbs_x$(x_nearest)_$(datetime).mp4"),
                    video_perturb)
   end

   #Plot static images (final frame, by default)
   
   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total = Axis(fig_total[2, 1];
                     title = "Total buoyancy (b)", axis_kwargs_yz...)
   ax_w_total = Axis(fig_total[2, 3];
                     title = "Total vertical velocity (w)", axis_kwargs_yz...)
   ax_u_total = Axis(fig_total[3, 1];
                     title = "Total zonal velocity (u)", axis_kwargs_yz...)
   ax_v_total = Axis(fig_total[3, 3];
                     title = "Total meridional velocity (v)", axis_kwargs_yz...)

   ax_b_perturb = Axis(fig_perturb[2, 1]; 
		       title = "Buoyancy perturbation (b')", axis_kwargs_yz...)
   ax_w_perturb = Axis(fig_perturb[2, 3];
                       title = "Vertical velocity perturbation (w')",
                       axis_kwargs_yz...)
   ax_u_perturb = Axis(fig_perturb[3, 1];
                       title = "Zonal velocity perturbation (u')",
                       axis_kwargs_yz...)
   ax_v_perturb = Axis(fig_perturb[3, 3];
                       title = "Meridional velocity perturbation (v')",
                       axis_kwargs_yz...)

   hm_b_total = heatmap!(ax_b_total, y, z[z_plt:end], b_total_f_yz,
                         colorrange = lims_b_total, colormap = :balance)
   hm_w_total = heatmap!(ax_w_total, y, z[z_plt:end], w_total_f_yz,
                         colorrange = lims_w_total, colormap = :balance)
   hm_u_total = heatmap!(ax_u_total, y, z[z_plt:end], u_total_f_yz,
                         colorrange = lims_u_total, colormap = :balance)
   hm_v_total = heatmap!(ax_v_total, y, z[z_plt:end], v_total_f_yz,
                         colorrange = lims_v_total, colormap = :balance)

   hm_b_perturb = heatmap!(ax_b_perturb, y, z[z_plt:end], Δb_f_yz,
                           colorrange = lims_Δb, colormap = :balance)
   hm_w_perturb = heatmap!(ax_w_perturb, y, z[z_plt:end], Δw_f_yz,
                           colorrange = lims_Δw, colormap = :balance)
   hm_u_perturb = heatmap!(ax_u_perturb, y, z[z_plt:end], Δu_f_yz,
                           colorrange = lims_Δu, colormap = :balance)
   hm_v_perturb = heatmap!(ax_v_perturb, y, z[z_plt:end], Δv_f_yz,
                           colorrange = lims_Δv, colormap = :balance)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_u_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_v_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_u_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_v_perturb, tickformat = "{:.1e}",
            label = "m/s")

   title_total   = @sprintf("Fields at x = %i km; t = %.2f days",
                            x_nearest, times[Nt]/(3600*24))
   title_perturb = @sprintf(
                          "Perturbation fields at x = %i km; t = %.2f days",
                            x_nearest, times[Nt]/(3600*24))

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "fields_x$(x_nearest)_tf_$(datetime).png"),
                 fig_total)
   save(joinpath("./Plots", "perturbs_x$(x_nearest)_tf_$(datetime).png"),
                 fig_perturb)
   close(ds)
end

function visualize_fields_const_y(datetime, y_idx; 
		                  plot_animation = false, t_idx_skip = 1)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)
   bb, ub, vb, wb         = get_background_fields(ds)

   z_plt = div(length(z[:]), 2) #z-index to start plotting at

   b_total_f_xz = ds["b"][:, y_idx, z_plt:end, Nt]
   u_total_f_xz = ds["u"][:, y_idx, z_plt:end, Nt]
   v_total_f_xz = ds["v"][:, y_idx, z_plt:end, Nt]
   w_total_f_xz = ds["w"][:, y_idx, z_plt:end-1, Nt]

   Δb_f_xz = b_total_f_xz .- bb[:, y_idx, z_plt:end]
   Δu_f_xz = u_total_f_xz .- ub[:, y_idx, z_plt:end]
   Δv_f_xz = v_total_f_xz .- vb[:, y_idx, z_plt:end]
   Δw_f_xz = w_total_f_xz .- wb[:, y_idx, z_plt:end]

   lims_b_total = get_range_lims(b_total_f_xz)
   lims_u_total = get_range_lims(u_total_f_xz; prescribed_max = 1e-16)
   lims_v_total = get_range_lims(v_total_f_xz; prescribed_max = 1e-16)
   lims_w_total = get_range_lims(w_total_f_xz; prescribed_max = 1e-16)

   lims_Δb = get_range_lims(Δb_f_xz; prescribed_max = 1e-16)
   lims_Δu = get_range_lims(Δu_f_xz; prescribed_max = 1e-16)
   lims_Δv = get_range_lims(Δv_f_xz; prescribed_max = 1e-16)
   lims_Δw = get_range_lims(Δw_f_xz; prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   y_nearest, axis_kwargs_xz = get_2D_spatial_axis_kwargs(x, y, z; y_idx = y_idx)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      n = Observable(1)

      b_total_xz = @lift ds["b"][:, y_idx, z_plt:end, $n]
      u_total_xz = @lift ds["u"][:, y_idx, z_plt:end, $n]
      v_total_xz = @lift ds["v"][:, y_idx, z_plt:end, $n]
      w_total_xz = @lift ds["w"][:, y_idx, z_plt:end-1, $n]

      Δb_xz = @lift $b_total_xz .- bb[:, y_idx, z_plt:end]
      Δu_xz = @lift $u_total_xz .- ub[:, y_idx, z_plt:end]
      Δv_xz = @lift $v_total_xz .- vb[:, y_idx, z_plt:end]
      Δw_xz = @lift $w_total_xz .- wb[:, y_idx, z_plt:end]
   
      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total = Axis(fig_total[2, 1];
                        title = "Total buoyancy (b)", axis_kwargs_xz...)
      ax_w_total = Axis(fig_total[2, 3];
                        title = "Total vertical velocity (w)", 
			axis_kwargs_xz...)
      ax_u_total = Axis(fig_total[3, 1];
                        title = "Total zonal velocity (u)", axis_kwargs_xz...)
      ax_v_total = Axis(fig_total[3, 3];
                        title = "Total meridional velocity (v)", 
			axis_kwargs_xz...)

      ax_b_perturb = Axis(fig_perturb[2, 1];
                          title = "Buoyancy perturbation (b')", 
			  axis_kwargs_xz...)
      ax_w_perturb = Axis(fig_perturb[2, 3];
                          title = "Vertical velocity perturbation (w')",
                          axis_kwargs_xz...)
      ax_u_perturb = Axis(fig_perturb[3, 1];
                          title = "Zonal velocity perturbation (u')",
                          axis_kwargs_xz...)
      ax_v_perturb = Axis(fig_perturb[3, 3];
                          title = "Meridional velocity perturbation (v')",
                          axis_kwargs_xz...)

      hm_b_total = heatmap!(ax_b_total, x, z[z_plt:end], b_total_xz,
                            colorrange = lims_b_total, colormap = :balance)
      hm_w_total = heatmap!(ax_w_total, x, z[z_plt:end], w_total_xz,
                            colorrange = lims_w_total, colormap = :balance)
      hm_u_total = heatmap!(ax_u_total, x, z[z_plt:end], u_total_xz,
                            colorrange = lims_u_total, colormap = :balance)
      hm_v_total = heatmap!(ax_v_total, x, z[z_plt:end], v_total_xz,
                            colorrange = lims_v_total, colormap = :balance)

      hm_b_perturb = heatmap!(ax_b_perturb, x, z[z_plt:end], Δb_xz,
                              colorrange = lims_Δb, colormap = :balance)
      hm_w_perturb = heatmap!(ax_w_perturb, x, z[z_plt:end], Δw_xz,
                              colorrange = lims_Δw, colormap = :balance)
      hm_u_perturb = heatmap!(ax_u_perturb, x, z[z_plt:end], Δu_xz,
                              colorrange = lims_Δu, colormap = :balance)
      hm_v_perturb = heatmap!(ax_v_perturb, x, z[z_plt:end], Δv_xz,
                              colorrange = lims_Δv, colormap = :balance)

      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", 
               label = "m/s²")
      Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", 
	       label = "m/s")
      Colorbar(fig_total[3, 2], hm_u_total, tickformat = "{:.1e}", 
	       label = "m/s")
      Colorbar(fig_total[3, 4], hm_v_total, tickformat = "{:.1e}", 
	       label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
               label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}",
               label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_u_perturb, tickformat = "{:.1e}",
               label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_v_perturb, tickformat = "{:.1e}",
               label = "m/s")

      title_total = @lift @sprintf("Fields at y = %i km; t = %.2f days",
                                   y_nearest, times[$n]/(3600*24))
      title_perturb = @lift @sprintf(
                            "Perturbation fields at y = %i km; t = %.2f days",
                             y_nearest, times[$n]/(3600*24))

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                                  tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
				  tellwidth = false)

      frames = 1:Nt
  
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
         recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end

      save(joinpath("./Plots", "fields_y$(y_nearest)_$(datetime).mp4"), 
                    video_total)
      save(joinpath("./Plots", "perturbs_y$(y_nearest)_$(datetime).mp4"),
                    video_perturb)
   end

   #Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total = Axis(fig_total[2, 1];
                     title = "Total buoyancy (b)", axis_kwargs_xz...)
   ax_w_total = Axis(fig_total[2, 3];
                     title = "Total vertical velocity (w)", axis_kwargs_xz...)
   ax_u_total = Axis(fig_total[3, 1];
                     title = "Total zonal velocity (u)", axis_kwargs_xz...)
   ax_v_total = Axis(fig_total[3, 3];
                     title = "Total meridional velocity (v)", axis_kwargs_xz...)

   ax_b_perturb = Axis(fig_perturb[2, 1];
                       title = "Buoyancy perturbation (b')",
                       axis_kwargs_xz...)
   ax_w_perturb = Axis(fig_perturb[2, 3];
                       title = "Vertical velocity perturbation (w')",
                       axis_kwargs_xz...)
   ax_u_perturb = Axis(fig_perturb[3, 1];
                       title = "Zonal velocity perturbation (u')",
                       axis_kwargs_xz...)
   ax_v_perturb = Axis(fig_perturb[3, 3];
                       title = "Meridional velocity perturbation (v')",
                       axis_kwargs_xz...)

   hm_b_total = heatmap!(ax_b_total, x, z[z_plt:end], b_total_f_xz,
                         colorrange = lims_b_total, colormap = :balance)
   hm_w_total = heatmap!(ax_w_total, x, z[z_plt:end], w_total_f_xz,
                         colorrange = lims_w_total, colormap = :balance)
   hm_u_total = heatmap!(ax_u_total, x, z[z_plt:end], u_total_f_xz,
                         colorrange = lims_u_total, colormap = :balance)
   hm_v_total = heatmap!(ax_v_total, x, z[z_plt:end], v_total_f_xz,
                         colorrange = lims_v_total, colormap = :balance)

   hm_b_perturb = heatmap!(ax_b_perturb, x, z[z_plt:end], Δb_f_xz,
                           colorrange = lims_Δb, colormap = :balance)
   hm_w_perturb = heatmap!(ax_w_perturb, x, z[z_plt:end], Δw_f_xz,
                           colorrange = lims_Δw, colormap = :balance)
   hm_u_perturb = heatmap!(ax_u_perturb, x, z[z_plt:end], Δu_f_xz,
                           colorrange = lims_Δu, colormap = :balance)
   hm_v_perturb = heatmap!(ax_v_perturb, x, z[z_plt:end], Δv_f_xz,
                           colorrange = lims_Δv, colormap = :balance)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_u_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_v_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_u_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_v_perturb, tickformat = "{:.1e}",
            label = "m/s")

   title_total   = @sprintf("Fields at y = %i km; t = %.2f days",
                            y_nearest, times[Nt]/(3600*24))
   title_perturb = @sprintf(
                          "Perturbation fields at y = %i km; t = %.2f days",
                            y_nearest, times[Nt]/(3600*24))

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "fields_y$(y_nearest)_tf_$(datetime).png"),
                 fig_total)
   save(joinpath("./Plots", "perturbs_y$(y_nearest)_tf_$(datetime).png"),
                 fig_perturb)
   close(ds)
end

function visualize_fields_const_z(datetime, z_idx; 
		                  plot_animation = true, t_idx_skip = 1)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)
   bb, ub, vb, wb         = get_background_fields(ds)
   
   b_total_f_xy = ds["b"][:, :, z_idx, Nt]
   u_total_f_xy = ds["u"][:, :, z_idx, Nt]
   v_total_f_xy = ds["v"][:, :, z_idx, Nt]
   w_total_f_xy = ds["w"][:, :, z_idx, Nt]

   Δb_f_xy = b_total_f_xy .- bb[:, :, z_idx]
   Δu_f_xy = u_total_f_xy .- ub[:, :, z_idx]
   Δv_f_xy = v_total_f_xy .- vb[:, :, z_idx]
   Δw_f_xy = w_total_f_xy .- wb[:, :, z_idx]

   lims_b_total = get_range_lims(b_total_f_xy)
   lims_u_total = get_range_lims(u_total_f_xy; prescribed_max = 1e-16)
   lims_v_total = get_range_lims(v_total_f_xy; prescribed_max = 1e-16)
   lims_w_total = get_range_lims(w_total_f_xy; prescribed_max = 1e-16)

   lims_Δb = get_range_lims(Δb_f_xy; prescribed_max = 1e-16)
   lims_Δu = get_range_lims(Δu_f_xy; prescribed_max = 1e-16)
   lims_Δv = get_range_lims(Δv_f_xy; prescribed_max = 1e-16)
   lims_Δw = get_range_lims(Δw_f_xy; prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   depth_nearest, axis_kwargs_xy = get_2D_spatial_axis_kwargs(x, y, z; z_idx = z_idx)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      n = Observable(1)

      b_total_xy = @lift ds["b"][:, :, z_idx, $n]
      u_total_xy = @lift ds["u"][:, :, z_idx, $n]
      v_total_xy = @lift ds["v"][:, :, z_idx, $n]
      w_total_xy = @lift ds["w"][:, :, z_idx, $n]

      Δb_xy = @lift $b_total_xy .- bb[:, :, z_idx]
      Δu_xy = @lift $u_total_xy .- ub[:, :, z_idx]
      Δv_xy = @lift $v_total_xy .- vb[:, :, z_idx]
      Δw_xy = @lift $w_total_xy .- wb[:, :, z_idx]
   
      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total = Axis(fig_total[2, 1];
                        title = "Total buoyancy (b)", axis_kwargs_xy...)
      ax_w_total = Axis(fig_total[2, 3];
                        title = "Total vertical velocity (w)", 
			axis_kwargs_xy...)
      ax_u_total = Axis(fig_total[3, 1];
                        title = "Total zonal velocity (u)", axis_kwargs_xy...)
      ax_v_total = Axis(fig_total[3, 3];
                        title = "Total meridional velocity (v)", 
			axis_kwargs_xy...)

      ax_b_perturb = Axis(fig_perturb[2, 1];
                          title = "Buoyancy perturbation (b')", 
			  axis_kwargs_xy...)
      ax_w_perturb = Axis(fig_perturb[2, 3];
                          title = "Vertical velocity perturbation (w')",
                          axis_kwargs_xy...)
      ax_u_perturb = Axis(fig_perturb[3, 1];
                          title = "Zonal velocity perturbation (u')",
                          axis_kwargs_xy...)
      ax_v_perturb = Axis(fig_perturb[3, 3];
                          title = "Meridional velocity perturbation (v')",
                          axis_kwargs_xy...)

      hm_b_total = heatmap!(ax_b_total, x, y, b_total_xy,
                            colorrange = lims_b_total, colormap = :balance)
      hm_w_total = heatmap!(ax_w_total, x, y, w_total_xy,
                            colorrange = lims_w_total, colormap = :balance)
      hm_u_total = heatmap!(ax_u_total, x, y, u_total_xy,
                            colorrange = lims_u_total, colormap = :balance)
      hm_v_total = heatmap!(ax_v_total, x, y, v_total_xy,
                            colorrange = lims_v_total, colormap = :balance)

      hm_b_perturb = heatmap!(ax_b_perturb, x, y, Δb_xy,
                              colorrange = lims_Δb, colormap = :balance)
      hm_w_perturb = heatmap!(ax_w_perturb, x, y, Δw_xy,
                              colorrange = lims_Δw, colormap = :balance)
      hm_u_perturb = heatmap!(ax_u_perturb, x, y, Δu_xy,
                              colorrange = lims_Δu, colormap = :balance)
      hm_v_perturb = heatmap!(ax_v_perturb, x, y, Δv_xy,
                              colorrange = lims_Δv, colormap = :balance)
   
      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 2], hm_u_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 4], hm_v_total, tickformat = "{:.1e}", label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
               label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}",
               label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_u_perturb, tickformat = "{:.1e}",
               label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_v_perturb, tickformat = "{:.1e}",
               label = "m/s")

      title_total = @lift @sprintf("Fields at %i-m depth; t = %.2f days",
                                   depth_nearest, times[$n]/(3600*24))
      title_perturb = @lift @sprintf(
                            "Perturbation fields at %i-m depth; t = %.2f days",
                                     depth_nearest, times[$n]/(3600*24))

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                                  tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                                  tellwidth = false)
   
      frames = 1:Nt
   
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
         recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end
   
      save(joinpath("./Plots", "fields_z-$(depth_nearest)_$(datetime).mp4"), 
	            video_total)
      save(joinpath("./Plots", "perturbs_z-$(depth_nearest)_$(datetime).mp4"),
	            video_perturb)
   end
   
   #Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total = Axis(fig_total[2, 1];
                     title = "Total buoyancy (b)", axis_kwargs_xy...)
   ax_w_total = Axis(fig_total[2, 3];
                     title = "Total vertical velocity (w)", 
		     axis_kwargs_xy...)
   ax_u_total = Axis(fig_total[3, 1];
                     title = "Total zonal velocity (u)", axis_kwargs_xy...)
   ax_v_total = Axis(fig_total[3, 3];
                     title = "Total meridional velocity (v)",
                     axis_kwargs_xy...)

   ax_b_perturb = Axis(fig_perturb[2, 1];
                       title = "Buoyancy perturbation (b')",
                       axis_kwargs_xy...)
   ax_w_perturb = Axis(fig_perturb[2, 3];
                       title = "Vertical velocity perturbation (w')",
                       axis_kwargs_xy...)
   ax_u_perturb = Axis(fig_perturb[3, 1];
                       title = "Zonal velocity perturbation (u')",
                       axis_kwargs_xy...)
   ax_v_perturb = Axis(fig_perturb[3, 3];
                       title = "Meridional velocity perturbation (v')",
                       axis_kwargs_xy...)

   hm_b_total = heatmap!(ax_b_total, x, y, b_total_f_xy,
                         colorrange = lims_b_total, colormap = :balance)
   hm_w_total = heatmap!(ax_w_total, x, y, w_total_f_xy,
                         colorrange = lims_w_total, colormap = :balance)
   hm_u_total = heatmap!(ax_u_total, x, y, u_total_f_xy,
                         colorrange = lims_u_total, colormap = :balance)
   hm_v_total = heatmap!(ax_v_total, x, y, v_total_f_xy,
                         colorrange = lims_v_total, colormap = :balance)

   hm_b_perturb = heatmap!(ax_b_perturb, x, y, Δb_f_xy,
                           colorrange = lims_Δb, colormap = :balance)
   hm_w_perturb = heatmap!(ax_w_perturb, x, y, Δw_f_xy,
                           colorrange = lims_Δw, colormap = :balance)
   hm_u_perturb = heatmap!(ax_u_perturb, x, y, Δu_f_xy,
                           colorrange = lims_Δu, colormap = :balance)
   hm_v_perturb = heatmap!(ax_v_perturb, x, y, Δv_f_xy,
                           colorrange = lims_Δv, colormap = :balance)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_w_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_u_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_v_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_w_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_u_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_v_perturb, tickformat = "{:.1e}",
            label = "m/s")

   title_total   = @sprintf("Fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt]/(3600*24))
   title_perturb = @sprintf(
                          "Perturbation fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt]/(3600*24))

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "fields_z-$(depth_nearest)_tf_$(datetime).png"),
        fig_total)
   save(joinpath("./Plots", "perturbs_z-$(depth_nearest)_tf_$(datetime).png"),
        fig_perturb)
   close(ds)
end

function open_computed_dataset(datetime, Δx, Δy, Δz, f)

   computed_file = joinpath("./Output", "computed_$(datetime).nc")

   #Only do computations if file does not already exist
   if !isfile(computed_file)

      ds, x, y, z, times, Nt = open_dataset(datetime)

      frames = 1:Nt
      x_idcs = 2:length(x)-2
      y_idcs = 2:length(y)-2
      z_idcs = 2:length(z)-2

      function update_data_array!(data_array, i, j, k, n, value)
         data_array[i, j, k, n] = value
         return data_array
      end
      
      NCDataset(computed_file, "c") do comp_ds
         
         i, j, k = Observable(2), Observable(2), Observable(2)
         n       = Observable(1)

         b  = @lift ds["b"][:, :, :, $n]
         u  = @lift ds["u"][1:end-1, :, :, $n]
         v  = @lift ds["v"][:, 1:end-1, :, $n]
         w  = @lift ds["w"][:, :, 1:end-1, $n]
         qn = @lift q($u, $v, $w, $b, f, $i, $j, $k, Δx, Δy, Δz)

	 defDim(comp_ds, "x", length(x)-2)
	 defDim(comp_ds, "y", length(y)-2)
	 defDim(comp_ds, "z", length(z)-2)
	 defDim(comp_ds, "time", length(times))

	 q_data  = Array{Float64, 4}(undef, 
				     length(x)-2, 
				     length(y)-2, 
				     length(z)-2, 
				     Nt)

         for t = 1:frames[end]
	    for z_idx = 2:z_idcs[end]
       	       for x_idx = 2:x_idcs[end]
                  for y_idx = 2:y_idcs[end]
	             update_data_array!(q_data, 
					x_idx, y_idx, z_idx, t, 
					to_value(qn))
		     yield()
		     j[] = y_idx
	          end
	          i[] = x_idx
               end
	       k[] = z_idx
	    end
	    print("Computing q for time $(t) of $(Nt)" * " \r")
            n[] = t
	 end
	 defVar(comp_ds, "q", q_data, ("x", "y", "z", "time"))
      end #comp_ds gets closed automatically
   end
   return NCDataset(computed_file, "r")
end

function visualize_q_const_x(datetime, Δx, Δy, Δz, f, x_idx)

   ds, x, y, z, times, Nt = open_dataset(datetime)
   
   comp_ds = open_computed_dataset(datetime, Δx, Δy, Δz, f)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n    = Observable(1)
   b    = @lift ds["b"][:, :, z_plt:end, $n]
   u    = @lift ds["u"][:, :, z_plt:end, $n]
   v    = @lift ds["v"][:, :, z_plt:end, $n]
   w    = @lift ds["w"][:, :, z_plt:end-1, $n]
   q_yz = @lift comp_ds["q"][x_idx, :, z_plt:end, $n] 
   
   q_yz_f = comp_ds["q"][x_idx, :, z_plt:end, Nt]

   lims_q = get_range_lims(q_yz_f)

   x_nearest, axis_kwargs_yz = get_2D_spatial_axis_kwargs(x, y, z; x_idx = x_idx)

   fig_q = Figure(size = (600, 600))
   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_yz...)
   hm_q  = heatmap!(ax_q, y[2:end-1], z[z_plt:end], q_yz, colorrange = lims_q,
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title_q       = @lift @sprintf("q at x = %i km; t = %.2f days",
			           x_nearest, times[$n]/(3600*24))
   fig_q[1, 1:2] = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   frames  = 1:Nt
   video_q = VideoStream(fig_q, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_x$(x_nearest)_$(datetime).mp4"), video_q)
   close(ds)
end

function visualize_q_const_y(datetime, Δx, Δy, Δz, f, y_idx)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)
   
   comp_ds = open_computed_dataset(datetime, Δx, Δy, Δz, f)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n    = Observable(1)
   b    = @lift ds["b"][:, :, z_plt:end, $n]
   u    = @lift ds["u"][:, :, z_plt:end, $n]
   v    = @lift ds["v"][:, :, z_plt:end, $n]
   w    = @lift ds["w"][:, :, z_plt:end-1, $n]
   q_xz = @lift comp_ds["q"][:, y_idx, z_plt:end, $n]

   q_xz_f = comp_ds["q"][:, y_idx, z_plt:end, Nt]

   lims_q = get_range_lims(q_xz_f)

   y_nearest, axis_kwargs_xz = get_2D_spatial_axis_kwargs(x, y, z; y_idx = y_idx)

   fig_q = Figure(size = (600, 600))
   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_xz...)
   hm_q  = heatmap!(ax_q, x[2:end-1], z[z_plt:end], q_xz, colorrange = lims_q, 
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title_q        = @lift @sprintf("q at y = %i km; t = %.2f days",
			            y_nearest, times[$n]/(3600*24))
   fig_q[1, 1:2]  = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   frames  = 1:Nt
   video_q = VideoStream(fig_q, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_y$(y_nearest)_$(datetime).mp4"), video_q)
   close(ds)
end

function visualize_q_const_z(datetime, Δx, Δy, Δz, f, z_idx)

   ds, x, y, z, times, Nt = open_dataset(datetime)
   
   comp_ds = open_computed_dataset(datetime, Δx, Δy, Δz, f)

   n    = Observable(1)
   b    = @lift ds["b"][:, :, :, $n]
   u    = @lift ds["u"][:, :, :, $n]
   v    = @lift ds["v"][:, :, :, $n]
   w    = @lift ds["w"][:, :, 1:end-1, $n]
   q_xy = @lift comp_ds["q"][:, :, z_idx, $n]

   q_xy_f = comp_ds["q"][:, :, z_idx, Nt]

   lims_q = get_range_lims(q_xy_f)

   depth_nearest, axis_kwargs_xy = get_2D_spatial_axis_kwargs(x, y, z; z_idx = z_idx)

   fig_q = Figure(size = (600, 600))
   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_xy...)
   hm_q  = heatmap!(ax_q, x[2:end-1], y[2:end-1], q_xy, colorrange = lims_q,
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title_q       = @lift @sprintf("q at %i-m depth; t = %.2f days",
			     depth_nearest, times[$n]/(3600*24))
   fig_q[1, 1:2] = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   frames   = 1:Nt
   video_q  = VideoStream(fig_q, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_z-$(depth_nearest)_$(datetime).mp4"), video_q)
   close(ds)
end

function plot_background_ζa(datetime, U, f, σr, σz; 
		            x_idx = nothing, y_idx = nothing)

   ds, x, y, z, times, Nt = open_dataset(datetime)

   nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z; 
					             x_idx = x_idx, y_idx = y_idx)

   fig = Figure(size = (800, 400))
   ax  = Axis(fig[2, 1]; axis_kwargs...)

   if !isnothing(x_idx)
      ζa_b_yz = ζa_b(U, f, σr, σz, x[x_idx], y, z) 
      lims_ζa = get_range_lims(ζa_b_yz)
      hm      = heatmap!(ax, y, z, ζa_b_yz, colorrange = lims_ζa, 
			 colormap = :balance)
      title   = ("Absolute vorticity of background state at x = $(nearest) km")
      fname   = "bkgd_zeta_abs_x$(nearest)_$(datetime).png"
   elseif !isnothing(y_idx)
      ζa_b_xz = ζa_b(U, f, σr, σz, x, y[y_idx], z)
      lims_ζa = get_range_lims(ζa_b_xz)
      hm      = heatmap!(ax, x, z, ζa_b_xz, colorrange = lims_ζa,
                         colormap = :balance)
      title   = ("Absolute vorticity of background state at y = $(nearest) km")
      fname   = "bkgd_zeta_abs_y$(nearest)_$(datetime).png"
   end

   Colorbar(fig[2, 2], hm, tickformat = "{:.1e}", label = "1/s")

   fig[1, 1:2] = Label(fig, title, fontsize = 24, tellwidth = false)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", fname), fig)
   close(ds)
end
