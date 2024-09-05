include("Library.jl")

using Oceananigans
using CairoMakie, DataStructures, NCDatasets, Printf
using .ComputeSecondaries

function open_dataset(datetime)

   outfilepath = joinpath("./Output", "output_$(datetime).nc")
   
   ds = NCDataset(outfilepath, "r")
   x  = ds["xC"][:]
   y  = ds["yC"][:]
   z  = ds["zC"][:]
   t  = ds["time"][:]
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
   field_lims = (3/4) * [-field_max, field_max]
end

function get_axis_kwargs(x, y, z; 
		         x_idx = nothing, y_idx = nothing, z_idx = nothing)
   if !isnothing(x_idx)
      nearest_m   = round(Int, x[x_idx])
      axis_kwargs = (xlabel = "y [m]", ylabel = "z [m]")
   elseif !isnothing(y_idx)
      nearest_m   = round(Int, y[y_idx])
      axis_kwargs = (xlabel = "x [m]", ylabel = "z [m]")
   elseif !isnothing(z_idx)
      nearest_m   = round(Int, -z[z_idx])
      axis_kwargs = (xlabel = "x [m]", ylabel = "y [m]")
   end
   return nearest_m, axis_kwargs
end

function visualize_perturbs_const_x(datetime, x_idx)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n = Observable(1)

   bb, ub, vb, wb = get_background_fields(ds)

   b_yz = @lift ds["b"][x_idx, :, z_plt:end, $n] .- bb[x_idx, :, z_plt:end]
   u_yz = @lift ds["u"][x_idx, :, z_plt:end, $n] .- ub[x_idx, :, z_plt:end]
   v_yz = @lift ds["v"][x_idx, :, z_plt:end, $n] .- vb[x_idx, :, z_plt:end]
   w_yz = @lift ds["w"][x_idx, :, z_plt:end-1, $n] .- wb[x_idx, :, z_plt:end]

   bf_yz = ds["b"][x_idx, :, z_plt:end, Nt] .- bb[x_idx, :, z_plt:end]
   uf_yz = ds["u"][x_idx, :, z_plt:end, Nt] .- ub[x_idx, :, z_plt:end]
   vf_yz = ds["v"][x_idx, :, z_plt:end, Nt] .- vb[x_idx, :, z_plt:end]
   wf_yz = ds["w"][x_idx, :, z_plt:end-1, Nt] .- wb[x_idx, :, z_plt:end]

   lims_b = get_range_lims(bf_yz)
   lims_u = get_range_lims(uf_yz)
   lims_v = get_range_lims(vf_yz)
   lims_w = get_range_lims(wf_yz)

   x_nearest_m, axis_kwargs_yz = get_axis_kwargs(x, y, z; x_idx = x_idx)

   fig  = Figure(size = (1200, 800))
   
   ax_b = Axis(fig[2, 1];
               title = "Buoyancy perturbation (b')", axis_kwargs_yz...)
   ax_w = Axis(fig[2, 3];
               title = "Vertical velocity perturbation (w')", axis_kwargs_yz...)
   ax_u = Axis(fig[3, 1];
               title = "Zonal velocity perturbation (u')", axis_kwargs_yz...)
   ax_v = Axis(fig[3, 3];
               title = "Meridional velocity perturbation (v')", 
	       axis_kwargs_yz...)

   hm_b = heatmap!(ax_b, y, z[z_plt:end], b_yz, 
		   colorrange = lims_b, colormap = :balance)
   hm_w = heatmap!(ax_w, y, z[z_plt:end], w_yz, 
                   colorrange = lims_w, colormap = :balance)
   hm_u = heatmap!(ax_u, y, z[z_plt:end], u_yz, 
                   colorrange = lims_u, colormap = :balance)
   hm_v = heatmap!(ax_v, y, z[z_plt:end], v_yz, 
                   colorrange = lims_v, colormap = :balance)

   Colorbar(fig[2, 2], hm_b, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig[2, 4], hm_w, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig[3, 2], hm_u, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig[3, 4], hm_v, tickformat = "{:.1e}", label = "m/s")

   title = @lift @sprintf("Perturbation fields at x = %i m; t = %.2f days",
                          x_nearest_m, times[$n]/(3600*24))

   fig[1, 1:4] = Label(fig, title, fontsize = 24, tellwidth = false)

   frames = 1:Nt
   video  = VideoStream(fig, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "bwuv_x$(x_nearest_m)_$(datetime).mp4"), video)
   close(ds)
end

function visualize_perturbs_const_y(datetime, y_idx)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n = Observable(1)

   bb, ub, vb, wb = get_background_fields(ds)

   b_xz = @lift ds["b"][:, y_idx, z_plt:end, $n] .- bb[:, y_idx, z_plt:end]
   u_xz = @lift ds["u"][:, y_idx, z_plt:end, $n] .- ub[:, y_idx, z_plt:end]
   v_xz = @lift ds["v"][:, y_idx, z_plt:end, $n] .- vb[:, y_idx, z_plt:end]
   w_xz = @lift ds["w"][:, y_idx, z_plt:end-1, $n] .- wb[:, y_idx, z_plt:end]

   bf_xz = ds["b"][:, y_idx, z_plt:end, Nt] .- bb[:, y_idx, z_plt:end]
   uf_xz = ds["u"][:, y_idx, z_plt:end, Nt] .- ub[:, y_idx, z_plt:end]
   vf_xz = ds["v"][:, y_idx, z_plt:end, Nt] .- vb[:, y_idx, z_plt:end]
   wf_xz = ds["w"][:, y_idx, z_plt:end-1, Nt] .- wb[:, y_idx, z_plt:end]

   lims_b = get_range_lims(bf_xz)
   lims_u = get_range_lims(uf_xz)
   lims_v = get_range_lims(vf_xz)
   lims_w = get_range_lims(wf_xz)

   y_nearest_m, axis_kwargs_xz = get_axis_kwargs(x, y, z; y_idx = y_idx)

   fig  = Figure(size = (1200, 800))
   
   ax_b = Axis(fig[2, 1];
               title = "Buoyancy perturbation (b')", axis_kwargs_xz...)
   ax_w = Axis(fig[2, 3];
               title = "Vertical velocity perturbation (w')", axis_kwargs_xz...)
   ax_u = Axis(fig[3, 1];
               title = "Zonal velocity perturbation (u')", axis_kwargs_xz...)
   ax_v = Axis(fig[3, 3];
               title = "Meridional velocity perturbation (v')", 
	       axis_kwargs_xz...)

   hm_b = heatmap!(ax_b, x, z[z_plt:end], b_xz, 
		   colorrange = lims_b, colormap = :balance)
   hm_w = heatmap!(ax_w, x, z[z_plt:end], w_xz, 
                   colorrange = lims_w, colormap = :balance)
   hm_u = heatmap!(ax_u, x, z[z_plt:end], u_xz, 
                   colorrange = lims_u, colormap = :balance)
   hm_v = heatmap!(ax_v, x, z[z_plt:end], v_xz, 
                   colorrange = lims_v, colormap = :balance)

   Colorbar(fig[2, 2], hm_b, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig[2, 4], hm_w, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig[3, 2], hm_u, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig[3, 4], hm_v, tickformat = "{:.1e}", label = "m/s")

   title = @lift @sprintf("Perturbation fields at y = %i m; t = %.2f days",
                          y_nearest_m, times[$n]/(3600*24))

   fig[1, 1:4] = Label(fig, title, fontsize = 24, tellwidth = false)

   frames = 1:Nt
   video  = VideoStream(fig, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "bwuv_y$(y_nearest_m)_$(datetime).mp4"), video)
   close(ds)
end

function visualize_perturbs_const_z(datetime, z_idx)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)

   n = Observable(1)

   bb, ub, vb, wb = get_background_fields(ds)

   b_xy = @lift ds["b"][:, :, z_idx, $n] .- bb[:, :, z_idx]
   u_xy = @lift ds["u"][:, :, z_idx, $n] .- ub[:, :, z_idx]
   v_xy = @lift ds["v"][:, :, z_idx, $n] .- vb[:, :, z_idx]
   w_xy = @lift ds["w"][:, :, z_idx, $n] .- wb[:, :, z_idx]

   bf_xy = ds["b"][:, :, z_idx, Nt] .- bb[:, :, z_idx]
   uf_xy = ds["u"][:, :, z_idx, Nt] .- ub[:, :, z_idx]
   wf_xy = ds["w"][:, :, z_idx, Nt] .- wb[:, :, z_idx]

   lims_b  = get_range_lims(bf_xy)
   lims_uv = get_range_lims(uf_xy)
   lims_w  = get_range_lims(wf_xy)

   depth_nearest_m, axis_kwargs_xy = get_axis_kwargs(x, y, z; z_idx = z_idx)

   fig  = Figure(size = (1200, 800))
   
   ax_b = Axis(fig[2, 1];
               title = "Buoyancy perturbation (b')", axis_kwargs_xy...)
   ax_w = Axis(fig[2, 3];
               title = "Vertical velocity perturbation (w')", axis_kwargs_xy...)
   ax_u = Axis(fig[3, 1];
               title = "Zonal velocity perturbation (u')", axis_kwargs_xy...)
   ax_v = Axis(fig[3, 3];
               title = "Meridional velocity perturbation (v')", 
	       axis_kwargs_xy...)

   hm_b = heatmap!(ax_b, x, y, b_xy, colorrange = lims_b, colormap = :balance)
   hm_w = heatmap!(ax_w, x, y, w_xy, colorrange = lims_w, colormap = :balance)
   hm_u = heatmap!(ax_u, x, y, u_xy, colorrange = lims_uv, colormap = :balance)
   hm_v = heatmap!(ax_v, x, y, v_xy, colorrange = lims_uv, colormap = :balance)

   Colorbar(fig[2, 2], hm_b, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig[2, 4], hm_w, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig[3, 2], hm_u, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig[3, 4], hm_v, tickformat = "{:.1e}", label = "m/s")

   title = @lift @sprintf("Perturbation fields at depth %i m; t = %.2f days",
                          depth_nearest_m, times[$n]/(3600*24))

   fig[1, 1:4] = Label(fig, title, fontsize = 24, tellwidth = false)

   frames = 1:Nt
   video  = VideoStream(fig, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "bwuv_z$(depth_nearest_m)_$(datetime).mp4"), video)
   close(ds)
end

function open_computed_dataset(datetime, Δx, Δy, Δz, f)

   computed_file = joinpath("./Output", "computed_$(datetime).nc")

   #Only do computations if file does not already exist
   if !isfile(computed_file)

      ds, x, y, z, times, Nt = open_dataset(datetime)

      n = Observable(1)

      b  = @lift ds["b"][:, :, :, $n]
      u  = @lift ds["u"][:, :, :, $n]
      v  = @lift ds["v"][:, :, :, $n]
      w  = @lift ds["w"][:, :, 1:end-1, $n]
      q  = @lift ertelQ($u, $v, $w, $b, f, Δx, Δy, Δz)
      qr = @lift ∂r_ertelQ($q, Δx, Δy, x[2:end-1], y[2:end-1])

      comp_ds = NCDataset(computed_file, "c")
      defVar(comp_ds, "xG", x[1:end-1], ("xG",))
      defVar(comp_ds, "yG", y[1:end-1], ("yG",))
      defVar(comp_ds, "zG", z[1:end-1], ("zG",))
      defVar(comp_ds, "xCC", x[2:end-1], ("xCC",))
      defVar(comp_ds, "yCC", y[2:end-1], ("yCC",))
      defVar(comp_ds, "time", times[:], ("time",))

      q_data  = defVar(comp_ds, "q", Float64, 
		       ("xG", "yG", "zG", "time"))
      qr_data = defVar(comp_ds, "qr", Float64, 
		       ("xCC", "yCC", "zG", "time"))

      frames = 1:Nt

      for i = 1:frames[end]
	 q_data[:, :, :, i]  = to_value(q)
	 qr_data[:, :, :, i] = to_value(qr)
         msg = string("Computing frame(s) ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
      end

      close(ds)
      close(comp_ds)
   end

   comp_ds = NCDataset(computed_file, "r")
   xG  = comp_ds["xG"][:]
   yG  = comp_ds["yG"][:]
   zG  = comp_ds["zG"][:]
   xCC = comp_ds["xCC"][:]
   yCC = comp_ds["yCC"][:]

   return comp_ds, xG, yG, zG, xCC, yCC
end

function visualize_q_const_x(datetime, Δx, Δy, Δz, f, x_idx)

   ds, x, y, z, times, Nt = open_dataset(datetime)
   comp_ds, xG, yG, zG, xCC, yCC = open_computed_dataset(datetime, 
							 Δx, Δy, Δz, f)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n = Observable(1)

   b     = @lift ds["b"][:, :, z_plt:end, $n]
   u     = @lift ds["u"][:, :, z_plt:end, $n]
   v     = @lift ds["v"][:, :, z_plt:end, $n]
   w     = @lift ds["w"][:, :, z_plt:end-1, $n]
   q_yz  = @lift comp_ds["q"][x_idx, :, z_plt:end, $n] 
   qr_yz = @lift comp_ds["qr"][x_idx, :, z_plt:end, $n]
   
   q_yz_f  = comp_ds["q"][x_idx, :, z_plt:end, Nt]
   qr_yz_f = comp_ds["qr"][x_idx, :, z_plt:end, Nt]

   lims_q  = get_range_lims(q_yz_f)
   lims_qr = get_range_lims(qr_yz_f)

   x_nearest_m, axis_kwargs_yz = get_axis_kwargs(x, y, z; x_idx = x_idx)

   fig_q  = Figure(size = (600, 600))
   fig_qr = Figure(size = (600, 600))

   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_yz...)
   ax_qr = Axis(fig_qr[2, 1]; axis_kwargs_yz...)

   hm_q  = heatmap!(ax_q, yG, zG[z_plt:end], q_yz, colorrange = lims_q,
		    colormap = :balance)
   hm_qr = heatmap!(ax_qr, yCC, zG[z_plt:end], qr_yz, colorrange = lims_qr,
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")
   Colorbar(fig_qr[2, 2], hm_qr, tickformat = "{:.1e}", label = "1/s³m")

   title_q  = @lift @sprintf("q at x = %i m; t = %.2f days",
			     x_nearest_m, times[$n]/(3600*24))
   title_qr = @lift @sprintf("∂q/∂r at x = %i m; t = %.2f days",
			     x_nearest_m, times[$n]/(3600*24))

   fig_q[1, 1:2]  = Label(fig_q, title_q, fontsize = 24, tellwidth = false)
   fig_qr[1, 1:2] = Label(fig_qr, title_qr, fontsize = 24, tellwidth = false)

   frames   = 1:Nt
   video_q  = VideoStream(fig_q, format = "mp4", framerate = 6)
   video_qr = VideoStream(fig_qr, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      recordframe!(video_qr)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_x$(x_nearest_m)_$(datetime).mp4"), video_q)
   save(joinpath("./Plots", "qr_x$(x_nearest_m)_$(datetime).mp4"), video_qr)
   close(ds)
end

function visualize_q_const_y(datetime, Δx, Δy, Δz, f, y_idx)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)
   comp_ds, xG, yG, zG, xCC, yCC = open_computed_dataset(datetime,
                                                         Δx, Δy, Δz, f)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n = Observable(1)

   b     = @lift ds["b"][:, :, z_plt:end, $n]
   u     = @lift ds["u"][:, :, z_plt:end, $n]
   v     = @lift ds["v"][:, :, z_plt:end, $n]
   w     = @lift ds["w"][:, :, z_plt:end-1, $n]
   q_xz  = @lift comp_ds["q"][:, y_idx, z_plt:end, $n]
   qr_xz = @lift comp_ds["qr"][:, y_idx, z_plt:end, $n]

   q_xz_f  = comp_ds["q"][:, y_idx, z_plt:end, Nt]
   qr_xz_f = comp_ds["qr"][:, y_idx, z_plt:end, Nt]

   lims_q  = get_range_lims(q_xz_f)
   lims_qr = get_range_lims(qr_xz_f)

   y_nearest_m, axis_kwargs_xz = get_axis_kwargs(x, y, z; y_idx = y_idx)

   fig_q  = Figure(size = (600, 600))
   fig_qr = Figure(size = (600, 600))

   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_xz...)
   ax_qr = Axis(fig_qr[2, 1]; axis_kwargs_xz...)

   hm_q  = heatmap!(ax_q, xG, zG[z_plt:end], q_xz, colorrange = lims_q, 
		    colormap = :balance)
   hm_qr = heatmap!(ax_qr, xCC, zG[z_plt:end], qr_xz, colorrange = lims_qr,
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")
   Colorbar(fig_qr[2, 2], hm_qr, tickformat = "{:.1e}", label = "1/s³m")

   title_q  = @lift @sprintf("q at y = %i m; t = %.2f days",
			     y_nearest_m, times[$n]/(3600*24))
   title_qr = @lift @sprintf("∂q/∂r at y = %i m; t = %.2f days",
			     y_nearest_m, times[$n]/(3600*24))

   fig_q[1, 1:2]  = Label(fig_q, title_q, fontsize = 24, tellwidth = false)
   fig_qr[1, 1:2] = Label(fig_qr, title_qr, fontsize = 24, tellwidth = false)

   frames   = 1:Nt
   video_q  = VideoStream(fig_q, format = "mp4", framerate = 6)
   video_qr = VideoStream(fig_qr, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      recordframe!(video_qr)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_y$(y_nearest_m)_$(datetime).mp4"), video_q)
   save(joinpath("./Plots", "qr_y$(y_nearest_m)_$(datetime).mp4"), video_qr)
   close(ds)
end

function visualize_q_const_z(datetime, Δx, Δy, Δz, f, z_idx)

   ds, x, y, z, times, Nt = open_dataset(datetime)
   comp_ds, xG, yG, zG, xCC, yCC = open_computed_dataset(datetime,
                                                         Δx, Δy, Δz, f)

   n = Observable(1)

   b     = @lift ds["b"][:, :, :, $n]
   u     = @lift ds["u"][:, :, :, $n]
   v     = @lift ds["v"][:, :, :, $n]
   w     = @lift ds["w"][:, :, 1:end-1, $n]
   q_xy  = @lift comp_ds["q"][:, :, z_idx, $n]
   qr_xy = @lift comp_ds["qr"][:, :, z_idx, $n]

   q_xy_f  = comp_ds["q"][:, :, z_idx, Nt]
   qr_xy_f = comp_ds["qr"][:, :, z_idx, Nt]

   lims_q = get_range_lims(q_xy_f)
   lims_qr = get_range_lims(qr_xy_f)

   depth_nearest_m, axis_kwargs_xy = get_axis_kwargs(x, y, z; z_idx = z_idx)

   fig_q  = Figure(size = (600, 600))
   fig_qr = Figure(size = (600, 600))

   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_xy...)
   ax_qr = Axis(fig_qr[2, 1]; axis_kwargs_xy...)
   
   hm_q  = heatmap!(ax_q, xG, yG, q_xy, colorrange = lims_q, 
		    colormap = :balance)
   hm_qr = heatmap!(ax_qr, xCC, yCC, qr_xy, colorrange = lims_qr, 
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")
   Colorbar(fig_qr[2, 2], hm_qr, tickformat = "{:.1e}", label = "1/s³m")

   title_q  = @lift @sprintf("q at depth %i m; t = %.2f days",
			     depth_nearest_m, times[$n]/(3600*24))
   title_qr = @lift @sprintf("∂q/∂r at depth %i m; t = %.2f days",
                             depth_nearest_m, times[$n]/(3600*24))

   fig_q[1, 1:2]  = Label(fig_q, title_q, fontsize = 24, tellwidth = false)
   fig_qr[1, 1:2] = Label(fig_qr, title_qr, fontsize = 24, tellwidth = false) 

   frames   = 1:Nt
   video_q  = VideoStream(fig_q, format = "mp4", framerate = 6)
   video_qr = VideoStream(fig_qr, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      recordframe!(video_qr)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_z$(depth_nearest_m)_$(datetime).mp4"), video_q)
   save(joinpath("./Plots", "qr_z$(depth_nearest_m)_$(datetime).mp4"), video_qr)
   close(ds)
end
