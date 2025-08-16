include("LibraryVisualization.jl")

using Adapt, CairoMakie 
using CommonDataModel, CUDA, DataStructures, Glob, LaTeXStrings, NCDatasets

update_theme!(fontsize = 16)

using Oceananigans
using Oceananigans.Fields
using Oceananigans.OutputReaders
using OffsetArrays: no_offset_view
using Printf

####################

function visualize_norms(datetime)

   scalars_ds = open_scalars_dataset("scalars_$(datetime).nc")

   b′_norm  = scalars_ds[:b′_norm][:]
   print(b′_norm)
   ux′_norm = scalars_ds[:ux′_norm][:]
   uy′_norm = scalars_ds[:uy′_norm][:]
   ur′_norm = scalars_ds[:ur′_norm][:]
   uφ′_norm = scalars_ds[:uφ′_norm][:]
   uz′_norm = scalars_ds[:uz′_norm][:]
   
   times = scalars_ds[:time][:] ./ 86400 #Convert to days for readability

   fig_cyl   = Figure(size = (1200, 700))
   ax_b_cyl  = Axis(fig_cyl[2, 1]; title = L"Norm of $b'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]",
                    yscale = log10)
   ax_ur     = Axis(fig_cyl[2, 2]; title = L"Norm of $u_r'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_r'||$ [m/s]",
                    yscale = log10)
   ax_uφ     = Axis(fig_cyl[3, 1]; title = L"Norm of $u_{\phi}'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_{\phi}'||$ [m/s]",
                    yscale = log10)
   ax_uz_cyl = Axis(fig_cyl[3, 2]; title = L"Norm of $u_z'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_z'||$ [m/s]",
                    yscale = log10)

   fig_Cart   = Figure(size = (1200, 700))
   ax_b_Cart  = Axis(fig_Cart[2, 1]; title = L"Norm of $b'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]",
	    	     yscale = log10)
   ax_ux      = Axis(fig_Cart[2, 2]; title = L"Norm of $u_x'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||u_x'||$ [m/s]",
                     yscale = log10)
   ax_uy      = Axis(fig_Cart[3, 1]; title = L"Norm of $u_y'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||u_y'||$ [m/s]",
                     yscale = log10)
   ax_uz_Cart = Axis(fig_Cart[3, 2]; title = L"Norm of $u_z'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||u_z'||$ [m/s]",
                     yscale = log10)

   scatter!(ax_b_cyl, times, b′_norm, color = :black)
   scatter!(ax_ur, times, ur′_norm, color = :black)
   scatter!(ax_uφ, times, uφ′_norm, color = :black)
   scatter!(ax_uz_cyl, times, uz′_norm, color = :black)
   
   scatter!(ax_b_Cart, times, b′_norm, color = :black)
   scatter!(ax_ux, times, ur′_norm, color = :black)
   scatter!(ax_uy, times, uφ′_norm, color = :black)
   scatter!(ax_uz_Cart, times, uz′_norm, color = :black)

   mkpath("./Plots") #Make visualization directory if nonexistent

   fig_cyl[1, 1:2]  = Label(fig_cyl, "Norms of perturbation fields",
                           fontsize = 24, tellwidth = false)
   fig_Cart[1, 1:2] = Label(fig_Cart, "Norms of perturbation fields",
                           fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "norm_fields_$(datetime).png"), fig_cyl)
   save(joinpath("./Plots", "norm_Cart_fields_$(datetime).png"), fig_Cart)
   close(scalars_ds)
end

function compute_and_visualize_norms(datetime, grid; 
		                     bkgd_datetime = nothing, 
				     do_Cartesian = false)
   
   #`grid` is a required argument for now
   #but I would like to update the open_dataset function to get it automatically from the output file

   if isnothing(bkgd_datetime)
      bkgd_datetime = datetime
   end

   bkgd_ds = open_bkgd_dataset(bkgd_datetime)
   B       = bkgd_ds[:B][:, :, :, 1]
   Uφ      = bkgd_ds[:Uφ][:, :, :, 1]

   fig_cyl   = Figure(size = (1200, 700))
   ax_b_cyl  = Axis(fig_cyl[2, 1]; title = L"Norm of $b'$",
		    xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]", 
		    yscale = log10)
   ax_ur     = Axis(fig_cyl[2, 2]; title = L"Norm of $u_r'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_r'||$ [m/s]", 
		    yscale = log10)
   ax_uφ     = Axis(fig_cyl[3, 1]; title = L"Norm of $u_{\phi}'$",
		    xlabel = L"$t$ [days]", ylabel = L"$||u_{\phi}'||$ [m/s]", 
		    yscale = log10)
   ax_uz_cyl = Axis(fig_cyl[3, 2]; title = L"Norm of $u_z'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_z'||$ [m/s]",
                    yscale = log10)

   if do_Cartesian

      Ux, Uy = bkgd_ds[:Ux][:, :, :, 1], bkgd_ds[:Uy][:, :, :, 1]

      fig_Cart   = Figure(size = (1200, 700))
      ax_b_Cart  = Axis(fig_Cart[2, 1]; title = L"Norm of $b'$",
                        xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]",
                        yscale = log10)
      ax_ux      = Axis(fig_Cart[2, 2]; title = L"Norm of $u_x'$",
                        xlabel = L"$t$ [days]", ylabel = L"$||u_x'||$ [m/s]",
                        yscale = log10)
      ax_uy      = Axis(fig_Cart[3, 1]; title = L"Norm of $u_y'$",
                        xlabel = L"$t$ [days]", ylabel = L"$||u_y'||$ [m/s]",
                        yscale = log10)
      ax_uz_Cart = Axis(fig_Cart[3, 2]; title = L"Norm of $u_z'$",
                        xlabel = L"$t$ [days]", ylabel = L"$||u_z'||$ [m/s]",
                        yscale = log10)
   end

   outfile_list = glob("output_$(datetime)*", "./Output")

   file_idx = 1

   while file_idx < length(outfile_list)

      outfilename = outfile_list[file_idx]
      print(outfilename)
      ds, x, y, z, times, Nt = open_dataset(outfile_list)

      b      = ds[:b][:, :, :, :]
      ur, uφ = ds[:ur][:, :, :, :], ds[:uφ][:, :, :, :]
      uz     = ds[:uz][:, :, :, :]

      if do_Cartesian
         ux, uy = ds[:ux][:, :, :, :], ds[:uy][:, :, :, :]
      end

      n = Observable(2)

      b_norm  = @lift field_norm(b, $n; ψ_bkgd = B)
      ur_norm = @lift field_norm(ur, $n)
      uφ_norm = @lift field_norm(uφ, $n; ψ_bkgd = Uφ)
      uz_norm = @lift field_norm(uz, $n)

      if do_Cartesian
         ux_norm = @lift field_norm(ux, $n, ψ_bkgd = Ux)
	 uy_norm = @lift field_norm(uy, $n, ψ_bkgd = Uy)
      end

      for i = 2:Nt
      
         @lift scatter!(ax_b_cyl, times[$n], $b_norm, color = :black)
         @lift scatter!(ax_ur, times[$n], $ur_norm, color = :black)
         @lift scatter!(ax_uφ, times[$n], $uφ_norm, color = :black)
         @lift scatter!(ax_uz_cyl, times[$n], $uz_norm, color = :black)
   
         if do_Cartesian
            @lift scatter!(ax_b_Cart, times[$n], $b_norm, color = :black)
            @lift scatter!(ax_ux, times[$n], $ux_norm, color = :black)
            @lift scatter!(ax_uy, times[$n], $uy_norm, color = :black)
            @lift scatter!(ax_uz_Cart, times[$n], $uz_norm, color = :black)
         end
      
         yield()
         n[] = i
      end

      close(ds)
      file_idx += 1
   end 

   mkpath("./Plots") #Make visualization directory if nonexistent

   fig_cyl[1, 1:2] = Label(fig_cyl, "Norms of perturbation fields",
			   fontsize = 24, tellwidth = false)
   save(joinpath("./Plots", "norm_fields_$(datetime).png"), fig_cyl)

   if do_Cartesian
      fig_Cart[1, 1:2] = Label(fig_Cart, "Norms of perturbation fields",
			       fontsize = 24, tellwidth = false)
      save(joinpath("./Plots", "norms_Cartesian_$(datetime).png"), fig_Cart)
   end

   close(bkgd_ds)
end

function visualize_norms_poster(datetime, grid;
		bkgd_datetime = nothing)

   if isnothing(bkgd_datetime)
      bkgd_datetime = datetime
   end

   bkgd_ds = open_bkgd_dataset(bkgd_datetime)
   Uφ      = bkgd_ds[:Uφ][:, :, :, 1]

   fig = Figure(size = (500, 300))
   ax = Axis(fig[1, 1]; title = L"Norms of $u_r', u_{{\phi}}'$",
	                xlabel = L"$t$ [days]", ylabel = L"$\ell^2$-norm [m/s]",
                        yscale = log10)
   
   outfile_list = glob("output_$(datetime)*", "./Output")
   
   file_idx = 1

   while file_idx < 7
      
      outfilename = outfile_list[file_idx]

      ds, x, y, z, times, Nt = open_dataset(outfile_list)
      
      ur, uφ = ds[:ur][:, :, :, :], ds[:uφ][:, :, :, :]

      n = Observable(2)

      ur_norm = @lift field_norm(ur, $n)
      uφ_norm = @lift field_norm(uφ, $n; ψ_bkgd = Uφ)

      for i = 2:Nt

         @lift scatter!(ax, times[$n], $ur_norm, label = L"$u_r'$", color = colorant"#DB3E3E")
	 @lift scatter!(ax, times[$n], $uφ_norm, label = L"$u_{{\phi}}'$", color = colorant"#FAB146")

         yield()
         n[] = i
      end

      close(ds)
      file_idx += 1
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   axislegend(position = :rb, unique = true)
   save(joinpath("./Plots", "poster_version_norm_fields_$(datetime).png"), fig)
   close(bkgd_ds)
end

function visualize_b_and_ωz(datetime, Δx, Δy; 
		            x_idx = nothing, y_idx = nothing, z_idx = nothing,
		            bkgd_datetime = nothing, plot_animation = false, 
			    t_idx_skip = 1)

   if isnothing(bkgd_datetime)
      bkgd_datetime = datetime
   end

   bkgd_ds = open_bkgd_dataset(bkgd_datetime)
   B       = bkgd_ds[:B][:, :, :, 1]
   Ux, Uy  = bkgd_ds[:Ux][:, :, :, 1], bkgd_ds[:Uy][:, :, :, 1]

   outfile_list             = glob("./Output/output_$(datetime)*")
   ds_f, x, y, z, times, Nt = open_dataset(outfile_list[length(outfile_list)])

   if !isnothing(x_idx)

      B_slice  = @views B[x_idx, :, :]
      ωb_slice = @views ωz(Ux, Uy, Δx, Δy; x_idx = x_idx)

      b_total_f_slice = @views adapt(Array, ds_f[:b])[x_idx, :, :, Nt]
      ω_total_f_slice = @views ωz(ds_f[:ux][:, :, :, Nt], ds_f[:uy][:, :, :, Nt], Δx, Δy;
                                  x_idx = x_idx)

      idx_kwargs           = (x_idx = x_idx,)
      nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
                                                        x_idx = x_idx)

      h_dim, v_dim, const_dim, units = y, z, "x", "km"

   elseif !isnothing(y_idx)

      B_slice  = @views B[:, y_idx, :]
      ωb_slice = @views ωz(Ux, Uy, Δx, Δy; y_idx = y_idx)

      b_total_f_slice = @views adapt(Array, ds_f[:b])[:, y_idx, :, Nt]
      ω_total_f_slice = @views ωz(ds_f[:ux][:, :, :, Nt], ds_f[:uy][:, :, :, Nt], Δx, Δy;
                                  y_idx = y_idx)

      idx_kwargs           = (y_idx = y_idx,)
      nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
                                                        y_idx = y_idx)

      h_dim, v_dim, const_dim, units = x, z, "y", "km"

   elseif !isnothing(z_idx)

      B_slice  = @views B[:, :, z_idx]
      ωb_slice = @views ωz(Ux, Uy, Δx, Δy; z_idx = z_idx)

      b_total_f_slice = @views adapt(Array, ds_f[:b])[:, :, z_idx, Nt]
      ω_total_f_slice = ωz(ds_f[:ux][:, :, :, Nt], ds_f[:uy][:, :, :, Nt], Δx, Δy;
                           z_idx = z_idx)

      idx_kwargs           = (z_idx = z_idx,)
      nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
                                                        z_idx = z_idx)

      h_dim, v_dim, const_dim, units = x, y, "z", "m"
   end

   Δb_f_slice = b_total_f_slice .- B_slice
   Δω_f_slice = ω_total_f_slice .- ωb_slice

   lims_b_total = get_range_lims(b_total_f_slice; max_fraction = 0.75)
   lims_ω_total = get_range_lims(ω_total_f_slice; max_fraction = 0.75)

   lims_Δb = get_range_lims(Δb_f_slice;
                            max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δω = get_range_lims(Δω_f_slice;
                            max_fraction = 0.75, prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   #Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 500))
   fig_perturb = Figure(size = (1200, 500))

   ax_b_total = Axis(fig_total[2, 1];
                     title = "Total buoyancy (b)", axis_kwargs...)
   ax_ω_total = Axis(fig_total[2, 3];
                     title = "Total vertical vorticity (ζ)", axis_kwargs...)

   ax_b_perturb = Axis(fig_perturb[2, 1];
                       title = "Buoyancy perturbation (b')", axis_kwargs...)
   ax_ω_perturb = Axis(fig_perturb[2, 3];
                       title = "Vertical vorticity perturbation (ζ')",
                       axis_kwargs...)

   hm_b_total = heatmap!(ax_b_total, h_dim, v_dim, b_total_f_slice,
                         colorrange = lims_b_total,
                         colormap = Reverse(:RdBu_5),
                         highclip = :red, lowclip = :blue)
   hm_ω_total = heatmap!(ax_ω_total, h_dim, v_dim, ω_total_f_slice,
                         colorrange = lims_ω_total,
                         colormap = Reverse(:RdBu_5),
                         highclip = :red, lowclip = :blue)

   hm_b_perturb = heatmap!(ax_b_perturb, h_dim, v_dim, Δb_f_slice,
                           colorrange = lims_Δb,
                           colormap = Reverse(:RdBu_5),
                           highclip = :red, lowclip = :blue)
   hm_ω_perturb = heatmap!(ax_ω_perturb, h_dim, v_dim, Δω_f_slice,
                           colorrange = lims_Δω,
                           colormap = Reverse(:RdBu_5),
                           highclip = :red, lowclip = :blue)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ω_total, tickformat = "{:.1e}", label = "1/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ω_perturb, tickformat = "{:.1e}",
            label = "1/s")
   
   title_total   = @sprintf("Fields at %s = %i %s; t = %.2f days",
                            const_dim, nearest, units, times[Nt])
   title_perturb = @sprintf("Perturbation fields at %s = %i %s; t = %.2f days",
                            const_dim, nearest, units, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots",
                 "bzeta_total_$(const_dim)$(nearest)_tf_$(datetime).png"),
        fig_total)
   save(joinpath("./Plots",
                 "bzeta_perturbs_$(const_dim)$(nearest)_tf_$(datetime).png"),
        fig_perturb)
 
   close(ds_f)
   
   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      fig_total   = Figure(size = (1200, 500))
      fig_perturb = Figure(size = (1200, 500))

      ax_b_total = Axis(fig_total[2, 1];
                        title = "Total buoyancy (b)", axis_kwargs...)
      ax_ω_total = Axis(fig_total[2, 3];
                        title = "Total vertical vorticity (ζ)", axis_kwargs...)

      ax_b_perturb = Axis(fig_perturb[2, 1];
                          title = "Buoyancy perturbation (b')", axis_kwargs...)
      ax_ω_perturb = Axis(fig_perturb[2, 3];
                          title = "Vertical vorticity perturbation (ζ')",
                          axis_kwargs...)

      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      file_idx = 1

      while file_idx < length(outfile_list)
        
         outfilename                 = outfile_list[file_idx] 
	 ds, x, y, z, times, Nt_file = open_dataset(outfilename)

         ux = @views adapt(Array, ds[:ux])[:, :, :, :]
         uy = @views adapt(Array, ds[:uy])[:, :, :, :]

         if !isnothing(x_idx)
            b_total_tseries_slice = @views adapt(Array, ds[:b])[x_idx, :, :, :]
         elseif !isnothing(y_idx)
            b_total_tseries_slice = @views adapt(Array, ds[:b])[:, y_idx, :, :]
         elseif !isnothing(z_idx)
            b_total_tseries_slice = @views adapt(Array, ds[:b])[:, :, z_idx, :]
         end

         n = Observable(1)

         b_total_n_slice = @lift b_total_tseries_slice[:, :, $n]
         ω_total_n_slice = @lift ωz(ds[:ux][:, :, :, $n], ds[:uy][:, :, :, $n],
                                 Δx, Δy; idx_kwargs...)

         Δb_n_slice = @lift $b_total_n_slice .- B_slice
         Δω_n_slice = @lift $ω_total_n_slice .- ωb_slice

         hm_b_total = heatmap!(ax_b_total, h_dim, v_dim, b_total_n_slice,
                            colorrange = lims_b_total,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
         hm_ω_total = heatmap!(ax_ω_total, h_dim, v_dim, ω_total_n_slice,
                            colorrange = lims_ω_total,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)

         hm_b_perturb = heatmap!(ax_b_perturb, h_dim, v_dim, Δb_n_slice,
                              colorrange = lims_Δb,
                              colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)
         hm_ω_perturb = heatmap!(ax_ω_perturb, h_dim, v_dim, Δω_n_slice,
                              colorrange = lims_Δω,
                              colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)

         Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}",
               label = "m/s²")
         Colorbar(fig_total[2, 4], hm_ω_total, tickformat = "{:.1e}",
               label = "1/s")

         Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
               label = "m/s²")
         Colorbar(fig_perturb[2, 4], hm_ω_perturb, tickformat = "{:.1e}",
               label = "1/s")

         title_total = @lift @sprintf("Fields at %s = %i %s; t = %.2f days",
                                   const_dim, nearest, units, times[$n])
         title_perturb = @lift @sprintf(
                            "Perturbation fields at %s = %i %s; t = %.2f days",
                            const_dim, nearest, units, times[$n])

         fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                                  tellwidth = false)
         fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                                  tellwidth = false)
 
	 for i = 1:t_idx_skip:Nt_file
            recordframe!(video_total)
	    recordframe!(video_perturb)
	    yield()
	    n[] = i
         end

	 close(ds)
	 file_idx += 1
      end

      save(joinpath("./Plots",
            "bzeta_total_$(const_dim)$(nearest)_$(datetime).mp4"),
           video_total)
      save(joinpath("./Plots",
            "bzeta_perturbs_$(const_dim)$(nearest)_$(datetime).mp4"),
           video_perturb)
   end
end

function visualize_z_grid(datetime, grid, zmin; zmax = 0.0)

   mkpath("./Plots") #Make visualization directory if nonexistent

   zc = znodes(grid, Center())
   zf = znodes(grid, Face())
   Δz = zspacings(grid, Center())
   
   fig  = Figure(size=(1200, 600))
   axz  = Axis(fig[1, 1], title = "z-grid")
   axΔz = Axis(fig[2, 1]; xlabel = "z (m)", ylabel = "z-spacing (m)")

   lines!(axz, [zmin, zmax], [0, 0], color = :gray)
   scatter!(axz, zf, 0 * zf, marker = :vline, color = :gray, markersize = 20)
   scatter!(axz, zc, 0 * zc)
   hidedecorations!(axz)
   hidespines!(axz)

   scatter!(axΔz, zc, Δz)
   hidespines!(axΔz, :t, :r)

   rowsize!(fig.layout, 1, Relative(0.1))

   save(joinpath("./Plots", "zgrid_$(datetime).png"), fig)
end

function visualize_fields_const_x(datetime, x_idx; 
		                  bkgd_datetime = nothing, 
				  plot_animation = false, t_idx_skip = 1)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)

   if isnothing(bkgd_datetime)
      bkgd_datetime = datetime
   end

   bkgd_ds = open_bkgd_dataset(bkgd_datetime)
   B       = bkgd_ds[:B][:, :, :, 1]
   Uφ      = bkgd_ds[:Uφ][:, :, :, 1]

   z_plt = 1 #div(length(z[:]), 2) #z-index to start plotting at

   b_total_f_yz  = ds[:b][x_idx, :, z_plt:end, Nt]
   ur_total_f_yz = ds[:ur][x_idx, :, z_plt:end, Nt]
   uφ_total_f_yz = ds[:uφ][x_idx, :, z_plt:end, Nt]
   uz_total_f_yz = ds[:uz][x_idx, :, z_plt:end-1, Nt]

   Δb_f_yz  = b_total_f_yz .- B[x_idx, :, z_plt:end]
   Δuφ_f_yz = uφ_total_f_yz .- Uφ[x_idx, :, z_plt:end]

   lims_b_total  = get_range_lims(b_total_f_yz; max_fraction = 0.75)
   lims_ur       = get_range_lims(ur_total_f_yz; 
				  max_fraction = 0.75, prescribed_max = 1e-16)
   lims_uφ_total = get_range_lims(uφ_total_f_yz; 
				  max_fraction = 0.75, prescribed_max = 1e-16)
   lims_uz       = get_range_lims(uz_total_f_yz; 
				  max_fraction = 0.75, prescribed_max = 1e-16)

   lims_Δb  = get_range_lims(Δb_f_yz; 
			     max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δuφ = get_range_lims(Δuφ_f_yz; 
			     max_fraction = 0.75, prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   x_nearest, axis_kwargs_yz = get_2D_spatial_axis_kwargs(x, y, z; x_idx = x_idx)
   
   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      n = Observable(1)

      b_total_yz  = @lift ds[:b][x_idx, :, z_plt:end, $n]
      ur_total_yz = @lift ds[:ur][x_idx, :, z_plt:end, $n]
      uφ_total_yz = @lift ds[:uφ][x_idx, :, z_plt:end, $n]
      uz_total_yz = @lift ds[:uz][x_idx, :, z_plt:end-1, $n]

      Δb_yz  = @lift $b_total_yz .- B[x_idx, :, z_plt:end]
      Δuφ_yz = @lift $uφ_total_yz .- Uφ[x_idx, :, z_plt:end]
      
      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total  = Axis(fig_total[2, 1];
                        title = L"Total buoyancy ($b$)", axis_kwargs_yz...)
      ax_ur_total = Axis(fig_total[2, 3];
                        title = L"Total radial velocity ($u_r$)", 
			axis_kwargs_yz...)
      ax_uφ_total = Axis(fig_total[3, 1];
			 title = L"Total azimuthal velocity ($u_{\phi}$)", axis_kwargs_yz...)
      ax_uz_total = Axis(fig_total[3, 3];
                        title = L"Total vertical velocity ($u_z$)", 
			axis_kwargs_yz...)

      ax_b_perturb  = Axis(fig_perturb[2, 1];
                          title = L"Buoyancy perturbation ($b'$)", 
			  axis_kwargs_yz...)
      ax_ur_perturb = Axis(fig_perturb[2, 3];
                          title = L"Radial velocity perturbation ($u_r'$)", 
		          axis_kwargs_yz...)
      ax_uφ_perturb = Axis(fig_perturb[3, 1];
			   title = L"Azimuthal velocity perturbation ($u_{\phi}$')", 
		          axis_kwargs_yz...)
      ax_uz_perturb = Axis(fig_perturb[3, 3];
                          title = L"Vertical velocity perturbation ($u_z'$)", 
	                  axis_kwargs_yz...)

      hm_b_total  = heatmap!(ax_b_total, y, z[z_plt:end], b_total_yz, 
			    colorrange = lims_b_total, 
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
      hm_ur_total = heatmap!(ax_ur_total, y, z[z_plt:end], ur_total_yz,
                             colorrange = lims_ur, colormap = Reverse(:RdBu_5),
                             highclip = :red, lowclip = :blue)
      hm_uφ_total = heatmap!(ax_uφ_total, y, z[z_plt:end], uφ_total_yz,
                             colorrange = lims_uφ_total,
			     colormap = Reverse(:RdBu_5),
                             highclip = :red, lowclip = :blue)
      hm_uz_total = heatmap!(ax_uz_total, y, z[z_plt:end], uz_total_yz,
                             colorrange = lims_uz, colormap = Reverse(:RdBu_5),
                             highclip = :red, lowclip = :blue)

      hm_b_perturb  = heatmap!(ax_b_perturb, y, z[z_plt:end], Δb_yz,
                               colorrange = lims_Δb, 
		               colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
      hm_ur_perturb = heatmap!(ax_ur_perturb, y, z[z_plt:end], ur_total_yz,
                               colorrange = lims_ur, 
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
      hm_uφ_perturb = heatmap!(ax_uφ_perturb, y, z[z_plt:end], Δuφ_yz,
                               colorrange = lims_Δuφ, 
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
      hm_uz_perturb = heatmap!(ax_uz_perturb, y, z[z_plt:end], uz_total_yz,
                               colorrange = lims_uz, 
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)

      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", 
	       label = "m/s²")
      Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", 
	       label = "m/s")
      Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", 
	       label = "m/s")
      Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", 
	       label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}", 
               label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}", 
               label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}", 
               label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}", 
	       label = "m/s")

      title_total   = @lift @sprintf("Fields at x = %i km; t = %.2f days",
                                     x_nearest, times[$n])
      title_perturb = @lift @sprintf(
			    "Perturbation fields at x = %i km; t = %.2f days",
                            x_nearest, times[$n])

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

   ax_b_total  = Axis(fig_total[2, 1];
                      title = L"Total buoyancy ($b$)", axis_kwargs_yz...)
   ax_ur_total = Axis(fig_total[2, 3];
                     title = L"Total radial velocity ($u_r$)", axis_kwargs_yz...)
   ax_uφ_total = Axis(fig_total[3, 1];
		     title = L"Total azimuthal velocity ($u_{\phi}$)", axis_kwargs_yz...)
   ax_uz_total = Axis(fig_total[3, 3];
                     title = L"Total vertical velocity ($u_z$)", axis_kwargs_yz...)

   ax_b_perturb  = Axis(fig_perturb[2, 1]; 
		        title = L"Buoyancy perturbation ($b'$)", axis_kwargs_yz...)
   ax_ur_perturb = Axis(fig_perturb[2, 3];
                       title = L"Radial velocity perturbation ($u_r'$)",
                       axis_kwargs_yz...)
   ax_uφ_perturb = Axis(fig_perturb[3, 1];
			title = L"Azimuthal velocity perturbation ($u_{\phi}$')",
                       axis_kwargs_yz...)
   ax_uz_perturb = Axis(fig_perturb[3, 3];
                       title = L"Vertical velocity perturbation ($u_z'$)",
                       axis_kwargs_yz...)

   hm_b_total  = heatmap!(ax_b_total, y, z[z_plt:end], b_total_f_yz,
                          colorrange = lims_b_total, colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_ur_total = heatmap!(ax_ur_total, y, z[z_plt:end], ur_total_f_yz,
                          colorrange = lims_ur, colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_uφ_total = heatmap!(ax_uφ_total, y, z[z_plt:end], uφ_total_f_yz,
                          colorrange = lims_uφ_total, 
			  colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_uz_total = heatmap!(ax_uz_total, y, z[z_plt:end], uz_total_f_yz,
                          colorrange = lims_uz, colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)

   hm_b_perturb  = heatmap!(ax_b_perturb, y, z[z_plt:end], Δb_f_yz,
                            colorrange = lims_Δb, colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_ur_perturb = heatmap!(ax_ur_perturb, y, z[z_plt:end], ur_total_f_yz,
                            colorrange = lims_ur, colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_uφ_perturb = heatmap!(ax_uφ_perturb, y, z[z_plt:end], Δuφ_f_yz,
                            colorrange = lims_Δuφ, colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_uz_perturb = heatmap!(ax_uz_perturb, y, z[z_plt:end], uz_total_f_yz,
                            colorrange = lims_uz, colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}",
            label = "m/s")

   title_total   = @sprintf("Fields at x = %i km; t = %.2f days",
                            x_nearest, times[Nt])
   title_perturb = @sprintf(
                          "Perturbation fields at x = %i km; t = %.2f days",
                            x_nearest, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "fields_x$(x_nearest)_tf_$(datetime).png"),
                 fig_total)
   save(joinpath("./Plots", "perturbs_x$(x_nearest)_tf_$(datetime).png"),
                 fig_perturb)
   close(bkgd_ds)
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

   lims_b_total = get_range_lims(b_total_f_xz; max_fraction = 0.75)
   lims_u_total = get_range_lims(u_total_f_xz; 
				 max_fraction = 0.75, prescribed_max = 1e-16)
   lims_v_total = get_range_lims(v_total_f_xz; 
				 max_fraction = 0.75, prescribed_max = 1e-16)
   lims_w_total = get_range_lims(w_total_f_xz; 
				 max_fraction = 0.75, prescribed_max = 1e-16)

   lims_Δb = get_range_lims(Δb_f_xz; 
			    max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δu = get_range_lims(Δu_f_xz; 
			    max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δv = get_range_lims(Δv_f_xz; 
			    max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δw = get_range_lims(Δw_f_xz; 
			    max_fraction = 0.75, prescribed_max = 1e-16)

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
                            colorrange = lims_b_total, 
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
      hm_w_total = heatmap!(ax_w_total, x, z[z_plt:end], w_total_xz,
                            colorrange = lims_w_total,
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
      hm_u_total = heatmap!(ax_u_total, x, z[z_plt:end], u_total_xz,
                            colorrange = lims_u_total, 
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
      hm_v_total = heatmap!(ax_v_total, x, z[z_plt:end], v_total_xz,
                            colorrange = lims_v_total,
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)

      hm_b_perturb = heatmap!(ax_b_perturb, x, z[z_plt:end], Δb_xz,
                              colorrange = lims_Δb, 
			      colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)
      hm_w_perturb = heatmap!(ax_w_perturb, x, z[z_plt:end], Δw_xz,
                              colorrange = lims_Δw, 
			      colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)
      hm_u_perturb = heatmap!(ax_u_perturb, x, z[z_plt:end], Δu_xz,
                              colorrange = lims_Δu, 
			      colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)
      hm_v_perturb = heatmap!(ax_v_perturb, x, z[z_plt:end], Δv_xz,
                              colorrange = lims_Δv,
			      colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)

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

      title_total   = @lift @sprintf("Fields at y = %i km; t = %.2f days",
                                   y_nearest, times[$n])
      title_perturb = @lift @sprintf(
                            "Perturbation fields at y = %i km; t = %.2f days",
                             y_nearest, times[$n])

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
                            y_nearest, times[Nt])
   title_perturb = @sprintf(
                          "Perturbation fields at y = %i km; t = %.2f days",
                            y_nearest, times[Nt])

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
		                  bkgd_datetime = nothing,
				  plot_animation = true, t_idx_skip = 1)
   
   if isnothing(bkgd_datetime)
      bkgd_datetime = datetime
   end

   bkgd_ds = open_bkgd_dataset(bkgd_datetime)
   B       = bkgd_ds[:B][:, :, :, 1]
   Uφ      = bkgd_ds[:Uφ][:, :, :, 1]

   outfile_list             = glob("./Output/output_$(datetime)*")
   ds_f, x, y, z, times, Nt = open_dataset(outfile_list[length(outfile_list)])

   b_total_f_xy  = ds_f[:b][:, :, z_idx, Nt]
   ur_total_f_xy = ds_f[:ur][:, :, z_idx, Nt]
   uφ_total_f_xy = ds_f[:uφ][:, :, z_idx, Nt]
   uz_total_f_xy = ds_f[:uz][:, :, z_idx, Nt]

   Δb_f_xy  = b_total_f_xy .- B[:, :, z_idx]
   Δuφ_f_xy = uφ_total_f_xy .- Uφ[:, :, z_idx]

   lims_b_total  = get_range_lims(b_total_f_xy; max_fraction = 0.75)
   lims_ur       = get_range_lims(ur_total_f_xy; 
			          max_fraction = 0.75, prescribed_max = 1e-16)
   lims_uφ_total = get_range_lims(uφ_total_f_xy; 
				  max_fraction = 0.75, prescribed_max = 1e-16)
   lims_uz       = get_range_lims(uz_total_f_xy; 
				  max_fraction = 0.75, prescribed_max = 1e-16)

   lims_Δb  = get_range_lims(Δb_f_xy; 
			     max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δuφ = get_range_lims(Δuφ_f_xy; 
			     max_fraction = 0.75, prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   depth_nearest, axis_kwargs_xy = get_2D_spatial_axis_kwargs(x, y, z; 
							      z_idx = z_idx)

   #Plot static images (final frame, by default)
   
   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total  = Axis(fig_total[2, 1];
                      title = L"Total buoyancy ($b$)", axis_kwargs_xy...)
   ax_ur_total = Axis(fig_total[2, 3];
                      title = L"Total radial velocity ($u_r$)",
                      axis_kwargs_xy...)
   ax_uφ_total = Axis(fig_total[3, 1];
                      title = L"Total azimuthal velocity ($u_{\phi}$)", axis_kwargs_xy...)
   ax_uz_total = Axis(fig_total[3, 3];
                      title = L"Total vertical velocity ($u_z$)",
                      axis_kwargs_xy...)

   ax_b_perturb  = Axis(fig_perturb[2, 1];
                        title = L"Buoyancy perturbation ($b'$)",
                        axis_kwargs_xy...)
   ax_ur_perturb = Axis(fig_perturb[2, 3];
                        title = L"Radial velocity perturbation ($u_r'$)",
                        axis_kwargs_xy...)
   ax_uφ_perturb = Axis(fig_perturb[3, 1];
                        title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                        axis_kwargs_xy...)
   ax_uz_perturb = Axis(fig_perturb[3, 3];
                        title = L"Vertical velocity perturbation ($u_z'$)",
                        axis_kwargs_xy...)

   hm_b_total  = heatmap!(ax_b_total, x, y, b_total_f_xy,
                          colorrange = lims_b_total,
                          colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_ur_total = heatmap!(ax_ur_total, x, y, ur_total_f_xy,
                          colorrange = lims_ur,
                          colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_uφ_total = heatmap!(ax_uφ_total, x, y, uφ_total_f_xy,
                          colorrange = lims_uφ_total,
                          colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_uz_total = heatmap!(ax_uz_total, x, y, uz_total_f_xy,
                          colorrange = lims_uz,
                          colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)

   hm_b_perturb  = heatmap!(ax_b_perturb, x, y, Δb_f_xy,
                            colorrange = lims_Δb,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_ur_perturb = heatmap!(ax_ur_perturb, x, y, ur_total_f_xy,
                            colorrange = lims_ur,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, y, Δuφ_f_xy,
                            colorrange = lims_Δuφ,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_uz_perturb = heatmap!(ax_uz_perturb, x, y, uz_total_f_xy,
                            colorrange = lims_uz,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}",
            label = "m/s")

   title_total   = @sprintf("Fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt])
   title_perturb = @sprintf(
                          "Perturbation fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "fields_z$(depth_nearest)_tf_$(datetime).png"),
        fig_total)
   save(joinpath("./Plots", "perturbs_z$(depth_nearest)_tf_$(datetime).png"),
        fig_perturb)

   close(ds_f)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total  = Axis(fig_total[2, 1];
                         title = L"Total buoyancy ($b$)", axis_kwargs_xy...)
      ax_ur_total = Axis(fig_total[2, 3];
                         title = L"Total radial velocity ($u_r$)", 
		         axis_kwargs_xy...)
      ax_uφ_total = Axis(fig_total[3, 1];
			 title = L"Total azimuthal velocity ($u_{\phi}$)", 
			 axis_kwargs_xy...)
      ax_uz_total = Axis(fig_total[3, 3];
                         title = L"Total vertical velocity ($u_z$)", 
			 axis_kwargs_xy...)

      ax_b_perturb  = Axis(fig_perturb[2, 1];
                           title = L"Buoyancy perturbation ($b'$)", 
		           axis_kwargs_xy...)
      ax_ur_perturb = Axis(fig_perturb[2, 3];
                           title = L"Radial velocity perturbation ($u_r'$)",
                           axis_kwargs_xy...)
      ax_uφ_perturb = Axis(fig_perturb[3, 1];
			   title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                           axis_kwargs_xy...)
      ax_uz_perturb = Axis(fig_perturb[3, 3];
                           title = L"Vertical velocity perturbation ($u_z'$)",
                           axis_kwargs_xy...)

      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      
      ds, x, y, z, times, Nt = open_dataset(outfile_list)

      n = Observable(1)

      b_total_xy  = @lift ds[:b][:, :, z_idx, $n]
      ur_total_xy = @lift ds[:ur][:, :, z_idx, $n]
      uφ_total_xy = @lift ds[:uφ][:, :, z_idx, $n]
      uz_total_xy = @lift ds[:uz][:, :, z_idx, $n]

      Δb_xy  = @lift $b_total_xy .- B[:, :, z_idx]
      Δuφ_xy = @lift $uφ_total_xy .- Uφ[:, :, z_idx]
      
      hm_b_total  = heatmap!(ax_b_total, x, y, b_total_xy,
                             colorrange = lims_b_total,
			     colormap = Reverse(:RdBu_5),
			     highclip = :red, lowclip = :blue)
      hm_ur_total = heatmap!(ax_ur_total, x, y, ur_total_xy,
                             colorrange = lims_ur, 
			     colormap = Reverse(:RdBu_5),
                             highclip = :red, lowclip = :blue)
      hm_uφ_total = heatmap!(ax_uφ_total, x, y, uφ_total_xy,
                             colorrange = lims_uφ_total,
			     colormap = Reverse(:RdBu_5),
                             highclip = :red, lowclip = :blue)
      hm_uz_total = heatmap!(ax_uz_total, x, y, uz_total_xy,
                             colorrange = lims_uz,
			     colormap = Reverse(:RdBu_5),
                             highclip = :red, lowclip = :blue)

      hm_b_perturb  = heatmap!(ax_b_perturb, x, y, Δb_xy,
                               colorrange = lims_Δb, 
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
      hm_ur_perturb = heatmap!(ax_ur_perturb, x, y, ur_total_xy,
                               colorrange = lims_ur,
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
      hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, y, Δuφ_xy,
                               colorrange = lims_Δuφ,
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
      hm_uz_perturb = heatmap!(ax_uz_perturb, x, y, uz_total_xy,
                               colorrange = lims_uz,
			       colormap = Reverse(:RdBu_5),
                               highclip = :red, lowclip = :blue)
   
      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
               label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}",
               label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}",
               label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}",
               label = "m/s")

      title_total   = @lift @sprintf("Fields at %i-m depth; t = %.2f days",
					depth_nearest, times[$n])
      title_perturb = @lift @sprintf(
                            "Perturbation fields at %i-m depth; t = %.2f days",
                                     depth_nearest, times[$n])

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                                  tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                                  tellwidth = false)
   
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)
         
      for i = 1:t_idx_skip:Nt
         print(i, " of ", Nt)
	 recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end

      close(ds)
      save(joinpath("./Plots", "fields_z$(depth_nearest)_$(datetime).mp4"), 
	            video_total)
      save(joinpath("./Plots", "perturbs_z$(depth_nearest)_$(datetime).mp4"),
	            video_perturb)
   end
   
   #=Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total  = Axis(fig_total[2, 1];
                      title = L"Total buoyancy ($b$)", axis_kwargs_xy...)
   ax_ur_total = Axis(fig_total[2, 3];
                      title = L"Total radial velocity ($u_r$)", 
		      axis_kwargs_xy...)
   ax_uφ_total = Axis(fig_total[3, 1];
		      title = L"Total azimuthal velocity ($u_{\phi}$)", axis_kwargs_xy...)
   ax_uz_total = Axis(fig_total[3, 3];
		      title = L"Total vertical velocity ($u_z$)",
                      axis_kwargs_xy...)

   ax_b_perturb  = Axis(fig_perturb[2, 1];
                        title = L"Buoyancy perturbation ($b'$)",
                        axis_kwargs_xy...)
   ax_ur_perturb = Axis(fig_perturb[2, 3];
                        title = L"Radial velocity perturbation ($u_r'$)",
                        axis_kwargs_xy...)
   ax_uφ_perturb = Axis(fig_perturb[3, 1];
		 	title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                        axis_kwargs_xy...)
   ax_uz_perturb = Axis(fig_perturb[3, 3];
                        title = L"Vertical velocity perturbation ($u_z'$)",
                        axis_kwargs_xy...)

   hm_b_total  = heatmap!(ax_b_total, x, y, b_total_f_xy,
                          colorrange = lims_b_total,
			  colormap = Reverse(:RdBu_5),
			  highclip = :red, lowclip = :blue)
   hm_ur_total = heatmap!(ax_ur_total, x, y, ur_total_f_xy,
                          colorrange = lims_ur,
			  colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_uφ_total = heatmap!(ax_uφ_total, x, y, uφ_total_f_xy,
                          colorrange = lims_uφ_total,
			  colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)
   hm_uz_total = heatmap!(ax_uz_total, x, y, uz_total_f_xy,
                          colorrange = lims_uz,
			  colormap = Reverse(:RdBu_5),
                          highclip = :red, lowclip = :blue)

   hm_b_perturb  = heatmap!(ax_b_perturb, x, y, Δb_f_xy,
                            colorrange = lims_Δb, 
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_ur_perturb = heatmap!(ax_ur_perturb, x, y, ur_total_f_xy,
                            colorrange = lims_ur,
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, y, Δuφ_f_xy,
                            colorrange = lims_Δuφ,
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
   hm_uz_perturb = heatmap!(ax_uz_perturb, x, y, uz_total_f_xy,
                            colorrange = lims_uz,
			    colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}",
            label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}",
            label = "m/s")

   title_total   = @sprintf("Fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt])
   title_perturb = @sprintf(
                          "Perturbation fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots", "fields_z$(depth_nearest)_tf_$(datetime).png"),
        fig_total)
   save(joinpath("./Plots", "perturbs_z$(depth_nearest)_tf_$(datetime).png"),
        fig_perturb) =#
   close(bkgd_ds)
   #close(ds)
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
			           x_nearest, times[$n])
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
			            y_nearest, times[$n])
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
			     depth_nearest, times[$n])
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
