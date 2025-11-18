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

   scalars_ds, times = open_scalars_dataset("scalars_$(datetime).nc")

   b′_norm  = scalars_ds[:b′_norm][2:end]
   ux′_norm = scalars_ds[:ux′_norm][2:end]
   uy′_norm = scalars_ds[:uy′_norm][2:end]
   ur′_norm = scalars_ds[:ur′_norm][2:end]
   uφ′_norm = scalars_ds[:uφ′_norm][2:end]
   uz′_norm = scalars_ds[:uz′_norm][2:end]
   
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

function visualize_energetics(datetime, grid)

   outfile_list         = glob("./Output/energetics_$(datetime)*")
   energetics_ds, t, Nt = open_energetics_dataset(outfile_list)

   pKE_data = energetics_ds[:pKE][:, :, :, :]
   
   fig = Figure(size = (1200, 700))
   ax  = Axis(fig[2, 1]; xlabel = "Time [days]", ylabel = "Energy [m^5/s^2]",
                         yscale = log10)

   n = Observable(1)

   pKE_Field_n = CenterField(grid)
   @lift set!(pKE_Field_n, pKE_data[:, :, :, $n])

   for i = 1:Nt

      integrated_pKE_n = Field(Integral(pKE_Field_n))
      
      compute!(integrated_pKE_n)
      scatter!(ax, t[i], integrated_pKE_n[1], color = :black)
      yield()
      
      n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent

   fig[1, 1] = Label(fig, "Volume-integrated perturbation kinetic energy",
                     fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "pKE_$(datetime).png"), fig)
   close(energetics_ds)
end

function visualize_b_and_ωz(datetime, Δx, Δy; 
		            x_idx = nothing, y_idx = nothing, z_idx = nothing,
		            plot_animation = true, t_idx_skip = 1)

   B  = no_offset_view(adapt(Array, B))[]
   Ux = no_offset_view(adapt(Array, Ux))[Hx+1:length(x)+Hx, 1, Hz+1:length(zC)+Hz]
   Uy = no_offset_view(adapt(Array, Uy))

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

function visualize_fields_const_y(datetime, y_idx, B, Uφ, Hx, Hz;
		                  plot_animation = true, t_idx_skip = 1)
  
   outfile_list                  = glob("./Output/output_$(datetime)*")
   ds_f, x, y, zC, zF, times, Nt = open_dataset(outfile_list[length(outfile_list)]; 
						Hx = Hx, Hz = Hz)

   xCzC_idcs, xCzF_idcs = get_2D_spatial_axis_idcs("y";
                                  Hx = Hx, Hy = Hy, Hz = Hz,
                                  y_idx = y_idx,
                                  xC = x, yC = y, zC = zC, zF = zF)

   B  = no_offset_view(adapt(Array, B))[xCzC_idcs...]
   Uφ = no_offset_view(adapt(Array, Uφ))[xCzC_idcs...]

   b_total_f_xz  = adapt(Array, ds_f[:b])[xCzC_idcs..., Nt]
   ur_total_f_xz = adapt(Array, ds_f[:ur])[xCzC_idcs..., Nt]
   uφ_total_f_xz = adapt(Array, ds_f[:uφ])[xCzC_idcs..., Nt]
   uz_total_f_xz = adapt(Array, ds_f[:uz])[xCzF_idcs..., Nt]

   Δb_f_xz  = b_total_f_xz .- B
   Δuφ_f_xz = uφ_total_f_xz .- Uφ

   lims_b_total  = get_range_lims(b_total_f_xz)
   lims_ur       = get_range_lims(ur_total_f_xz)
   lims_uφ_total = get_range_lims(uφ_total_f_xz)
   lims_uz       = get_range_lims(uz_total_f_xz)

   lims_Δb  = get_range_lims(Δb_f_xz)
   lims_Δuφ = get_range_lims(Δuφ_f_xz)

   mkpath("./Plots") #Make visualization directory if nonexistent

   y_nearest, ax_kwargs_xz = get_2D_spatial_axis_kwargs(x, y, zC; 
							y_idx = y_idx, 
							const_y = true)

   #Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total = Axis(fig_total[2, 1];
                      title = L"Total buoyancy ($b$)", ax_kwargs_xz...)
   ax_ur_total = Axis(fig_total[2, 3];
                      title = L"Total radial velocity ($u_r$)", ax_kwargs_xz...)
   ax_uφ_total = Axis(fig_total[3, 1];
		      title = L"Total azimuthal velocity ($u_{\phi}$)", 
		      ax_kwargs_xz...)
   ax_uz_total = Axis(fig_total[3, 3];
                      title = L"Total vertical velocity ($u_z$)", 
		      ax_kwargs_xz...)

   ax_b_perturb  = Axis(fig_perturb[2, 1];
                        title = L"Buoyancy perturbation ($b'$)", 
			ax_kwargs_xz...)
   ax_ur_perturb = Axis(fig_perturb[2, 3];
                        title = L"Radial velocity perturbation ($u_r'$)",
                        ax_kwargs_xz...)
   ax_uφ_perturb = Axis(fig_perturb[3, 1];
		        title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                        ax_kwargs_xz...)
   ax_uz_perturb = Axis(fig_perturb[3, 3];
                        title = L"Vertical velocity perturbation ($u_z'$)",
                        ax_kwargs_xz...)

   hm_b_total = heatmap!(ax_b_total, x, zC, b_total_f_xz, colorrange = lims_b_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_ur_total = heatmap!(ax_ur_total, x, zC, ur_total_f_xz, colorrange = lims_ur, colormap = Reverse(:RdBu_5))
   hm_uφ_total = heatmap!(ax_uφ_total, x, zC, uφ_total_f_xz, colorrange = lims_uφ_total, colormap = Reverse(:RdBu_5))
   hm_uz_total = heatmap!(ax_uz_total, x, zF, uz_total_f_xz, colorrange = lims_uz, colormap = Reverse(:RdBu_5))

   hm_b_perturb = heatmap!(ax_b_perturb, x, zC, Δb_f_xz, colorrange = lims_Δb, colormap = Reverse(:RdBu_5))
   hm_ur_perturb = heatmap!(ax_ur_perturb, x, zC, ur_total_f_xz, colorrange = lims_ur, colormap = Reverse(:RdBu_5))
   hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, zC, Δuφ_f_xz, colorrange = lims_Δuφ, colormap = Reverse(:RdBu_5))
   hm_uz_perturb = heatmap!(ax_uz_perturb, x, zF, uz_total_f_xz, colorrange = lims_uz, colormap = Reverse(:RdBu_5))

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}", label = "m/s")

   title_total   = @sprintf("Fields at y = %i km; t = %.2f days",
                            y_nearest, times[Nt])
   title_perturb = @sprintf(
                          "Perturbation fields at y = %i km; t = %.2f days",
                            y_nearest, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24, tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "fields_y$(y_nearest)_tf_$(datetime).png"), fig_total)
   save(joinpath("./Plots", "perturbs_y$(y_nearest)_tf_$(datetime).png"), fig_perturb)
   close(ds_f)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip
   
      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total  = Axis(fig_total[2, 1];
                         title = L"Total buoyancy ($b$)", ax_kwargs_xz...)
      ax_ur_total = Axis(fig_total[2, 3];
                         title = L"Total radial velocity ($u_r$)", 
			 ax_kwargs_xz...)
      ax_uφ_total = Axis(fig_total[3, 1];
			 title = L"Total azimuthal velocity ($u_{\phi}$)",
			 ax_kwargs_xz...)
      ax_uz_total = Axis(fig_total[3, 3];
                         title = L"Total vertical velocity ($u_z$)", 
			 ax_kwargs_xz...)

      ax_b_perturb  = Axis(fig_perturb[2, 1];
                           title = L"Buoyancy perturbation ($b'$)", 
		 	   ax_kwargs_xz...)
      ax_ur_perturb = Axis(fig_perturb[2, 3];
                           title = L"Radial velocity perturbation ($u_r'$)",
                           ax_kwargs_xz...)
      ax_uφ_perturb = Axis(fig_perturb[3, 1];
			   title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                           ax_kwargs_xz...)
      ax_uz_perturb = Axis(fig_perturb[3, 3];
                           title = L"Vertical velocity perturbation ($u_z'$)",
                           ax_kwargs_xz...)

      ds, x, y, zC, zF, times, Nt = open_dataset(outfile_list)

      n = Observable(1)

      b_total_xz  = @lift ds[:b][xCzC_idcs..., $n]
      ur_total_xz = @lift ds[:ur][xCzC_idcs..., $n]
      uφ_total_xz = @lift ds[:uφ][xCzC_idcs..., $n]
      uz_total_xz = @lift ds[:uz][xCzF_idcs..., $n]

      Δb_xz  = @lift $b_total_xz .- B
      Δuφ_xz = @lift $uφ_total_xz .- Uφ
      
      hm_b_total  = heatmap!(ax_b_total, x, zC, b_total_xz, colorrange = lims_b_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_ur_total = heatmap!(ax_ur_total, x, zC, ur_total_xz, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uφ_total = heatmap!(ax_uφ_total, x, zC, uφ_total_xz, colorrange = lims_uφ_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uz_total = heatmap!(ax_uz_total, x, zF, uz_total_xz, colorrange = lims_uz, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)

      hm_b_perturb  = heatmap!(ax_b_perturb, x, zC, Δb_xz, colorrange = lims_Δb, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_ur_perturb = heatmap!(ax_ur_perturb, x, zC, ur_total_xz, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, zC, Δuφ_xz, colorrange = lims_Δuφ, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uz_perturb = heatmap!(ax_uz_perturb, x, zF, uz_total_xz, colorrange = lims_uz, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)

      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}", label = "m/s")

      title_total   = @lift @sprintf("Fields at y = %i km; t = %.2f days",
                                     y_nearest, times[$n])
      title_perturb = @lift @sprintf(
                            "Perturbation fields at y = %i km; t = %.2f days",
                                     y_nearest, times[$n])

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24, tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24, tellwidth = false)

      frames = 1:Nt
  
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
	 print(i, " of ", Nt, "\n")
         recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end

      save(joinpath("./Plots", "fields_y$(y_nearest)_$(datetime).mp4"), video_total)
      save(joinpath("./Plots", "perturbs_y$(y_nearest)_$(datetime).mp4"), video_perturb)
   end
   close(ds)
end

function visualize_fields_const_z(datetime, z_idx, B, Uφ; Hx = 3, Hy = 3, 
		                  plot_animation = true, t_idx_skip = 1)

   outfile_list                  = glob("./Output/output_$(datetime)*")
   ds_f, x, y, zC, zF, times, Nt = open_dataset(outfile_list[length(outfile_list)];
						Hx = Hx, Hy = Hy)

   xyz_idcs = get_2D_spatial_axis_idcs("z";
				       Hx = Hx, Hy = Hy,
                                       z_idx = z_idx,
                                       xC = x, yC = y)

   B  = no_offset_view(adapt(Array, B))[xyz_idcs...]
   Uφ = no_offset_view(adapt(Array, Uφ))[xyz_idcs...]

   b_total_f_xy  = adapt(Array, ds_f[:b])[xyz_idcs..., Nt]
   ur_total_f_xy = adapt(Array, ds_f[:ur])[xyz_idcs..., Nt]
   uφ_total_f_xy = adapt(Array, ds_f[:uφ])[xyz_idcs..., Nt]
   uz_total_f_xy = adapt(Array, ds_f[:uz])[xyz_idcs..., Nt]
   
   Δb_f_xy  = b_total_f_xy .- B
   Δuφ_f_xy = uφ_total_f_xy .- Uφ
   
   lims_b_total  = get_range_lims(b_total_f_xy; max_fraction = 1.0)
   lims_ur       = get_range_lims(ur_total_f_xy; 
   			          max_fraction = 1.0, prescribed_max = 1e-16)
   lims_uφ_total = get_range_lims(uφ_total_f_xy; 
				  max_fraction = 1.0, prescribed_max = 1e-16)
   lims_uz       = get_range_lims(uz_total_f_xy; 
				  max_fraction = 1.0, prescribed_max = 1e-16)

   lims_Δb  = get_range_lims(Δb_f_xy; 
   			     max_fraction = 1.0, prescribed_max = 1e-16)
   lims_Δuφ = get_range_lims(Δuφ_f_xy; 
			     max_fraction = 1.0, prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   depth_nearest, ax_kwargs_xy = get_2D_spatial_axis_kwargs(x, y, zC; z_idx = z_idx)

   #Plot static images (final frame, by default)
   
   fig_total   = Figure(size = (1200, 800))
   fig_perturb = Figure(size = (1200, 800))

   ax_b_total  = Axis(fig_total[2, 1];
                      title = L"Total buoyancy ($b$)", ax_kwargs_xy...)
   ax_ur_total = Axis(fig_total[2, 3];
                      title = L"Total radial velocity ($u_r$)",
                      ax_kwargs_xy...)
   ax_uφ_total = Axis(fig_total[3, 1];
                      title = L"Total azimuthal velocity ($u_{\phi}$)", 
		      ax_kwargs_xy...)
   ax_uz_total = Axis(fig_total[3, 3];
                      title = L"Total vertical velocity ($u_z$)",
                      ax_kwargs_xy...)

   ax_b_perturb  = Axis(fig_perturb[2, 1];
                        title = L"Buoyancy perturbation ($b'$)",
                        ax_kwargs_xy...)
   ax_ur_perturb = Axis(fig_perturb[2, 3];
                        title = L"Radial velocity perturbation ($u_r'$)",
                        ax_kwargs_xy...)
   ax_uφ_perturb = Axis(fig_perturb[3, 1];
                        title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                        ax_kwargs_xy...)
   ax_uz_perturb = Axis(fig_perturb[3, 3];
                        title = L"Vertical velocity perturbation ($u_z'$)",
                        ax_kwargs_xy...)

   hm_b_total  = heatmap!(ax_b_total, x, y, b_total_f_xy, colorrange = lims_b_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_ur_total = heatmap!(ax_ur_total, x, y, ur_total_f_xy, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_uφ_total = heatmap!(ax_uφ_total, x, y, uφ_total_f_xy, colorrange = lims_uφ_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_uz_total = heatmap!(ax_uz_total, x, y, uz_total_f_xy, colorrange = lims_uz, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)

   hm_b_perturb  = heatmap!(ax_b_perturb, x, y, Δb_f_xy, colorrange = lims_Δb, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_ur_perturb = heatmap!(ax_ur_perturb, x, y, ur_total_f_xy, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, y, Δuφ_f_xy, colorrange = lims_Δuφ, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_uz_perturb = heatmap!(ax_uz_perturb, x, y, uz_total_f_xy, colorrange = lims_uz, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}", label = "m/s")

   title_total   = @sprintf("Fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt])
   title_perturb = @sprintf(
                          "Perturbation fields at %i-m depth; t = %.2f days",
                            depth_nearest, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24, tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "fields_z$(depth_nearest)_tf_$(datetime).png"), fig_total)
   save(joinpath("./Plots", "perturbs_z$(depth_nearest)_tf_$(datetime).png"), fig_perturb)
   close(ds_f)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      fig_total   = Figure(size = (1200, 800))
      fig_perturb = Figure(size = (1200, 800))

      ax_b_total  = Axis(fig_total[2, 1];
                         title = L"Total buoyancy ($b$)", ax_kwargs_xy...)
      ax_ur_total = Axis(fig_total[2, 3];
                         title = L"Total radial velocity ($u_r$)", 
		         ax_kwargs_xy...)
      ax_uφ_total = Axis(fig_total[3, 1];
			 title = L"Total azimuthal velocity ($u_{\phi}$)", 
			 ax_kwargs_xy...)
      ax_uz_total = Axis(fig_total[3, 3];
                         title = L"Total vertical velocity ($u_z$)", 
			 ax_kwargs_xy...)

      ax_b_perturb  = Axis(fig_perturb[2, 1];
                           title = L"Buoyancy perturbation ($b'$)", 
		           ax_kwargs_xy...)
      ax_ur_perturb = Axis(fig_perturb[2, 3];
                           title = L"Radial velocity perturbation ($u_r'$)",
                           ax_kwargs_xy...)
      ax_uφ_perturb = Axis(fig_perturb[3, 1];
			   title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                           ax_kwargs_xy...)
      ax_uz_perturb = Axis(fig_perturb[3, 3];
                           title = L"Vertical velocity perturbation ($u_z'$)",
                           ax_kwargs_xy...)

      ds, x, y, zC, zF, times, Nt = open_dataset(outfile_list)

      n = Observable(1)

      b_total_xy  = @lift ds[:b][xyz_idcs..., $n]
      ur_total_xy = @lift ds[:ur][xyz_idcs..., $n]
      uφ_total_xy = @lift ds[:uφ][xyz_idcs..., $n]
      uz_total_xy = @lift ds[:uz][xyz_idcs..., $n]

      Δb_xy  = @lift $b_total_xy .- B
      Δuφ_xy = @lift $uφ_total_xy .- Uφ
      
      hm_b_total  = heatmap!(ax_b_total, x, y, b_total_xy, colorrange = lims_b_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_ur_total = heatmap!(ax_ur_total, x, y, ur_total_xy, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uφ_total = heatmap!(ax_uφ_total, x, y, uφ_total_xy, colorrange = lims_uφ_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uz_total = heatmap!(ax_uz_total, x, y, uz_total_xy, colorrange = lims_uz, colormap = Reverse(:RdBu_5), highclip = :red, lowclip = :blue)

      hm_b_perturb  = heatmap!(ax_b_perturb, x, y, Δb_xy, colorrange = lims_Δb, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_ur_perturb = heatmap!(ax_ur_perturb, x, y, ur_total_xy, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uφ_perturb = heatmap!(ax_uφ_perturb, x, y, Δuφ_xy, colorrange = lims_Δuφ, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uz_perturb = heatmap!(ax_uz_perturb, x, y, uz_total_xy, colorrange = lims_uz, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   
      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

      Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_perturb[2, 4], hm_ur_perturb, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_perturb[3, 2], hm_uφ_perturb, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_perturb[3, 4], hm_uz_perturb, tickformat = "{:.1e}", label = "m/s")

      title_total   = @lift @sprintf("Fields at %i-m depth; t = %.2f days",
				     depth_nearest, times[$n])
      title_perturb = @lift @sprintf(
                            "Perturbation fields at %i-m depth; t = %.2f days",
                                     depth_nearest, times[$n])

      fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24, tellwidth = false)
      fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24, tellwidth = false)
   
      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)
         
      for i = 1:t_idx_skip:Nt
         print(i, " of ", Nt, "\n")
	 recordframe!(video_total)
         recordframe!(video_perturb)
         yield()
         n[] = i
      end

      save(joinpath("./Plots", "fields_z$(depth_nearest)_$(datetime).mp4"), video_total)
      save(joinpath("./Plots", "perturbs_z$(depth_nearest)_$(datetime).mp4"), video_perturb)
   end
   close(ds)
end

function open_computed_dataset(datetime, Δx, Δy, Δz, f)

   computed_file = joinpath("./Output", "computed_$(datetime).nc")

   if !isfile(computed_file) #Only compute if file does not already exist

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
