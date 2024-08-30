include("Library.jl")

using Oceananigans
using CairoMakie, NCDatasets, Printf
using .ComputeSecondaries

function open_dataset(datetime)

   outfilepath = joinpath("./Output", "output_$(datetime).nc")
   ds          = NCDataset(outfilepath, "r")

   x = ds["xC"][:]
   y = ds["yC"][:]
   z = ds["zC"][:]

   times = ds["time"][:]
   Nt    = length(times)

   return ds, x, y, z, times, Nt

end

function get_background_fields(ds)
   bb = ds["b"][:, :, :, 1]
   ub = ds["u"][:, :, :, 1]
   vb = ds["v"][:, :, :, 1]
   wb = ds["w"][:, :, 1:end-1, 1]
   return bb, ub, vb, wb
end

function get_range_lims(final_field; prescribed_max = 0)
   field_max = max(maximum(abs.(final_field)), prescribed_max)
   field_lims = (3/4) * [-field_max, field_max]
end

function get_axis_kwargs(x, y, z; 
		x_idx = nothing, y_idx = nothing, z_idx = nothing)
   if !isnothing(z_idx)
      nearest_m   = round(Int, -z[z_idx])
      axis_kwargs = (xlabel = "x [m]", ylabel = "y [m]")
   end
   return nearest_m, axis_kwargs
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
         n[]=i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "bwuv_z$(depth_nearest_m)_$(datetime).mp4"), video)
   close(ds)

end

function visualize_q_const_z(datetime, Δx, Δy, Δz, f, z_idx)

   ds, x, y, z, times, Nt = open_dataset(datetime)

   n = Observable(1)

   b    = @lift ds["b"][:, :, :, $n]
   u    = @lift ds["u"][:, :, :, $n]
   v    = @lift ds["v"][:, :, :, $n]
   w    = @lift ds["w"][:, :, :, $n]
   q_xy = @lift ertelQ_2D($u, $v, $w, $b, f, Δx, Δy, Δz; z_idx = z_idx)

   qf_xy = ertelQ_2D(ds["u"][:, :, :, Nt], ds["v"][:, :, :, Nt], 
		     ds["w"][:, :, :, Nt], ds["b"][:, :, :, Nt],
   		     f, Δx, Δy, Δz; z_idx = z_idx)

   lims_q = get_range_lims(qf_xy)

   depth_nearest_m, axis_kwargs_xy = get_axis_kwargs(x, y, z; z_idx = z_idx)

   fig  = Figure(size = (600, 600))
   ax_q = Axis(fig[2, 1]; axis_kwargs_xy...)
   hm_q = heatmap!(ax_q, x, y, q_xy, colorrange = lims_q, colormap = :balance)

   Colorbar(fig[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title = @lift @sprintf("Potential vorticity at depth %i m; t = %.2f days",
			  depth_nearest_m, times[$n]/(3600*24))
   
   fig[1, 1:2] = Label(fig, title, fontsize = 24, tellwidth = false)

   frames = 1:Nt
   video  = VideoStream(fig, format = "mp4", framerate = 6)

   for i=1:frames[end]
      recordframe!(video)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[]=i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_z$(depth_nearest_m)_$(datetime).mp4"), video)
   close(ds)

end
#=
## Next, plot at constant x

x     = ds["xC"][x_idx]
z     = ds["zC"][4:end-2]
const z_plt = div(length(z[:]),2) #z-index to start plot at

#Use initial frame to define background fields
bb = ds["b"][x_idx, :, 4:end-2, 1]
ub = ds["u"][x_idx, :, :, 1]
vb = ds["v"][x_idx, :, :, 1]
wb = ds["w"][x_idx, :, :, 1]

n = Observable(1)

b    = @lift ds["b"][x_idx, :, 4:end-2, $n] .- bb
btot = @lift ds["b"][x_idx, :, 4:end-2, $n]
u    = @lift ds["u"][:, :, 4:end-2, $n]
v    = @lift ds["v"][:, :, 4:end-2, $n]
w    = @lift ds["w"][:, :, 4:end-3, $n]

ωtotal = @lift ζ_2D($u, $v, $w, Δx, Δy, Δz, x_idx, nothing, nothing)
∇_b    = @lift ∇b_2D($btot, Δx, Δy, Δz, 'x')
fq     = @lift @. f*($ωtotal + f) * $∇_b
fq_r   = @lift ∂r_fq($fq, Δx, Δy, x, y[2:end-1])

#Restrict fields to plotting domain
u_yz  = @lift $u[x_idx, :, z_plt:length(z[:])-1]
v_yz  = @lift $v[x_idx, :, z_plt:length(z[:])-1]
w_yz  = @lift $w[x_idx, :, z_plt:length(z[:])-1]
b_yz  = @lift $b[:, z_plt:length(z[:])-1]
fq_yz = @lift $fq[:, z_plt:length(z[:])-1]
fqr_yz = @lift $fq_r[:, z_plt:length(z[:])-2]

#Use final frames to define maximum values
bf    = ds["b"][x_idx, :, z_plt:length(z[:])-1, length(times)] .- ds["b"][x_idx, :, z_plt:length(z[:])-1, 1]
uf_yz = ds["u"][x_idx, :, z_plt:length(z[:])-1, length(times)]
wf_yz = ds["w"][x_idx, :, z_plt:length(z[:])-1, length(times)]

max_b  =  max(maximum(bf), 5e-9)
lim_b  =  (3/4) * [-max_b, max_b]
max_u  =  max(maximum(uf_yz), 1e-10)
lim_u  =  (3/4) * [-max_u, max_u]
max_w  =  max(maximum(wf_yz), 1e-10)
lim_w  =  (3/4) * [-max_w, max_w]
lim_fq =  1 .* [-1e-11, 1e-11]
lim_fqr = 1 .* [-1e-15, 1e-15]

fig3 = Figure(size = (1200, 800))
axis_kwargs_yz = (xlabel = "y [m]", ylabel = "z [m]")
ax_b    = Axis(fig3[2, 1]; title = "Buoyancy perturbation (b')", axis_kwargs_yz...)
ax_w    = Axis(fig3[2, 3]; title = "Vertical velocity (w)", axis_kwargs_yz...)
ax_u    = Axis(fig3[3, 1]; title = "Zonal velocity (u)", axis_kwargs_yz...)
ax_v    = Axis(fig3[3, 3]; title = "Meridional velocity (v)", axis_kwargs_yz...)

hm_b = heatmap!(ax_b, y, z[z_plt:length(z[:])-1], b_yz, colorrange = lim_b, colormap = :balance)
Colorbar(fig3[2, 2], hm_b, tickformat = "{:.1e}", label = "m/s²")
hm_w = heatmap!(ax_w, y, z[z_plt:length(z[:])-1], w_yz, colorrange = lim_w, colormap = :balance)
Colorbar(fig3[2, 4], hm_w, tickformat = "{:.1e}", label = "m/s")
hm_u = heatmap!(ax_u, y, z[z_plt:length(z[:])-1], u_yz, colorrange = lim_u, colormap = :balance)
Colorbar(fig3[3, 2], hm_u, tickformat = "{:.1e}", label = "m/s")
hm_v = heatmap!(ax_v, y, z[z_plt:length(z[:])-1], v_yz, colorrange = lim_u, colormap = :balance)
Colorbar(fig3[3, 4], hm_v, tickformat = "{:.1e}", label = "m/s")

fig4 = Figure(size = (800, 500))
ax_fq = Axis(fig4[2, 1]; axis_kwargs_yz...)
hm_fq = heatmap!(ax_fq, y, z[z_plt:length(z[:])-1], fq_yz, colorrange = lim_fq, colormap = :balance)# cm)
Colorbar(fig4[2, 2], hm_fq, tickformat = "{:.1e}", label = "1/s³")

fig5 = Figure(size = (800, 500))
ax_fqr = Axis(fig5[2, 1]; axis_kwargs_yz...)
hm_fqr = heatmap!(ax_fqr, y, z[z_plt:length(z[:])-1], fqr_yz, colorrange = lim_fqr, colormap = :balance)
Colorbar(fig5[2, 2], hm_fqr, tickformat = "{:.1e}", label = "1/s³m")

title3 = @lift @sprintf("Fields at x = 0; t = %.2f days", times[$n]/(3600*24))
fig3[1, 1:4] = Label(fig3, title3, fontsize = 24, tellwidth = false)

title4 = @lift @sprintf("Potential vorticity at x = 0; t = %.2f days", times[$n]/(3600*24))
fig4[1, 1:2] = Label(fig4, title4, fontsize = 24, tellwidth = false)

title5 = @lift @sprintf("Radial derivative of potential vorticity at x = 0; t = %.2f days", times[$n]/(3600*24))
fig5[1, 1:2] = Label(fig5, title5, fontsize = 24, tellwidth = false)

frames = 1:length(times)

video3 = VideoStream(fig3, format = "mp4", framerate = 6)
video4 = VideoStream(fig4, format = "mp4", framerate = 6)
video5 = VideoStream(fig5, format = "mp4", framerate = 6)

for i=1:frames[end]
    recordframe!(video3)
    recordframe!(video4)
    recordframe!(video5)
    msg = string("Plotting frame(s) ", i, " of ", frames[end])
        print(msg * " \r")
        n[]=i
end

save(joinpath("./Plots", "bwuv_x0_$(datetime).mp4"), video3)
save(joinpath("./Plots", "fq_x0_$(datetime).mp4"), video4)
save(joinpath("./Plots", "fqr_x0_$(datetime).mp4"), video5)

close(ds)
=#
