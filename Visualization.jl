include("Library.jl")

using Oceananigans
using CairoMakie, NCDatasets, Printf
using .ComputeSecondaries

function visualize_const_z(outfilepath, z_idx, Δx, Δy, Δz, f, datetime)
   
   ds = NCDataset(outfilepath, "r")

   x = ds["xC"][:]
   y = ds["yC"][:]
   z = ds["zC"][:]

   depth_nearest_m = round(Int, -(ds["zC"][z_idx]))

   times = ds["time"][:]
   Nt    = length(times)

   n     = Observable(1)

   #Use initial frame to define background fields
   bb = ds["b"][:, :, :, 1]
   #ub = ds["u"][:, :, z_idx, 1]
   #vb = ds["v"][:, :, z_idx, 1]
   #wb = ds["w"][:, :, z_idx, 1]

   b    = @lift ds["b"][:, :, :, $n] .- bb
   btot = @lift ds["b"][:, :, :, $n]
   u    = @lift ds["u"][:, :, :, $n]
   v    = @lift ds["v"][:, :, :, $n]
   w    = @lift ds["w"][:, :, 1:end-1, $n]

   ωtot = @lift ζ_2D($u, $v, $w, Δx, Δy, Δz, nothing, nothing, z_idx)
   ∇_b  = @lift ∇b_2D($btot, Δx, Δy, Δz, nothing, nothing, z_idx)
   q    = @lift ertelQ_2D($u, $v, $w, $b, f, Δx, Δy, Δz, nothing, nothing, z_idx)

   b_xy = @lift $b[:, :, z_idx]
   u_xy = @lift $u[:, :, z_idx]
   v_xy = @lift $v[:, :, z_idx]
   w_xy = @lift $w[:, :, z_idx]

   #Use final frames to estimate maximum values
   bf    = ds["b"][:, :, z_idx, Nt] .- bb[:, :, z_idx]
   uf_xy = ds["u"][:, :, z_idx, Nt]
   wf_xy = ds["w"][:, :, z_idx, Nt]
   qf    = ertelQ_2D(ds["u"][:,:,:,Nt], ds["v"][:,:,:,Nt], ds["w"][:,:,:,Nt], ds["b"][:,:,:,Nt],
		     f, Δx, Δy, Δz, nothing, nothing, z_idx)

   max_b =  max(maximum(bf), 5e-9)
   lim_b =  (3/4) * [-max_b, max_b]
   max_u =  max(maximum(uf_xy), 1e-10)
   lim_u =  (3/4) * [-max_u, max_u]
   max_w =  max(maximum(wf_xy), 1e-10)
   lim_w =  (3/4) * [-max_w, max_w] 
   max_q =  max(maximum(qf), 1e-15)
   lim_q =  (3/4) * [-max_q, max_q]

   fig1 = Figure(size = (1200, 800))
   axis_kwargs_xy = (xlabel = "x [m]", ylabel = "y [m]")
   ax_b = Axis(fig1[2, 1]; 
	       title = "Buoyancy perturbation (b')", axis_kwargs_xy...)
   ax_w = Axis(fig1[2, 3]; 
	       title = "Vertical velocity (w)", axis_kwargs_xy...)
   ax_u = Axis(fig1[3, 1]; 
	       title = "Zonal velocity (u)", axis_kwargs_xy...)
   ax_v = Axis(fig1[3, 3]; 
	       title = "Meridional velocity (v)", axis_kwargs_xy...)

   hm_b = heatmap!(ax_b, x, y, b_xy, colorrange = lim_b, colormap = :balance)
   Colorbar(fig1[2, 2], hm_b, tickformat = "{:.1e}", label = "m/s²")
   hm_w = heatmap!(ax_w, x, y, w_xy, colorrange = lim_w, colormap = :balance)
   Colorbar(fig1[2, 4], hm_w, tickformat = "{:.1e}", label = "m/s")
   hm_u = heatmap!(ax_u, x, y, u_xy, colorrange = lim_u, colormap = :balance)
   Colorbar(fig1[3, 2], hm_u, tickformat = "{:.1e}", label = "m/s")
   hm_v = heatmap!(ax_v, x, y, v_xy, colorrange = lim_u, colormap = :balance)
   Colorbar(fig1[3, 4], hm_v, tickformat = "{:.1e}", label = "m/s")

   title1 = @lift @sprintf("Fields at depth %i m; t = %.2f days", 
			   depth_nearest_m, times[$n]/(3600*24))
   fig1[1, 1:4] = Label(fig1, title1, fontsize = 24, tellwidth = false)

   fig2 = Figure(size = (600, 600))
   ax_q = Axis(fig2[2, 1]; axis_kwargs_xy...)
    
   hm_q = heatmap!(ax_q, x, y, q, colorrange = lim_q, colormap = :balance)
   Colorbar(fig2[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title2 = @lift @sprintf("Potential vorticity at depth %i m; t = %.2f days",
			   depth_nearest_m, times[$n]/(3600*24))
   fig2[1, 1:2] = Label(fig2, title2, fontsize = 24, tellwidth = false)

   frames = 1:Nt

   video1 = VideoStream(fig1, format = "mp4", framerate = 6)
   video2 = VideoStream(fig2, format = "mp4", framerate = 6)

   for i=1:frames[end]
      recordframe!(video1)
      recordframe!(video2)
      msg = string("Plotting frame(s) ", i, " of ", frames[end])
         print(msg * " \r")
         n[]=i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "bwuv_z$(depth_nearest_m)_$(datetime).mp4"), 
	video1)
   save(joinpath("./Plots", "q_z$(depth_nearest_m)_$(datetime).mp4"),
        video2)

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
