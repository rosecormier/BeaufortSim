include("Library.jl")

using Oceananigans
using CairoMakie, NCDatasets, Printf
using .VisualizationFunctions

const f = 2 * (7.2921 * 1e-5) * sin(74 * pi/180) #s^-1

const Lx = 2000 * 1e3 #m
const Ly = 2000 * 1e3 #m
const Lz = 1000       #m

const datetime = ARGS[1]
output_filename = joinpath("./Output", "output_$datetime.nc")
ds = NCDataset(output_filename, "r")

const depth_idx = parse(Int, ARGS[2])

x  =   ds["xC"][:]
y  =   ds["yC"][:]
z  =   ds["zC"][depth_idx]
depth_nearest_m = round(Int, -z)

Nx = length(ds["xC"][:])
Ny = length(ds["yC"][:])
Nz = length(ds["zC"][:])

Δx = Lx / Nx
Δy = Ly / Ny
Δz = Lz / Nz

times = ds["time"][:]

n = Observable(1)

#Use initial frame to define background fields
bb = ds["b"][:, :, depth_idx, 1]
ub = ds["u"][:, :, depth_idx, 1]
vb = ds["v"][:, :, depth_idx, 1]
wb = ds["w"][:, :, depth_idx, 1]

b    = @lift ds["b"][:, :, depth_idx, $n] .- bb
btot = @lift ds["b"][:, :, depth_idx, $n]
u    = @lift ds["u"][:, :, :, $n]
v    = @lift ds["v"][:, :, :, $n]
w    = @lift ds["w"][:, :, :, $n]

ωtotal = @lift ζ_2D($u, $v, $w, Δx, Δy, Δz, depth_idx)
∇_b    = @lift ∇b_2D($btot, Δx, Δy, Δz)
fq     = @lift @. f*($ωtotal + f) * $∇_b

u_xy = @lift $u[:, :, depth_idx]
v_xy = @lift $v[:, :, depth_idx]
w_xy = @lift $w[:, :, depth_idx]

#Use final frames to define maximum values
bf    = ds["b"][:, :, depth_idx, length(times)]
uf_xy = ds["u"][:, :, depth_idx, length(times)]
wf_xy = ds["w"][:, :, depth_idx, length(times)]

max_b  =  max(maximum(bf), 5e-9)
lim_b  =  (3/4) * [-max_b, max_b]
max_u  =  max(maximum(uf_xy), 1e-10)
lim_u  =  (3/4) * [-max_u, max_u]
max_w  =  (1/4) * max(maximum(wf_xy), 1e-10)
lim_w  =  (3/4) * [-max_w, max_w] 
lim_fq =  1 .* [-0.6e-13, 1e-13/5]

cm = [Makie.to_colormap(Reverse(:roma))[1:1:128];ones(2,1).*RGBAf(1.0,1.0,1.0,1.0); Makie.to_colormap(Reverse(:roma))[214:254]]
cm = cm[:,1];

fig1 = Figure(size = (1200, 1200))
axis_kwargs_xy = (xlabel = "x [m]", ylabel = "y [m]")
ax_b    = Axis(fig1[2, 1]; title = "Buoyancy perturbation", axis_kwargs_xy...)
ax_w    = Axis(fig1[2, 3]; title = "w", axis_kwargs_xy...)
ax_u    = Axis(fig1[3, 1]; title = "u", axis_kwargs_xy...)
ax_v    = Axis(fig1[3, 3]; title = "v", axis_kwargs_xy...)

hm_b = heatmap!(ax_b, x, y, b, colorrange = lim_b, colormap = :balance)
Colorbar(fig1[2, 2], hm_b, tickformat = "{:.1e}", label = "m/s²")
hm_w = heatmap!(ax_w, x, y, w_xy, colorrange = lim_w, colormap = :balance)
Colorbar(fig1[2, 4], hm_w, tickformat = "{:.1e}", label = "m/s")
hm_u = heatmap!(ax_u, x, y, u_xy, colorrange = lim_u, colormap = :balance)
Colorbar(fig1[3, 2], hm_u, tickformat = "{:.1e}", label = "m/s")
hm_v = heatmap!(ax_v, x, y, v_xy, colorrange = lim_u, colormap = :balance)
Colorbar(fig1[3, 4], hm_v, tickformat = "{:.1e}", label = "m/s")

fig2 = Figure(size = (600, 600))
ax_fq = Axis(fig2[2, 1]; title = "fq", axis_kwargs_xy...)
hm_fq = heatmap!(ax_fq, x, y, fq, colorrange = lim_fq, colormap = cm)
Colorbar(fig2[2, 2], hm_fq, tickformat = "{:.1e}", label = "1/s³kg")

title1 = @lift @sprintf("Fields at depth %i m; t = %.2f days", depth_nearest_m, times[$n]/(3600*24))
fig1[1, 1:4] = Label(fig1, title1, fontsize = 24, tellwidth = false)

title2 = @lift @sprintf("Potential vorticity at depth %i m; t = %.2f days", depth_nearest_m, times[$n]/(3600*24))
fig2[1, 1:2] = Label(fig2, title2, fontsize = 24, tellwidth = false)

frames = 1:length(times)

video1 = VideoStream(fig1, format = "mp4", framerate = 6)
video2 = VideoStream(fig2, format = "mp4", framerate = 6)

for i=1:frames[end]
    recordframe!(video1)
    recordframe!(video2)
    msg = string("Plotting frame(s) ", i, " of ", frames[end])
        print(msg * " \r")
        n[]=i
end
save("bwuv_z$(depth_nearest_m)_$(datetime).mp4", video1)
save("fq_z$(depth_nearest_m)_$(datetime).mp4", video2)

close(ds)
