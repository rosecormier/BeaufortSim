include("Library.jl")

using Oceananigans
using CairoMakie, NCDatasets, Printf
using .VisualizationFunctions

#const f     = 0.864e-4
#const fₕ    = 0.0
#const N²    = (3.7e-3)^2
#const ν     = 0.36*1.5/2.2e4
#const Umax  = 0.36*1.5
const Lx     = 2000 * 1e3
const Ly     = 2000 * 1e3
#const Lⱼ    = 3000    #1.5x to compensate for Umax.
const Lz     = 1000
#const D     = 200
#const z0    = -Lz/2
#const y0    = 0 # center of the jet in y

const datetime = ARGS[1]
output_filename = joinpath("./Output", "output_$datetime.nc")
ds = NCDataset(output_filename, "r")

const depth_idx = parse(Int, ARGS[2])

x  =   ds["xC"][:]
y  =   ds["yC"][:]
z  =   ds["zC"][depth_idx]

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

ωtotal = @lift ζ_2D($u, $v, $w, Δx, Δy, Δz)
#∇_b    = @lift ∇b_2D($btot, Δy, Δz)
#fq     = @lift @. f*($ωtotal + f) .* $∇_b

u_xy = @lift $u[:, :, depth_idx]
v_xy = @lift $v[:, :, depth_idx]
w_xy = @lift $w[:, :, depth_idx]

max_b  =  @lift max(maximum($b), 5e-9)
lim_b  =  @lift 3/4* [-$max_b,$max_b]
max_v  = @lift max(maximum($v_xy), 1e-10)
lim_v  =  @lift 3/4* [-$max_v,$max_v]
max_w  = @lift 1/4* max(maximum($w_xy), 1e-10)
lim_w  =  @lift 3/4* [-$max_w,$max_w]
lim_fq = 1 .* [-0.6e-13, 1e-13/5]

cm = [Makie.to_colormap(Reverse(:roma))[1:1:128];ones(2,1).*RGBAf(1.0,1.0,1.0,1.0); Makie.to_colormap(Reverse(:roma))[214:254]]
cm = cm[:,1];

fig1=Figure(size = (1200, 1200))
axis_kwargs_xy = (xlabel = "x", ylabel = "y")
ax_b    = Axis(fig1[2, 1]; title = "b'", axis_kwargs_xy...)
ax_u    = Axis(fig1[2, 3]; title = "u", axis_kwargs_xy...)
ax_v    = Axis(fig1[3, 1]; title = "v", axis_kwargs_xy...)
ax_w    = Axis(fig1[3, 3]; title = "w", axis_kwargs_xy...)
#ax_fq   = Axis(fig1[3, 3]; title = "fq", axis_kwargs_xy...)

hm_b = heatmap!(ax_b, x, y, b, colorrange=lim_b, colormap = :balance)
Colorbar(fig1[2, 2], hm_b, tickformat= "{:.1e}")
hm_u = heatmap!(ax_u, x, y, u_xy, colorrange=lim_v, colormap = :balance)
Colorbar(fig1[2, 4], hm_u, tickformat= "{:.1e}")
hm_v = heatmap!(ax_v, x, y, v_xy, colorrange=lim_v, colormap = :balance)
Colorbar(fig1[3, 2], hm_v, tickformat= "{:.1e}")
hm_w = heatmap!(ax_w, x, y, w_xy, colorrange=lim_w, colormap = :balance)
Colorbar(fig1[3, 4], hm_w, tickformat= "{:.1e}") #2], hm_w, tickformat= "{:.1e}")
#hm_fq = heatmap!(ax_fq, y, z, fq, colorrange=lim_fq ,colormap = cm)
#Colorbar(fig1[3, 4], hm_fq, tickformat= "{:.1e}")

title =  @lift @sprintf("t = %.2f days", times[$n]/(3600*24))
fig1[1, 1:4] = Label(fig1, title, fontsize=24, tellwidth=false)
frames = 1:length(times)
video1=VideoStream(fig1,format="mp4", framerate=12)
for i=1:frames[end]
    recordframe!(video1)
    msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[]=i
end
close(ds)
save("zoomPlots_fq.mp4", video1)
