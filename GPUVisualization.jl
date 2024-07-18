include("Library.jl")

using Oceananigans
using NCDatasets, Printf, CairoMakie
using .VisualizationFunctions

const f     = 0.864e-4
const fₕ    = 0.0
const N²    = (3.7e-3)^2
const ν     = 0.36*1.5/2.2e4
const Umax  = 0.36*1.5
const Lx     = 20000
const Ly     = 30000
const Lⱼ    = 3000    #1.5x to compensate for Umax.
const Lz    = 1000
const D     = 200
const z0    = -Lz/2
const y0    = 0 # center of the jet in y

const datetime = ARGS[1]
print(datetime)
ds = NCDataset("output_$datetime.nc", "r")
x  =   ds["xC"][:]
y  =   ds["yC"][:]
z  =   ds["zC"][:]
Nx = length(ds["xC"][:])
Ny = length(ds["yC"][:])
Nz = length(ds["zC"][:])
times = ds["time"][:]

# imports the parameters directly from the model if running inside the model,
# assigns a value to them if this script is running by itself.

   Δx = 2*Lx/Nx
   Δy = 2*Ly/Ny
   Δz = Lz/Nz
# pinpoint jet position
jp_c = (y0 + Ly) ./ Δy
jp_ini = floor(Int,  jp_c - Lⱼ ./ Δy) - 1125
jp_end = floor(Int, jp_c + Lⱼ ./ Δy) -125

z_jet = -z0 ./ Δz
z_jet_ini = floor(Int, z_jet - D/ Δz)
z_jet_end = floor(Int, z_jet + D/ Δz)

n=Observable(1)
# get variables

bb = ds["b"][1 , jp_ini : jp_end, z_jet_ini : z_jet_end, 1]
ub = ds["u"][: , jp_ini : jp_end, z_jet_ini : z_jet_end, 1]
vb = ds["v"][: , jp_ini : jp_end, z_jet_ini : z_jet_end, 1]
wb = ds["w"][: , jp_ini : jp_end, z_jet_ini : z_jet_end, 1]

b = @lift ds["b"][1 , jp_ini : jp_end, z_jet_ini : z_jet_end, $n] .- bb
bt = @lift ds["b"][1 , jp_ini : jp_end, z_jet_ini : z_jet_end, $n]
u = @lift ds["u"][: , jp_ini : jp_end, z_jet_ini : z_jet_end, $n]
v = @lift ds["v"][: , jp_ini : jp_end, z_jet_ini : z_jet_end, $n]
w = @lift ds["w"][: , jp_ini : jp_end, z_jet_ini : z_jet_end, $n]

ωtotal = @lift ζ_2D($u, $v, $w, Δy, Δz)
∇_b= @lift ∇b_2D($bt, Δy, Δz)
fq = @lift @. f*($ωtotal + f) .* $∇_b

y  =   y[jp_ini: jp_end]
z  =   z[z_jet_ini : z_jet_end]

u_yz = @lift $u[1 , :, :]
v_yz = @lift $v[1 , :, :]
w_yz = @lift $w[1 , :, :]

# plot data

max_b  =  @lift max(maximum($b), 5e-9)
lim_b  =  @lift 3/4* [-$max_b,$max_b]
max_v  =  @lift max(maximum($v_yz), 1e-10)
lim_v  =  @lift 3/4* [-$max_v,$max_v]
max_w  =  @lift 1/4* max(maximum($w_yz), 1e-10)
lim_w  =  @lift 3/4* [-$max_w,$max_w]
lim_fq = 1 .* [-0.6e-13, 1e-13/5]


cm = [Makie.to_colormap(Reverse(:roma))[1:1:128];ones(2,1).*RGBAf(1.0,1.0,1.0,1.0); Makie.to_colormap(Reverse(:roma))[214:254]]
cm = cm[:,1];


fig1=Figure(resolution = (1200, 1200))
axis_kwargs_yz = (xlabel = "y", ylabel = "z")
ax_b    = Axis(fig1[2, 1]; title = "b'", axis_kwargs_yz...)
ax_v    = Axis(fig1[2, 3]; title = "v", axis_kwargs_yz...)
ax_w    = Axis(fig1[3, 1]; title = "w", axis_kwargs_yz...)
ax_fq   = Axis(fig1[3, 3]; title = "fq", axis_kwargs_yz...)

hm_b = heatmap!(ax_b, y, z, b, colorrange=lim_b ,colormap = :balance)
Colorbar(fig1[2, 2], hm_b, tickformat= "{:.1e}")
hm_v = heatmap!(ax_v, y, z, v_yz, colorrange=lim_v ,colormap = :balance)
Colorbar(fig1[2, 4], hm_v, tickformat= "{:.1e}")
hm_w = heatmap!(ax_w, y, z, w_yz, colorrange=lim_w ,colormap = :balance)
Colorbar(fig1[3, 2], hm_w, tickformat= "{:.1e}")
hm_fq = heatmap!(ax_fq, y, z, fq, colorrange=lim_fq ,colormap = cm)
Colorbar(fig1[3, 4], hm_fq, tickformat= "{:.1e}")

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
