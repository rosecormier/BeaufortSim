"""
Visualization of fields from Beaufort-Gyre simulation:
    -Vertical vorticity component
    -Buoyancy
    -Buoyancy perturbation
    -Velocity components
    
R. Cormier, 2024
"""

import argparse
import matplotlib
import numpy as np
import os
import pandas as pd
import xarray as xr

matplotlib.use("Agg")

from math import e
from os.path import join

import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib.pyplot as plt

####################

plt.rcParams["font.size"] = 16
plt.rcParams["axes.titlesize"] = 18
plt.rcParams["text.usetex"] = True

####################

#SET UP PARSER AND GET COMMAND-LINE ARGUMENTS

parser = argparse.ArgumentParser()
parser.add_argument("label", type=str,
                    help="Label (yymmdd-HHMMSS) of specific output file")
parser.add_argument("-slice", type=int, default=1,
                    help="Timeseries index slice (>1 speeds up plotting)")
parser.add_argument("-depth", type=int, default=1,
                    help="Depth (i.e., zC) index at which to plot data")
args = parser.parse_args()

file_label, slice_len, depth_idx = args.label, args.slice, args.depth

####################

#SET UP DIRECTORY INFORMATION

output_dir = "Output"
vis_dir = "Plots"

if not os.path.exists(join(vis_dir)):
    os.makedirs(join(vis_dir)) #Make vis directory if nonexistent
    
output_filepath = join(output_dir, "output_{}.nc".format(file_label))

####################

#FUNCTION DEFINITIONS

def load_data(filename, slice_len):
    sim_ds = xr.open_dataset(filename)
    xvort_data, yvort_data = sim_ds["ωx"], sim_ds["ωy"]
    zvort_data = sim_ds["ωz"]
    buoyancy_data, b_perturb_data = sim_ds["b"], sim_ds["b_perturb"]
    u_data, v_data, w_data = sim_ds["u"], sim_ds["v"], sim_ds["w"]
    C_grid = xr.Dataset(coords={"xC": ("xC", sim_ds["xC"].data), 
                                "yC": ("yC", sim_ds["yC"].data), 
                                "zC": ("zC", sim_ds["zC"].data), 
                                "time": ("time", sim_ds["time"].data)})
    C_grid = C_grid.assign(xvorticity=xvort_data.interp(yF=C_grid["yC"], 
                                            zF=C_grid["zC"]).squeeze(), 
                        yvorticity=yvort_data.interp(xF=C_grid["xC"], 
                                            zF=C_grid["zC"]).squeeze(),
                        zvorticity=zvort_data.interp(xF=C_grid["xC"], 
                                            yF=C_grid["yC"]).squeeze(), 
                        buoyancy=buoyancy_data.squeeze(),
                        b_perturb=b_perturb_data.squeeze(),
                        u=u_data.interp(xF=C_grid["xC"]).squeeze(), 
                        v=v_data.interp(yF=C_grid["yC"]).squeeze(), 
                        w=w_data.interp(zF=C_grid["zC"]).squeeze())
    times = C_grid["time"]
    time_iter = np.arange(len(times)-2) #Index -1 may contain NaN
    time_iter = time_iter[::slice_len] #Slice iterable to speed up visualization
    return C_grid, time_iter

def get_title_time(np_timedelta_obj):
    pd_timedelta_obj = pd.to_timedelta(np_timedelta_obj)
    day = pd_timedelta_obj.days + 1
    total_seconds = pd_timedelta_obj.seconds
    hours, remaining_seconds = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remaining_seconds, 60)
    title_time = "Day {}, t={:02}:{:02}:{:02}".format(day, 
                                                      hours, minutes, seconds)
    return title_time

def animate_zvorticity(time, C_grid, vmax, depth_str=""): 
    time_title = get_title_time(C_grid.time[time].values)
    fig.suptitle("z-Vorticity{}; {}".format(depth_str, time_title))
    frame_vorticity = C_grid["zvorticity"].isel(zC=depth_idx, time=time)
    frame_vorticity.drop_sel(xC=C_grid.xC[-1], yC=C_grid.yC[-1]) #Remove NaNs
    pcm = ax.pcolormesh(C_grid["xC"]*1e-3, C_grid["yC"]*1e-3, 
                        frame_vorticity.values, 
                        cmap="PuOr_r", vmin=-vmax, vmax=vmax)
    return pcm

def animate_velocity_uv_comps(time, C_grid, vmax, depth_str=""):
    time_title = get_title_time(C_grid.time[time].values)
    fig.suptitle("Horizontal Velocity Components{}; {}".format(depth_str,
                                                                 time_title))
    frame_u = C_grid["u"].isel(zC=depth_idx, time=time)
    frame_v = C_grid["v"].isel(zC=depth_idx, time=time)
    frame_u.drop_sel(xC=C_grid.xC[-1], yC=C_grid.yC[-1]) #Remove NaNs
    frame_v.drop_sel(xC=C_grid.xC[-1], yC=C_grid.yC[-1]) #Remove NaNs
    pcm1 = ax1.pcolormesh(C_grid["xC"]*1e-3, C_grid["yC"]*1e-3, 
                          frame_u.values, 
                          cmap="seismic",
                          vmin=-vmax, vmax=vmax)
    pcm2 = ax2.pcolormesh(C_grid["xC"]*1e-3, C_grid["yC"]*1e-3, 
                          frame_v.values, 
                          cmap="seismic", vmin=-vmax, vmax=vmax)
    return pcm1, pcm2

def animate_velocity_w_comp(time, C_grid, vmax, depth_str=""):
    time_title = get_title_time(C_grid.time[time].values)
    fig.suptitle("Vertical Velocity Component{}; {}".format(depth_str, 
                                                              time_title))
    frame_w = C_grid["w"].isel(zC=depth_idx, time=time)
    frame_w.drop_sel(xC=C_grid.xC[-1], yC=C_grid.yC[-1]) #Remove NaNs
    pcm = ax.pcolormesh(C_grid["xC"]*1e-3, C_grid["yC"]*1e-3, 
                        frame_w.values, 
                        cmap="seismic", vmin=-vmax, vmax=vmax)
    return pcm

def animate_buoyancy(time, C_grid, vmax, depth_str=""):
    time_title = get_title_time(C_grid.time[time].values)
    fig.suptitle("Buoyancy{}; {}".format(depth_str, time_title))
    frame_buoyancy = C_grid["buoyancy"].isel(zC=depth_idx, time=time)
    frame_buoyancy.drop_sel(xC=C_grid.xC[-1], yC=C_grid.yC[-1]) #Remove NaNs
    pcm = ax.pcolormesh(C_grid["xC"]*1e-3, C_grid["yC"]*1e-3, 
                        frame_buoyancy.values, 
                        cmap="Oranges", vmin=0, vmax=vmax)
    return pcm

def animate_b_perturbation(time, C_grid, vmax, depth_str=""):
    time_title = get_title_time(C_grid.time[time].values)
    fig.suptitle("Perturbation in buoyancy field{}; {}".format(depth_str, 
                                                                 time_title))
    frame_perturbation = C_grid["b_perturb"].isel(zC=depth_idx, time=time)
    frame_perturbation.drop_sel(xC=C_grid.xC[-1], yC=C_grid.yC[-1]) #Remove NaNs
    pcm = ax.pcolormesh(C_grid["xC"]*1e-3, C_grid["yC"]*1e-3,
                        frame_perturbation.values,
                        cmap="BrBG_r", vmin=-vmax, vmax=vmax)
    return pcm

####################

#LOAD SIMULATION DATA

C_grid, time_iter = load_data(output_filepath, slice_len)

#For plot titles
depth_title_str = " at {}m depth".format(-C_grid.zC[depth_idx].values) 

####################

#CREATE ANIMATIONS

#Plot z-vorticity
fig, ax = plt.subplots(figsize=(10,8))
ax.set_xlabel(r"x ($km$)")
ax.set_ylabel(r"y ($km$)")
vmax = np.max(abs(C_grid["zvorticity"].isel(zC=depth_idx, time=0)))
fig.colorbar(animate_zvorticity(0, C_grid, vmax), extend="both", 
             label=r"$s^{-1}$")
anim = animation.FuncAnimation(fig, animate_zvorticity, 
                               fargs=(C_grid, vmax, depth_title_str), 
                               frames=time_iter)
anim.save(join(vis_dir, "zvorticity_{}.gif".format(file_label)), 
          progress_callback=lambda i, n: print(f'saving frame {i} of {n}'))
plt.close()

#Plot buoyancy
fig, ax = plt.subplots(figsize=(10,8))
ax.set_xlabel(r"x ($km$)")
ax.set_ylabel(r"y ($km$)")
vmax = np.max(abs(C_grid["buoyancy"].isel(zC=depth_idx, time=0)))
fig.colorbar(animate_buoyancy(0, C_grid, vmax), extend="max", label=r"$m/s^2$")
anim = animation.FuncAnimation(fig, animate_buoyancy, 
                               fargs=(C_grid, vmax, depth_title_str), 
                               frames=time_iter)
anim.save(join(vis_dir, "buoyancy_{}.gif".format(file_label)), 
          progress_callback=lambda i, n: print(f'saving frame {i} of {n}'))
plt.close()

#Plot buoyancy perturbation
fig, ax = plt.subplots(figsize=(10,8))
ax.set_xlabel(r"x ($km$)")
ax.set_ylabel(r"y ($km$)")
vmax = np.max(abs(C_grid["b_perturb"].isel(zC=depth_idx, time=-1)))
fig.colorbar(animate_b_perturbation(0, C_grid, vmax), extend="both", 
             label=r"$m/s^2$")
anim = animation.FuncAnimation(fig, animate_b_perturbation, 
                               fargs=(C_grid, vmax, depth_title_str), 
                               frames=time_iter)
anim.save(join(vis_dir, "b_perturb_{}.gif".format(file_label)), 
          progress_callback=lambda i, n: print(f'saving frame {i} of {n}'))
plt.close()

#Plot horizontal velocity components
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14,8))
ax1.set_xlabel(r"x ($km$)")
ax2.set_xlabel(r"x ($km$)")
ax1.set_ylabel(r"y ($km$)")
ax1.set_title(r"$u$-Component")
ax2.set_title(r"$v$-Component")
vmax_u = np.max(abs(C_grid["u"].isel(zC=depth_idx, time=0)))
vmax_v = np.max(abs(C_grid["v"].isel(zC=depth_idx, time=0)))
vmax = max(vmax_u, vmax_v)
fig.colorbar(animate_velocity_uv_comps(0, C_grid, vmax)[0], ax=[ax1, ax2], 
             extend="both", label=r"$m/s$", location="bottom", shrink=0.5)
anim = animation.FuncAnimation(fig, animate_velocity_uv_comps, 
                               fargs=(C_grid, vmax, depth_title_str), 
                               frames=time_iter)
anim.save(join(vis_dir, "velocity_uv_{}.gif".format(file_label)), 
          progress_callback=lambda i, n: print(f'saving frame {i} of {n}'))
plt.close()

#Plot vertical velocity component
fig, ax = plt.subplots(figsize=(10,8))
ax.set_xlabel(r"x ($km$)")
ax.set_ylabel(r"y ($km$)")
vmax = np.max(abs(C_grid["w"].isel(zC=depth_idx, time=0)))
fig.colorbar(animate_velocity_w_comp(0, C_grid, vmax, depth_title_str), 
             extend="both", label=r"$m/s$")
anim = animation.FuncAnimation(fig, animate_velocity_w_comp, 
                               fargs=(C_grid, vmax, depth_title_str), 
                               frames=time_iter)
anim.save(join(vis_dir, "velocity_w_{}.gif".format(file_label)), 
          progress_callback=lambda i, n: print(f'saving frame {i} of {n}'))
plt.close()
