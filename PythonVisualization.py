from glob import glob
from netCDF4 import Dataset
from numpy.linalg import norm
from os.path import join

import numpy as np
import matplotlib.pyplot as plt

datetime = "250706-143810"

bkgd_file = fr"./Output/bkgd_{datetime}.nc"
openbkgd  = Dataset(bkgd_file, "r")
uφ_bkgd   = openbkgd.variables["Uφ"][0, :, :, :]

openbkgd.close()

all_times      = []
ur_norms_log   = []
uphi_norms_log = []
uz_norms_log   = []

ur_data   = []
uphi_data = []
uz_data   = []

file_list = glob(fr"./Output/output_{datetime}*.nc")

for filepath in file_list:

    openfile = Dataset(filepath, "r")
     
    x     = openfile.variables["xC"][:]
    y     = openfile.variables["yC"][:]
    z     = openfile.variables["zC"][:]
    times = openfile.variables["time"][:]
    ur    = openfile.variables["ur"][:, :, :, :]
    uphi  = openfile.variables["uφ"][:, :, :, :]
    uz    = openfile.variables["uz"][:, :, :, :]

    for t in range(len(times)):
    
       ur_norm_t   = norm(ur[t, :, :, :])
       uphi_norm_t = norm(uphi[t, :, :, :] - uφ_bkgd)
       uz_norm_t   = norm(uz[t, :, :, :])

       all_times.append(times[t])
       ur_norms_log.append(np.log(ur_norm_t))
       uphi_norms_log.append(np.log(uphi_norm_t))
       uz_norms_log.append(np.log(uz_norm_t))

       ur_data.append(ur[t, :, :, :])
       uphi_data.append(uphi[t, :, :, :])
       uz_data.append(uz[t, :, :, :])

    openfile.close()

all_times      = np.array(all_times)

ur_norms_log   = np.array(ur_norms_log)
uphi_norms_log = np.array(uphi_norms_log)
uz_norms_log   = np.array(uz_norms_log)

ur_data   = np.array(ur_data)
uphi_data = np.array(uphi_data)
uz_data   = np.array(uz_data)

time_indexing = np.argsort(all_times)

fit_idx_start, fit_idx_stop = 55, 140

times              = (all_times[time_indexing])[fit_idx_start:fit_idx_stop]
ur_norms_log_fit   = (ur_norms_log[time_indexing])[fit_idx_start:fit_idx_stop]
uphi_norms_log_fit = (uphi_norms_log[time_indexing])[fit_idx_start:fit_idx_stop]
uz_norms_log_fit   = (uz_norms_log[time_indexing])[fit_idx_start:fit_idx_stop]

ur_fit_params   = np.polyfit(times, ur_norms_log_fit, 1)
uphi_fit_params = np.polyfit(times, uphi_norms_log_fit, 1)
uz_fit_params   = np.polyfit(times, uz_norms_log_fit, 1)

print("Best-fit parameters for log|u_r|: [slope, intercept] = ", ur_fit_params, "\n")
print("Best-fit parameters for log|u_phi|: [slope, intercept] = ", uphi_fit_params, "\n")
print("Best-fit parameters for log|u_z|: [slope, intercept] = ", uz_fit_params)

fig, ax = plt.subplots(1, 1)
ax.scatter(all_times / 86400, ur_norms_log, color = "red", 
           label = "$||u_r'||$")
ax.scatter(all_times / 86400, uphi_norms_log, color = "orange", 
           label = "$||u_{{\phi}}'||$")
ax.scatter(all_times / 86400, uz_norms_log, color = "green",
           label = "$||u_z'||$")
ax.plot(np.array(times) / 86400, 
        ur_fit_params[0] * np.array(times) + ur_fit_params[1], 
        color = "mediumblue", label = "Best fit to $||u_r'||$")
ax.plot(np.array(times) / 86400, 
        uphi_fit_params[0] * np.array(times) + uphi_fit_params[1], 
        color = "mediumpurple", label = "Best fit to $||u_{{\phi}}'||$")
ax.plot(np.array(times) / 86400,
        uz_fit_params[0] * np.array(times) + uz_fit_params[1],
        color = "black", label = "Best fit to $||u_z'||$")
ax.set(xlabel = "Time [days]", ylabel = "Log of $\ell^2$-norm [m/s]")
ax.legend()
plt.grid(True)
#plt.show()
fig.savefig(join("Plots",
            f"norms_{datetime}_fit_{fit_idx_start}-{fit_idx_stop}.png"))
plt.close(fig)

ur_data   = ur_data[time_indexing]
uphi_data = uphi_data[time_indexing]
uz_data   = uz_data[time_indexing]

for k in range(10, 11, 1): #len(z)-3, 1):

    for t in range(fit_idx_start, fit_idx_stop, 3):
        
        time = (all_times[time_indexing])[t]
        zval = z[k]
        
        fig, axs = plt.subplots(1, 2, sharey = True)
        pcm_ur = axs[0].pcolormesh(x[200:606] / 1e3, y[200:606] / 1e3,
                                     ur_data[t, k, 200:606, 200:606],
                                     cmap = "seismic")
        axs[0].set(xlabel = "x [km]", ylabel = "y [km]")
        plt.colorbar(pcm_ur, ax = axs[0])
        pcm_uphi = axs[1].pcolormesh(x[200:606] / 1e3, y[200:606] / 1e3,
                            (uphi_data[t, k, 200:606, 200:606]
                             - uφ_bkgd[k, 200:606, 200:606]),
                            cmap = "seismic")
        axs[1].set(xlabel = "x [km]", ylabel = "y [km]")
        plt.colorbar(pcm_uphi, ax = axs[1])
        fig.suptitle(f"Perturbations in horizontal velocity components; \n t = {time/86400:.2f} days, z = {zval} m")
        #plt.show()
        fig.savefig(join("Plots",
                    f"ur_uphi_perturb_{datetime}_t{str(t).zfill(5)}_z{zval}.png"))
        plt.close(fig)

