from glob import glob
from netCDF4 import Dataset
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt

bkgd_file = r"./Output/bkgd_250706-143810.nc"
openbkgd  = Dataset(bkgd_file, "r")
uφ_bkgd   = openbkgd.variables["Uφ"][0, :, :, :]

openbkgd.close()

all_times      = []
ur_norms_log   = []
uphi_norms_log = []

file_list = glob(r"./Output/output_250706-143810*.nc")

for filepath in file_list:

    openfile = Dataset(filepath, "r")
  
    times = openfile.variables["time"][:]
    ur    = openfile.variables["ur"][:, :, :, :]
    uphi  = openfile.variables["uφ"][:, :, :, :]

    for t in range(len(times)):
       ur_norm_t   = norm(ur[t, :, :, :])
       uphi_norm_t = norm(uphi[t, :, :, :] - uφ_bkgd)

       all_times.append(times[t])
       ur_norms_log.append(np.log(ur_norm_t))
       uphi_norms_log.append(np.log(uphi_norm_t))

    openfile.close()

all_times      = np.array(all_times)
ur_norms_log   = np.array(ur_norms_log)
uphi_norms_log = np.array(uphi_norms_log)

time_indexing = np.argsort(all_times)

times              = (all_times[time_indexing])[60:-160]
ur_norms_log_fit   = (ur_norms_log[time_indexing])[60:-160]
uphi_norms_log_fit = (uphi_norms_log[time_indexing])[60:-160]

ur_fit_params   = np.polyfit(times, ur_norms_log_fit, 1)
uphi_fit_params = np.polyfit(times, uphi_norms_log_fit, 1)

print("Best-fit parameters for log|u_r|: [slope, intercept] = ", ur_fit_params, "\n")
print("Best-fit parameters for log|u_phi|: [slope, intercept] = ", uphi_fit_params, "\n")

plt.figure()
plt.scatter(all_times, ur_norms_log)
plt.scatter(all_times, uphi_norms_log, color="red")
plt.plot(np.array(times), ur_fit_params[0] * np.array(times) + ur_fit_params[1], color = "orange")
plt.plot(np.array(times), uphi_fit_params[0] * np.array(times) + uphi_fit_params[1], color = "red")
plt.show()
