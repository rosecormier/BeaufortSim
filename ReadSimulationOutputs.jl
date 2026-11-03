using Glob, NCDatasets
using Oceananigans.OutputReaders

#######################################
# FUNCTIONS TO HANDLE NETCDF DATASETS #
#######################################

function pad_filenames(datetime; prefix = "output")
   #=
   Prefix filenames with zeros until alphanumeric order of filenames is the same
    as chronological order of data intervals contained by the respective files.
   =#

   outfile_paths = glob("./Output/$(prefix)_$(datetime)_part*")

   if length(outfile_paths) > 9
      for file_path in outfile_paths
         if length(file_path) == 38
	          mv(file_path, replace(file_path, "_part" => "_part0"))
         end
      end
   end
   
   return glob("./Output/$(prefix)_$(datetime)*")
end

function sort_times(unsortedTimes)
   #=
   Ordered t-values and indices to sort other data chronologically (essential if
    simulation was ever picked up from a checkpoint).
   =#
   
   sorting_indices = sortperm(unsortedTimes)
   t               = unsortedTimes[sorting_indices]
   
   return t, sorting_indices
end

function load_times_in_days(dataset)

   unsorted_t               = dataset[:time][:] ./ 86400
   t, chronological_indices = sort_times(unsorted_t)
   Nt                       = length(t)
   
   return t, Nt, chronological_indices
end

function open_nc_dataset(outfilename; Hx = 3, Hy = 3, Hz = 3)
   #=
   Open a generic NetCDF dataset output by Oceananigans simulation.
   Return the dataset itself as well as gridpoints and timeseries info.   
   =#

   ds = NCDataset(outfilename)

   x, y, z = nothing, nothing, nothing #Defaults in case any dimension is Flat
   
   #Load coords of non-Flat dimensions; convert them to km for readability

   if length(ds[:x_caa]) > 1
      x = ds[:x_caa][:] ./ 1000
   end

   if length(ds[:y_aca]) > 1
      y = ds[:y_aca][:] ./1000
   end
   
   if length(ds[:z_aac]) > 1
      zC = ds[:z_aac][:] ./ 1000
      zF = ds[:z_aaf][:] ./ 1000
   end

   t, Nt, chronological_indices = load_times_in_days(ds)

   return ds, x, y, zC, zF, t, Nt, chronological_indices
end

function open_energetics_dataset(energeticsfilename)
   
   energetics_ds                = NCDataset(energeticsfilename)
   t, Nt, chronological_indices = load_times_in_days(energetics_ds)

   return energetics_ds, t, Nt, chronological_indices
end

function open_scalars_dataset(scalarfilename)

   scalars_ds                   = NCDataset(scalarfilename)
   t, Nt, chronological_indices = load_times_in_days(scalars_ds)

   return scalars_ds, t, Nt, chronological_indices
end

#####################################
# FUNCTIONS TO HANDLE JLD2 DATASETS #
#####################################

function LoadFieldTimeseries(field::String, runLabel::String)
   TS = FieldTimeSeries("./Output/output_$(runLabel)", field; backend = OnDisk())
end

function OpenJLD2Output(outfilename)

   ds = FieldDataset(outfilename; backend = OnDisk())

   ur_Timeseries = ds[:ur]
   uφ_Timeseries = ds[:uφ]
	 ux_Timeseries = ds[:ux]
	 uy_Timeseries = ds[:uy]
	 uz_Timeseries = ds[:uz]
	 b_Timeseries  = ds[:b]
   ζ_Timeseries  = ds[:ζ]
   
   t = ur_Timeseries.times

   return t, ur_Timeseries, uφ_Timeseries, ux_Timeseries, uy_Timeseries, uz_Timeseries, b_Timeseries, ζ_Timeseries
end