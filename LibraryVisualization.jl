using Adapt, CUDA, Glob, LinearAlgebra, NCDatasets
using Oceananigans.AbstractOperations, Oceananigans.Fields

###########################################################################
# HELPER FUNCTIONS TO PLOT 2D SLICES OF PRIMITIVE OR DIAGNOSTIC VARIABLES #
###########################################################################

function get_symm_range_lims(scalar_field; 
                             max_fraction = 1, prescribed_max = 1e-16)
   #=
   Given a scalar field, extract appropriate limits for a colorbar that is
    symmetric about 0.
   Scale the absolute maximum of the data by 'max_fraction'.
   If the absolute maximum of the data is less than 'prescribed_max', treat
    'prescribed_max' as the absolute maximum of the data instead.
   =#
   
   field_max  = max(maximum(abs.(scalar_field)), prescribed_max)
   field_lims = [-(max_fraction * field_max), (max_fraction * field_max)]
end

function get_2D_spatial_axis_idcs(const_dimension;
                                  Hx = 3, Hy = 3, Hz = 3,
		                              x_idx = nothing, 
                                  y_idx = nothing, 
                                  z_idx = nothing,
				                          xC = nothing, yC = nothing, zC = nothing,
                      				    zF = nothing)         
   #=
   Return two arrays of coordinate indices partitioning a 2D slice (constant
    along 'const_dimension').
   The first array can be used to project a field located at
    (Center, Center, Center)-points onto the 2D slice; the second array can be
    used to project a field located at (Center, Center, Face)-points onto the 2D
    slice.
   =#

   if const_dimension == "x"

      if isnothing(x_idx) #Grid is 2D with only y and z axes
         yCzC_idcs = (1, 
                      (Hy + 1):(length(yC) + Hy - 1), 
                      (Hz + 1):(length(zC) + Hz)
                     )
         yCzF_idcs = (1, 
                      (Hy + 1):(length(yC) + Hy), 
                      (Hz + 1):(length(zF) + Hz)
                     ) 
      else #Grid is 3D
         yCzC_idcs = (x_idx, 
                      (Hy + 1):(length(yC) + Hy), 
                      (Hz + 1):(length(zC) + Hz)
                     )
         yCzF_idcs = (x_idx, 
                      (Hy + 1):(length(yC) + Hy), 
                      (Hz + 1):(length(zF) + Hz)
                     )
      end

      return yCzC_idcs, yCzF_idcs

   elseif const_dimension == "y"

      if isnothing(y_idx) #Grid is 2D with only x and z axes
         xCzC_idcs = ((Hx + 1):(length(xC) + Hx), 
                      1, 
                      (Hz + 1):(length(zC) + Hz)
                     )
         xCzF_idcs = ((Hx + 1):(length(xC) + Hx), 
                      1, 
                      (Hz + 1):(length(zF) + Hz)
                     )
      else #Grid is 3D
         xCzC_idcs = ((Hx + 1):(length(xC) + Hx), 
                      y_idx, 
                      (Hz + 1):(length(zC) + Hz)
                     )
         xCzF_idcs = ((Hx + 1):(length(xC) + Hx), 
                      y_idx, 
                      (Hz + 1):(length(zF) + Hz)
                     )
      end

      return xCzC_idcs, xCzF_idcs

   elseif const_dimension == "z"

      if isnothing(z_idx) #Grid is 2D with only x and y axes
         xCyC_idcs = ((Hx + 1):(length(xC) + Hx), 
                      (Hy + 1):(length(yC) + Hy),
                      1
                     )
      else #Grid is 3D
         xCyC_idcs = ((Hx + 1):(length(xC) - (2 * Hx)), 
                      (Hy + 1):(length(yC) - (2 * Hy)), 
                      (z_idx - Hz)
                     )
      end

      return xCyC_idcs, xCyC_idcs #zC or zF is irrelevant on a constant-z slice
   end
end

#function get_2D_spatial_axis_kwargs(x, y, z, const_dimension_idx)#;
#                                    #x_idx = nothing,#
#				                            #y_idx = nothing,
#                                   #z_idx = nothing)

function get_2D_spatial_axis_kwargs(const_dimension, const_dimension_idx, const_dimension_coords)
   #=
   Return the rounded-off value of the constant dimension at the specified
    index.
   Also return a tuple of axis labels for plots on a 2D slice.
   =#

   if const_dimension == "x"
      nearest     = round(Int, const_dimension_coords[const_dimension_idx])
      axis_kwargs = (xlabel = "y [km]", ylabel = "z [km]")
   elseif const_dimension == "y"
      nearest     = round(Int, const_dimension_coords[const_dimension_idx])
      axis_kwargs = (xlabel = "x [km]", ylabel = "z [km]")
   elseif const_dimension == "z"
      nearest     = round(const_dimension_coords[const_dimension_idx],
                          digits = 2)
      axis_kwargs = (xlabel = "x [km]", ylabel = "y [km]")
   end

   #nearest = 0

   #=
   if !isnothing(x_idx)
      const_dim = "x"
      nearest   = round(Int, x[x_idx])
   elseif !isnothing(y_idx)
      const_dim = "y"
      nearest   = round(Int, y[y_idx])
   elseif !isnothing(z_idx)
      const_dim = "z"
      nearest   = round(z[z_idx], digits = 2)
   end

   if const_dim == "x"
      axis_kwargs = (xlabel = "y [km]", ylabel = "z [km]")
   elseif const_dim == "y"
      axis_kwargs = (xlabel = "x [km]", ylabel = "z [km]")
   elseif const_dim == "z"
      axis_kwargs = (xlabel = "x [km]", ylabel = "y [km]")
   end
   =#
   
   return nearest, axis_kwargs
end

#################################
# FINITE-DIFFERENCING FUNCTIONS #
#################################

function order1_forward_difference(t, u)
   #=
   Compute du/dt by first-order forward-differencing.
   =#
   
   return @. (u[2:end] - u[1:end-1]) / (t[2:end] - t[1:end-1])
end

function centered_difference(t, u)
   #=
   Compute du/dt by centered-differencing.
   =#

   u_i         = u[2:end-1]
   u_i_minus_1 = u[1:end-2]
   u_i_plus_1  = u[3:end]
   
   if length(t) > 3
      Delta_t_minus = t[2:end-1] .- t[1:end-2] # t[2:end-1] .- t[1:end-2]
      Delta_t_plus  = t[3:end] .- t[2:end-1] #t[3:end] .- t[2:end-1]
   elseif length(t) == 3
      Delta_t_minus = t[2] - t[1]
      Delta_t_plus  = t[3] - t[2]
   end
   
   A = -Delta_t_plus ./ Delta_t_minus
   B = (Delta_t_plus ./ Delta_t_minus) .- (Delta_t_minus ./ Delta_t_plus)
   C = Delta_t_minus ./ Delta_t_plus
   
   return ((A .* u_i_minus_1 .+ B .* u_i .+ C .* u_i_plus_1) 
           ./ (Delta_t_minus .+ Delta_t_plus))
end

#########################################################
# DIAGNOSTICS INVOLVING DERIVATIVES OF PRIMITIVE FIELDS #
#########################################################

function pointwise_∇b(i, j, k, b, xC, yC, zC)
   #=
   Provided buoyancy as a 3D array, compute Cartesian buoyancy-gradient
    components via centered differencing at the specified [i, j, k] coordinate
    triple.
   Useful for constructing KernelFunctionOperations or computing gradient from
    buoyancy data saved to NetCDF file.
   =#

   ∂x_b = centered_difference(xC[(i - 1):(i + 1)], b[(i - 1):(i + 1), j, k])
   ∂y_b = centered_difference(yC[(j - 1):(j + 1)], b[i, (j - 1):(j + 1), k])
   ∂z_b = centered_difference(zC[(k - 1):(k + 1)], b[i, j, (k - 1):(k + 1)])
   
   return ∂x_b, ∂y_b, ∂z_b
end

function pointwise_ω(i, j, k, ux, uy, uz, xC, yC, zC)
   #=
   Provided Cartesian velocity components as 3D arrays, compute Cartesian 
    vorticity components via centered differencing at the specified [i, j, k] 
    coordinate triple.
   Useful for constructing KernelFunctionOperations or computing vorticity from
    velocity data saved to NetCDF file.
   =#
   
   ∂y_ux = centered_difference(yC[(j - 1):(j + 1)], uy[i, (j - 1):(j + 1), k])
   ∂z_ux = centered_difference(zC[(k - 1):(k + 1)], uz[i, j, (k - 1):(k + 1)])
   ∂x_uy = centered_difference(xC[(i - 1):(i + 1)], uy[(i - 1):(i + 1), j, k])
   ∂z_uy = centered_difference(zC[(k - 1):(k + 1)], uy[i, j, (k - 1):(k + 1)])
   ∂x_uz = centered_difference(xC[(i - 1):(i + 1)], uz[(i - 1):(i + 1), j, k])
   ∂y_uz = centered_difference(yC[(j - 1):(j + 1)], uz[i, (j - 1):(j + 1), k]) 
   
   ωx = ∂y_uz - ∂z_uy
   ωy = ∂z_ux - ∂x_uz
   ωz = ∂x_uy - ∂y_ux
   
   return ωx, ωy, ωz
end

function CenterFields_ω(grid, ux_Field, uy_Field, uz_Field)
   #=
   Provided Cartesian velocity-component Fields (XFaceField, etc.) on 'grid', 
    compute and return Cartesian vorticity components as CenterFields on
    'grid'.
   Useful for computing vorticity directly on Oceananigans grid, including at
    simulation runtime.
   =#

   ωx = CenterField(grid)
   ωy = CenterField(grid)
   ωz = CenterField(grid)

   set!(ωx, ∂y(uz_Field) - ∂z(uy_Field))
   set!(ωy, ∂z(ux_Field) - ∂x(uz_Field))
   set!(ωz, ∂x(uy_Field) - ∂y(ux_Field))
   
   return ωx, ωy, ωz
end

function CenterField_ωr(grid, uφ_Field, uz_Field, rCenterField, φCenterField)
   #=
   Provided Cartesian velocity-component Fields (XFaceField, etc.) on 'grid', 
    compute and return radial vorticity component as a CenterField on 'grid'.
   Useful for computing vorticity directly on Oceananigans grid, including at
    simulation runtime.
   =#
   
   @inline S(i, j, k, g) = @inbounds sin(φCenterField)[i, j, k]
   @inline C(i, j, k, g) = @inbounds cos(φCenterField)[i, j, k]

   @inline ωr_ccc(i, j, k, g) = @inbounds ((C(i, j, k, g) 
                                            * ∂y(uz_Field)[i, j, k] 
                                            - S(i, j, k, g) 
                                            * ∂x(uz_Field)[i, j, k]
                                           ) - (∂z(uφ_Field)[i, j, k] 
                                                / rCenterField[i, j, k]
                                               )
                                          )
   
   ωr_op = KernelFunctionOperation{Center, Center, Center}(ωr_ccc, grid)
   
   @compute ωr = Field(ωr_op)
   
   return ωr
end

function ωz(ux, uy, x, y; x_idx = nothing, y_idx = nothing, z_idx = nothing)
   #=
   Compute only the vertical component of vorticity, via centered differencing, 
    on a 2D slice.
   =#

   if !isnothing(x_idx)
      ΔuyΔx = centered_difference(x[(x_idx - 1):(x_idx + 1)], 
                                  uy[(x_idx - 1):(x_idx + 1), :, :])
      ΔuxΔy = centered_difference(y, ux[x_idx, :, :])
   elseif !isnothing(y_idx)
      ΔuyΔx = centered_difference(x, uy[:, y_idx, :])
      ΔuxΔy = centered_difference(y[(y_idx - 1):(y_idx + 1)], 
                                  ux[:, (y_idx -1):(y_idx + 1), :])
   elseif !isnothing(z_idx)
      ΔuyΔx = centered_difference(x, uy[:, :, z_idx])
      ΔuxΔy = centered_difference(y, ux[:, :, z_idx])
   end
   
   ωz = ΔuyΔx .- ΔuxΔy
end

function pointwise_q_Ertel(i, j, k, b, ux, uy, uz, f, xC, yC, zC)

   ω  = pointwise_ω(i, j, k, ux, uy, uz, xC, yC, zC)
   ∇b = pointwise_∇b(i, j, k, b, xC, yC, zC)
   
   ωx, ωy, ωz       = ω[1][1], ω[2][1], ω[3][1]
   ∂x_b, ∂y_b, ∂z_b = ∇b[1][1], ∇b[2][1], ∇b[3][1]

   q = (ωx * ∂x_b) + (ωy * ∂y_b) + ((f + ωz) * ∂z_b)
end

function CenterFields_q_Ertel(b_Field, ux_Field, uy_Field, uz_Field, grid, f; 
                              returnAsArray = false)
   
   ∂x_b_Field, ∂y_b_Field, ∂z_b_Field = ∂x(b_Field), ∂y(b_Field), ∂z(b_Field)
   ωx_Field, ωy_Field, ωz_Field       = CenterFields_ω(grid, ux_Field, uy_Field,
                                                       uz_Field)
   
   q = (ωx_Field * ∂x_b_Field) + (ωy_Field * ∂y_b_Field) 
        + ((f + ωz_Field) * ∂z_b_Field)

   if returnAsArray
      return adapt(Array, q)
   elseif !returnAsArray
      return q
   end
end

function normalStrainInHorizontalFlow(ux, uy, Δx, Δy; 
                                      x_idx = nothing, 
                                      y_idx = nothing, 
                                      z_idx = nothing)
   #=
   On a 2D slice, compute normal strain in horizontal (i.e., to xy-plane) 
    projection of flow.
   =#

   if !isnothing(x_idx)
      ΔuxΔx = centered_difference(x[(x_idx - 1):(x_idx + 1)],
                                  ux[(x_idx - 1):(x_idx + 1), :, :])
      ΔuyΔy = centered_difference(y, uy[x_idx, :, :])
   elseif !isnothing(y_idx)
      ΔuxΔx = centered_difference(x, ux[:, y_idx, :])
      ΔuyΔy = centered_difference(y[(y_idx - 1):(y_idx + 1)], 
                                  uy[:, (y_idx -1):(y_idx + 1), :])
   elseif !isnothing(z_idx)
      ΔuxΔx = centered_difference(x, ux[:, :, z_idx])
      ΔuyΔy = centered_difference(y, uy[:, :, z_idx])
   end
   
   Sn = ΔuxΔx .- ΔuyΔy
end

function shearStrainInHorizontalFlow(ux, uy, Δx, Δy; 
                                     x_idx = nothing, 
                                     y_idx = nothing,
                                     z_idx = nothing)
   #=
   On a 2D slice, compute shear strain in horizontal (i.e., to xy-plane)
    projection of flow.
   =#
                                     
   if !isnothing(x_idx)
      ΔuyΔx = centered_difference(x[(x_idx - 1):(x_idx + 1)],
                                  uy[(x_idx - 1):(x_idx + 1), :, :])
      centered_difference(y, ux[x_idx, :, :])
   elseif !isnothing(y_idx)
      ΔuyΔx = centered_difference(x, uy[:, y_idx, :])
      ΔuxΔy = centered_difference(y[(y_idx - 1):(y_idx + 1)], 
                                  ux[:, (y_idx -1):(y_idx + 1), :])
   elseif !isnothing(z_idx)
      ΔuyΔx = centered_difference(x, uy[:, :, z_idx])
      ΔuxΔy = centered_difference(y, ux[:, :, z_idx])
   end
   
   Ss = ΔuyΔx .+ ΔuxΔy
end

function OkuboWeiss(u, v, Δx, Δy; 
                    x_idx = nothing, y_idx = nothing, z_idx = nothing)
   #=
   Compute W pointwise on a 2D slice.
   =#

   Sn = normalStrainInHorizontalFlow(u, v, Δx, Δy; 
                                     x_idx = x_idx, y_idx = y_idx, z_idx = z_idx)
   Ss = shearStrainInHorizontalFlow(u, v, Δx, Δy; 
                                    x_idx = x_idx, y_idx = y_idx, z_idx = z_idx)
   ζ  = ωz(u, v, Δx, Δy; x_idx = x_idx, y_idx = y_idx, z_idx = z_idx)

   return @. Sn^2 + Ss^2 - ζ^2
end

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
   Sort data chronologically (essential if simulation was ever picked up from a
    checkpoint).
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

function open_dataset(outfilename; Hx = 3, Hy = 3, Hz = 3)
   #=
   Open a generic NetCDF dataset output by Oceananigans simulation.
   Return the dataset itself as well as gridpoints and timeseries info.   
   =#

   ds = NCDataset(outfilename)

   x, y, z = nothing, nothing, nothing #Defaults in case any dimension is Flat
   
   #Load coords of non-Flat dimensions; convert them to km for readability

   if length(ds[:x_caa][:]) > 1
      x = ds[:x_caa][:] ./ 1000 #ds[:x_caa][(Hx + 1):(length(ds[:x_caa][:]) - Hx)] ./ 1000
   end

   if length(ds[:y_aca][:]) > 1
      y = ds[:y_aca][:] ./ 1000 #ds[:y_aca][(Hy + 1):(length(ds[:y_aca][:]) - Hy)] ./ 1000
   end
   
   if length(ds[:z_aac][:]) > 1
      zC = ds[:z_aac][:] ./ 1000 # ds[:z_aac][(Hz + 1):(length(ds[:z_aac][:]) - Hz)] ./ 1000
      zF = ds[:z_aaf][:] ./ 1000 # ds[:z_aaf][(Hz + 1):(length(ds[:z_aaf][:]) - Hz)] ./ 1000
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

###################################

function ζa_b(U, f, σr, σz, x, y, z)
   r2_arr = @. x^2 + y^2
   z2_arr = transpose(z.^2 .* ones(Float64, (1, length(r2_arr))))
   ζa_b   = @. (f + (2 * U / σr) * (r2_arr / (σr^2) - 1)
	              * exp(1 - r2_arr / (σr^2) - z2_arr / (σz^2)))
   return ζa_b
end

function ζa(f, u, v, w, Δx, Δy, Δz)
   ωx, ωy, ωz = ω(u, v, w, Δx, Δy, Δz)
   ζa         = f + ωz
end

function field_norm(ψ, n; ψ_bkgd = 0)
   #=
   Compute L2-norm of a perturbation field.
   =#

   ψ_n          = ψ[:, :, :, n]
   ψ_perturb_n  = ψ_n .- ψ_bkgd
   perturb_norm = norm(ψ_perturb_n)
end

function ∂r_q(q, x, y, zk, all_z)
   
   k = findfirst(==(zk), all_z)
   
   ∂x_qk = 0 .* q[2:end-1, 2:end-1, k]
   ∂y_qk = 0 .* q[2:end-1, 2:end-1, k]
   
   for j in 1:(length(y)):1
      ∂x_qk[:, j] += centered_difference(x, q[:, j, k])
   end
   
   for i in 1:(length(x)):1
      ∂y_qk[i, :] += centered_difference(y, q[i, :, k])
   end

   x = x[2:end-1]
   y = y[2:end-1]
   
   r     = sqrt.(x.^2 + y.^2) 
   ∂r_qk = @. (x * ∂x_qk + y * ∂y_qk) / r
end

function empirical_growth_rate(t, perturb_norm; differencing = "forward")
   #=
   Diagnose growth-rate timeseries of a field with L2-norm perturb_norm.
   =#
   
   if differencing == "forward"
      growth_rate = order1_forward_difference(t, log.(perturb_norm))
      new_t       = t[1:end-1]
   elseif differencing == "centered"
      growth_rate = centered_difference(t, log.(perturb_norm))
      new_t       = t[2:end-1]
   end
   return new_t, growth_rate
end

#=
function buoyancy_L(x,y,z,N²,f,Umax,D,Lⱼ,z0,y0)
    Bba = @. 2*f*Umax*Lⱼ/D^2
    Bbb = @. (tanh((y-y0)/Lⱼ)+1)
    Bbc = @. (z-z0) * exp(-(z-z0)^2/D^2)
    b_L = N².*z.+ Bba.* Bbc .* transpose(Bbb)   #Left Boundary = N²*z
return b_L
end

function buoyancy_R(x,y,z,N²,f,Umax,D,Lⱼ,z0,y0)
    Bba = @. 2*f*Umax*Lⱼ/D^2
    Bbb = @. (tanh((y-y0)/Lⱼ) - 1)
    Bbc = @. (z-z0) * exp(-(z-z0)^2/D^2)
    b_R = N².*z .+ Bba.* Bbc .* transpose(Bbb)   #Right Boundary = N²*z
return b_R
end

function buoyancy_C(x,y,z,N²,f,Umax,D,Lⱼ,z0,y0)
    Bba = @. 2*f*Umax*Lⱼ/D^2
    Bbb = @. (tanh((y-y0)/Lⱼ))
    Bbc = @. (z-z0) * exp(-(z-z0)^2/D^2)
    b_C = N².*z.+ Bba.* Bbc .* transpose(Bbb)   #Center = N²*z
return b_C
end

function vel_B(x,y,z,N²,f,Umax,D,Lⱼ,z0,y0)
    Uba = @. Umax/cosh((y-y0)/Lⱼ)^2
    Ubb = @. exp(-(z-z0)^2/D^2)
    Ub= Ubb*transpose(Uba)   #background velocity field
return Ub
end

function density_from_buoyancy(b, ρ₀)
    # Compute the density (ρ) from the buoyancy (b) and reference density (ρ₀)
    ρ = @. ρ₀ * (1 - b / 9.81)
    return ρ
end

function ωt(x,y,z,Umax,D,Lⱼ,z0,y0)
    ωba = @. sech((y-y0)/Lⱼ) * tanh((y-y0)/Lⱼ)
    ωb = @. 2 * Umax/Lⱼ * exp(-(z-z0)^2/D^2)
    ωt = ωb .* transpose(ωba)   #background vorticity field
end

function BestFit(degree,interval, abscissa,ordenate)
    linear_fit_polynomial = fit(abscissa[interval], log.(ordenate[interval]), degree, var = :abscissa)
    constant, slope = linear_fit_polynomial[0], linear_fit_polynomial[1]
    best_fit = @. exp(constant + slope * abscissa)
    @sprintf "The growth rate is approximately %5.1e" slope
    return best_fit, constant, slope
end
=#