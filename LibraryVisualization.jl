using Adapt, CUDA, Glob, LinearAlgebra, NCDatasets
using Oceananigans.AbstractOperations, Oceananigans.Fields

####################

function compute_polar_coords(grid)

   function r_coord(i, j, k, grid)
      xCi = @views adapt(CuArray, xnodes(grid, Center()))[i]
      yCj = @views adapt(CuArray, ynodes(grid, Center()))[j]
      rij = sqrt(xCi^2 + yCj^2)
   end

   function φ_coord(i, j, k, grid)
      xCi = @views adapt(CuArray, xnodes(grid, Center()))[i]
      yCj = @views adapt(CuArray, ynodes(grid, Center()))[j]
      φij = atan(yCj, xCi)
   end

   r_KernOp = KernelFunctionOperation{Center, Center, Center}(r_coord, grid)
   r        = Field(r_KernOp)
   φ_KernOp = KernelFunctionOperation{Center, Center, Center}(φ_coord, grid)
   φ        = Field(φ_KernOp)

   compute!(r)
   compute!(φ)
   return(r, φ)
end

function xy_vector_to_rφ(vx, vy, grid, useGPU)

   r, φ = compute_polar_coords(grid)

   if useGPU
   
      function interpolate_vx(i, j, k, grid)
         vxCi = @views interpolate((adapt(CuArray, xnodes(grid, Center()))[i],
                                    adapt(CuArray, ynodes(grid, Center()))[j],
                                    adapt(CuArray, znodes(grid, Center()))[k]),
                                   vx, (Face(), Center(), Center()), grid)
      end

      function interpolate_vy(i, j, k, grid)
         vyCj = @views interpolate((adapt(CuArray, xnodes(grid, Center()))[i],
                                    adapt(CuArray, ynodes(grid, Center()))[j],
                                    adapt(CuArray, znodes(grid, Center()))[k]),
                                   vy, (Center(), Face(), Center()), grid)
      end
   
      vxC_KernOp = KernelFunctionOperation{Center, Center, Center}(
                                        interpolate_vx, grid)
      vxC        = Field(vxC_KernOp)
      vyC_KernOp = KernelFunctionOperation{Center, Center, Center}(
                                        interpolate_vy, grid)
      vyC        = Field(vyC_KernOp)
      
      compute!(vxC)
      compute!(vyC)

      vr_BinaryOp = (vxC * cos(φ)) + (vyC * sin(φ))
      vr          = Field(vr_BinaryOp)
      vφ_BinaryOp = (vyC * cos(φ)) - (vxC * sin(φ))
      vφ          = Field(vφ_BinaryOp)
      
   elseif !useGPU
   
      vr = @at (Center, Center, Center) (vx * cos(φ)) + (vy * sin(φ))
      vφ = @at (Center, Center, Center) (vy * cos(φ)) - (vx * sin(φ))
   end

   compute!(vr)
   compute!(vφ)
   return vr, vφ
end

function ω(u, v, w, i, j, k, Δx, Δy, Δz)
   
   ωx = @. ((w[i, j:j+1, k] - w[i, j-1:j, k]) / Δy 
	          - (v[i, j, k:k+1] - v[i, j, k-1:k]) / Δz)
   
   ωy = @. ((u[i, j, k:k+1] - u[i, j, k-1:k]) / Δz
	          - (w[i:i+1, j, k] - w[i-1:i, j, k]) / Δx)
   
   ωz = @. ((v[i:i+1, j, k] - v[i-1:i, j, k]) / Δx
	          - (u[i, j:j+1, k] - u[i, j-1:j, k]) / Δy)
   
   return (ωx[1] + ωx[2]) / 2, (ωy[1] + ωy[2]) / 2, (ωz[1] + ωz[2]) / 2
end

function ωz(u, v, Δx, Δy; 
            x_idx = nothing, y_idx = nothing, z_idx = nothing)

   if !isnothing(x_idx)
      ωz = @. (((v[x_idx+1, 2:end-1, :] - v[x_idx, 2:end-1, :]) 
                 + v[x_idx, 2:end-1, :] - v[x_idx-1, 2:end-1, :]) / (2 * Δx)
                - (u[x_idx, 2:end, :] - u[x_idx, 1:end-1, :]) / Δy)
   
   elseif !isnothing(y_idx)
      ωz = @. ((v[2:end, y_idx, :] - v[1:end-1, y_idx, :]) / Δx
	             - ((u[2:end-1, y_idx+1, :] - u[2:end-1, y_idx, :])
		               + u[2:end-1, y_idx, :] - u[2:end-1, y_idx-1, :]) / (2*Δy))

   elseif !isnothing(z_idx)
      ωz = @. ((v[2:end, 2:end-1, z_idx] - v[1:end-1, 2:end-1, z_idx]) / Δx
               - (u[2:end-1, 2:end, z_idx] - u[2:end-1, 1:end-1, z_idx]) / Δy)
   end
   return ωz
end

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

function ∇b(b, i, j, k, Δx, Δy, Δz)
   ∂x_b = @. (b[i:i+1, j, k] - b[i-1:i, j, k]) / Δx
   ∂y_b = @. (b[i, j:j+1, k] - b[i, j-1:j, k]) / Δy
   ∂z_b = @. (b[i, j, k:k+1] - b[i, j, k-1:k]) / Δz
   
   return ((∂x_b[1] + ∂x_b[2]) / 2, 
           (∂y_b[1] + ∂y_b[2]) / 2, 
           (∂z_b[1] + ∂z_b[2]) / 2)
end

function q(u, v, w, b, f, x_idx, y_idx, z_idx, Δx, Δy, Δz)
   ωx, ωy, ωz       = ω(u, v, w, x_idx, y_idx, z_idx, Δx, Δy, Δz)
   ∂x_b, ∂y_b, ∂z_b = ∇b(b, x_idx, y_idx, z_idx, Δx, Δy, Δz)
   q                = (ωx * ∂x_b) + (ωy * ∂y_b) + ((f + ωz) * ∂z_b)
end

function ∂r_q(q, x, y, i, j, k, Δx, Δy)
   
   ∂x_q = @. (q[i:i+1, j, k] - q[i-1:i, j, k]) / Δx
   ∂y_q = @. (q[i, j:j+1, k] - q[i, j-1:j, k]) / Δy
   
   r    = sqrt(x^2 + y^2) 
   ∂r_q = @. (x * ∂x_q + y * ∂y_q) / r
   
   return (∂r_q[1] + ∂r_q[2]) / 2
end

function field_norm(ψ, n; ψ_bkgd = 0)
   ψ_n          = ψ[:, :, :, n]
   ψ_perturb_n  = ψ_n .- ψ_bkgd
   perturb_norm = norm(ψ_perturb_n)
end

function pad_filenames(datetime; prefix = "output")

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

function load_times_in_days(dataset)

   t  = dataset[:time][:] ./ 86400
   Nt = length(t)
   
   return t, Nt
end

function open_dataset(outfilename; Hx = 3, Hy = 3, Hz = 3)

   ds = NCDataset(outfilename)

   x, y, z = nothing, nothing, nothing #Defaults in case any dimension is Flat
   
   #Load coords of non-Flat dimensions; convert them to km for readability

   if length(ds[:x_caa][:]) > 1
      x = ds[:x_caa][Hx+1:length(ds[:x_caa][:])-Hx] ./ 1000
   end

   if length(ds[:y_aca][:]) > 1
      y = ds[:y_aca][Hy+1:length(ds[:y_aca][:])-Hy] ./ 1000
   end
   
   if length(ds[:z_aac][:]) > 1
      zC = ds[:z_aac][Hz+1:length(ds[:z_aac][:])-Hz] ./ 1000
      zF = ds[:z_aaf][Hz+1:length(ds[:z_aaf][:])-Hz] ./ 1000
   end

   t, Nt = load_times_in_days(ds)

   return ds, x, y, zC, zF, t, Nt
end

function open_energetics_dataset(energeticsfilename)
   
   energetics_ds = NCDataset(energeticsfilename)
   t, Nt         = load_times_in_days(energetics_ds)

   return energetics_ds, t, Nt
end

function open_scalars_dataset(scalarfilename)

   scalars_ds = NCDataset(scalarfilename)
   t, Nt      = load_times_in_days(scalars_ds)

   return scalars_ds, t, Nt
end

function order1_forward_difference(t, u)
   return @. (u[2:end] - u[1:end-1]) / (t[2:end] - t[1:end-1])
end

function get_range_lims(final_field; max_fraction = 1, prescribed_max = 1e-16)
   field_max  = max(maximum(abs.(final_field)), prescribed_max)
   field_lims = [-(max_fraction * field_max), (max_fraction * field_max)]
end

function get_2D_spatial_axis_idcs(const_dim;
                                  Hx = 3, Hy = 3, Hz = 3,
		                              x_idx = nothing, y_idx = nothing, z_idx = nothing,
				                          xC = nothing, yC = nothing, zC = nothing,
                      				    zF = nothing)

   if const_dim == "x"

      if isnothing(x_idx) #Grid is 2D with only y and z axes
         yCzC_idcs = (1, Hy+1:length(yC)+Hy, Hz+1:length(zC)+Hz)
         yCzF_idcs = (1, Hy+1:length(yC)+Hy, Hz+1:length(zF)+Hz) 
      else #Grid is 3D
         yCzC_idcs = (x_idx, Hy+1:length(yC)+Hy, Hz+1:length(zC)+Hz)
         yCzF_idcs = (x_idx, Hy+1:length(yC)+Hy, Hz+1:length(zF)+Hz)
      end

      return yCzC_idcs, yCzF_idcs

   elseif const_dim == "y"

      if isnothing(y_idx) #Grid is 2D with only x and z axes
         xCzC_idcs = (Hx+1:length(xC)+Hx, 1, Hz+1:length(zC)+Hz)
         xCzF_idcs = (Hx+1:length(xC)+Hx, 1, Hz+1:length(zF)+Hz)
      else #Grid is 3D
         xCzC_idcs = (Hx+1:length(xC)+Hx, y_idx, Hz+1:length(zC)+Hz)
         xCzF_idcs = (Hx+1:length(xC)+Hx, y_idx, Hz+1:length(zF)+Hz)
      end

      return xCzC_idcs, xCzF_idcs

   elseif const_dim == "z"

      if isnothing(z_idx) #Grid is 2D with only x and y axes
         xCyC_idcs = (Hx+1:length(xC)+Hx, Hy+1:length(yC)+Hy, 1)
      else #Grid is 3D
         xCyC_idcs = (Hx+1:length(xC)+Hx, Hy+1:length(yC)+Hy, z_idx)
      end

      return xCyC_idcs, xCyC_idcs
   end
end

function get_2D_spatial_axis_kwargs(x, y, z, const_dim;
                                    x_idx = nothing,
				                            y_idx = nothing,
                                    z_idx = nothing)

   nearest = 0

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
   
   return nearest, axis_kwargs
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