using LinearAlgebra, Printf
using Oceananigans.Fields

####################

module ComputeSecondaries
   export ω, ωz, ζa_b, ζa, ∇b, q, ∂r_q, field_norm, φ
end

####################

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
	        + v[x_idx, 2:end-1, :] - v[x_idx-1, 2:end-1, :]) / (2*Δx)
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
   ζa_b = @. (f + (2*U/σr) * (r2_arr/(σr^2) - 1)
	      * exp(1 - r2_arr/(σr^2) - z2_arr/(σz^2)))
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
   ∂r_q = @. (x*∂x_q + y*∂y_q) / r
   
   return (∂r_q[1] + ∂r_q[2]) / 2
end

function field_norm(φ, n)

   φ_initial, φ_n = φ[:, :, :, 1], φ[:, :, :, n]
   φ_perturb_n    = φ_n .- φ_initial

   perturb_norm  = norm(φ_perturb_n)
   relative_norm = perturb_norm / norm(φ_initial)
   
   return perturb_norm, relative_norm
end

function φ(i, j, k, grid)
   xc_i = xnodes(grid, Center())[i]
   yc_j = ynodes(grid, Center())[j]
   
   return atan(yc_j, xc_i)
end

#####################

using Adapt, NCDatasets

####################

module VisFunctions
   export open_dataset, get_background_fields, get_range_lims, 
   get_2D_spatial_axis_kwargs
end

####################

function open_dataset(datetime)

   outfilepath = joinpath("./Output", "output_$(datetime).nc")

   ds = NCDataset(outfilepath, "r")
   x  = ds["xC"][:] ./ 1000 #Convert to km for readability
   y  = ds["yC"][:] ./ 1000 #Convert to km for readability
   z  = ds["zC"][:]
   t  = ds["time"][:]
   Nt = length(t)

   return ds, x, y, z, t, Nt
end

function get_background_fields(ds)
   bb = @views adapt(Array, ds["b"])[:, :, :, 1]
   ub = @views adapt(Array, ds["u"])[:, :, :, 1]
   vb = @views adapt(Array, ds["v"])[:, :, :, 1]
   wb = @views adapt(Array, ds["w"])[:, :, 1:end-1, 1]
   return bb, ub, vb, wb
end

function get_range_lims(final_field; prescribed_max = 0)
   field_max  = max(maximum(abs.(final_field)), prescribed_max)
   field_lims = [-field_max, field_max]
end

function get_2D_spatial_axis_kwargs(x, y, z;
                                    x_idx = nothing, y_idx = nothing,
                                    z_idx = nothing)
   if !isnothing(x_idx)
      nearest     = round(Int, x[x_idx])
      axis_kwargs = (xlabel = "y [km]", ylabel = "z [m]")
   elseif !isnothing(y_idx)
      nearest     = round(Int, y[y_idx])
      axis_kwargs = (xlabel = "x [km]", ylabel = "z [m]")
   elseif !isnothing(z_idx)
      nearest     = round(Int, z[z_idx])
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
