using Printf

module ComputeSecondaries
   export ω, ω_2D, ζa_b, ζa, ∇b, ∇b_2D, ertelQ, ertelQ_2D, ∂r_ertelQ
end

function ω(u, v, w, i, j, k, Δx, Δy, Δz)
   ωx = @. ((w[i, j:j+1, k] - w[i, j-1:j, k]) / Δy 
            - (v[i, j, k:k+1] - v[i, j, k-1:k]) / Δz)
   ωy = @. ((u[i, j, k:k+1]- u[i, j, k-1:k]) / Δz
	    - (w[i:i+1, j, k] - w[i-1:i, j, k]) / Δx)
   ωz = @. ((v[i:i+1, j, k] - v[i-1:i, j, k]) / Δx
	    - (u[i, j:j+1, k] - u[i, j-1:j, k]) / Δy)
   return ωx, ωy, ωz
end

function ω_2D(u, v, w, Δx, Δy, Δz; 
	      x_idx = nothing, y_idx = nothing, z_idx = nothing)
   if !isnothing(x_idx)
      ω_2D = @. ((w[x_idx,2:end,2:end] - w[x_idx,1:end-1,2:end]) / Δy 
		- (v[x_idx,2:end,2:end] - v[x_idx,2:end,1:end-1]) / Δz)
   elseif !isnothing(y_idx)
      ω_2D = @. ((u[2:end,y_idx,2:end] - u[2:end,y_idx,1:end-1]) / Δz 
		- (w[2:end,y_idx,2:end] - w[1:end-1,y_idx,2:end]) / Δx)
   elseif !isnothing(z_idx)
      ω_2D = @. ((v[2:end,2:end,z_idx] - v[1:end-1,2:end,z_idx]) / Δx 
		- (u[2:end,2:end,z_idx] - u[2:end,1:end-1,z_idx]) / Δy)
   end
   return ω_2D
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
   ∂x_b = @. (b[i:i+1, j, k] - b[i-1:i, j, k]) / (2*Δx) #(b[2:end,2:end,2:end] - b[1:end-1,2:end,2:end]) / Δx
   ∂y_b = @. (b[i, j:j+1, k] - b[i, j-1:j, k]) / (2*Δy) # (b[2:end,2:end,2:end] - b[2:end,1:end-1,2:end]) / Δy
   ∂z_b = @. (b[i, j, k:k+1] - b[i, j, k-1:k]) / (2*Δz) #(b[2:end,2:end,2:end] - b[2:end,2:end,1:end-1]) / Δz
   return ∂x_b, ∂y_b, ∂z_b
end

function ∇b_2D(b, Δx, Δy, Δz; 
	       x_idx = nothing, y_idx = nothing, z_idx = nothing)
   if !isnothing(x_idx)
      ∇b_2D = @. ((b[x_idx,2:end,2:end] - b[x_idx,1:end-1,2:end]) / Δy 
		 + (b[x_idx,2:end,2:end] - b[x_idx,2:end,1:end-1]) / Δz)
   elseif !isnothing(y_idx)
      ∇b_2D = @. ((b[2:end,y_idx,2:end] - b[1:end-1,y_idx,2:end]) / Δx 
	         + (b[2:end,y_idx,2:end] - b[2:end,y_idx,1:end-1]) / Δz)
   elseif !isnothing(z_idx)
      ∇b_2D = @. ((b[2:end,2:end,z_idx] - b[1:end-1,2:end,z_idx]) / Δx 
		 + (b[2:end,2:end,z_idx] - b[2:end,1:end-1,z_idx]) / Δy)
   end
   return ∇b_2D
end

function ertelQ(u, v, w, b, f, x_idx, y_idx, z_idx, Δx, Δy, Δz)
   ωx, ωy, ωz = ω(u, v, w, x_idx, y_idx, z_idx, Δx, Δy, Δz)
   ∂x_b, ∂y_b, ∂z_b = ∇b(b, x_idx, y_idx, z_idx, Δx, Δy, Δz)
   Q = @. (ωx * ∂x_b) + (ωy * ∂y_b) + ((f + ωz) * ∂z_b)
end

function ertelQ_2D(u, v, w, b, f, Δx, Δy, Δz; 
		   x_idx = nothing, y_idx = nothing, z_idx = nothing)
   if !isnothing(x_idx)
      ωx, ωy, ωz = ω(u[x_idx-1:x_idx+1,:,:], v[x_idx-1:x_idx+1,:,:],
		     w[x_idx-1:x_idx+1,:,:], Δx, Δy, Δz)
      ∂x_b, ∂y_b, ∂z_b = ∇b(b[x_idx-1:x_idx+1,:,:], Δx, Δy, Δz)
      Q_2D = @. (ωy * ∂y_b) + ((f + ωz) * ∂z_b)
      avg_Q_2D = @. (Q_2D[1,:,:] + Q_2D[2,:,:]) / 2
   elseif !isnothing(y_idx)
      ωx, ωy, ωz = ω(u[:,y_idx-1:y_idx+1,:], v[:,y_idx-1:y_idx+1,:],
                     w[:,y_idx-1:y_idx+1,:], Δx, Δy, Δz)
      ∂x_b, ∂y_b, ∂z_b = ∇b(b[:,y_idx-1:y_idx+1,:], Δx, Δy, Δz)
      Q_2D = @. (ωx * ∂x_b) + ((f + ωz) * ∂z_b)
      avg_Q_2D = @. (Q_2D[:,1,:] + Q_2D[:,2,:]) / 2
   elseif !isnothing(z_idx)
      ωx, ωy, ωz = ω(u[:,:,z_idx-1:z_idx+1], v[:,:,z_idx-1:z_idx+1], 
		     w[:,:,z_idx-1:z_idx+1], Δx, Δy, Δz)
      ∂x_b, ∂y_b, ∂z_b = ∇b(b[:,:,z_idx-1:z_idx+1], Δx, Δy, Δz)
      Q_2D = @. (ωx * ∂x_b) + (ωy * ∂y_b)
      avg_Q_2D = @. (Q_2D[:,:,1] + Q_2D[:,:,2]) / 2
   end
   return avg_Q_2D
end

function ∂r_ertelQ(Q, Δx, Δy, xC, yC)
   ∂x_Q = @. (Q[2:end,2:end,:] - Q[1:end-1,2:end,:]) / Δx
   ∂y_Q = @. (Q[2:end,2:end,:] - Q[2:end,1:end-1,:]) / Δy
   r    = @. (xC^2 + yC^2)^0.5
   ∂r_Q = @. (xC * ∂x_Q + yC * ∂y_Q) / r
   return ∂r_Q
end

#=function u_background(x,y,z,Umax,D,Lⱼ,z0,y0)
    Uba = @. Umax/cosh((y-y0)/Lⱼ)^2
    Ubb = @. exp(-(z-z0)^2/D^2)
    Ub= Ubb*transpose(Uba) 
end

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
# function ω_b(x,y,z,Umax,D,Lⱼ,z0)
#     ωa = @. 2 / Lⱼ * Umax
#     ωb = @. exp(-(z-z0)^2 / D^2)
#     ωc = transpose(tanh.(y./Lⱼ) ./(cosh.(y./Lⱼ).^2))
#     ω_background =  ωa .* ωb .* ωc   #background vorticity field
# end

function BestFit(degree,interval, abscissa,ordenate)
    linear_fit_polynomial = fit(abscissa[interval], log.(ordenate[interval]), degree, var = :abscissa)
    constant, slope = linear_fit_polynomial[0], linear_fit_polynomial[1]
    best_fit = @. exp(constant + slope * abscissa)
    @sprintf "The growth rate is approximately %5.1e" slope
    return best_fit, constant, slope
end
=#
