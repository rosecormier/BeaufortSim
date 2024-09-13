using Printf

module ComputeSecondaries
   export ω, ω_2D, ∇b, ∇b_2D, ertelQ, ertelQ_2D, ∂r_ertelQ
end

function ω(u, v, w, Δx, Δy, Δz)
   ωx = @. ((w[2:end,2:end,2:end] - w[2:end,1:end-1,2:end]) / Δy 
           - (v[2:end,2:end,2:end] - v[2:end,2:end,1:end-1]) / Δz)
   ωy = @. ((u[2:end,2:end,2:end] - u[2:end,2:end,1:end-1]) / Δz 
	   - (w[2:end,2:end,2:end] - w[1:end-1,2:end,2:end]) / Δx)
   ωz = @. ((v[2:end,2:end,2:end] - v[1:end-1,2:end,2:end]) / Δx 
	   - (u[2:end,2:end,2:end] - u[2:end,1:end-1,2:end]) / Δy)
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

function ∇b(b, Δx, Δy, Δz)
   ∂x_b = (b[2:end,2:end,2:end] - b[1:end-1,2:end,2:end]) / Δx
   ∂y_b = (b[2:end,2:end,2:end] - b[2:end,1:end-1,2:end]) / Δy
   ∂z_b = (b[2:end,2:end,2:end] - b[2:end,2:end,1:end-1]) / Δz
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

function ertelQ(u, v, w, b, f, Δx, Δy, Δz)
   ωx, ωy, ωz = ω(u, v, w, Δx, Δy, Δz)
   ∂x_b, ∂y_b, ∂z_b = ∇b(b, Δx, Δy, Δz)
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
#=
function randomSine(x,y,z, Nx,Ny,Nz)
#Random coefficients (gaussian-distributed) for phase (ϕ) and amplitude (A)

amplitude = FileIO.load("cosineVariables.jld2","amplitude")
phase = FileIO.load("cosineVariables.jld2","phase")

ϕx = phase[:,1]
Ax = amplitude[:,1]
ϕy = phase[:,2]
Ay = amplitude[:,2]
ϕz = phase[:,3]
Az = amplitude[:,3]

rs = 0
    for k = 1:400
        rs = rs *(Ax[k].*sin.(2*π .* x .* (k) / Nx .+ ϕx[k]) .* Ay[k].*sin.(2*π.*y .* (k) / Ny .+ ϕy[k]) .* 
        Az[k].*sin.(2*π.*z .* (k) / Nz .+ ϕz[k]))
    end
    return rs
end

function randomSineGPU(x,y,z, Nx,Ny,Nz)
#Random coefficients (gaussian-distributed) for phase (ϕ) and amplitude (A)

amplitude = FileIO.load("cosineVariables.jld2","amplitude")
phase = FileIO.load("cosineVariables.jld2","phase")

ϕx = phase[:,1]
Ax = amplitude[:,1]
ϕy = phase[:,2]
Ay = amplitude[:,2]
ϕz = phase[:,3]
Az = amplitude[:,3]

rs = 0
Threads.@threads for k = 1:400
        rs = rs + (Ax[k].*sin.(2*π .* x .* (k) / Nx .+ ϕx[k]) .* Ay[k].*sin.(2*π.*y .* (k) / Ny .+ ϕy[k]) .*
        Az[k].*sin.(2*π.*z .* (k) / Nz .+ ϕz[k]))
    end
    return rs
end

function randomSineGPU_2D(y,z, Ny,Nz)
#Random coefficients (gaussian-distributed) for phase (ϕ) and amplitude (A)

amplitude = FileIO.load("cosineVariables.jld2","amplitude")
phase = FileIO.load("cosineVariables.jld2","phase")

#ϕx = phase[:,1]
#Ax = amplitude[:,1]
ϕy = phase[:,2]
Ay = amplitude[:,2]
ϕz = phase[:,3]
Az = amplitude[:,3]

rs = 0
Threads.@threads for k = 1:400
        rs = rs + 1/400 .*(Ay[k].*sin.(2*π.*y .* (k) / Ny .+ ϕy[k]) .*
        Az[k].*sin.(2*π.*z .* (k) / Nz .+ ϕz[k]))
    end
    return rs
end
=#
