using Printf
module VisualizationFunctions
    #export buoyancy_C, buoyancy_L, buoyancy_R, u_background, density_from_buoyancy, ωt, ω_b
    export ζ, ζ_2D, ∇b, ∇b_2D, ∂r_fq #, ertel_q, ertel_q_2D, BestFit, randomSine,randomSineGPU 
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
=#
function ζ(u, v, w, Δx, Δy, Δz)
    ζx = (w[2:end,2:end,3:end]-w[2:end,1:end-1,3:end])./Δy .- (v[2:end,3:end,2:end]-v[2:end,3:end,1:end-1])./Δz
    ζy = (u[2:end,2:end,2:end]-u[2:end,2:end,1:end-1])./Δz - (w[2:end,2:end,3:end]-w[1:end-1,2:end,2:end-1])./Δx
    ζz = (v[2:end,3:end,2:end]-v[1:end-1,3:end,2:end])./Δx - (u[2:end,2:end,2:end]-u[2:end,1:end-1,2:end])./Δy
    return ζz #, ζx, ζy, ζz
end

function ζ_2D(u, v, w, Δx, Δy, Δz, x_idx, y_idx, z_idx)
    if !isnothing(x_idx)
        ζx = (w[x_idx,2:end,2:end]-w[x_idx,1:end-1,2:end])./Δy -(v[x_idx,2:end,2:end]-v[x_idx,2:end,1:end-1])./Δz
        return ζx
    elseif !isnothing(y_idx)
        ζy = (u[2:end,y_idx,2:end]-u[2:end,y_idx,1:end-1])./Δz -(w[2:end,y_idx,2:end]-w[1:end-1,y_idx,2:end])./Δx
	return ζy
    elseif !isnothing(z_idx)
        ζz = (v[2:end,2:end,z_idx]-v[1:end-1,2:end,z_idx])./Δx -(u[2:end,2:end,z_idx]-u[2:end,1:end-1,z_idx])./Δy
        return ζz
    end
end

function ∇b(b, Δx, Δy, Δz)
    ∂x_b =  (b[2:end,2:end,2:end]-b[1:end-1,2:end,2:end])/Δx
    ∂y_b =  (b[2:end,2:end,2:end]-b[2:end,1:end-1,2:end])/Δy
    ∂z_b =  (b[2:end,2:end,2:end]-b[2:end,2:end,1:end-1])/Δz
    db_dX = @. ∂x_b + ∂y_b + ∂z_b
    return ∂z_b
end

function ∇b_2D(b, Δx, Δy, Δz, excluded_dim)
    if excluded_dim == 'x'
        ∂x_b =  0.0
        ∂y_b =  (b[2:end,2:end]-b[1:end-1,2:end])/Δy
        ∂z_b =  (b[2:end,2:end]-b[2:end,1:end-1])/Δz
    elseif excluded_dim == 'y'
        ∂x_b =  (b[2:end,2:end]-b[1:end-1,2:end])/Δx
        ∂y_b =  0.0
        ∂z_b =  (b[2:end,2:end]-b[2:end,1:end-1])/Δz
    elseif excluded_dim == 'z'
        ∂x_b =  (b[2:end,2:end]-b[1:end-1,2:end])/Δx
        ∂y_b =  (b[2:end,2:end]-b[2:end,1:end-1])/Δy
        ∂z_b =  0.0
    end
    db_dX = @. ∂x_b + ∂y_b + ∂z_b
    return db_dX
end

function ∂r_fq(fq, Δx, Δy, xC, yC)
    ∂x_fq = (fq[2:end,2:end,:].-fq[1:end-1,2:end,:])./Δx
    ∂y_fq = (fq[2:end,2:end,:].-fq[2:end,1:end-1,:])./Δy
    r     = (xC.^2 .+ yC.^2).^0.5
    ∂r_fq = (xC.*∂x_fq .+ yC.*∂y_fq)./r
    return ∂r_fq
end

#=
function ertel_q(u, v, w, b, Δx, Δy, Δz)
    ζt,ζx,ζy,ζz = ζ(u, v, w, Δx, Δy, Δz)
    ∇b, ∂x_b, ∂y_b, ∂z_b = ∇b(b, Δx, Δy, Δz)
    qvert = (f .+ ζz).* ∂z_b
    q_bc = @. ζx * ∂x_b + ζy * ∂y_b
    qt = @. q_bc .+ qvert
    return qt, qvert, q_bc
end

function ertel_q_2D(u, v, w, b, Δy, Δz)
    ζt, ζx, ζy, ζz = ζ(u, v, w, Δy, Δz)
    db_dX, ∂x_b, ∂y_b, ∂z_b = ∇b(b, Δy, Δz)
    qvert = (f .+ ζz).*∂z_b
    q_bc = @. ζx * ∂x_b + ζy * ∂y_b
    qt = @. q_bc .+ qvert
    return qt, qvert, q_bc
end

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
