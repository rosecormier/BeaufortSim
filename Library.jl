function ReadInputFile(filename)
    
    input_params = Dict()
    
    for (i, line) in enumerate(eachline(filename))
        
        var_len, line_len = 0, length(line)
        if (line_len > 1 && line[1] != '#')

            for char in line
                if char == '='
                    break
                else
                    var_len += 1
                end
            end

            var = strip(line[1:var_len])
            line = strip(line[var_len:end])

            val_len = 0
            for char in line
                if char == '#' 
                    break
                else
                    val_len += 1
                end
            end

            val_str = strip(line[2:val_len]) #Index 1 is "="

            val = tryparse(Int, val_str)
            if isnothing(val)
                val = parse(Float64, val_str)
            end

            input_params[var] = val

        end    
    end
    
    input_params
    
end
    
#= For visualization =#
    
using Printf
#using FileIO, JLD2
module VisualizationFunctions
    export buoyancy_C, buoyancy_L, buoyancy_R, u_background, density_from_buoyancy, ωt, ω_b
    export ζ, ζ_2D, ∇b, ∇b_2D, ertel_q, ertel_q_2D, BestFit, randomSine,randomSineGPU 
end

function u_background(x,y,z,Umax,D,Lⱼ,z0,y0)
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

function ζ(u, v, w, Δx, Δy, Δz)
    ζx = (w[2:end,2:end,3:end]-w[2:end,1:end-1,3:end])./Δy .- (v[2:end,3:end,2:end]-v[2:end,3:end,1:end-1])./Δz
    ζy = (u[2:end,2:end,2:end]-u[2:end,2:end,1:end-1])./Δz - (w[2:end,2:end,3:end]-w[1:end-1,2:end,2:end-1])./Δx
    ζz = (v[2:end,3:end,2:end]-v[1:end-1,3:end,2:end])./Δx - (u[2:end,2:end,2:end]-u[2:end,1:end-1,2:end])./Δy
    return ζz #, ζx, ζy, ζz
end

function ζ_2D(u, v, w, Δx, Δy, Δz)
#    ζx = (w[1,2:end,3:end]-w[1,1:end-1,3:end])./Δy .- (v[1,3:end,2:end]-v[1,3:end,1:end-1])./Δz
#    ζy = (u[1,2:end,2:end]-u[1,2:end,1:end-1])./Δz
    #ζz = - (u[1,2:end,2:end]-u[1,1:end-1,2:end])./Δy
    ζz = (v[1,2:end,2:end]-v[1,1:end-1,2:end])./Δx -(u[1,2:end,2:end]-u[1,1:end-1,2:end])./Δy
    return ζz #, ζx, ζy, ζz
end


function ∇b(b, Δx, Δy, Δz)
    ∂x_b =  (b[2:end,2:end,2:end]-b[1:end-1,2:end,2:end])/Δx
    ∂y_b =  (b[2:end,2:end,2:end]-b[2:end,1:end-1,2:end])/Δy
    ∂z_b =  (b[2:end,2:end,2:end]-b[2:end,2:end,1:end-1])/Δz
    db_dX = @. ∂x_b + ∂y_b + ∂z_b
    return ∂z_b
end

function ∇b_2D(b, Δy, Δz)
    ∂x_b =  0.0
    ∂y_b =  (b[2:end,2:end]-b[1:end-1,2:end])/Δy
    ∂z_b =  (b[2:end,2:end]-b[2:end,1:end-1])/Δz
    db_dX = @. ∂x_b + ∂y_b + ∂z_b
    return db_dX#, ∂x_b, ∂y_b, ∂z_b
end

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
