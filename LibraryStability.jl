module Stability
   export check_inert_stability, check_grav_stability, compute_Bu
end

function check_inert_stability(σr, σz, f, U, x, y, z)
   r2   = @. x^2 + y^2
   z2   = transpose(z.^2 .* ones(Float64, (1, length(r2))))
   ζa_b = @. f - ((2 * sqrt(2) * U / σr) 
		  * exp((1/2) - r2/(σr^2) - z2/(σz^2))
		  * (1 - r2/(σr^2)))
   if !all(ζa_b .> 0)
      print("Warning: system is inertially unstable.")
   end
end

function check_grav_stability(σr, σz, f, U, N2, x, y, z)
   r2   = @. x^2 + y^2
   z2   = transpose(z.^2 .* ones(Float64, (1, length(r2))))
   ∂z_b = @. N2 + ((sqrt(2)*f*U*σr/(σz^2)) * exp((1/2) - z2/(σz^2)) 
		* (1 - exp(-r2/(σr^2))) * (1 - 2*z2/(σz^2)))
   if !all(∂z_b .> 0)
      print("Warning: system is gravitationally unstable.")
   end
end

function compute_Bu(σr, σz, f, N2)
   Bu = N2 * (σz / (f * σr))^2
end
