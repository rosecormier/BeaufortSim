module Stability
   export U_upper_bound, check_inert_stability, 
          N2_lower_bound, check_grav_stability,
	  compute_Bu
end

function U_upper_bound(σr, f)
   upper_bound = f * σr / (2 * sqrt(2 * exp(1)))
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

function N2_lower_bound(σr, σz, f, U)
   lower_bound = 2 * sqrt(2) * f * U * σr / (exp(1) * (σz^2))
end

function check_grav_stability(σr, σz, f, U, N2)
   lower_bound = N2_lower_bound(σr, σz, f, U)
   if N2 <= lower_bound
      print("Warning: system is gravitationally unstable.")
   end
end

function compute_Bu(σr, σz, f, N2)
   Bu = N2 * (σz / (f * σr))^2
end
