module Stability
   export U_upper_bound, check_inert_stability, 
          N2_lower_bound, check_grav_stability
end

function U_upper_bound(σr, f)
   upper_bound = f * σr / (2 * exp(1))
end

function check_inert_stability(σr, σz, f, U, x, y, z)
   r2_arr = @. x^2 + y^2
   z2_arr = transpose(z.^2 .* ones(Float64, (1, length(r2_arr))))
   ζa_b = @. (f + (2*U/σr) * (r2_arr/(σr^2) - 1)
              * exp(1 - r2_arr/(σr^2) - z2_arr/(σz^2)))
   if !all(ζa_b .> 0)
      print("Warning: system is inertially unstable.")
   end
end

function N2_lower_bound(σr, σz, f, U)
   lower_bound = exp(1) * U * f * σr / (σz^2)
end

function check_grav_stability(σr, σz, f, U, N2)
   lower_bound = N2_lower_bound(σr, σz, f, U)
   if N2 < lower_bound
      print("Warning: system is gravitationally unstable.")
   end
end
