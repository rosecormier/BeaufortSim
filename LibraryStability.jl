module Stratification
   export N2_lower_bound, grav_stable
end

function N2_lower_bound(σr, σz, f, U)
   lower_bound = exp(1) * U * f * σr / (σz^2)
end

function grav_stable(σr, σz, f, U, N2)
   lower_bound = N2_lower_bound(σr, σz, f, U)
   if N2 < lower_bound
      print("Warning: system is gravitationally unstable.")
   end
end
