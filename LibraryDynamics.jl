using SpecialFunctions

function lognormal_strat(N²₀, N²_max, d_ML, z; σ = 0.5)

   if N²_max == N²₀
      N² = N²₀ #Uniform-stratification case
      b  = @. N² * z

   elseif N²_max > N²₀

      μ  = log(-d_ML) + (σ^2)/2
      A  = -d_ML * σ * sqrt(2*pi) * (N²_max - N²₀) * exp(-(σ^2))

      if z == 0.0
         N² = 0.0
	 b  = 0.0
      elseif z < 0.0
         N² = @. N²₀ - (A / (z*σ*sqrt(2*pi))) * exp(-(log(-z) - μ)^2 / (2*σ^2))
         b  = @. N²₀*z - (A/2) * (1 + erf((log(z/d_ML) - (σ^2/2)) / (σ*sqrt(2))))
      end
   end

   return N², b
end

function chebyshev_spaced_faces(i, ξ_min, ξ_max, Nξ; ξ_centre = 0.0)

   Lξ = ξ_max - ξ_min
   
   N_above_ξ_centre = (Nξ/pi) * acos(1 + (2*ξ_centre)/Lξ)

   if i <= Nξ - N_above_ξ_centre
      i_face = ξ_centre - (Lξ/2) * (1 + cos((pi/Nξ)*(i+N_above_ξ_centre)))
   elseif i > Nξ - N_above_ξ_centre
      i_face = ξ_centre + (Lξ/2) * (1 + cos((pi/Nξ)*(i+N_above_ξ_centre)))
   end

   return i_face
end
