using Oceananigans.BoundaryConditions
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

function chebyshev_spaced_faces(i, ξ_min, Nξ; ξ_max = 0.0, ξ_centre = 0.0)

   Lξ = ξ_max - ξ_min
   
   arg_shift = asin(1 + (ξ_centre/Lξ))

   N_below_ξ_centre = (Nξ*pi) / (2*(pi-arg_shift)) 

   if i <= N_below_ξ_centre
      i_face = ξ_centre + Lξ * (sin((pi-arg_shift)*i/Nξ) - 1)
   elseif i > N_below_ξ_centre
      i_face = ξ_centre - Lξ * (sin((pi-arg_shift)*i/Nξ) - 1)
   end

   return i_face
end

function bkgd_fields(f, σr, σz, U, bkgd_N²_top, bkgd_N²_bot)
   
   #this will all be cleaner if we convert to polar coords upfront; i plan to change this
   
   if σz == "infinity" #Barotropic case
  
      #b̄ = (x, y, z) -> lognormal_strat(N²₀, N²_max, d_ML, z)[2]
      #ū = (x, y, z) -> ((sqrt(2)*U*y/σr)
      #                  * exp((1/2) - (x^2 + y^2)/(σr^2)))
      #v̄ = (x, y, z) -> -((sqrt(2)*U*x/σr)
      #                   * exp((1/2) - (x^2 + y^2)/(σr^2)))

      b̄ = (x, z) -> lognormal_strat(N²₀, N²_max, d_ML, z)[2]
      ū = (x, z) -> 0
      v̄ = (x, z) -> -((sqrt(2)*U*x/σr)
                         * exp((1/2) - (x^2)/(σr^2)))
	
      b̄z_top = (x, t) -> bkgd_N²_top
      b̄z_bot = (x, t) -> bkgd_N²_bot
   
      #b̄z_top = (x, y, t) -> bkgd_N²_top
      #b̄z_bot = (x, y, t) -> bkgd_N²_bot

   else #Baroclinic case
      
      b̄ = (x, y, z) -> (lognormal_strat(N²₀, N²_max, d_ML, z)[2]
                 + ((sqrt(2)*f*U*σr*z/(σz^2))
                    * exp((1/2) - (z/σz)^2)
                    * (1 - exp(-(x^2 + y^2)/(σr^2)))
                    * (1 - ((sqrt(2)*U/(f*σr)) * exp((1/2) - (z/σz)^2)
                             * (1 + exp(-(x^2 + y^2)/(σr^2)))
                           )
                      )
                   )
                )
      ū = (x, y, z) -> ((sqrt(2)*U*y/σr)
                        * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))
      v̄ = (x, y, z) -> -((sqrt(2)*U*x/σr)
                         * exp((1/2) - (x^2 + y^2)/(σr^2) - (z/σz)^2))

      b̄z_top = (x, y, t) -> (bkgd_N²_top
                                .+ (sqrt(2)*f*U*σr/(σz^2)
                                   * exp(1/2)
                                   * (1 - exp(-(x^2 + y^2)/(σr^2)))))
      b̄z_bot = (x, y, t) -> (bkgd_N²_bot
                                .+ (sqrt(2)*f*U*σr/(σz^2)
                                   * exp((1/2) - (Lz/σz)^2)
                                   * (1 - exp(-(x^2 + y^2)/(σr^2)))
                                   * (1 - 2 * (Lz/σz)^2)))
   end

   b̄_BCs = FieldBoundaryConditions(top    = GradientBoundaryCondition(b̄z_top),
				   bottom = GradientBoundaryCondition(b̄z_bot))

   return b̄, ū, v̄, b̄_BCs
end
