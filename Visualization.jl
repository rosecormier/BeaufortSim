include("LibraryDynamics.jl")
include("LibraryVisualization.jl")

using Adapt, CairoMakie, CommonDataModel, CUDA, DataStructures
using Glob, LaTeXStrings, NCDatasets
using Oceananigans
using Oceananigans.Fields
using Oceananigans.OutputReaders
using OffsetArrays: no_offset_view
using Polynomials: fit
using Printf

update_theme!(theme_latexfonts(), fontsize = 16)

####################

function visualize_B_U_Q_Ψ_vs_r_and_z(U, grid, f, σr, σz, N²_far, 
                                      doubleTanhParams, ambientStrat, Nr, Nz, 
                                      Lr, Lz)

   r = range(0, stop = Lr, length = Nr ÷ 2)
   z = range(-Lz, stop = 0, length = Nz)
   
   B_function  = bkgd_B_cylindrical_coords(f, U, σr/1000, σz, N²_far, 
                                           doubleTanhParams, ambientStrat)
   Q_function  = bkgd_Q_cylindrical_coords(f, U, σr/1000, σz, N²_far, 
                                           doubleTanhParams, ambientStrat)
   Uφ_function = bkgd_Uφ_cylindrical_coords(σr/1000, σz, U)
   Ψ_function  = bkgd_Ψ_cylindrical_coords(σr/1000, σz, U)

   fig  = Figure(size = (1200, 600))
   ax_B = Axis(fig[1, 1], xlabel = L"$r$ [km]", ylabel = L"$z$ [m]", 
               title = "Background buoyancy with QG-PV contours")
   ax_U = Axis(fig[1, 2], xlabel = L"$r$ [km]", ylabel = L"$z$ [m]", 
               title = "Background velocity with QG-streamfunction contours")

   hm_B = heatmap!(ax_B, r/1000, z, B_function, colormap = :balance)
   hm_U = heatmap!(ax_U, r/1000, z, Uφ_function, colormap = Reverse(:Blues), colorrange = (-U, 0))
   
   contour!(ax_B, r[2:end]/1000, z, Q_function, color = :yellow, levels = 10)
   contour!(ax_U, r/1000, z, Ψ_function, color = :yellow, levels = 10)
   
   Colorbar(fig[2, 1], hm_B, tickformat = "{:.1e}", label = "m/s²", vertical = false, width = Relative(3/4))
   Colorbar(fig[2, 2], hm_U, tickformat = "{:.1e}", label = "m/s", vertical = false, width = Relative(3/4))

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "background_fields.png"), fig)
end

function visualize_B_and_N²_vs_z(B, grid, x_idx, y_idx, doubleTanhParams, f, 
                                 σr, σz, U, N²_far; Hz = 3)

   B_total    = no_offset_view(adapt(Array, B)
                              )[x_idx, y_idx, Hz:length(grid.z.cᵃᵃᶜ) - Hz]
   b_TWB      = no_offset_view(TWB_b_field(grid, f, σr, σz, U)
                              )[x_idx, y_idx, Hz:length(grid.z.cᵃᵃᶜ) - Hz]
   ∂B∂z_total = no_offset_view(∂b∂z_field(B, grid))[x_idx, y_idx, 
                                                    Hz:(length(grid.z.cᵃᵃᶜ)
                                                        - Hz)
                                                   ]
   ∂b∂z_TWB   = no_offset_view(TWB_∂b∂z_field(grid, N²_far, f, σr, σz, U)
                              )[x_idx, y_idx, Hz:length(grid.z.cᵃᵃᶜ) - Hz]
   
   x = no_offset_view(adapt(Array, grid.xᶜᵃᵃ))[x_idx]
   y = no_offset_view(adapt(Array, grid.yᵃᶜᵃ))[y_idx]
   z = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ))[Hz:(length(grid.z.cᵃᵃᶜ) - Hz)]
   
   fig   = Figure(size = (1400, 700))
   ax_B  = Axis(fig[2, 1], xlabel = L"$B$ [m/s$^{2}$]", ylabel = L"$z$ [m]")
   ax_N2 = Axis(fig[2, 2], xlabel = L"$N^2$ [s$^{-2}$]", ylabel = L"$z$ [m]")
   
   lines!(ax_B, B_total, z, label = "Total")
   scatter!(ax_B, B_total, z, label = "Total (at gridpoints)")
   lines!(ax_B, buoyancyDoubleTanh(z, doubleTanhParams), z, 
          label = "From double-tanh function")
   lines!(ax_B, b_TWB, z, label = "Thermal-wind contribution")
   lines!(ax_B, N²_far .* z, z, label = "Linear term")
   
   lines!(ax_N2, ∂B∂z_total, z, label = "Total")
   scatter!(ax_N2, ∂B∂z_total, z, label = "Total (at gridpoints)")
   lines!(ax_N2, N²DoubleTanh(z, doubleTanhParams), z, 
          label = "From double-tanh function")
   lines!(ax_N2, ∂b∂z_TWB, z, label = "Thermal-wind contribution")
   lines!(ax_N2, N²_far .+ 0*z, z, label = "Linear term")
   
   fig[2, 3] = Legend(fig, ax_B)
   fig[1, 1] = Label(fig, "Background buoyancy", fontsize = 24, 
                      tellwidth = false)
   fig[1, 2] = Label(fig, L"Background $N^2$", fontsize = 24, tellwidth = false)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "b_and_N2_profiles.png"), fig)
end

function visualize_norms(datetime; 
		idxStartLinGrowth_b = 2, idxEndLinGrowth_b = -1,
    idxStartLinGrowth_ur = nothing, idxEndLinGrowth_ur = nothing,
		idxStartLinGrowth_uφ = nothing, idxEndLinGrowth_uφ = nothing,
		idxStartLinGrowth_ux = nothing, idxEndLinGrowth_ux = nothing,
		idxStartLinGrowth_uy = nothing, idxEndLinGrowth_uy = nothing,
		idxStartLinGrowth_uz = nothing, idxEndLinGrowth_uz = nothing,
		idxStartPlot = 2, idxEndPlot = -1)

   scalars_ds, times, Nt, chron_idcs = open_scalars_dataset(
                      glob("./Output/scalars_$(datetime)*"))
   
   if idxEndPlot < 0
      idxEndPlot = Nt + idxEndPlot #Make idxEndPlot a positive integer
   end
   
   #Retain only necessary portion of 'times' data and update Nt accordingly
   times = times[idxStartPlot:idxEndPlot]
   Nt    = length(times)

   #Load data, sorted chronologically and then restricted to plot interval
   b′_norm  = scalars_ds[:b′_norm][chron_idcs][idxStartPlot:idxEndPlot]
   ux′_norm = scalars_ds[:ux′_norm][chron_idcs][idxStartPlot:idxEndPlot]
   uy′_norm = scalars_ds[:uy′_norm][chron_idcs][idxStartPlot:idxEndPlot]
   ur′_norm = scalars_ds[:ur′_norm][chron_idcs][idxStartPlot:idxEndPlot]
   uφ′_norm = scalars_ds[:uφ′_norm][chron_idcs][idxStartPlot:idxEndPlot]
   uz′_norm = scalars_ds[:uz′_norm][chron_idcs][idxStartPlot:idxEndPlot]
   
   fig_cyl   = Figure(size = (1200, 700))
   ax_b_cyl  = Axis(fig_cyl[2, 1]; title = L"Norm of $b'$", 
                    xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]",
                    yscale = log10)
   ax_ur     = Axis(fig_cyl[2, 2]; title = L"Norm of $u_r'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_r'||$ [m/s]",
                    yscale = log10)
   ax_uφ     = Axis(fig_cyl[3, 1]; title = L"Norm of $u_{\phi}'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_{\phi}'||$ [m/s]",
                    yscale = log10)
   ax_uz_cyl = Axis(fig_cyl[3, 2]; title = L"Norm of $u_z'$",
                    xlabel = L"$t$ [days]", ylabel = L"$||u_z'||$ [m/s]",
                    yscale = log10)

   fig_Cart   = Figure(size = (1200, 700))
   ax_b_Cart  = Axis(fig_Cart[2, 1]; title = L"Norm of $b'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||b'||$ [m/s^2]",
                     yscale = log10)
   ax_ux      = Axis(fig_Cart[2, 2]; title = L"Norm of $u_x'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||u_x'||$ [m/s]",
                     yscale = log10)
   ax_uy      = Axis(fig_Cart[3, 1]; title = L"Norm of $u_y'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||u_y'||$ [m/s]",
                     yscale = log10)
   ax_uz_Cart = Axis(fig_Cart[3, 2]; title = L"Norm of $u_z'$",
                     xlabel = L"$t$ [days]", ylabel = L"$||u_z'||$ [m/s]",
                     yscale = log10)

   scatter!(ax_b_cyl, times, b′_norm, color = :black)
   scatter!(ax_ur, times, ur′_norm, color = :black)
   scatter!(ax_uφ, times, uφ′_norm, color = :black)
   scatter!(ax_uz_cyl, times, uz′_norm, color = :black)

   scatter!(ax_b_Cart, times, b′_norm, color = :black)
   scatter!(ax_ux, times, ux′_norm, color = :black)
   scatter!(ax_uy, times, uy′_norm, color = :black)
   scatter!(ax_uz_Cart, times, uz′_norm, color = :black)

   if idxEndLinGrowth_b > 0
      growthIdcs_b = (idxStartLinGrowth_b, idxEndLinGrowth_b)
   elseif idxEndLinGrowth_b < 0
      growthIdcs_b = (idxStartLinGrowth_b, Nt + idxEndLinGrowth_b)
   end

   function updateGrowthIdcs!(idxStart, idxEnd)
      
      if isnothing(idxStart) #Default to the start index used for b
         idxStart = growthIdcs_b[1]
      end

      if isnothing(idxEnd) #Default to the end index used for b
         idxEnd = growthIdcs_b[2]
      end

      return idxStart, idxEnd
   end

   idxStartLinGrowth_ur, idxEndLinGrowth_ur = updateGrowthIdcs!(idxStartLinGrowth_ur, 
                                                                idxEndLinGrowth_ur)
   idxStartLinGrowth_uφ, idxEndLinGrowth_uφ = updateGrowthIdcs!(idxStartLinGrowth_uφ, 
                                                                idxEndLinGrowth_uφ)
   idxStartLinGrowth_ux, idxEndLinGrowth_ux = updateGrowthIdcs!(idxStartLinGrowth_ux,
                                                                idxEndLinGrowth_ux)
   idxStartLinGrowth_uy, idxEndLinGrowth_uy = updateGrowthIdcs!(idxStartLinGrowth_uy, 
                                                                idxEndLinGrowth_uy)
   idxStartLinGrowth_uz, idxEndLinGrowth_uz = updateGrowthIdcs!(idxStartLinGrowth_uz,
                                                                idxEndLinGrowth_uz)

   b′NormFitInterval  = b′_norm[growthIdcs_b[1]:growthIdcs_b[2]]
   ur′NormFitInterval = ur′_norm[idxStartLinGrowth_ur:idxEndLinGrowth_ur]
   uφ′NormFitInterval = uφ′_norm[idxStartLinGrowth_uφ:idxEndLinGrowth_uφ]
   ux′NormFitInterval = ux′_norm[idxStartLinGrowth_ux:idxEndLinGrowth_ux]
   uy′NormFitInterval = uy′_norm[idxStartLinGrowth_uy:idxEndLinGrowth_uy]
   uz′NormFitInterval = uz′_norm[idxStartLinGrowth_uz:idxEndLinGrowth_uz]

   b′NormLinearFitParams  = fit(times[growthIdcs_b[1]:growthIdcs_b[2]], 
                                log.(b′NormFitInterval), 1, var = :times)
   ur′NormLinearFitParams = fit(times[idxStartLinGrowth_ur:idxEndLinGrowth_ur], 
                                log.(ur′NormFitInterval), 1, var = :times)
   uφ′NormLinearFitParams = fit(times[idxStartLinGrowth_uφ:idxEndLinGrowth_uφ], 
                                log.(uφ′NormFitInterval), 1, var = :times)
   ux′NormLinearFitParams = fit(times[idxStartLinGrowth_ux:idxEndLinGrowth_ux], 
                                log.(ux′NormFitInterval), 1, var = :times)
   uy′NormLinearFitParams = fit(times[idxStartLinGrowth_uy:idxEndLinGrowth_uy], 
                                log.(uy′NormFitInterval), 1, var = :times)
   uz′NormLinearFitParams = fit(times[idxStartLinGrowth_uz:idxEndLinGrowth_uz], 
                                log.(uz′NormFitInterval), 1, var = :times)
   
   @printf("Empirical growth rate:\n From b′-norm: %.5f per day\n From ur′-norm: %.5f per day\n From uφ′-norm: %.5f per day\n From ux′-norm: %.5f per day\n From uy′-norm: %.5f per day\n From uz′-norm: %.5f per day\n",
	    b′NormLinearFitParams[1], ur′NormLinearFitParams[1],
	    uφ′NormLinearFitParams[1], ux′NormLinearFitParams[1],
	    uy′NormLinearFitParams[1], uz′NormLinearFitParams[1])

   @inline linearFunction(fitParams, tFitInterval; offset = 2) = @. offset * 
   				exp(fitParams[0] + fitParams[1] * tFitInterval)

   lines!(ax_b_cyl, times[growthIdcs_b[1]:growthIdcs_b[2]],
          linearFunction(b′NormLinearFitParams, 
                         times[growthIdcs_b[1]:growthIdcs_b[2]])
         )
   lines!(ax_ur, times[idxStartLinGrowth_ur:idxEndLinGrowth_ur], 
          linearFunction(ur′NormLinearFitParams, 
                         times[idxStartLinGrowth_ur:idxEndLinGrowth_ur])
         )
   lines!(ax_uφ, times[idxStartLinGrowth_uφ:idxEndLinGrowth_uφ], 
          linearFunction(uφ′NormLinearFitParams, 
                         times[idxStartLinGrowth_uφ:idxEndLinGrowth_uφ])
         )
   lines!(ax_uz_cyl, times[idxStartLinGrowth_uz:idxEndLinGrowth_uz], 
          linearFunction(uz′NormLinearFitParams, 
                         times[idxStartLinGrowth_uz:idxEndLinGrowth_uz])
         )

   lines!(ax_b_Cart, times[growthIdcs_b[1]:growthIdcs_b[2]], 
          linearFunction(b′NormLinearFitParams, 
                         times[growthIdcs_b[1]:growthIdcs_b[2]])
         )
   lines!(ax_ux, times[idxStartLinGrowth_ux:idxEndLinGrowth_ux], 
          linearFunction(ux′NormLinearFitParams, 
                         times[idxStartLinGrowth_ux:idxEndLinGrowth_ux])
         )
   lines!(ax_uy, times[idxStartLinGrowth_uy:idxEndLinGrowth_uy], 
          linearFunction(uy′NormLinearFitParams, 
                         times[idxStartLinGrowth_uy:idxEndLinGrowth_uy])
         )
   lines!(ax_uz_Cart, times[idxStartLinGrowth_uz:idxEndLinGrowth_uz], 
          linearFunction(uz′NormLinearFitParams, 
                         times[idxStartLinGrowth_uz:idxEndLinGrowth_uz])
         )

   mkpath("./Plots") #Make visualization directory if nonexistent

   fig_cyl[1, 1:2]  = Label(fig_cyl, 
			   "Norms of perturbation fields with best linear fits",
                           fontsize = 24, tellwidth = false)
   fig_Cart[1, 1:2] = Label(fig_Cart, 
			   "Norms of perturbation fields with best linear fits",
                           fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "norm_fields_$(datetime).png"), fig_cyl)
   save(joinpath("./Plots", "norm_Cart_fields_$(datetime).png"), fig_Cart)
   close(scalars_ds)
end

function visualize_total_energy_budgets(datetime, grid)

   outfile_list          = glob("./Output/energetics_$(datetime)*")
   ds, t, Nt, chron_idcs = open_energetics_dataset(outfile_list)

   total_KE            = ds[:total_KE][chron_idcs]
   total_KE_adv_flux   = ds[:total_KE_adv_flux][chron_idcs]
   total_KE_production = ds[:total_KE_production][chron_idcs] #0.5
   total_pressure_work = ds[:total_pressure_work][chron_idcs] #0.5
   total_PE            = ds[:total_PE][chron_idcs]
   total_b_adv_flux    = ds[:total_b_adv_flux][chron_idcs]
   total_gravity_work  = ds[:total_gravity_work][chron_idcs]
   
   total_ME = total_KE + total_PE #Mechanical energy
   
   fig_ME = Figure(size = (1200, 700))
   ax_ME  = Axis(fig_ME[2, 1]; xlabel = "Time [days]",
                 ylabel = L"[m$^5$ s$^{-2}$]")
                 
   scatter!(ax_ME, t, 0.5 * total_ME, label = "50% of total ME", color = :black)
   scatter!(ax_ME, t, total_KE, label = "Total KE", color = :green)
   scatter!(ax_ME, t, total_PE, label = "Total PE", color = :purple)
   
   fig_ME[1, 1] = Label(fig_ME, "Total-ME budget", fontsize = 24, 
                        tellwidth = false)

   axislegend(ax_ME)
   save(joinpath("./Plots", "MEbudget_$(datetime).png"), fig_ME)
   
   #Convert time to s to compute time-derivative of KE
   dt_KE_total = order1_forward_difference(86400 * t, total_KE)
   
   fig_KE_budget = Figure(size = (1200, 700))
   ax_KE_budget  = Axis(fig_KE_budget[2, 1]; xlabel = "Time [days]",
                        ylabel = L"[m$^5$ s$^{-3}$]")
                     
   scatter!(ax_KE_budget, t[1:end-1], dt_KE_total, 
            label = "Time derivative of total KE", color = :yellowgreen)
   scatter!(ax_KE_budget, t, total_KE_adv_flux, label = "Advective flux", color = :royalblue)
   scatter!(ax_KE_budget, t, total_pressure_work, label = "Pressure work", color = :sienna)
   scatter!(ax_KE_budget, t, total_KE_production, label = "Buoyant production", color = :orange)
   scatter!(ax_KE_budget, t[1:end-1], 
            (dt_KE_total - total_KE_adv_flux[1:end-1] 
             - total_pressure_work[1:end-1] - total_KE_production[1:end-1]), 
            label = "Residual", color = :red)

   fig_KE_budget[1, 1] = Label(fig_KE_budget, "Terms in total-KE budget",
                               fontsize = 24, tellwidth = false)

   axislegend(ax_KE_budget)
   save(joinpath("./Plots", "KEbudget_$(datetime).png"), fig_KE_budget)

   #Convert time to s to compute time-derivative of PE
   dt_PE_total = order1_forward_difference(86400 * t, total_PE)

   fig_PE_budget = Figure(size = (1200, 700))
   ax_PE_budget  = Axis(fig_PE_budget[2, 1]; xlabel = "Time [days]",
                        ylabel = L"[m$^5$ s$^{-3}$]")
                     
   scatter!(ax_PE_budget, t[1:end-1], dt_PE_total, 
            label = "Time derivative of total PE", color = :yellowgreen)
   scatter!(ax_PE_budget, t, total_b_adv_flux, label = "Internal sources", color = :royalblue)
   scatter!(ax_PE_budget, t, total_gravity_work, label = "Gravity work", color = :sienna)
   scatter!(ax_PE_budget, t, -total_KE_production, label = "Buoyant production (of KE)", color = :orange)
   scatter!(ax_PE_budget, t[1:end-1], 
            (dt_PE_total - total_b_adv_flux[1:end-1] 
             - total_gravity_work[1:end-1] + total_KE_production[1:end-1]), 
            label = "Residual", color = :red)

   fig_PE_budget[1, 1] = Label(fig_PE_budget, "Terms in total-PE budget",
                               fontsize = 24, tellwidth = false)

   axislegend(ax_PE_budget, position = :lb)
   save(joinpath("./Plots", "PEbudget_$(datetime).png"), fig_PE_budget)

   close(ds)
end

function visualize_PKE(datetime, grid)

   outfile_list          = glob("./Output/energetics_$(datetime)*")
   ds, t, Nt, chron_idcs = open_energetics_dataset(outfile_list)
   
   initialKE = ds[:total_KE][chron_idcs][1]

   PKE          = ds[:total_PKE][chron_idcs] / initialKE
   PAPE_to_PKE  = ds[:total_PAPE_to_PKE][chron_idcs] / initialKE
   BTI_transfer = ds[:total_BTI_transfer][chron_idcs] / initialKE
   BCI_transfer = ds[:total_BCI_transfer][chron_idcs] / initialKE
   
   PKE          = PKE .- PKE[1]
   PAPE_to_PKE  = PAPE_to_PKE .- PAPE_to_PKE[1]
   BTI_transfer = BTI_transfer .- BTI_transfer[1]
   BCI_transfer = BCI_transfer .- BCI_transfer[1]
   
   gyre_PKE          = ds[:gyre_PKE][chron_idcs]
   gyre_PAPE_to_PKE  = ds[:gyre_PAPE_to_PKE][chron_idcs]
   gyre_BTI_transfer = ds[:gyre_BTI_transfer][chron_idcs]
   gyre_BCI_transfer = ds[:gyre_BCI_transfer][chron_idcs]
   
   gyre_PKE          = gyre_PKE .- gyre_PKE[1]
   gyre_PAPE_to_PKE  = gyre_PAPE_to_PKE .- gyre_PAPE_to_PKE[1]
   gyre_BTI_transfer = gyre_BTI_transfer .- gyre_BTI_transfer[1]
   gyre_BCI_transfer = gyre_BCI_transfer .- gyre_BCI_transfer[1]

   fig_PKEtotal = Figure(size = (1200, 700))
   ax_PKEtotal  = Axis(fig_PKEtotal[2, 1]; xlabel = "Time [days]",
                       ylabel = "pKE per initial KE", yscale = log10)

   scatter!(ax_PKEtotal, t[2:end], PKE[2:end], color = :black)

   fig_PKEtotal[1, 1] = Label(fig_PKEtotal,
	                            "Ratio of volume-integrated perturbation kinetic energy to volume-integrated initial kinetic energy",
                     	        fontsize = 24, tellwidth = false)
 
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "PKEtotal_$(datetime).png"), fig_PKEtotal)
   
   fig_PKEgyre = Figure(size = (1200, 700))
   ax_PKEgyre  = Axis(fig_PKEgyre[2, 1]; xlabel = "Time [days]",
                      ylabel = "PKE in gyre region", yscale = log10)

   scatter!(ax_PKEgyre, t[2:end], gyre_PKE[2:end], color = :black)

   fig_PKEgyre[1, 1] = Label(fig_PKEgyre,
                             "Perturbation kinetic energy integrated over gyre region",
                     	       fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "PKEgyre_$(datetime).png"), fig_PKEgyre)
   
   #Convert time to s to compute time-derivative of PKE
   dt_PKE_total = order1_forward_difference(86400 * t, PKE)
   
   fig_budget = Figure(size = (1200, 700))
   ax_budget  = Axis(fig_budget[2, 1]; xlabel = "Time [days]",
                     ylabel = L"[s$^{-1}$]")

   scatter!(ax_budget, t[1:end-1], dt_PKE_total, 
            label = "Time derivative of total PKE", color = :yellowgreen)
   scatter!(ax_budget, t, PAPE_to_PKE, label = "PAPE to PKE", color = :navy)
   scatter!(ax_budget, t, BTI_transfer, label = "BTI transfer", 
            color = :hotpink)
   scatter!(ax_budget, t, BCI_transfer, label = "BCI transfer", 
            color = :darkgreen)
   scatter!(ax_budget, t[1:end-1], 
           (dt_PKE_total - PAPE_to_PKE[1:end-1] - BTI_transfer[1:end-1]
            - BCI_transfer[1:end-1]), 
           label = "Residual", color = :red)
   axislegend(ax_budget)

   fig_budget[1, 1] = Label(fig_budget, "Terms in PKE budget", fontsize = 24,
                            tellwidth = false)

   save(joinpath("./Plots", "PKEbudget_$(datetime).png"), fig_budget)
   
   #Convert time to s to compute time-derivative of PKE
   dt_PKE_gyre = order1_forward_difference(86400 * t, gyre_PKE)

   fig_gyreBudget = Figure(size = (1200, 700))
   ax_gyreBudget  = Axis(fig_gyreBudget[2, 1]; xlabel = "Time [days]",
                         ylabel = L"[s$^{-1}$]")

   lines!(ax_gyreBudget, t[1:end-1], dt_PKE_gyre, 
          label = "Time derivative of total PKE", color = :yellowgreen)
   lines!(ax_gyreBudget, t, gyre_PAPE_to_PKE, label = "PAPE to PKE",
          color = :navy)
   lines!(ax_gyreBudget, t, gyre_BTI_transfer, label = "BTI transfer", 
          color = :hotpink)
   lines!(ax_gyreBudget, t, gyre_BCI_transfer, label = "BCI transfer", 
          color = :darkgreen)
   lines!(ax_gyreBudget, t[1:end-1], 
          (dt_PKE_gyre - gyre_PAPE_to_PKE[1:end-1] - gyre_BTI_transfer[1:end-1]
           - gyre_BCI_transfer[1:end-1]), 
          label = "Residual", color = :red)
   axislegend(ax_gyreBudget)

   fig_gyreBudget[1, 1] = Label(fig_gyreBudget,
			    "Terms in PKE budget, integrated over gyre region",
			    fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "PKEgyreBudget_$(datetime).png"), fig_gyreBudget)
   close(ds)
end

function visualize_b_and_ωz(datetime, Δx, Δy; 
		            x_idx = nothing, y_idx = nothing, z_idx = nothing,
		            plot_animation = true, t_idx_skip = 1)

   B  = no_offset_view(adapt(Array, B))[]
   Ux = no_offset_view(adapt(Array, Ux))[Hx+1:length(x)+Hx, 1, Hz+1:length(zC)+Hz]
   Uy = no_offset_view(adapt(Array, Uy))

   bkgd_ds = open_bkgd_dataset(bkgd_datetime)
   B       = bkgd_ds[:B][:, :, :, 1]
   Ux, Uy  = bkgd_ds[:Ux][:, :, :, 1], bkgd_ds[:Uy][:, :, :, 1]

   outfile_list             = glob("./Output/output_$(datetime)*")
   ds_f, x, y, z, times, Nt = open_dataset(outfile_list[length(outfile_list)])

   if !isnothing(x_idx)

      B_slice  = @views B[x_idx, :, :]
      ωb_slice = @views ωz(Ux, Uy, Δx, Δy; x_idx = x_idx)

      b_total_f_slice = @views adapt(Array, ds_f[:b])[x_idx, :, :, Nt]
      ω_total_f_slice = @views ωz(ds_f[:ux][:, :, :, Nt], ds_f[:uy][:, :, :, Nt], Δx, Δy;
                                  x_idx = x_idx)

      idx_kwargs           = (x_idx = x_idx,)
      nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
                                                        x_idx = x_idx)

      h_dim, v_dim, const_dim, units = y, z, "x", "km"

   elseif !isnothing(y_idx)

      B_slice  = @views B[:, y_idx, :]
      ωb_slice = @views ωz(Ux, Uy, Δx, Δy; y_idx = y_idx)

      b_total_f_slice = @views adapt(Array, ds_f[:b])[:, y_idx, :, Nt]
      ω_total_f_slice = @views ωz(ds_f[:ux][:, :, :, Nt], ds_f[:uy][:, :, :, Nt], Δx, Δy;
                                  y_idx = y_idx)

      idx_kwargs           = (y_idx = y_idx,)
      nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
                                                        y_idx = y_idx)

      h_dim, v_dim, const_dim, units = x, z, "y", "km"

   elseif !isnothing(z_idx)

      B_slice  = @views B[:, :, z_idx]
      ωb_slice = @views ωz(Ux, Uy, Δx, Δy; z_idx = z_idx)

      b_total_f_slice = @views adapt(Array, ds_f[:b])[:, :, z_idx, Nt]
      ω_total_f_slice = ωz(ds_f[:ux][:, :, :, Nt], ds_f[:uy][:, :, :, Nt], Δx, Δy;
                           z_idx = z_idx)

      idx_kwargs           = (z_idx = z_idx,)
      nearest, axis_kwargs = get_2D_spatial_axis_kwargs(x, y, z;
                                                        z_idx = z_idx)

      h_dim, v_dim, const_dim, units = x, y, "z", "m"
   end

   Δb_f_slice = b_total_f_slice .- B_slice
   Δω_f_slice = ω_total_f_slice .- ωb_slice

   lims_b_total = get_range_lims(b_total_f_slice; max_fraction = 0.75)
   lims_ω_total = get_range_lims(ω_total_f_slice; max_fraction = 0.75)

   lims_Δb = get_range_lims(Δb_f_slice;
                            max_fraction = 0.75, prescribed_max = 1e-16)
   lims_Δω = get_range_lims(Δω_f_slice;
                            max_fraction = 0.75, prescribed_max = 1e-16)

   mkpath("./Plots") #Make visualization directory if nonexistent

   #Plot static images (final frame, by default)

   fig_total   = Figure(size = (1200, 500))
   fig_perturb = Figure(size = (1200, 500))

   ax_b_total = Axis(fig_total[2, 1];
                     title = "Total buoyancy (b)", axis_kwargs...)
   ax_ω_total = Axis(fig_total[2, 3];
                     title = "Total vertical vorticity (ζ)", axis_kwargs...)

   ax_b_perturb = Axis(fig_perturb[2, 1];
                       title = "Buoyancy perturbation (b')", axis_kwargs...)
   ax_ω_perturb = Axis(fig_perturb[2, 3];
                       title = "Vertical vorticity perturbation (ζ')",
                       axis_kwargs...)

   hm_b_total = heatmap!(ax_b_total, h_dim, v_dim, b_total_f_slice,
                         colorrange = lims_b_total,
                         colormap = Reverse(:RdBu_5),
                         highclip = :red, lowclip = :blue)
   hm_ω_total = heatmap!(ax_ω_total, h_dim, v_dim, ω_total_f_slice,
                         colorrange = lims_ω_total,
                         colormap = Reverse(:RdBu_5),
                         highclip = :red, lowclip = :blue)

   hm_b_perturb = heatmap!(ax_b_perturb, h_dim, v_dim, Δb_f_slice,
                           colorrange = lims_Δb,
                           colormap = Reverse(:RdBu_5),
                           highclip = :red, lowclip = :blue)
   hm_ω_perturb = heatmap!(ax_ω_perturb, h_dim, v_dim, Δω_f_slice,
                           colorrange = lims_Δω,
                           colormap = Reverse(:RdBu_5),
                           highclip = :red, lowclip = :blue)

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ω_total, tickformat = "{:.1e}", label = "1/s")

   Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
            label = "m/s²")
   Colorbar(fig_perturb[2, 4], hm_ω_perturb, tickformat = "{:.1e}",
            label = "1/s")
   
   title_total   = @sprintf("Fields at %s = %i %s; t = %.2f days",
                            const_dim, nearest, units, times[Nt])
   title_perturb = @sprintf("Perturbation fields at %s = %i %s; t = %.2f days",
                            const_dim, nearest, units, times[Nt])

   fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                               tellwidth = false)
   fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                               tellwidth = false)

   save(joinpath("./Plots",
                 "bzeta_total_$(const_dim)$(nearest)_tf_$(datetime).png"),
        fig_total)
   save(joinpath("./Plots",
                 "bzeta_perturbs_$(const_dim)$(nearest)_tf_$(datetime).png"),
        fig_perturb)
 
   close(ds_f)
   
   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      fig_total   = Figure(size = (1200, 500))
      fig_perturb = Figure(size = (1200, 500))

      ax_b_total = Axis(fig_total[2, 1];
                        title = "Total buoyancy (b)", axis_kwargs...)
      ax_ω_total = Axis(fig_total[2, 3];
                        title = "Total vertical vorticity (ζ)", axis_kwargs...)

      ax_b_perturb = Axis(fig_perturb[2, 1];
                          title = "Buoyancy perturbation (b')", axis_kwargs...)
      ax_ω_perturb = Axis(fig_perturb[2, 3];
                          title = "Vertical vorticity perturbation (ζ')",
                          axis_kwargs...)

      video_total   = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_perturb = VideoStream(fig_perturb, format = "mp4", framerate = 6)

      file_idx = 1

      while file_idx < length(outfile_list)
        
         outfilename                 = outfile_list[file_idx] 
         ds, x, y, z, times, Nt_file = open_dataset(outfilename)

         ux = @views adapt(Array, ds[:ux])[:, :, :, :]
         uy = @views adapt(Array, ds[:uy])[:, :, :, :]

         if !isnothing(x_idx)
            b_total_tseries_slice = @views adapt(Array, ds[:b])[x_idx, :, :, :]
         elseif !isnothing(y_idx)
            b_total_tseries_slice = @views adapt(Array, ds[:b])[:, y_idx, :, :]
         elseif !isnothing(z_idx)
            b_total_tseries_slice = @views adapt(Array, ds[:b])[:, :, z_idx, :]
         end

         n = Observable(1)

         b_total_n_slice = @lift b_total_tseries_slice[:, :, $n]
         ω_total_n_slice = @lift ωz(ds[:ux][:, :, :, $n], ds[:uy][:, :, :, $n],
                                 Δx, Δy; idx_kwargs...)

         Δb_n_slice = @lift $b_total_n_slice .- B_slice
         Δω_n_slice = @lift $ω_total_n_slice .- ωb_slice

         hm_b_total = heatmap!(ax_b_total, h_dim, v_dim, b_total_n_slice,
                            colorrange = lims_b_total,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)
         hm_ω_total = heatmap!(ax_ω_total, h_dim, v_dim, ω_total_n_slice,
                            colorrange = lims_ω_total,
                            colormap = Reverse(:RdBu_5),
                            highclip = :red, lowclip = :blue)

         hm_b_perturb = heatmap!(ax_b_perturb, h_dim, v_dim, Δb_n_slice,
                              colorrange = lims_Δb,
                              colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)
         hm_ω_perturb = heatmap!(ax_ω_perturb, h_dim, v_dim, Δω_n_slice,
                              colorrange = lims_Δω,
                              colormap = Reverse(:RdBu_5),
                              highclip = :red, lowclip = :blue)

         Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}",
               label = "m/s²")
         Colorbar(fig_total[2, 4], hm_ω_total, tickformat = "{:.1e}",
               label = "1/s")

         Colorbar(fig_perturb[2, 2], hm_b_perturb, tickformat = "{:.1e}",
               label = "m/s²")
         Colorbar(fig_perturb[2, 4], hm_ω_perturb, tickformat = "{:.1e}",
               label = "1/s")

         title_total = @lift @sprintf("Fields at %s = %i %s; t = %.2f days",
                                   const_dim, nearest, units, times[$n])
         title_perturb = @lift @sprintf(
                            "Perturbation fields at %s = %i %s; t = %.2f days",
                            const_dim, nearest, units, times[$n])

         fig_total[1, 1:4]   = Label(fig_total, title_total, fontsize = 24,
                                  tellwidth = false)
         fig_perturb[1, 1:4] = Label(fig_perturb, title_perturb, fontsize = 24,
                                  tellwidth = false)
 
         for i = 1:t_idx_skip:Nt_file
            recordframe!(video_total)
            recordframe!(video_perturb)
            yield()
            n[] = i
         end

         close(ds)
         file_idx += 1
      end

      save(joinpath("./Plots",
            "bzeta_total_$(const_dim)$(nearest)_$(datetime).mp4"),
           video_total)
      save(joinpath("./Plots",
            "bzeta_perturbs_$(const_dim)$(nearest)_$(datetime).mp4"),
           video_perturb)
   end
end

function visualize_z_grid(datetime, grid, zmin; zmax = 0.0)

   mkpath("./Plots") #Make visualization directory if nonexistent

   zc = znodes(grid, Center())
   zf = znodes(grid, Face())
   Δz = zspacings(grid, Center())
   
   fig  = Figure(size=(1200, 600))
   axz  = Axis(fig[1, 1], title = "z-grid")
   axΔz = Axis(fig[2, 1]; xlabel = "z (m)", ylabel = "z-spacing (m)")

   lines!(axz, [zmin, zmax], [0, 0], color = :gray)
   scatter!(axz, zf, 0 * zf, marker = :vline, color = :gray, markersize = 20)
   scatter!(axz, zc, 0 * zc)
   hidedecorations!(axz)
   hidespines!(axz)

   scatter!(axΔz, zc, Δz)
   hidespines!(axΔz, :t, :r)

   rowsize!(fig.layout, 1, Relative(0.1))

   save(joinpath("./Plots", "zgrid_$(datetime).png"), fig)
end

function visualize_fields_2D_slice(datetime, const_dim, const_idx, B, Uφ;
                                   Hx = 0, Hy = 0, Hz = 0, 
                                   plot_animation = true, t_idx_skip = 1,
                                   plot_speed_animation = false)
   #=
   Plot 2D slices of prognostic fields and, optionally, horizontal speed. 
    By default, data are assumed to exclude halos.
   =#

   outfile_list = glob("./Output/output_$(datetime)*")
   
   ds_f, x, y, zC, zF, times, Nt, chron_idcs = open_dataset(
					                                 outfile_list[length(outfile_list)];
                                           Hx = Hx, Hy = Hy, Hz = Hz)

   if const_dim == "x"
      x_idx, y_idx, z_idx       = const_idx, nothing, nothing
      axis1, axis2_zC, axis2_zF = x, zC, zF
   elseif const_dim == "y"
      x_idx, y_idx, z_idx       = nothing, const_idx, nothing
      axis1, axis2_zC, axis2_zF = y, zC, zF
   elseif const_dim == "z"
      x_idx, y_idx, z_idx       = nothing, nothing, const_idx
      axis1, axis2_zC, axis2_zF = x, y, y
   end

   xyzC_idcs, xyzF_idcs = get_2D_spatial_axis_idcs(const_dim;
                                  Hx = Hx, Hy = Hy, Hz = Hz,
                                  x_idx = x_idx, y_idx = y_idx, z_idx = z_idx,
                                  xC = x, yC = y, zC = zC, zF = zF)

   B  = adapt(Array, B)[xyzC_idcs...]
   Uφ = no_offset_view(adapt(Array, Uφ))[xyzC_idcs...]

   b_total_f  = adapt(Array, ds_f[:b])[xyzC_idcs..., chron_idcs[Nt]]
   ur_total_f = adapt(Array, ds_f[:ur])[xyzC_idcs..., chron_idcs[Nt]]
   uφ_total_f = adapt(Array, ds_f[:uφ])[xyzC_idcs..., chron_idcs[Nt]]
   uz_total_f = adapt(Array, ds_f[:uz])[xyzF_idcs..., chron_idcs[Nt]]

   Δb_f  = b_total_f .- B
   Δuφ_f = uφ_total_f .- Uφ

   lims_b_total  = get_range_lims(b_total_f)
   lims_ur       = get_range_lims(ur_total_f)
   lims_uφ_total = get_range_lims(uφ_total_f)
   lims_uz       = get_range_lims(uz_total_f)
   lims_Δb       = get_range_lims(Δb_f)
   lims_Δuφ      = get_range_lims(Δuφ_f)

   mkpath("./Plots") #Make visualization directory if nonexistent

   nearest, ax_kwargs = get_2D_spatial_axis_kwargs(x, y, zC, const_dim;
                                                   x_idx = x_idx, 
                                                   y_idx = y_idx,
                                                   z_idx = z_idx)

   #Plot static images (final frame, by default)

   fig_total = Figure(size = (1200, 800))
   fig_pert  = Figure(size = (1200, 800))

   ax_b_total  = Axis(fig_total[2, 1];
                      title = L"Total buoyancy ($b$)", ax_kwargs...)
   ax_ur_total = Axis(fig_total[2, 3];
                      title = L"Total radial velocity ($u_r$)", ax_kwargs...)
   ax_uφ_total = Axis(fig_total[3, 1];
                      title = L"Total azimuthal velocity ($u_{\phi}$)", 
		      ax_kwargs...)
   ax_uz_total = Axis(fig_total[3, 3];
                      title = L"Total vertical velocity ($u_z$)", ax_kwargs...)

   ax_b_pert  = Axis(fig_pert[2, 1];
                     title = L"Buoyancy perturbation ($b'$)", ax_kwargs...)
   ax_ur_pert = Axis(fig_pert[2, 3];
                     title = L"Radial velocity perturbation ($u_r'$)", 
	             ax_kwargs...)
   ax_uφ_pert = Axis(fig_pert[3, 1];
                     title = L"Azimuthal velocity perturbation ($u_{\phi}'$)", 
	             ax_kwargs...)
   ax_uz_pert = Axis(fig_pert[3, 3];
                     title = L"Vertical velocity perturbation ($u_z'$)", 
		     ax_kwargs...)

   hm_b_total  = heatmap!(ax_b_total, axis1, axis2_zC, b_total_f, colorrange = lims_b_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
   hm_ur_total = heatmap!(ax_ur_total, axis1, axis2_zC, ur_total_f, colorrange = lims_ur, colormap = Reverse(:RdBu_5))
   hm_uφ_total = heatmap!(ax_uφ_total, axis1, axis2_zC, uφ_total_f, colorrange = lims_uφ_total, colormap = Reverse(:RdBu_5))
   hm_uz_total = heatmap!(ax_uz_total, axis1, axis2_zF, uz_total_f, colorrange = lims_uz, colormap = Reverse(:RdBu_5))

   hm_b_pert  = heatmap!(ax_b_pert, axis1, axis2_zC, Δb_f, colorrange = lims_Δb, colormap = Reverse(:RdBu_5))
   hm_ur_pert = heatmap!(ax_ur_pert, axis1, axis2_zC, ur_total_f, colorrange = lims_ur, colormap = Reverse(:RdBu_5))
   hm_uφ_pert = heatmap!(ax_uφ_pert, axis1, axis2_zC, Δuφ_f, colorrange = lims_Δuφ, colormap = Reverse(:RdBu_5))
   hm_uz_pert = heatmap!(ax_uz_pert, axis1, axis2_zF, uz_total_f, colorrange = lims_uz, colormap = Reverse(:RdBu_5))

   Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

   Colorbar(fig_pert[2, 2], hm_b_pert, tickformat = "{:.1e}", label = "m/s²")
   Colorbar(fig_pert[2, 4], hm_ur_pert, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_pert[3, 2], hm_uφ_pert, tickformat = "{:.1e}", label = "m/s")
   Colorbar(fig_pert[3, 4], hm_uz_pert, tickformat = "{:.1e}", label = "m/s")

   title_total = @sprintf("Fields at %s = %.2f km; t = %.2f days",
                          const_dim, nearest, times[Nt])
   title_pert  = @sprintf("Perturbation fields at %s = %.2f km; t = %.2f days",
                          const_dim, nearest, times[Nt])

   fig_total[1, 1:4] = Label(fig_total, title_total, fontsize = 24, tellwidth = false)
   fig_pert[1, 1:4]  = Label(fig_pert, title_pert, fontsize = 24, tellwidth = false)

   save(joinpath("./Plots", "fields_$(const_dim)$(nearest)_tf_$(datetime).png"), fig_total)
   save(joinpath("./Plots", "perturbs_$(const_dim)$(nearest)_tf_$(datetime).png"), fig_pert)
   close(ds_f)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      fig_total = Figure(size = (1200, 800))
      fig_pert  = Figure(size = (1200, 800))

      ax_b_total  = Axis(fig_total[2, 1];
                         title = L"Total buoyancy ($b$)", ax_kwargs...)
      ax_ur_total = Axis(fig_total[2, 3];
                         title = L"Total radial velocity ($u_r$)",
                         ax_kwargs...)
      ax_uφ_total = Axis(fig_total[3, 1];
                         title = L"Total azimuthal velocity ($u_{\phi}$)",
                         ax_kwargs...)
      ax_uz_total = Axis(fig_total[3, 3];
                         title = L"Total vertical velocity ($u_z$)",
                         ax_kwargs...)

      ax_b_pert  = Axis(fig_pert[2, 1];
                        title = L"Buoyancy perturbation ($b'$)", ax_kwargs...)
      ax_ur_pert = Axis(fig_pert[2, 3];
                        title = L"Radial velocity perturbation ($u_r'$)",
                        ax_kwargs...)
      ax_uφ_pert = Axis(fig_pert[3, 1];
                        title = L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                        ax_kwargs...)
      ax_uz_pert = Axis(fig_pert[3, 3];
                        title = L"Vertical velocity perturbation ($u_z'$)",
                        ax_kwargs...)
                        
      if plot_speed_animation #Plot animation of (horizontal) speed
         fig_speed = Figure(size = (800, 500))
         ax_speed  = Axis(fig_speed[2, 1]; ax_kwargs...)
      end

      ds, x, y, zC, zF, times, Nt, chron_idcs = open_dataset(outfile_list,
                                                             Hx = Hx, Hy = Hy, 
                                                             Hz = Hz)

      n = Observable(1)

      b_total  = @lift ds[:b][xyzC_idcs..., chron_idcs[$n]]
      ur_total = @lift ds[:ur][xyzC_idcs..., chron_idcs[$n]]
      uφ_total = @lift ds[:uφ][xyzC_idcs..., chron_idcs[$n]]
      uz_total = @lift ds[:uz][xyzF_idcs..., chron_idcs[$n]]

      Δb  = @lift $b_total .- B
      Δuφ = @lift $uφ_total .- Uφ

      hm_b_total  = heatmap!(ax_b_total, axis1, axis2_zC, b_total, colorrange = lims_b_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_ur_total = heatmap!(ax_ur_total, axis1, axis2_zC, ur_total, colorrange = lims_ur, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uφ_total = heatmap!(ax_uφ_total, axis1, axis2_zC, uφ_total, colorrange = lims_uφ_total, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)
      hm_uz_total = heatmap!(ax_uz_total, axis1, axis2_zF, uz_total, colorrange = lims_uz, colormap = Reverse(:RdBu_5))#, highclip = :red, lowclip = :blue)

      hm_b_pert  = heatmap!(ax_b_pert, axis1, axis2_zC, Δb, colormap = Reverse(:RdBu_5))#, colorrange = lims_Δb, highclip = :red, lowclip = :blue)
      hm_ur_pert = heatmap!(ax_ur_pert, axis1, axis2_zC, ur_total, colormap = Reverse(:RdBu_5))#, colorrange = lims_ur, highclip = :red, lowclip = :blue)
      hm_uφ_pert = heatmap!(ax_uφ_pert, axis1, axis2_zC, Δuφ, colormap = Reverse(:RdBu_5))#, colorrange = lims_Δuφ, highclip = :red, lowclip = :blue)
      hm_uz_pert = heatmap!(ax_uz_pert, axis1, axis2_zF, uz_total, colormap = Reverse(:RdBu_5))#, colorrange = lims_uz, highclip = :red, lowclip = :blue)

      Colorbar(fig_total[2, 2], hm_b_total, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_total[2, 4], hm_ur_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 2], hm_uφ_total, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_total[3, 4], hm_uz_total, tickformat = "{:.1e}", label = "m/s")

      Colorbar(fig_pert[2, 2], hm_b_pert, tickformat = "{:.1e}", label = "m/s²")
      Colorbar(fig_pert[2, 4], hm_ur_pert, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_pert[3, 2], hm_uφ_pert, tickformat = "{:.1e}", label = "m/s")
      Colorbar(fig_pert[3, 4], hm_uz_pert, tickformat = "{:.1e}", label = "m/s")

      title_total = @lift @sprintf("Fields at %s = %.2f km; t = %.2f days",
                                   const_dim, nearest, times[$n])
      title_pert  = @lift @sprintf(
                            "Perturbation fields at %s = %.2f km; t = %.2f days",
                                   const_dim, nearest, times[$n])

      fig_total[1, 1:4] = Label(fig_total, title_total, fontsize = 24, tellwidth = false)
      fig_pert[1, 1:4]  = Label(fig_pert, title_pert, fontsize = 24, tellwidth = false)

      if plot_speed_animation

         speed_total = @lift sqrt.($ur_total.^2 .+ $uφ_total.^2)
         
         hm_speed_total = heatmap!(ax_speed, axis1, axis2_zC, speed_total, colormap = :Reds)
         
         Colorbar(fig_speed[2, 2], hm_speed_total, tickformat = "{:.1e}", label = "m/s")
         
         title_speed = @lift @sprintf("Total horizontal speed at %s = %.2f km; t = %.2f days",
                                      const_dim, nearest, times[$n])
                                        
         fig_speed[1, 1:2] = Label(fig_speed, title_speed, fontsize = 24, tellwidth = false)
         video_speed = VideoStream(fig_speed, format = "mp4", framerate = 6)
      end

      frames = 1:Nt

      video_total = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_pert  = VideoStream(fig_pert, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
         
         print(i, " of ", Nt, "\n")
         
         recordframe!(video_total)
         recordframe!(video_pert)
         
         if plot_speed_animation
            recordframe!(video_speed)
         end
         
         yield()
         n[] = i
      end

      save(joinpath("./Plots", "fields_$(const_dim)$(nearest)_$(datetime).mp4"), video_total)
      save(joinpath("./Plots", "perturbs_$(const_dim)$(nearest)_$(datetime).mp4"), video_pert)
      
      if plot_speed_animation
         save(joinpath("./Plots", "horizontalSpeed_$(const_dim)$(nearest)_$(datetime).mp4"), video_speed)
      end
   end
   close(ds)
end

function open_computed_dataset(datetime, Δx, Δy, Δz, f)
   #=
   Produced nc file containing computed diagnostics, if it does not already
    exist, and return the opened file in read mode.
   =#

   computed_file = joinpath("./Output", "computed_$(datetime).nc")

   if !isfile(computed_file) #Only compute if file does not already exist

      ds, x, y, z, times, Nt = open_dataset(datetime)

      frames = 1:Nt
      x_idcs = 2:length(x)-2
      y_idcs = 2:length(y)-2
      z_idcs = 2:length(z)-2

      function update_data_array!(data_array, i, j, k, n, value)
         data_array[i, j, k, n] = value
         return data_array
      end
      
      NCDataset(computed_file, "c") do comp_ds
         
         i, j, k = Observable(2), Observable(2), Observable(2)
         n       = Observable(1)

         b  = @lift ds["b"][:, :, :, $n]
         u  = @lift ds["u"][1:end-1, :, :, $n]
         v  = @lift ds["v"][:, 1:end-1, :, $n]
         w  = @lift ds["w"][:, :, 1:end-1, $n]
         qn = @lift q($u, $v, $w, $b, f, $i, $j, $k, Δx, Δy, Δz)

         defDim(comp_ds, "x", length(x)-2)
         defDim(comp_ds, "y", length(y)-2)
         defDim(comp_ds, "z", length(z)-2)
         defDim(comp_ds, "time", length(times))

         q_data  = Array{Float64, 4}(undef, length(x)-2, length(y)-2, 
                                     length(z)-2, Nt)

         for t = 1:frames[end]
            for z_idx = 2:z_idcs[end]
               for x_idx = 2:x_idcs[end]
               
                  for y_idx = 2:y_idcs[end]
                     update_data_array!(q_data, x_idx, y_idx, z_idx, t, 
                                        to_value(qn))
                     yield()
                     j[] = y_idx
                  end
                  
                  i[] = x_idx
               end
               
               k[] = z_idx
            end
            
            print("Computing q for time $(t) of $(Nt)" * " \r")
            n[] = t
         end
         
         defVar(comp_ds, "q", q_data, ("x", "y", "z", "time"))
      end #comp_ds gets closed automatically
   end

   return NCDataset(computed_file, "r")
end

#=
function visualize_OkuboWeiss(datetime, )

end
=#

function visualize_q_2D_slice(datetime, const_dim, const_idx, f;
                              Hx = 0, Hy = 0, Hz = 0, 
                              plot_animation = true)
                              
   #=
   Plot 2D slices of potential vorticity. By default, data are assumed to exclude
    halos.
   =#

   outfile_list = glob("./Output/output_$(datetime)*")
   
   ds_f, x, y, zC, zF, times, Nt, chron_idcs = open_dataset(
					                                 outfile_list[length(outfile_list)];
                                           Hx = Hx, Hy = Hy, Hz = Hz)
   
   comp_ds = open_computed_dataset(datetime, Δx, Δy, Δz, f)
   
   if const_dim == "x"
      x_idx, y_idx, z_idx       = const_idx, nothing, nothing
      axis1, axis2_zC, axis2_zF = x, zC, zF
   elseif const_dim == "y"
      x_idx, y_idx, z_idx       = nothing, const_idx, nothing
      axis1, axis2_zC, axis2_zF = y, zC, zF
   elseif const_dim == "z"
      x_idx, y_idx, z_idx       = nothing, nothing, const_idx
      axis1, axis2_zC, axis2_zF = x, y, y
   end

   xyzC_idcs, xyzF_idcs = get_2D_spatial_axis_idcs(const_dim;
                                  Hx = Hx, Hy = Hy, Hz = Hz,
                                  x_idx = x_idx, y_idx = y_idx, z_idx = z_idx,
                                  xC = x, yC = y, zC = zC, zF = zF)
                                  
   q_f = adapt(Array, comp_ds[:q])[xyzC_idcs..., chron_idcs[Nt]]
   
   lims_q = get_range_lims(q_f)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   
   nearest, ax_kwargs = get_2D_spatial_axis_kwargs(x, y, zC, const_dim;
                                                   x_idx = x_idx, 
                                                   y_idx = y_idx,
                                                   z_idx = z_idx)

   #Plot static images (final frame, by default)
   
   
   fig_q = Figure(size = (1200, 800))
   ax_q  = Axis(fig_q[2, 1]; title = L"Potential vorticity ($q$)", ax_kwargs...)

   #=z_plt = div(length(z[:]), 2) #z-index to start plot at

   n    = Observable(1)
   b    = @lift ds[:b][:, :, z_plt:end, $n]
   u    = @lift ds[:u][:, :, z_plt:end, $n]
   v    = @lift ds[:v][:, :, z_plt:end, $n]
   w    = @lift ds[:w][:, :, z_plt:end-1, $n]
   q_yz = @lift comp_ds[:q][x_idx, :, z_plt:end, $n] 
   
   q_yz_f = comp_ds[:q][x_idx, :, z_plt:end, Nt]

   lims_q = get_range_lims(q_yz_f)

   x_nearest, axis_kwargs_yz = get_2D_spatial_axis_kwargs(x, y, z; x_idx = x_idx)

   fig_q = Figure(size = (600, 600))
   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_yz...)
   hm_q  = heatmap!(ax_q, y[2:end-1], z[z_plt:end], q_yz, colorrange = lims_q,
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title_q       = @lift @sprintf("q at x = %i km; t = %.2f days",
			           x_nearest, times[$n])
   fig_q[1, 1:2] = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   frames  = 1:Nt
   video_q = VideoStream(fig_q, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   save(joinpath("./Plots", "q_x$(x_nearest)_$(datetime).mp4"), video_q)
   =#
   close(ds)
end

function visualize_q_const_y(datetime, Δx, Δy, Δz, f, y_idx)
   
   ds, x, y, z, times, Nt = open_dataset(datetime)
   
   comp_ds = open_computed_dataset(datetime, Δx, Δy, Δz, f)

   z_plt = div(length(z[:]), 2) #z-index to start plot at

   n    = Observable(1)
   b    = @lift ds["b"][:, :, z_plt:end, $n]
   u    = @lift ds["u"][:, :, z_plt:end, $n]
   v    = @lift ds["v"][:, :, z_plt:end, $n]
   w    = @lift ds["w"][:, :, z_plt:end-1, $n]
   q_xz = @lift comp_ds["q"][:, y_idx, z_plt:end, $n]

   q_xz_f = comp_ds["q"][:, y_idx, z_plt:end, Nt]

   lims_q = get_range_lims(q_xz_f)

   y_nearest, axis_kwargs_xz = get_2D_spatial_axis_kwargs(x, y, z; y_idx = y_idx)

   fig_q = Figure(size = (600, 600))
   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_xz...)
   hm_q  = heatmap!(ax_q, x[2:end-1], z[z_plt:end], q_xz, colorrange = lims_q, 
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title_q        = @lift @sprintf("q at y = %i km; t = %.2f days",
			            y_nearest, times[$n])
   fig_q[1, 1:2]  = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   frames  = 1:Nt
   video_q = VideoStream(fig_q, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_y$(y_nearest)_$(datetime).mp4"), video_q)
   close(ds)
end

function visualize_q_const_z(datetime, Δx, Δy, Δz, f, z_idx)

   ds, x, y, z, times, Nt = open_dataset(datetime)
   
   comp_ds = open_computed_dataset(datetime, Δx, Δy, Δz, f)

   n    = Observable(1)
   b    = @lift ds["b"][:, :, :, $n]
   u    = @lift ds["u"][:, :, :, $n]
   v    = @lift ds["v"][:, :, :, $n]
   w    = @lift ds["w"][:, :, 1:end-1, $n]
   q_xy = @lift comp_ds["q"][:, :, z_idx, $n]

   q_xy_f = comp_ds["q"][:, :, z_idx, Nt]

   lims_q = get_range_lims(q_xy_f)

   depth_nearest, axis_kwargs_xy = get_2D_spatial_axis_kwargs(x, y, z; z_idx = z_idx)

   fig_q = Figure(size = (600, 600))
   ax_q  = Axis(fig_q[2, 1]; axis_kwargs_xy...)
   hm_q  = heatmap!(ax_q, x[2:end-1], y[2:end-1], q_xy, colorrange = lims_q,
		    colormap = :balance)

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")

   title_q       = @lift @sprintf("q at %i-m depth; t = %.2f days",
			     depth_nearest, times[$n])
   fig_q[1, 1:2] = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   frames   = 1:Nt
   video_q  = VideoStream(fig_q, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video_q)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_z-$(depth_nearest)_$(datetime).mp4"), video_q)
   close(ds)
end