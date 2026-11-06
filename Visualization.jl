include("FourierDiagnostics.jl")
include("LibraryDynamics.jl")
include("LibraryVisualization.jl")
include("ReadSimulationOutputs.jl")

using Adapt, CairoMakie, CommonDataModel, CUDA, DataStructures
using Glob, LaTeXStrings
using Oceananigans
using Oceananigans.AbstractOperations
using Oceananigans.Fields
using Oceananigans.OutputReaders
using OffsetArrays: no_offset_view
using Polynomials: fit
using Printf

update_theme!(theme_latexfonts(), fontsize = 16)

################################################################################

function Fouriertesting(datetimenow, model, gridParams, backgroundFields)

   outfile_paths = glob("./Output/output_$(datetimenow)*")
   
   b_TS = LoadFieldTimeseries("b", datetimenow) #FieldTimeSeries("./Output/output_$(datetimenow)", "b"; backend = OnDisk())
   ux_TS = FieldTimeSeries("./Output/output_$(datetimenow)", "ux"; backend = OnDisk())
   uy_TS = FieldTimeSeries("./Output/output_$(datetimenow)", "uy"; backend = OnDisk())
   uz_TS = FieldTimeSeries("./Output/output_$(datetimenow)", "uz"; backend = OnDisk())
   
   b_f = b_TS[end]
   ux_f = ux_TS[end]
   uy_f = uy_TS[end]
   uz_f = uz_TS[end]

   B = adapt(Array, interior(backgroundFields[:B]))
   
   print(b_TS)
   
   n = Observable(1)
      
   b′ = @lift b_TS[$n] .- B
   
   Nt = length(b_TS)
   frames = 1:Nt

   t_idx_skip = 1
   for i = 1:t_idx_skip:frames[end]
         
      print(i, " of ", Nt, "\n")
         
      yield()
      n[] = i
   end
   
   
   plot_rMean_φDFT_field(uy_TS[end], model.grid, gridParams, 10, 8, "uy")
   #=
   plot_rzMean_φDFT_field(b_TS[end], model.grid, gridParams, 500, 20, "buoyancy")
   plot_rMean_φDFT_field(b_TS[end], model.grid, gridParams, 500, 20, "buoyancy")
   plot_zMean_φDFT_field(b_TS[end], model.grid, gridParams, 500, 20, "buoyancy")

   plot_rzMean_φDFT_field(ux_TS[end], model.grid, gridParams, 500, 20, "ux")
   plot_rMean_φDFT_field(ux_TS[end], model.grid, gridParams, 500, 20, "ux")
   plot_zMean_φDFT_field(ux_TS[end], model.grid, gridParams, 500, 20, "ux")

   plot_rzMean_φDFT_field(uy_TS[end], model.grid, gridParams, 500, 20, "uy")
   plot_rMean_φDFT_field(uy_TS[end], model.grid, gridParams, 500, 20, "uy")
   plot_zMean_φDFT_field(uy_TS[end], model.grid, gridParams, 500, 20, "uy")
   
   plot_rzMean_φDFT_field(uz_TS[end], model.grid, gridParams, 500, 20, "uz")
   plot_rMean_φDFT_field(uz_TS[end], model.grid, gridParams, 500, 20, "uz")
   plot_zMean_φDFT_field(uz_TS[end], model.grid, gridParams, 500, 20, "uz")
   =#
end

################################################################################

function visualize_Ψ_and_z_derivs(simGrid, tall_Ψ_ccf_Field, 
                                  tall_∂Ψ∂z_ccc_Field, tall_∂2Ψ∂z2_ccf_Field)

   fig    = Figure(size = (1400, 700))
   ax_dz0 = Axis(fig[1, 1], xlabel = L"$\Psi$ [m^2/s]", ylabel = L"$z$ [m]")
   ax_dz1 = Axis(fig[1, 2], xlabel = L"$d\Psi /dz$ [m/s]", ylabel = L"$z$ [m]")
   ax_dz2 = Axis(fig[1, 3], xlabel = L"$d^2 \Psi /dz^2$ [1/s]", 
                 ylabel = L"$z$ [m]")
  
   lines!(ax_dz0, no_offset_view(adapt(Array, tall_Ψ_ccf_Field[28, 28, :])), 
          no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), label = "Ψ_ccf_Field")
   scatter!(ax_dz0, no_offset_view(adapt(Array, tall_Ψ_ccf_Field[28, 28, :])),
            no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), 
            label = "Ψ_ccf_Field (at gridpoints)")
     
   lines!(ax_dz1, no_offset_view(adapt(Array, tall_∂Ψ∂z_ccc_Field[28, 28, :])),
          no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶜ)), 
          label = "∂Ψ∂z_ccc_Field")
   scatter!(ax_dz1, 
            no_offset_view(adapt(Array, tall_∂Ψ∂z_ccc_Field[28, 28, :])), 
            no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶜ)), 
            label = "∂Ψ∂z_ccc_Field (at gridpoints)")
  
   lines!(ax_dz2, 
          no_offset_view(adapt(Array, tall_∂2Ψ∂z2_ccf_Field[28, 28, :])), 
          no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), 
          label = "∂2Ψ∂z2_ccf_Field")
   scatter!(ax_dz2, 
            no_offset_view(adapt(Array, tall_∂2Ψ∂z2_ccf_Field[28, 28, :])), 
            no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)), 
            label = "∂2Ψ∂z2_ccf_Field (at gridpoints)")
   lines!(ax_dz2, 
          no_offset_view(adapt(Array, 
                               Field(∂z(tall_∂Ψ∂z_ccc_Field))[28, 28, :])
                        ), no_offset_view(adapt(Array, simGrid.z.cᵃᵃᶠ)),
          label = "z-deriv of ∂Ψ∂z_ccc_Field")
     
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "initialConditionPsi.png"), fig)
end

################################################################################

function visualize_B_and_N²_vs_z(B_total, b_TWB, b_linear, grid, x_idx, y_idx,
                                 gyreParams, doubleTanhParams, Hz)

   B_total_slice = no_offset_view(adapt(Array, B_total))[x_idx, y_idx, :]
   b_TWB_slice  = no_offset_view(adapt(Array, b_TWB)
                                )[x_idx, y_idx,
                                  (Hz + 1):(length(grid.z.cᵃᵃᶜ) - Hz)]
   
   ∂B∂z_total_Field = ZFaceField(grid)
   ∂b∂z_TWB_Field   = ZFaceField(grid)
   
   set!(∂B∂z_total_Field, ∂z(B_total))
   set!(∂b∂z_TWB_Field, ∂z(b_TWB))
   
   ∂B∂z_total_slice = no_offset_view(adapt(Array, ∂B∂z_total_Field)
                                    )[x_idx, y_idx, :]
   ∂b∂z_TWB_slice   = no_offset_view(adapt(Array, ∂b∂z_TWB_Field)
                                    )[x_idx, y_idx, 
                                      (Hz + 1):(length(grid.z.cᵃᵃᶠ) - Hz)]
   
   x      = no_offset_view(adapt(Array, grid.xᶜᵃᵃ))[x_idx]
   y      = no_offset_view(adapt(Array, grid.yᵃᶜᵃ))[y_idx]
   zC_int = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ)
                          )[(Hz + 1):(length(grid.z.cᵃᵃᶜ) - Hz)]
   zC_all = no_offset_view(adapt(Array, grid.z.cᵃᵃᶜ))
   zF_int = no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ)
                          )[(Hz + 1):(length(grid.z.cᵃᵃᶠ) - Hz)]
   zF_all = no_offset_view(adapt(Array, grid.z.cᵃᵃᶠ))
   
   fig   = Figure(size = (1400, 700))
   ax_B  = Axis(fig[2, 1], xlabel = L"$B$ [m/s$^{2}$]", ylabel = L"$z$ [m]")
   ax_N2 = Axis(fig[2, 2], xlabel = L"$N^2$ [s$^{-2}$]", ylabel = L"$z$ [m]")
   
   lines!(ax_B, B_total_slice, zC_all, label = "Total")
   scatter!(ax_B, B_total_slice, zC_all, label = "Total (at gridpoints)")
   #lines!(ax_B, buoyancyDoubleTanh(zC_all, doubleTanhParams), zC_all, 
   #       label = "From double-tanh function")
   lines!(ax_B, b_TWB_slice, zC_int, label = "Thermal-wind contribution")
   lines!(ax_B, gyreParams.N²_far .* zC_int, zC_int, label = "Linear term")
   
   lines!(ax_N2, ∂B∂z_total_slice, zF_all, label = "Total")
   scatter!(ax_N2, ∂B∂z_total_slice, zF_all, label = "Total (at gridpoints)")
   #lines!(ax_N2, N²DoubleTanh(zF_all, doubleTanhParams), zF_all, 
   #       label = "From double-tanh function")
   lines!(ax_N2, ∂b∂z_TWB_slice, zF_int, label = "Thermal-wind contribution")
   lines!(ax_N2, gyreParams.N²_far .+ 0 * zF_int, zF_int, label = "Linear term")

   fig[2, 3] = Legend(fig, ax_B)
   fig[1, 1] = Label(fig, "Background buoyancy", fontsize = 24, 
                     tellwidth = false)
   fig[1, 2] = Label(fig, L"Background $N^2$", fontsize = 24, tellwidth = false)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "b_and_N2_profiles_i$(x_idx)_j$(y_idx)_Nz$(grid.Nz).png"), fig)
end

################################################################################

function visualize_norms(datetime; 
		idxStartLinGrowth_b = 2, idxEndLinGrowth_b = -1,
    idxStartLinGrowth_ur = nothing, idxEndLinGrowth_ur = nothing,
		idxStartLinGrowth_uφ = nothing, idxEndLinGrowth_uφ = nothing,
		idxStartLinGrowth_ux = nothing, idxEndLinGrowth_ux = nothing,
		idxStartLinGrowth_uy = nothing, idxEndLinGrowth_uy = nothing,
		idxStartLinGrowth_uz = nothing, idxEndLinGrowth_uz = nothing,
		idxStartPlot = 2, idxEndPlot = -1,
    growth_rate = "linear_best_fit", t_skip_idx = 1)

   scalars_ds, times, Nt, chron_idcs = open_scalars_dataset(glob("./Output/scalars_$(datetime)*"))
   
   if idxEndPlot < 0
      idxEndPlot = Nt + idxEndPlot #Make idxEndPlot a positive integer
   end
   
   #Define range of t-indices to be plotted (this is only for readability)
   tPlotRange = idxStartPlot:t_skip_idx:idxEndPlot
   
   #Retain only necessary portion of 'times' data and update Nt accordingly
   times = times[tPlotRange]
   Nt    = length(times)

   #Load data, sorted chronologically and then restricted to plot interval
   b′_norm  = scalars_ds[:b′_norm][chron_idcs][tPlotRange]
   ux′_norm = scalars_ds[:ux′_norm][chron_idcs][tPlotRange]
   uy′_norm = scalars_ds[:uy′_norm][chron_idcs][tPlotRange]
   ur′_norm = scalars_ds[:ur′_norm][chron_idcs][tPlotRange]
   uφ′_norm = scalars_ds[:uφ′_norm][chron_idcs][tPlotRange]
   uz′_norm = scalars_ds[:uz′_norm][chron_idcs][tPlotRange]
   
   fig_cyl, ax_b_cyl, ax_ur, ax_uφ, ax_uz_cyl = set_up_fig_perturb_norms(times, b′_norm, ur′_norm, uφ′_norm, uz′_norm; growth_rate = growth_rate)
   fig_Cart, ax_b_Cart, ax_ux, ax_uy, ax_uz_Cart = set_up_fig_perturb_norms(times, b′_norm, ux′_norm, uy′_norm, uz′_norm; Cartesian = true, growth_rate = growth_rate)
   
   if growth_rate == "linear_best_fit" #Plot linear fit on restricted t-interval

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

      idxStartLinGrowth_ur, idxEndLinGrowth_ur = updateGrowthIdcs!(idxStartLinGrowth_ur, idxEndLinGrowth_ur)
      idxStartLinGrowth_uφ, idxEndLinGrowth_uφ = updateGrowthIdcs!(idxStartLinGrowth_uφ, idxEndLinGrowth_uφ)
      idxStartLinGrowth_ux, idxEndLinGrowth_ux = updateGrowthIdcs!(idxStartLinGrowth_ux, idxEndLinGrowth_ux)
      idxStartLinGrowth_uy, idxEndLinGrowth_uy = updateGrowthIdcs!(idxStartLinGrowth_uy, idxEndLinGrowth_uy)
      idxStartLinGrowth_uz, idxEndLinGrowth_uz = updateGrowthIdcs!(idxStartLinGrowth_uz, idxEndLinGrowth_uz)

      #=
      @printf("Empirical growth rate:\n From b′-norm: %.5f per day\n From ur′-norm: %.5f per day\n From uφ′-norm: %.5f per day\n From ux′-norm: %.5f per day\n From uy′-norm: %.5f per day\n From uz′-norm: %.5f per day\n",
	       b′NormLinearFitParams[1], ur′NormLinearFitParams[1],
	       uφ′NormLinearFitParams[1], ux′NormLinearFitParams[1],
	       uy′NormLinearFitParams[1], uz′NormLinearFitParams[1])
      =#
      
      plot_growth_fit_line!(ax_b_cyl, times, growthIdcs_b, b′_norm)
      plot_growth_fit_line!(ax_ur, times, (idxStartLinGrowth_ur, 
                                          idxEndLinGrowth_ur), 
                            ur′_norm)
      plot_growth_fit_line!(ax_uφ, times, (idxStartLinGrowth_uφ, 
                                           idxEndLinGrowth_uφ), 
                            uφ′_norm)
      plot_growth_fit_line!(ax_uz_cyl, times, (idxStartLinGrowth_uz, 
                                               idxEndLinGrowth_uz), 
                            uz′_norm)
      
      plot_growth_fit_line!(ax_b_Cart, times, growthIdcs_b, b′_norm)
      plot_growth_fit_line!(ax_ux, times, (idxStartLinGrowth_ux, 
                                           idxEndLinGrowth_ux),
                            ux′_norm)
      plot_growth_fit_line!(ax_uy, times, (idxStartLinGrowth_uy, 
                                           idxEndLinGrowth_uy), 
                            uy′_norm)
      plot_growth_fit_line!(ax_uz_Cart, times, (idxStartLinGrowth_uz, 
                                                idxEndLinGrowth_uz), 
                            uz′_norm)
      
   elseif growth_rate == "timeseries" #Plot time-derivative at all t
   
      superimpose_growth_plots!(fig_cyl, times, b′_norm, ur′_norm, uφ′_norm, 
                                uz′_norm)
      superimpose_growth_plots!(fig_Cart, times, b′_norm, ux′_norm, uy′_norm, 
                                uz′_norm)
   end

   mkpath("./Plots") #Make visualization directory if nonexistent
   
   save(joinpath("./Plots", "norm_cyl_fields_$(datetime).png"), fig_cyl)
   save(joinpath("./Plots", "norm_Cart_fields_$(datetime).png"), fig_Cart)
   
   close(scalars_ds)
end

################################################################################

function visualize_ME(datetime, grid, simParams::SimulationParams)

   outfile_list          = glob("./Output/energetics_$(datetime)*")
   ds, t, Nt, chron_idcs = open_energetics_dataset(outfile_list)

   total_KE = ds[:total_KE][chron_idcs]
   total_PE = ds[:total_PE][chron_idcs]
   total_ME = total_KE + total_PE
   
   #########################################
   # Quantify simulation nonconservativity #
   #########################################
   
   percentage_ME_loss = quantify_simulation_nonconservativity(total_KE, total_PE)
   
   @printf("Percentage of mechanical energy lost during simulation duration: %.5f", percentage_ME_loss)

   fig_ME = Figure(size = (1200, 700))
   ax_ME  = Axis(fig_ME[2, 1]; xlabel = "Time [days]",
                 ylabel = L"[m$^5$ s$^{-2}$]")
                 
   lines!(ax_ME, t, 0.5 * total_ME, label = "50% of total ME", color = :black,
          linewidth = 2)
   lines!(ax_ME, t, total_KE, label = "Total KE", color = :green, linewidth = 2)
   lines!(ax_ME, t, total_PE, label = "Total PE", color = :purple, 
          linewidth = 2)
   
   fig_ME[1, 1] = Label(fig_ME, "Total-ME budget", fontsize = 24, 
                        tellwidth = false)

   axislegend(ax_ME)
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "MEbudget_$(datetime).png"), fig_ME)
   
   if simParams.save_total_ME_budgets_terms
   
      total_KE_adv_flux   = ds[:total_KE_adv_flux][chron_idcs]
      total_KE_production = ds[:total_KE_production][chron_idcs]
      total_pressure_work = ds[:total_pressure_work][chron_idcs]
      total_gravity_work  = ds[:total_gravity_work][chron_idcs]
      
      #Convert time to s to compute time-derivative of KE
      dt_KE_total = order1_forward_difference(86400 * t, total_KE)
   
      fig_KE_budget = Figure(size = (1200, 700))
      ax_KE_budget  = Axis(fig_KE_budget[2, 1]; xlabel = "Time [days]",
                           ylabel = L"[m$^5$ s$^{-3}$]")
                     
      lines!(ax_KE_budget, t[1:end-1], dt_KE_total, 
             label = "Time derivative of total KE", color = :yellowgreen, 
             linewidth = 2)
      lines!(ax_KE_budget, t, total_KE_adv_flux, label = "Advective flux",
             color = :royalblue, linewidth = 2)
      lines!(ax_KE_budget, t, total_pressure_work, label = "Pressure work", 
             color = :sienna, linewidth = 2)
      lines!(ax_KE_budget, t, total_KE_production, label = "Buoyant production",
             color = :orange, linewidth = 2)
      lines!(ax_KE_budget, t[1:end-1], 
             (dt_KE_total - total_KE_adv_flux[1:end-1] 
              - total_pressure_work[1:end-1] - total_KE_production[1:end-1]),
             label = "Residual", color = :red)

      fig_KE_budget[2, 2] = Legend(fig_KE_budget, ax_KE_budget)
      fig_KE_budget[1, 1] = Label(fig_KE_budget, "Terms in total-KE budget",
                                  fontsize = 24, tellwidth = false)

      save(joinpath("./Plots", "KEbudget_$(datetime).png"), fig_KE_budget)
      
      #Convert time to s to compute time-derivative of PE
      dt_PE_total = order1_forward_difference(86400 * t, total_PE)

      fig_PE_budget = Figure(size = (1200, 700))
      ax_PE_budget  = Axis(fig_PE_budget[2, 1]; xlabel = "Time [days]",
                           ylabel = L"[m$^5$ s$^{-3}$]")
                     
      scatter!(ax_PE_budget, t[1:end-1], dt_PE_total, 
               label = "Time derivative of total PE", color = :yellowgreen)
      scatter!(ax_PE_budget, t, total_gravity_work, label = "Gravity work", 
               color = :royalblue)
      scatter!(ax_PE_budget, t, -total_KE_production, 
               label = "Buoyant production (of KE)", color = :orange)
      scatter!(ax_PE_budget, t[1:end-1], 
               (dt_PE_total - total_gravity_work[1:end-1] 
                + total_KE_production[1:end-1]), 
               label = "Residual", color = :red)

      fig_PE_budget[1, 1] = Label(fig_PE_budget, "Terms in total-PE budget",
                                  fontsize = 24, tellwidth = false)

      axislegend(ax_PE_budget, position = :lb)
      save(joinpath("./Plots", "PEbudget_$(datetime).png"), fig_PE_budget)
   end

   close(ds)
end

################################################################################

function visualize_PKE(datetime, grid, simParams::SimulationParams)

   outfile_list          = glob("./Output/energetics_$(datetime)*")
   ds, t, Nt, chron_idcs = open_energetics_dataset(outfile_list)

   if simParams.save_PKE_full_domain_budget_terms
   
      PKE           = ds[:total_PKE][chron_idcs][2:end]
      PAPE_to_PKE   = ds[:total_PAPE_to_PKE][chron_idcs][2:end]
      BTI_transfer  = ds[:total_BTI_transfer][chron_idcs][2:end]
      #BCI_transfer = ds[:total_BCI_transfer][chron_idcs][2:end]
      bkgd_transfer = ds[:total_bkgd_transfer][chron_idcs][2:end]
      
      total_bkgd_xxx = ds[:total_bkgd_xxx][chron_idcs][2:end]
      total_bkgd_xyx = ds[:total_bkgd_xyx][chron_idcs][2:end]
      total_bkgd_xzx = ds[:total_bkgd_xzx][chron_idcs][2:end]
      total_bkgd_yxy = ds[:total_bkgd_yxy][chron_idcs][2:end]
      total_bkgd_yyy = ds[:total_bkgd_yyy][chron_idcs][2:end]
      total_bkgd_yzy = ds[:total_bkgd_yzy][chron_idcs][2:end]
      total_bkgd_yzy = ds[:total_bkgd_yzy][chron_idcs][2:end]
      total_bkgd_zxz = ds[:total_bkgd_zxz][chron_idcs][2:end]
      total_bkgd_zyz = ds[:total_bkgd_zyz][chron_idcs][2:end]
      total_bkgd_zzz = ds[:total_bkgd_zzz][chron_idcs][2:end]
   
      ##################
      # Plot total PKE #
      ##################
   
      fig_PKEtotal = Figure(size = (1200, 700))
      ax_PKEtotal  = Axis(fig_PKEtotal[2, 1]; xlabel = "Time [days]",
                          ylabel = "Domain-integrated specific PKE", 
                          yscale = log10)

      scatter!(ax_PKEtotal, t[2:end], PKE, color = :black)

      fig_PKEtotal[1, 1] = Label(fig_PKEtotal,
	                               "Volume-integrated perturbation kinetic energy",
                                 fontsize = 24, tellwidth = false)
 
      mkpath("./Plots") #Make visualization directory if nonexistent
      save(joinpath("./Plots", "PKEtotal_$(datetime).png"), fig_PKEtotal)
      
      ################################
      # Plot whole-domain PKE budget #
      ################################
      
      #Convert time to s to compute time-derivative of PKE
      dt_PKE_total = centered_difference(86400 * t[2:end], PKE)
   
      fig_budget = Figure(size = (1200, 700))
      ax_budget  = Axis(fig_budget[2, 1]; xlabel = "Time [days]",
                        ylabel = L"[s$^{-1}$]")

      lines!(ax_budget, t[3:end-1], dt_PKE_total, 
             label = "Time derivative of total PKE", color = :yellowgreen,
             linewidth = 2)
      lines!(ax_budget, t[2:end], PAPE_to_PKE, label = "PAPE to PKE",
             color = :navy, linewidth = 2)
      scatter!(ax_budget, t[2:end], BTI_transfer, label = "BTI transfer", 
               color = :blue) #:hotpink)
      #scatter!(ax_budget, t[2:end], BCI_transfer, label = "BCI transfer", 
      #         color = :darkgreen)
      #scatter!(ax_budget, t[3:end-1], 
      #        (dt_PKE_total - PAPE_to_PKE[2:end-1] - BTI_transfer[2:end-1]
      #         - BCI_transfer[2:end-1]), 
      #       label = "Residual", color = :red)
      lines!(ax_budget, t[2:end], bkgd_transfer,
             label = "Total background transfer", color = :pink, linewidth = 2)
      lines!(ax_budget, t[2:end], total_bkgd_xyx[1:end] + total_bkgd_xxx[1:end] + total_bkgd_yxy[1:end] + total_bkgd_yyy[1:end], label = "From horizontal shear", color = :pink, linestyle = :dash, linewidth = 2)
      lines!(ax_budget, t[2:end], total_bkgd_xzx[1:end] + total_bkgd_yzy[1:end], label = "From vertical shear", color = :pink, linestyle = :dot, linewidth = 2)
      lines!(ax_budget, t[3:end-1], (dt_PKE_total - PAPE_to_PKE[2:end-1] - bkgd_transfer[2:end-1]), label = "Residual", color = :red, linewidth = 2)
   
      fig_budget[2, 2] = Legend(fig_budget, ax_budget)
      fig_budget[1, 1] = Label(fig_budget, "Terms in PKE budget", fontsize = 24,
                               tellwidth = false)

      save(joinpath("./Plots", "PKEbudget_$(datetime).png"), fig_budget)
      
      ##########################################
      # Plot PKE-budget terms on separate axes #
      ##########################################
      
      fig_separated = Figure(size = (1200, 1050))
      ax_sep_dt     = Axis(fig_separated[2, 1]; ylabel = L"[s$^{-1}$]", 
                           title = L"$\frac{d}{dt} \iiint$ PKE $dV$")
      ax_sep_PAPE   = Axis(fig_separated[2, 2], title = L"$\iiint b' u_z'\, dV$")
      ax_sep_BTI    = Axis(fig_separated[3, 1]; xlabel = "Time [days]", 
                           ylabel = L"[s$^{-1}$]", title = "Energy conversion due to background horizontal shear")
      ax_sep_BCI    = Axis(fig_separated[3, 2]; xlabel = "Time [days]", 
                           title = "Energy conversion due to background vertical shear")
   
      lines!(ax_sep_dt, t[3:end-1], dt_PKE_total, color = :yellowgreen,
             linewidth = 2)
      lines!(ax_sep_PAPE, t[2:end], PAPE_to_PKE, color = :navy, linewidth = 2)
      #=scatter!(ax_sep_BTI, t[2:end], BTI_transfer, color = :hotpink, linewidth = 2)
      scatter!(ax_sep_BCI, t[2:end], BCI_transfer, color = :darkgreen, linewidth = 2)
      =#
      lines!(ax_sep_BTI, t[2:end], total_bkgd_xyx[1:end], color = :hotpink, 
             linewidth = 2)
      lines!(ax_sep_BCI, t[2:end], total_bkgd_xzx[1:end], color = :darkgreen, 
             linewidth = 2)
   
      fig_separated[1, 1:2] = Label(fig_separated, 
                                    "Individual terms in PKE budget", 
                                    fontsize = 24, tellwidth = false)

      save(joinpath("./Plots", "PKEbudgetTerms_$(datetime).png"), fig_separated)
   end
   
   if simParams.save_individual_Cartesian_bkgd_energy_to_PKE_terms
   
      total_bkgd_xxx = ds[:total_bkgd_xxx][chron_idcs][2:end]
      total_bkgd_xyx = ds[:total_bkgd_xyx][chron_idcs][2:end]
      total_bkgd_xzx = ds[:total_bkgd_xzx][chron_idcs][2:end]
      total_bkgd_yxy = ds[:total_bkgd_yxy][chron_idcs][2:end]
      total_bkgd_yyy = ds[:total_bkgd_yyy][chron_idcs][2:end]
      total_bkgd_yzy = ds[:total_bkgd_yzy][chron_idcs][2:end]
      total_bkgd_yzy = ds[:total_bkgd_yzy][chron_idcs][2:end]
      total_bkgd_zxz = ds[:total_bkgd_zxz][chron_idcs][2:end]
      total_bkgd_zyz = ds[:total_bkgd_zyz][chron_idcs][2:end]
      total_bkgd_zzz = ds[:total_bkgd_zzz][chron_idcs][2:end]
      
      #############################################################
      # Plot (Cartesian) terms in total background-->PKE transfer #
      #############################################################
      
      fig_bkgd_terms = Figure(size = (1200, 700))
      ax_bkgd_terms  = Axis(fig_bkgd_terms[2, 1]; xlabel = "Time [days]", 
                            ylabel = L"[s$^{-1}$]")
                            
      beta_x = lines!(ax_bkgd_terms, t[2:end], total_bkgd_xxx, color = :blue)
      lines!(ax_bkgd_terms, t[2:end], total_bkgd_xyx, color = :blue, 
             linestyle = :dash)
      lines!(ax_bkgd_terms, t[2:end], total_bkgd_xzx, color = :blue, 
             linestyle = :dot)
   
      beta_y = lines!(ax_bkgd_terms, t[2:end], total_bkgd_yxy, color = :orange)
      lines!(ax_bkgd_terms, t[2:end], total_bkgd_yyy, color = :orange, 
             linestyle = :dash)
      lines!(ax_bkgd_terms, t[2:end], total_bkgd_yzy, color = :orange, 
             linestyle = :dot)
   
      beta_z = lines!(ax_bkgd_terms, t[2:end], total_bkgd_zxz, color = :red)
      lines!(ax_bkgd_terms, t[2:end], total_bkgd_zyz, color = :red, 
             linestyle = :dash)
      lines!(ax_bkgd_terms, t[2:end], total_bkgd_zzz, color = :red, 
             linestyle = :dot)
             
      #Plot elements for use in legend
      alpha_x, alpha_y, alpha_z = ([LineElement(color = :gray)],
                                   [LineElement(color = :gray, 
                                                linestyle = :dash)],
                                   [LineElement(color = :gray, 
                                                linestyle = :dot)])
                                
      legend_group_color     = [beta_x, beta_y, beta_z]
      legend_group_linestyle = [alpha_x, alpha_y, alpha_z]
      
      legend = Legend(fig_bkgd_terms[2, 2], 
                      [legend_group_linestyle, legend_group_color],
                      [[L"$x$ (i.e., along-background-flow shear)", 
                        L"$y$ (i.e., across-background-flow shear)", 
                        L"$z$ (i.e., vertical background-flow shear)"],
                       [L"$x$", L"$y$", L"$z$"]],
                      [L"$\alpha$", L"$\beta$"], nbanks = 1, 
                      titleposition = :left)
                   
      fig_bkgd_terms[1, 1:2] = Label(fig_bkgd_terms, 
                                     L"Background-energy-to-PKE terms $\iiint u_{\alpha}′ u_{\beta}′ \frac{\partial U_{\beta}}{\partial x_{\alpha}}\, dV$")

      save(joinpath("./Plots", "backgroundTransferToPKE_individualTerms_$(datetime).png"), fig_bkgd_terms)
   end

   if simParams.save_PKE_gyre_budget_terms
     
      gyre_PKE          = ds[:gyre_PKE][chron_idcs]
      gyre_PAPE_to_PKE  = ds[:gyre_PAPE_to_PKE][chron_idcs]
      gyre_BTI_transfer = ds[:gyre_BTI_transfer][chron_idcs]
      gyre_BCI_transfer = ds[:gyre_BCI_transfer][chron_idcs]
      
      #Convert time to s to compute time-derivative of PKE
      dt_PKE_gyre = order1_forward_difference(86400 * t, gyre_PKE)
      
      ###########################################
      # Plot PKE timeseries in gyre-core region #
      ###########################################
      
      fig_PKEgyre = Figure(size = (1200, 700))
      ax_PKEgyre  = Axis(fig_PKEgyre[2, 1]; xlabel = "Time [days]",
                         ylabel = "PKE in gyre region", yscale = log10)

      scatter!(ax_PKEgyre, t[2:end], gyre_PKE[2:end], color = :black)

      fig_PKEgyre[1, 1] = Label(fig_PKEgyre,
                                "Specific PKE integrated over gyre region",
                     	          fontsize = 24, tellwidth = false)

      save(joinpath("./Plots", "PKEgyre_$(datetime).png"), fig_PKEgyre)
      
      #######################################
      # Plot PKE budget in gyre-core region #
      #######################################

      fig_gyreBudget = Figure(size = (1200, 700))
      ax_gyreBudget  = Axis(fig_gyreBudget[2, 1]; xlabel = "Time [days]",
                            ylabel = L"[s$^{-1}$]")

      lines!(ax_gyreBudget, t[1:end-1], dt_PKE_gyre, 
             label = "Time derivative of total PKE", color = :yellowgreen, 
             linewidth = 2)
      lines!(ax_gyreBudget, t, gyre_PAPE_to_PKE, label = "PAPE to PKE",
             color = :navy, linewidth = 2)
      lines!(ax_gyreBudget, t, gyre_BTI_transfer, label = "BTI transfer", 
             color = :hotpink, linewidth = 2)
      lines!(ax_gyreBudget, t, gyre_BCI_transfer, label = "BCI transfer", 
             color = :darkgreen, linewidth = 2)
      lines!(ax_gyreBudget, t[1:end-1], 
             (dt_PKE_gyre - gyre_PAPE_to_PKE[1:end-1] 
              - gyre_BTI_transfer[1:end-1] - gyre_BCI_transfer[1:end-1]), 
             label = "Residual", color = :red, linewidth = 2)

      fig_gyreBudget[2, 2] = Legend(fig_gyreBudget, ax_gyreBudget)
      fig_gyreBudget[1, 1] = Label(fig_gyreBudget,
			                             "Terms in PKE budget, integrated over gyre region",
			                             fontsize = 24, tellwidth = false)

      save(joinpath("./Plots", "PKEgyreBudget_$(datetime).png"), fig_gyreBudget)
   end

   close(ds)
end

################################################################################

function visualize_b_and_ζ(datetime, const_dimension, const_dimension_idx, grid,
                           backgroundFields, Hx, Hy, Hz; 
		                       plot_animation = true, t_idx_skip = 1, 
                           colorbar_max_fraction = 1)
   #=
   Plot 2D slices of buoyancy and vertical component of vorticity. 
   By default, data are assumed to exclude halos.
   =#

   outfile_list = glob("./Output/output_$(datetime)*")
   
   ds_i = open_nc_dataset(outfile_list[1]; Hx = Hx, Hy = Hy, Hz = Hz)[1]
   ds_f, x, y, zC, zF, times, Nt, chron_idcs = open_nc_dataset(outfile_list[length(outfile_list)]; Hx = Hx, Hy = Hy, Hz = Hz)

   if const_dimension == "x"
      const_dimension_coords, axis1, axis2_zC, axis2_zF = x, y, zC, zF
      H_axis1, H_axis2_zC, H_axis2_zF                   = Hy, Hz, Hz
   elseif const_dimension == "y"
      const_dimension_coords, axis1, axis2_zC, axis2_zF = y, x, zC, zF
      H_axis1, H_axis2_zC, H_axis2_zF                   = Hx, Hz, Hz
   elseif const_dimension == "z"
      const_dimension_coords, axis1, axis2_zC, axis2_zF = zC, x, y, y
      H_axis1, H_axis2_zC, H_axis2_zF                   = Hx, Hy, Hy
   end

   #If needed for testing, can set halos to 0 here
   xyzC_idcs, xyzF_idcs = get_2D_spatial_axis_idcs(const_dimension, const_dimension_idx;
                                  Hx = Hx, Hy = Hy, Hz = Hz,
                                  xC = x, yC = y, zC = zC, zF = zF)

   B      = no_offset_view(adapt(Array, backgroundFields[:B]))[xyzC_idcs...]
   ζ_bkgd = no_offset_view(adapt(Array, backgroundFields[:ζ_bkgd])
                          )[xyzC_idcs...]

   b_total_i = ds_i[:b][xyzC_idcs..., chron_idcs[1]]
   b_total_f = ds_f[:b][xyzC_idcs..., chron_idcs[Nt]]
   ζ_total_i = ds_i[:ζ][xyzC_idcs..., chron_idcs[1]]
   ζ_total_f = ds_f[:ζ][xyzC_idcs..., chron_idcs[Nt]]

   Δb_i = b_total_i .- B
   Δb_f = b_total_f .- B
   Δζ_i = ζ_total_i .- ζ_bkgd
   Δζ_f = ζ_total_f .- ζ_bkgd

   nearest, ax_kwargs = get_2D_spatial_axis_kwargs(const_dimension,
                                                   const_dimension_idx, 
                                                   const_dimension_coords)

   #Plot static images (first and final frames, by default)
                               
   fig_total_i = plot_fig_heatmaps_b_ζ(Dict("b" => b_total_i, 
                                            "ζ" => ζ_total_i),
                                       axis1, H_axis1, axis2_zC, H_axis2_zC, 
                                       axis2_zF, H_axis2_zF,
                                       nearest, ax_kwargs, 
              "Fields at $(const_dimension) = $(nearest) km; initial condition",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Total buoyancy ($b$)",
                               "ζ" => L"Total vertical vorticity ($\zeta$)")
                                      )
   fig_total_f = plot_fig_heatmaps_b_ζ(Dict("b" => b_total_f, 
                                            "ζ" => ζ_total_f), 
                                       axis1, H_axis1, axis2_zC, H_axis2_zC,
                                       axis2_zF, H_axis2_zF, 
                                       nearest, ax_kwargs, 
    "Fields at $(const_dimension) = $(nearest) km; t = $(@sprintf("%.2f", times[Nt])) days",
              colorbar_max_fraction; 
              ax_titles = Dict("b" => L"Total buoyancy ($b$)",
                               "ζ" => L"Total vertical vorticity ($\zeta$)")
                                      )
   fig_perturb_i = plot_fig_heatmaps_b_ζ(Dict("b" => Δb_i, 
                                              "ζ" => Δζ_i),
                                       axis1, H_axis1, axis2_zC, H_axis2_zC, 
                                       axis2_zF, H_axis2_zF,
                                       nearest, ax_kwargs, 
              "Perturbations at $(const_dimension) = $(nearest) km; initial condition",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Buoyancy perturbation ($b'$)",
                               "ζ" => L"Vertical vorticity perturbation ($\zeta'$)")
                                      )
   fig_perturb_f = plot_fig_heatmaps_b_ζ(Dict("b" => Δb_f, 
                                              "ζ" => Δζ_f), 
                                       axis1, H_axis1, axis2_zC, H_axis2_zC,
                                       axis2_zF, H_axis2_zF, 
                                       nearest, ax_kwargs, 
    "Perturbations at $(const_dimension) = $(nearest) km; t = $(@sprintf("%.2f", times[Nt])) days",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Buoyancy perturbation ($b'$)",
                               "ζ" => L"Vertical vorticity perturbation ($\zeta'$)")
                                      )

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "bzeta_total_$(const_dimension)$(nearest)_ti_$(datetime).png"), fig_total_i)
   save(joinpath("./Plots", "bzeta_total_$(const_dimension)$(nearest)_tf_$(datetime).png"), fig_total_f)
   save(joinpath("./Plots", "bzeta_perturbs_$(const_dimension)$(nearest)_ti_$(datetime).png"), fig_perturb_i)
   save(joinpath("./Plots", "bzeta_perturbs_$(const_dimension)$(nearest)_tf_$(datetime).png"), fig_perturb_f)
 
   close(ds_i)
   close(ds_f)
   
   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip

      ds, x, y, zC, zF, times, Nt, chron_idcs = open_nc_dataset(outfile_list, Hx = Hx, Hy = Hy, Hz = Hz)

      n = Observable(1)
      
      b_total = @lift ds[:b][xyzC_idcs..., chron_idcs[$n]]
      ζ_total = @lift ds[:ζ][xyzC_idcs..., chron_idcs[$n]]

      Δb = @lift $b_total .- B
      Δζ = @lift $ζ_total .- ζ_bkgd
      
      fig_total = set_up_anim_heatmaps_b_ζ(
               Dict("b" => b_total, "ζ" => ζ_total), 
               axis1, H_axis1, axis2_zC, H_axis2_zC, axis2_zF, H_axis2_zF, nearest, ax_kwargs;
                                  colorranges = Dict("b" => get_symm_range_lims(b_total_f),
                                                     "ζ" => get_symm_range_lims(ζ_total_f)) 
                                  )
      fig_pert = set_up_anim_heatmaps_b_ζ(
               Dict("b" => Δb, "ζ" => Δζ), 
               axis1, H_axis1, axis2_zC, H_axis2_zC, axis2_zF, H_axis2_zF, nearest, ax_kwargs; 
                                  ax_titles = Dict("b" => L"$b'$", "ζ" => L"$\zeta'$"),
                                  colorranges = Dict("b" => get_symm_range_lims(Δb_f),
                                                     "ζ" => get_symm_range_lims(Δζ_f)) 
                                  )
                                  
      title_total = @lift "Fields at $(const_dimension) = $(@sprintf("%.2f", nearest)) km; t = $(@sprintf("%.2f", times[$n])) days"
      title_pert  = @lift "Perturbation fields at $(const_dimension) = $(@sprintf("%.2f", nearest)) km; t = $(@sprintf("%.2f", times[$n])) days"

      fig_total[1, 1:4] = Label(fig_total, title_total, fontsize = 24,
                                tellwidth = false)
      fig_pert[1, 1:4]  = Label(fig_pert, title_pert, fontsize = 24, 
                                tellwidth = false)
                                
      frames = 1:Nt

      video_total = VideoStream(fig_total, format = "mp4", framerate = 6)
      video_pert  = VideoStream(fig_pert, format = "mp4", framerate = 6)

      for i = 1:t_idx_skip:frames[end]
         
         print(i, " of ", Nt, "\n")
         
         recordframe!(video_total)
         recordframe!(video_pert)

         yield()
         n[] = i
      end

      save(joinpath("./Plots", "bzeta_total_$(const_dimension)$(nearest)_$(datetime).mp4"), video_total)
      save(joinpath("./Plots", "bzeta_perturbs_$(const_dimension)$(nearest)_$(datetime).mp4"), video_pert)

      close(ds)
   end
end

################################################################################

function visualize_z_grid(datetime, grid, zmin; zmax = 0.0)

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

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "zgrid_$(datetime).png"), fig)
end

################################################################################

function visualize_fields_2D_slice(datetime, const_dimension, 
                                   const_dimension_idx, backgroundFields,
                                   Hx, Hy, Hz;
                                   plot_animation = true, t_idx_skip = 1,
                                   colorbar_max_fraction = 1,
                                   plot_speed_animation = false)
   #=
   Plot 2D slices of prognostic fields and, optionally, horizontal speed. 
   By default, data are assumed to exclude halos.
   =#

   outfile_list = glob("./Output/output_$(datetime)*")
   
   #Open initial and final datasets
   ds_i = open_nc_dataset(outfile_list[1]; Hx = Hx, Hy = Hy, Hz = Hz)[1]
   ds_f, x, y, zC, zF, times, Nt, chron_idcs = open_nc_dataset(outfile_list[length(outfile_list)]; Hx = Hx, Hy = Hy, Hz = Hz)

   if const_dimension == "x"
      const_dimension_coords, axis1, axis2_zC, axis2_zF = x, y, zC, zF
      H_axis1, H_axis2_zC, H_axis2_zF                   = Hy, Hz, Hz
   elseif const_dimension == "y"
      const_dimension_coords, axis1, axis2_zC, axis2_zF = y, x, zC, zF
      H_axis1, H_axis2_zC, H_axis2_zF                   = Hx, Hz, Hz
   elseif const_dimension == "z"
      const_dimension_coords, axis1, axis2_zC, axis2_zF = zC, x, y, y
      H_axis1, H_axis2_zC, H_axis2_zF                   = Hx, Hy, Hy
   end

   #If needed for testing, can set halos to 0 here
   xyzC_idcs, xyzF_idcs = get_2D_spatial_axis_idcs(const_dimension, const_dimension_idx;
                                  Hx = Hx, Hy = Hy, Hz = Hz,
                                  xC = x, yC = y, zC = zC, zF = zF)

   B  = no_offset_view(adapt(Array, backgroundFields[:B]))[xyzC_idcs...]
   Ur = no_offset_view(adapt(Array, backgroundFields[:Ur]))[xyzC_idcs...]
   Uφ = no_offset_view(adapt(Array, backgroundFields[:Uφ]))[xyzC_idcs...]
   Uz = no_offset_view(adapt(Array, backgroundFields[:Uz]))[xyzF_idcs...]

   b_total_i  = ds_i[:b][xyzC_idcs..., chron_idcs[1]]
   b_total_f  = ds_f[:b][xyzC_idcs..., chron_idcs[Nt]]
   ur_total_i = ds_i[:ur][xyzC_idcs..., chron_idcs[1]]
   ur_total_f = ds_f[:ur][xyzC_idcs..., chron_idcs[Nt]]
   uφ_total_i = ds_i[:uφ][xyzC_idcs..., chron_idcs[1]]
   uφ_total_f = ds_f[:uφ][xyzC_idcs..., chron_idcs[Nt]]
   uz_total_i = ds_i[:uz][xyzF_idcs..., chron_idcs[1]]
   uz_total_f = ds_f[:uz][xyzF_idcs..., chron_idcs[Nt]]

   Δb_i  = b_total_i .- B
   Δb_f  = b_total_f .- B
   Δur_i = ur_total_i .- Ur
   Δur_f = ur_total_f .- Ur
   Δuφ_i = uφ_total_i .- Uφ
   Δuφ_f = uφ_total_f .- Uφ
   Δuz_i = uz_total_i .- Uz
   Δuz_f = uz_total_f .- Uz

   nearest, ax_kwargs = get_2D_spatial_axis_kwargs(const_dimension,
                                                   const_dimension_idx, 
                                                   const_dimension_coords)

   #Plot static images of initial and final frames
   
   fig_total_i = plot_fig_heatmaps_b_ur_uφ_uz(
              Dict("b" => b_total_i, "ur" => ur_total_i, "uφ" => uφ_total_i, 
                   "uz" => uz_total_i), 
              axis1, H_axis1, axis2_zC, H_axis2_zC, axis2_zF, H_axis2_zF,
              nearest, ax_kwargs, 
              "Fields at $(const_dimension) = $(nearest) km; initial condition",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Total buoyancy ($b$)",
                               "ur" => L"Total radial velocity ($u_r$)", 
                               "uφ" => L"Total azimuthal velocity ($u_{\phi}$)",
                               "uz" => L"Total vertical velocity ($u_z$)")
                                         )
   fig_total_f = plot_fig_heatmaps_b_ur_uφ_uz(
              Dict("b" => b_total_f, "ur" => ur_total_f, "uφ" => uφ_total_f, 
                   "uz" => uz_total_f), 
              axis1, H_axis1, axis2_zC, H_axis2_zC, axis2_zF, H_axis2_zF,
              nearest, ax_kwargs, 
              "Fields at $(const_dimension) = $(nearest) km; t = $(@sprintf("%.2f", times[Nt])) days",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Total buoyancy ($b$)",
                               "ur" => L"Total radial velocity ($u_r$)", 
                               "uφ" => L"Total azimuthal velocity ($u_{\phi}$)",
                               "uz" => L"Total vertical velocity ($u_z$)")
                                         )
   fig_pert_i  = plot_fig_heatmaps_b_ur_uφ_uz(
              Dict("b" => Δb_i, "ur" => Δur_i, "uφ" => Δuφ_i, "uz" => Δuz_i), 
              axis1, H_axis1, axis2_zC, H_axis2_zC, axis2_zF, H_axis2_zF,
              nearest, ax_kwargs, 
              "Perturbation fields at $(const_dimension) = $(nearest) km; initial condition",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Buoyancy perturbation ($b'$)",
                               "ur" => L"Radial velocity perturbation ($u_r'$)", 
                               "uφ" => L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                               "uz" => L"Vertical velocity perturbation ($u_z'$)")
                                         )
   fig_pert_f  = plot_fig_heatmaps_b_ur_uφ_uz(
              Dict("b" => Δb_f, "ur" => Δur_f, "uφ" => Δuφ_f, "uz" => Δuz_f), 
              axis1, H_axis1, axis2_zC, H_axis2_zC, axis2_zF, H_axis2_zF,
              nearest, ax_kwargs, 
              "Perturbation fields at $(const_dimension) = $(nearest) km; t = $(@sprintf("%.2f", times[Nt])) days",
              colorbar_max_fraction;
              ax_titles = Dict("b" => L"Buoyancy perturbation ($b'$)",
                               "ur" => L"Radial velocity perturbation ($u_r'$)", 
                               "uφ" => L"Azimuthal velocity perturbation ($u_{\phi}'$)",
                               "uz" => L"Vertical velocity perturbation ($u_z$')")
                                         )

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "fields_$(const_dimension)$(nearest)_ti_$(datetime).png"), fig_total_i)
   save(joinpath("./Plots", "fields_$(const_dimension)$(nearest)_tf_$(datetime).png"), fig_total_f)
   save(joinpath("./Plots", "perturbs_$(const_dimension)$(nearest)_ti_$(datetime).png"), fig_pert_i)
   save(joinpath("./Plots", "perturbs_$(const_dimension)$(nearest)_tf_$(datetime).png"), fig_pert_f)
   
   close(ds_i)
   close(ds_f)

   if plot_animation #Plot animated fields, slicing timeseries at t_idx_skip
                    
      if plot_speed_animation #Plot animation of (horizontal) speed
      
         fig_speed = Figure(size = (800, 500))
         ax_speed  = Axis(fig_speed[2, 1]; ax_kwargs...)
         
         speed_f = sqrt.(ur_total_f.^2 .+ uφ_total_f.^2)
         
         lims_speed    = get_symm_range_lims(speed_f)
         lims_speed[1] = 0
      end

      ds, x, y, zC, zF, times, Nt, chron_idcs = open_nc_dataset(outfile_list, Hx = Hx, Hy = Hy, Hz = Hz)

      n = Observable(1)

      b_total  = @lift ds[:b][xyzC_idcs..., chron_idcs[$n]]
      ur_total = @lift ds[:ur][xyzC_idcs..., chron_idcs[$n]]
      uφ_total = @lift ds[:uφ][xyzC_idcs..., chron_idcs[$n]]
      uz_total = @lift ds[:uz][xyzF_idcs..., chron_idcs[$n]]

      Δb  = @lift $b_total .- B
      Δur = @lift $ur_total .- Ur
      Δuφ = @lift $uφ_total .- Uφ
      Δuz = @lift $uz_total .- Uz
      
      fig_total = set_up_anim_heatmaps_b_ur_uφ_uz(
               Dict("b" => b_total, "ur" => ur_total, "uφ" => uφ_total, 
                    "uz" => uz_total), 
                                  axis1, H_axis1, axis2_zC, H_axis2_zC, 
                                  axis2_zF, H_axis2_zF, nearest, ax_kwargs,
                                  colorbar_max_fraction;
                                  colorranges = Dict("b" => get_symm_range_lims(b_total_f),
                                                     "ur" => get_symm_range_lims(ur_total_f), 
                                                     "uφ" => get_symm_range_lims(uφ_total_f), 
                                                     "uz" => get_symm_range_lims(uz_total_f)), 
                                  )
      fig_pert = set_up_anim_heatmaps_b_ur_uφ_uz(
               Dict("b" => Δb, "ur" => Δur, "uφ" => Δuφ, "uz" => Δuz), 
                                  axis1, H_axis1, axis2_zC, H_axis2_zC, 
                                  axis2_zF, H_axis2_zF, nearest, ax_kwargs,
                                  colorbar_max_fraction; 
                                  ax_titles = Dict("b" => L"$b'$", 
                                                   "ur" => L"$u_r'$", 
                                                   "uφ" => L"$u_{\phi}'$", 
                                                   "uz" => L"$u_z'$"),
                                  colorranges = Dict("b" => get_symm_range_lims(Δb_f),
                                                     "ur" => get_symm_range_lims(Δur_f), 
                                                     "uφ" => get_symm_range_lims(Δuφ_f), 
                                                     "uz" => get_symm_range_lims(Δuz_f)), 
                                  )

      title_total = @lift "Fields at $(const_dimension) = $(@sprintf("%.2f", nearest)) km; t = $(@sprintf("%.2f", times[$n])) days"
      title_pert  = @lift "Perturbation fields at $(const_dimension) = $(@sprintf("%.2f", nearest)) km; t = $(@sprintf("%.2f", times[$n])) days"

      fig_total[1, 1:4] = Label(fig_total, title_total, fontsize = 24,
                                tellwidth = false)
      fig_pert[1, 1:4]  = Label(fig_pert, title_pert, fontsize = 24, 
                                tellwidth = false)

      if plot_speed_animation

         speed_total = @lift sqrt.($ur_total.^2 .+ $uφ_total.^2)
         
         hm_speed_total = heatmap!(ax_speed, axis1, axis2_zC, speed_total,
                                   colormap = :Reds, colorrange = lims_speed)
         
         ux_total = @lift ds[:ux][xyzC_idcs..., chron_idcs[$n]]
         uy_total = @lift ds[:uy][xyzC_idcs..., chron_idcs[$n]]

         Nh, Nv = length(axis1), length(axis2_zC)
         
         #arrows2d!(ax_speed, axis1[1:(Nh ÷ 25):end], 
         #                axis2_zC[1:(Nv ÷ 25):end],
         #                ux_total[1:(Nh ÷ 25):end, 1:(Nv ÷ 25):end], 
         #                uy_total[1:(Nh ÷ 25):end, 1:(Nv ÷ 25):end],
         #                minshaftlength = 0, lengthscale = 500)
         
         Colorbar(fig_speed[2, 2], hm_speed_total, tickformat = "{:.1e}", label = "m/s")
         
         title_speed = @lift @sprintf("Total horizontal speed at %s = %.2f km; t = %.2f days",
                                      const_dim, nearest, times[$n])
                                        
         fig_speed[1, 1:2] = Label(fig_speed, title_speed, fontsize = 24, 
                                   tellwidth = false)
         video_speed       = VideoStream(fig_speed, format = "mp4", 
                                         framerate = 6)
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

      save(joinpath("./Plots", "fields_$(const_dimension)$(nearest)_$(datetime).mp4"), video_total)
      save(joinpath("./Plots", "perturbs_$(const_dimension)$(nearest)_$(datetime).mp4"), video_pert)
      
      if plot_speed_animation
         save(joinpath("./Plots", "horizontalSpeed_$(const_dim)$(nearest)_$(datetime).mp4"), video_speed)
      end
      
      close(ds)
   end
end

################################################################################

function plot_rzMean_φDFT_field(scalarField, grid, gridParams, Nr, Nφ, 
                                fieldName; 
                                Hx = nothing, Hy = nothing, Hz = nothing)
   #=
   Compute azimuthal Fourier transform of 'scalarField', which lives on an
    Oceananigans rectilinear 'grid', by first interpolating to a uniformly spaced
    Nr x Nφ x Nz cylindrical grid. 
   By default, 'scalarField' is assumed to include the same halos as 'grid'.
   Compute the combined radial + vertical mean, then plot the average azimuthal
    spectrum.
   =#
   
   kφs, rzMean_φFFTfield = compute_φFFT_of_sim_Field(scalarField, grid,
                                                     gridParams, Nr, Nφ;
                                                     Hx = Hx, Hy = Hy, Hz = Hz,
                                                     rAvg = true, zAvg = true)
   
   FFT_vec = rzMean_φFFTfield[1, :, 1] #Index to remove singleton dimensions
   Nkφ     = length(FFT_vec)
   
   fig = Figure(size = (1200, 700))
   ax  = Axis(fig[2, 1]; xlabel = "Azimuthal wavenumber",
                      ylabel = "Fourier component")
                      
   scatter!(ax, kφs, (abs2.(FFT_vec) ./ Nkφ), color = :black)
   
   fig[1, 1] = Label(fig, 
                     L"$rz$-averaged azimuthal Fourier components of %$(fieldName)",
                     fontsize = 24, tellwidth = false)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "rzMeanAzimuthalSpectrum_$(fieldName).png"), fig)
end

################################################################################

function plot_rMean_φDFT_field(scalarField, grid, gridParams, Nr, Nφ, fieldName;
                               Hx = nothing, Hy = nothing, Hz = nothing)
   #=
   Compute azimuthal Fourier transform of 'scalarField', which lives on an
    Oceananigans rectilinear 'grid', by first interpolating to a uniformly spaced
    Nr x Nφ x Nz cylindrical grid. 
   By default, 'scalarField' is assumed to include the same halos as 'grid'.
   Compute the radial mean, then plot the r-averaged azimuthal spectrum as a
    heatmap in {kφ} x {z} space.
   =#
   
   if isnothing(Hz)
      Hz = gridParams.Hz
   end

   kφs, rMean_φFFTfield = compute_φFFT_of_sim_Field(scalarField, grid,
                                                    gridParams, Nr, Nφ;
                                                    Hx = Hx, Hy = Hy, Hz = Hz,
                                                    rAvg = true, zAvg = false)
   
   FFT_array = rMean_φFFTfield[1, :, :] #Index to remove singleton dimension
   Nkφ       = length(FFT_array[:, 1])
   
   fig = Figure(size = (1200, 700))
   ax  = Axis(fig[2, 1]; xlabel = "Azimuthal wavenumber", ylabel = "z [m]")

   hm = heatmap!(ax, kφs, grid.z.cᵃᵃᶜ[1:(end - Hz)], (abs2.(FFT_array) ./ Nkφ),
                 colormap = :Reds)
   
   Colorbar(fig[2, 2], hm, tickformat = "{:.1e}", label = "Fourier component")
   
   fig[1, 1] = Label(fig, L"$r$-averaged Fourier azimuthal components of %$(fieldName)",
                     fontsize = 24, tellwidth = false)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "rMeanAzimuthalSpectrum_$(fieldName).png"), fig)
end

################################################################################

function plot_zMean_φDFT_field(scalarField, grid, gridParams::GridParams, Nr, 
                               Nφ, fieldName;
                               Hx = nothing, Hy = nothing, Hz = nothing)
   #=
   Compute azimuthal Fourier transform of 'scalarField', which lives on an
    Oceananigans rectilinear 'grid', by first interpolating to a uniformly spaced
    Nr x Nφ x Nz cylindrical grid. 
   By default, 'scalarField' is assumed to include the same halos as 'grid'.
   Compute the vertical mean, then plot the z-averaged azimuthal spectrum as a
    heatmap in {r} x {kφ} space.
   =#

   kφs, zMean_φFFTfield = compute_φFFT_of_sim_Field(scalarField, grid,
                                                    gridParams, Nr, Nφ;
                                                    Hx = Hx, Hy = Hy, Hz = Hz,
                                                    rAvg = false, zAvg = true)
   
   FFT_array = zMean_φFFTfield[:, :, 1] #Index to remove singleton dimension
   Nkφ       = length(FFT_array[1, :])
   
   r_coords = construct_uniform_polar_grid(gridParams.Lr, Nr, Nφ)[1]
   
   fig = Figure(size = (1200, 700))
   ax  = Axis(fig[2, 1]; xlabel = "r [m]", ylabel = "Azimuthal wavenumber")

   hm = heatmap!(ax, r_coords, kφs, (abs2.(FFT_array) ./ Nkφ), colormap = :Reds)
   
   Colorbar(fig[2, 2], hm, tickformat = "{:.1e}", label = "Fourier component")
   
   fig[1, 1] = Label(fig, L"$z$-averaged Fourier azimuthal components of %$(fieldName)",
                     fontsize = 24, tellwidth = false)
   
   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "zMeanAzimuthalSpectrum_$(fieldName).png"), fig)
end

################################################################################

function plot_NonhydrostaticModel_pressures(model; 
               x_idx = nothing, y_idx = nothing, z_idx = nothing, 
               Hx = 3, Hy = 3, Hz = 3)

   #Note: this function needs to be updated to generalize to plots on a constant-x or constant-z slice. Only handles constant-y case right now.

   @at (Center, Center, Center) @compute residual = ∂z(model.pressures.pHY′) .- model.tracers.b

   #Convert coordinates to km for readability
   x  = no_offset_view(model.grid.xᶜᵃᵃ) ./ 1000
   y  = no_offset_view(model.grid.yᵃᶜᵃ) ./ 1000
   zC = no_offset_view(model.grid.z.cᵃᵃᶜ) ./ 1000
   
   #Coordinates in grid interior
   x_int  = x[(Hx + 1):(end - Hx)]
   y_int  = y[(Hy + 1):(end - Hy)]
   zC_int = zC[(Hz + 1):(end - Hz)]

   fig = Figure(size = (1000, 700))
   ax  = Axis(fig[2, 1], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]")

   hm = heatmap!(ax, x_int, zC_int, no_offset_view(adapt(Array, residual))[:, 25, :], colormap = :balance)

   Colorbar(fig[3, 1], hm, tickformat = "{:.1e}", vertical = false)
   
   fig[1, 1] = Label(fig, L"$\partial p_{HY'} /\partial z - b$", 
                     fontsize = 24, tellwidth = false)

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "hydrostatic_pressure_residual.png"), fig)

   fig        = Figure(size = (2000, 700))
   ax_HY      = Axis(fig[1, 1], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]", 
                     title = "Hydrostatic pressure anomaly")
   ax_NHS     = Axis(fig[1, 2], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]", 
                     title = "Nonhydrostatic pressure")
   ax_b       = Axis(fig[1, 3], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]",
                     title = "Buoyancy")
   ax_HY_int  = Axis(fig[3, 1], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]", 
                     title = "Hydrostatic pressure anomaly, domain interior")
   ax_NHS_int = Axis(fig[3, 2], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]", 
                     title = "Nonhydrostatic pressure, domain interior")
   ax_b_int   = Axis(fig[3, 3], xlabel = L"$x$ [km]", ylabel = L"$z$ [km]", 
                     title = "Buoyancy, domain interior")
   
   hm_HY      = heatmap!(ax_HY, x, zC, no_offset_view(adapt(Array, model.pressures.pHY′))[:, 25, :], colormap = :balance)
   hm_NHS     = heatmap!(ax_NHS, x, zC, no_offset_view(adapt(Array, model.pressures.pNHS))[:, 25, :], colormap = :balance)
   hm_b       = heatmap!(ax_b, x, zC, no_offset_view(adapt(Array, model.tracers.b))[:, 25, :], colormap = :balance)
   hm_HY_int  = heatmap!(ax_HY_int, x_int, zC_int, no_offset_view(adapt(Array, interior(model.pressures.pHY′)))[:, 25, :], colormap = :balance)
   hm_NHS_int = heatmap!(ax_NHS_int, x_int, zC_int, no_offset_view(adapt(Array, interior(model.pressures.pNHS)))[:, 25, :], colormap = :balance)
   hm_b_int   = heatmap!(ax_b_int, x_int, zC_int, no_offset_view(adapt(Array, interior(model.tracers.b)))[:, 25, :], colormap = :balance)

   Colorbar(fig[2, 1], hm_HY, tickformat = "{:.1e}", label = "m²/s²", 
            vertical = false)
   Colorbar(fig[2, 2], hm_NHS, tickformat = "{:.1e}", label = "m²/s²", 
            vertical = false)
   Colorbar(fig[2, 3], hm_b, tickformat = "{:.1e}", label = "m/s²", 
            vertical = false)
   Colorbar(fig[4, 1], hm_HY_int, tickformat = "{:.1e}", label = "m²/s²", 
            vertical = false)
   Colorbar(fig[4, 2], hm_NHS_int, tickformat = "{:.1e}", label = "m²/s²",
            vertical = false)
   Colorbar(fig[4, 3], hm_b_int, tickformat = "{:.1e}", label = "m/s²",
            vertical = false)

   save(joinpath("./Plots", "model_pressures.png"), fig)
end

################################################################################

function open_computed_dataset(datetime, f)
   #=
   Produce nc file containing computed diagnostics, if it does not already
    exist, and return the opened file in read mode.
   =#

   computed_file = joinpath("./Output", "computed_$(datetime).nc")

   if !isfile(computed_file) #Only compute if file does not already exist

      outfile_list = glob("./Output/output_$(datetime)*")

      ds, x, y, zC, zF, times, Nt, chron_idcs = open_nc_dataset(outfile_list)

      if !("q_Ertel" in ds) #Only proceed if PV wasn't saved during simulation

         frames = 1:Nt
         x_idcs = 2:(length(x) - 2)
         y_idcs = 2:(length(y) - 2)
         z_idcs = 2:(length(zC) - 2)

         function update_data_array!(data_array, i, j, k, n, value)
            data_array[i, j, k, n] = value
            return data_array
         end
      
         NCDataset(computed_file, "c") do comp_ds
         
            i, j, k = Observable(2), Observable(2), Observable(2)
            n       = Observable(1)

            b  = @lift ds["b"][:, :, :, $n]
            ux = @lift ds["ux"][:, :, :, $n]
            uy = @lift ds["uy"][:, :, :, $n]
            uz = @lift ds["uz"][:, :, :, $n]
            qn = @lift pointwise_q_Ertel($i, $j, $k, $b, $ux, $uy, $uz, f, 
                                         x, y, zC)

            defDim(comp_ds, "x", length(x)-2)
            defDim(comp_ds, "y", length(y)-2)
            defDim(comp_ds, "z", length(zC)-2)
            defDim(comp_ds, "time", length(times))

            q_data = Array{Float64, 4}(undef, length(x)-2, length(y)-2, 
                                       length(zC)-2, Nt)

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
         
            defVar(comp_ds, "q_Ertel", q_data, ("x", "y", "zC", "time"))
         end #comp_ds gets closed automatically
      end
      
      close(ds)
   end

   return NCDataset(computed_file, "r")
end

function visualize_q_2D_slice(datetime, const_dim, const_idx, f;
                              Hx = 0, Hy = 0, Hz = 0, 
                              plot_animation = true)
                              
   #=
   Plot 2D slices of potential vorticity. By default, data are assumed to exclude
    halos.
   =#

   outfile_list = glob("./Output/output_$(datetime)*")
   
   ds_f, x, y, zC, zF, times, Nt, chron_idcs = open_nc_dataset(outfile_list[length(outfile_list)]; Hx = Hx, Hy = Hy, Hz = Hz)
   
   comp_ds = open_computed_dataset(datetime, f)
   
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
                                  Hx = 1, Hy = 1, Hz = 1,
                                  x_idx = x_idx, y_idx = y_idx, z_idx = z_idx,
                                  xC = x, yC = y, zC = zC, zF = zF)

   q_f    = adapt(Array, comp_ds[:q_Ertel])[xyzC_idcs..., chron_idcs[Nt]]
   lims_q = get_symm_range_lims(q_f)
   
   close(ds_f)
   
   nearest, ax_kwargs = get_2D_spatial_axis_kwargs(x, y, zC, const_dim;
                                                   x_idx = x_idx, 
                                                   y_idx = y_idx,
                                                   z_idx = z_idx)

   #Plot static images (final frame, by default)

   fig_q = Figure(size = (1200, 800))
   ax_q  = Axis(fig_q[2, 1]; ax_kwargs...)
   hm_q  = heatmap!(ax_q, axis1, axis2_zC, q_f, colorrange = lims_q, 
                    colormap = Reverse(:RdBu_5))

   Colorbar(fig_q[2, 2], hm_q, tickformat = "{:.1e}", label = "1/s³")
   
   title_q       = @sprintf("Ertel PV at %s = %.2f km; t = %.2f days",
                            const_dim, nearest, times[Nt])
   fig_q[1, 1:2] = Label(fig_q, title_q, fontsize = 24, tellwidth = false)

   mkpath("./Plots") #Make visualization directory if nonexistent
   save(joinpath("./Plots", "q_$(const_dim)$(nearest)_tf_$(datetime).png"), fig_q)
   
   #Plot animation, slicing timeseries at t_idx_skip

   fig_anim = Figure(size = (1200, 800))
   ax_anim  = Axis(fig_anim[2, 1]; ax_kwargs...)
   
   n = Observable(1)

   q = @lift comp_ds[:q][xyzC_idcs..., chron_idcs[$n]]
   
   hm_anim = heatmap!(ax_anim, axis1, axis2_zC, q, colorrange = lims_q, 
                      colormap = Reverse(:RdBu_5))
   
   Colorbar(fig_anim[2, 2], hm_anim, tickformat = "{:.1e}", label = "1/s³")
   
   title_anim       = @lift @sprintf("Ertel PV at %s = %.2f km; t = %.2f days",
                                     const_dim, nearest, times[$n])
   fig_anim[1, 1:2] = Label(fig_total, title_total, fontsize = 24, 
                            tellwidth = false)

   frames = 1:Nt
   video  = VideoStream(fig_anim, format = "mp4", framerate = 6)

   for i = 1:frames[end]
      recordframe!(video)
      yield()
      print("Plotting frame(s) $(i) of $(frames[end])" * " \r")
      n[] = i
   end

   save(joinpath("./Plots", "q_$(const_dim)$(nearest)_$(datetime).mp4"), video)
   close(comp_ds)
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

   lims_q = get_symm_range_lims(q_xz_f)

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

   lims_q = get_symm_range_lims(q_xy_f)

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