include("LibraryCoordinates.jl")
include("Parameters.jl")

using Oceananigans
using Oceananigans.Units

function Progress(s::Simulation, modelParams::ModelParams)

   u_absmax             = maximum(abs, s.model.velocities.u)
   w_absmax             = maximum(abs, s.model.velocities.w)
   b_absmax             = maximum(abs, s.model.tracers.b)
   
   if modelParams.useNonhydrostaticModel
      pHY′_absmax          = maximum(abs, s.model.pressures.pHY′)
      pNHS_absmax          = maximum(abs, s.model.pressures.pNHS)
      interior_pNHS_absmax = maximum(abs, interior(s.model.pressures.pNHS))
      interior_pNHS_absmin = minimum(abs, interior(s.model.pressures.pNHS))
   end
   
   interior_hdiv_absmax = maximum(abs, interior(∂x(s.model.velocities.u) 
                                                + ∂y(s.model.velocities.v)
                                               )
                                 )
   interior_vdiv_absmax = maximum(abs, interior(∂z(s.model.velocities.w)))
   
   @info @sprintf("Iter: %d; time: %.2e days; Δt: %s",
                  iteration(s), (time(s)/day),  prettytime(s.Δt))
   @info @sprintf("max|u|: %.2e; max|w|: %.2e; max|b|: %.2e", 
                  u_absmax, w_absmax, b_absmax)
                  
   if modelParams.useNonhydrostaticModel
      @info @sprintf("max|pHY′|: %.2e", pHY′_absmax)
      @info @sprintf("max|pNHS|: %.2e", pNHS_absmax)
      @info @sprintf("max|pNHS| in interior: %.2e", interior_pNHS_absmax)
      @info @sprintf("min|pNHS| in interior: %.2e", interior_pNHS_absmin)
   end
   
   @info @sprintf("max|horizontal divergence of velocity| in interior: %.2e",
                  interior_hdiv_absmax)
   @info @sprintf("max|vertical divergence of velocity| in interior: %.2e", 
                  interior_vdiv_absmax)
   
   return nothing
end

function StandardSimulationOutputFields(model, backgroundFields)
   #=
   Assemble named tuple of prognostic and secondary Fields, and a named tuple of
    scalars, to be passed to output writers for every simulation, by default.
   The former includes Cartesian and cylindrical velocity components, buoyancy,
    and ζ.
   The latter includes norms of the Cartesian and cylindrical velocity-
    perturbations and of buoyancy perturbation.
   =#
   
   ux, uy, uz = model.velocities

   ur = KernelFunctionOperation{Center, Center, Center}(vector_r_component_ccc, model.grid, ux, uy)
   uφ = KernelFunctionOperation{Center, Center, Center}(vector_φ_component_ccc, model.grid, ux, uy)

   FieldOutputs = (ur = ur,
                   uφ = uφ,
	                 ux = ux,
	                 uy = uy,
	                 uz = uz,
	                 b = model.tracers.b,
                   ζ = CenterField_ζ′(model.grid, ux, uy, backgroundFields[:Ux], backgroundFields[:Uy])
                  )

   ux′_norm(model) = perturbation_norm(ux, backgroundFields[:Ux])
   uy′_norm(model) = perturbation_norm(uy, backgroundFields[:Uy])
   ur′_norm(model) = perturbation_norm(ur, backgroundFields[:Ur])
   uφ′_norm(model) = perturbation_norm(uφ, backgroundFields[:Uφ])
   uz′_norm(model) = perturbation_norm(uz, backgroundFields[:Uz])
   b′_norm(model)  = perturbation_norm(model.tracers.b, backgroundFields[:B])

   scalarOutputs = (ux′_norm = ux′_norm,
		                uy′_norm = uy′_norm,
		                ur′_norm = ur′_norm,
		                uφ′_norm = uφ′_norm,
		                uz′_norm = uz′_norm,
		                b′_norm = b′_norm
                   )

   return FieldOutputs, scalarOutputs
end

function EnergyDiagnosticOutputs(s::Simulation, simParams::SimulationParams, gyreParams::GyreParams, backgroundFields)

   energy_diagnostics = (; )
   
   if simParams.save_PKE_full_domain_budget_terms
      energy_diagnostics = merge(energy_diagnostics,
                                 (total_PKE          = PKE(s, backgroundFields),
                                  total_PAPE_to_PKE  = PAPE_to_PKE(s, backgroundFields),
                                  total_BTI_transfer = BTI_transfer(s, backgroundFields),
                                  total_BCI_transfer = BCI_transfer(s, backgroundFields)
                                 )
                                )
   end
   if simParams.save_PKE_gyre_budget_terms
      energy_diagnostics = merge(energy_diagnostics,
                                 (gyre_PKE            = gyre_PKE(s, gyreParams),
                                  gyre_PAPE_to_PKE    = gyre_PAPE_to_PKE(s, gyreParams),
                                  gyre_BTI_transfer   = gyre_BTI_transfer(s, gyreParams, backgroundFields),
                                  gyre_BCI_transfer   = gyre_BCI_transfer(s, gyreParams, backgroundFields)
                                 )
                                )
   end
   if simParams.save_individual_Cartesian_bkgd_energy_to_PKE_terms
      energy_diagnostics = merge(energy_diagnostics,
                                 (total_bkgd_transfer = bkgd_to_PKE_transfer(s, backgroundFields),
                                  total_bkgd_xxx = bkgd_transfer_ux′ux′∂xUx(s, backgroundFields),
                                  total_bkgd_xyx = bkgd_transfer_ux′uy′∂yUx(s, backgroundFields),
                                  total_bkgd_xzx = bkgd_transfer_ux′uz′∂zUx(s, backgroundFields),
                                  total_bkgd_yxy = bkgd_transfer_uy′ux′∂xUy(s, backgroundFields),
                                  total_bkgd_yyy = bkgd_transfer_uy′uy′∂yUy(s, backgroundFields),
                                  total_bkgd_yzy = bkgd_transfer_uy′uz′∂zUy(s, backgroundFields),
                                  total_bkgd_zxz = bkgd_transfer_uz′ux′∂xUz(s, backgroundFields),
                                  total_bkgd_zyz = bkgd_transfer_uz′uy′∂yUz(s, backgroundFields),
                                  total_bkgd_zzz = bkgd_transfer_uz′uz′∂zUz(s, backgroundFields)
                                 )
                                )
   end
   if simParams.save_total_ME_budgets_terms
      energy_diagnostics = merge(energy_diagnostics,
                                 (total_KE            = total_KE(s),
                                  total_KE_adv_flux   = total_KE_advective_flux(s),
                                  total_pressure_work = total_pressure_work(s, modelParams),
                                  total_KE_production = total_KE_production(s, modelParams),
                                  total_PE            = total_PE(s),
                                  total_int_PE_source = total_internal_PE_sources(s),
                                  total_gravity_work  = total_gravity_work(s, GravitationalAcceleration)
                                 )
                                )
   end
end

function CreateSimulation(model, Δt, tf, Δt_save, simParams::SimulationParams,
                          modelParams::ModelParams, gyreParams::GyreParams, 
                          backgroundFields, runLabel::String)

   sim = Simulation(model;
                    Δt = Δt,
                    stop_time = tf, 
                    align_time_step = simParams.align_time_step, 
                    minimum_relative_step = simParams.minimum_relative_step
                   )
   
   #Callback to print information about simulation at regular intervals
   add_callback!(sim, Progress, TimeInterval(Δt_save); parameters = modelParams)
   
   
   #Diagnostics to output
   Field_outputs, scalar_outputs = StandardSimulationOutputFields(model, 
                                                               backgroundFields)
   energetics_outputs            = EnergyDiagnosticOutputs(sim, simParams,
                                                           gyreParams,
                                                           backgroundFields)
   
   outputFieldsFilePath  = joinpath("./Output", "output_$(runLabel).$(simParams.output_Fields_file_extension)")
   outputScalarsFilePath = joinpath("./Output", "scalars_$(runLabel).nc")
   
   if energetics_outputs != (; )
      energeticsFilePath = joinpath("./Output", "energetics_$(datetimenow).nc")
   end
   
   scalar_writer = NetCDFWriter(model, 
                                   scalar_outputs,
		                               filename = outputScalarsFilePath, 
				                           schedule = TimeInterval(Δt_save),
                                   array_type = Array{Float64},
                                   file_splitting = FileSizeLimit(30GiB),
		                               dimensions = (ux′_norm = (),
					                                       uy′_norm = (),
					                                       ur′_norm = (),
					                                       uφ′_norm = (),
					                                       uz′_norm = (),
					                                       b′_norm  = ()
                                                )
                                  )
                                  
      energy_writer = NetCDFWriter(model, energetics_outputs,
                             filename = energeticsFilePath,
                             schedule = TimeInterval(Δt_save),
                             array_type = Array{Float64},
                             file_splitting = FileSizeLimit(30GiB)
                            )
   
   if simParams.output_Fields_file_extension == "nc"
      
      Field_writer = NetCDFWriter(model, 
                                  Field_outputs,
                                  filename = outputFieldsFilePath, 
                                  schedule = TimeInterval(Δt_save),
                                  array_type = Array{Float64},
                                  with_halos = true,
                                  file_splitting = FileSizeLimit(30GiB)
                                 )

   elseif simParams.output_Fields_file_extension == "jld2"
   
      Field_writer = JLD2Writer(model, 
                                Field_outputs;
                                filename = outputFieldsFilePath,
                                schedule = TimeInterval(Δt_save),
                                array_type = Array{Float64},
                                with_halos = true,
                                file_splitting = FileSizeLimit(30GiB),
                                including = []
                               )

   else
      print("'output_Fields_file_extension' in 'simParams' must be 'nc' or 'jld2'")
   end
   
   sim.output_writers[:field_writer]  = Field_writer
   sim.output_writers[:scalar_writer] = scalar_writer
   
   if energetics_outputs != (; )
      sim.output_writers[:energy_writer] = energy_writer
   end
   
   return sim
end