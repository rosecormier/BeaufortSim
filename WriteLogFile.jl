function WriteLogFile(simulation, duration, logfilepath, outfilepath, gridParams::GridParams, gyreParams::GyreParams, stratParams::DoubleTanhStratParams, perturbParams::PerturbationParams, simParams::SimulationParams)

   open(logfilepath, "w") do file
   
      write(file, "Nx, Ny, Nz = $(gridParams.Nx), $(gridParams.Ny), $(gridParams.Nz) \n")
      write(file, "Lr, Lz = $(gridParams.Lr), $(gridParams.Lz) \n\n")
      
      write(file, "σr, σz = $(gyreParams.σr), $(gyreParams.σz) \n")
      write(file, "U = $(gyreParams.U) \n")
      write(file, "Far-field N² term = $(gyreParams.N²_far) \n\n")
      
      if gyreParams.doubleTanhAmbientStrat
         write(file, "Ambient stratification type: double-tanh \n")
         write(file, "Stratification parameters: g = $(stratParams.gravity),
               ρ₀ = $(stratParams.ρ₀), A_s = $(stratParams.A_s), 
               C_s = $(stratParams.C_s), z_s = $(stratParams.z_s), 
               A_d = $(stratParams.A_d), C_d = $(stratParams.C_d), 
               z_d = $(stratParams.z_d) \n\n")
      else
         write(file, "Ambient stratification type: constant \n")
      end   
      
      write(file, "Max. ψ′ = $(perturbParams.max_ψ′) \n")
      write(file, "Random-number seed = $(perturbParams.seed)\n\n")
      
      write(file, "Checkpoint interval = $(simParams.Δt_checkpt) \n")
      write(file, "align_time_step = $(simParams.align_time_step) \n")
      write(file, "minimum_relative_step = $(simParams.minimum_relative_step) \n\n")
      
      #write(file, "Δt, tf = $(Δt), $(tf) \n\n")
      write(file, "Total number of iterations = $(iteration(simulation)) \n")
      write(file, "Simulation runtime = $(duration) \n")
      write(file, "Output filesize = $(pretty_filesize(filesize(outfilepath)))")
   end

end