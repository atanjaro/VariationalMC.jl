@doc raw"""

    model_summary( simulation_info::SimulationInfo,
                   metadata )::Nothing

Writes model summary to file.

"""
function model_summary(
    simulation_info::SimulationInfo,
    metadata
)::Nothing

    # if process ID is 1
    if iszero(simulation_info.pID)

        # construct full filename, including filepath
        fn = joinpath(simulation_info.datafolder, "model_summary.toml")

        # open file to write to
        open(fn, "w") do fout
            # write pht
            @printf fout "Particle-hole transformed = %d\n\n" metadata["pht"]

            # write dt
            @printf fout "Optimization rate = %d\n\n" metadata["dt"]

            # write N_equil
            @printf fout "Number of equilibration/thermalization steps = %d\n\n" metadata["N_equil"]

            # write N_opt
            @printf fout "Number of optimization steps = %d\n\n" metadata["N_opt"]

            # write N_sim
            @printf fout "Number of simulation steps = %d\n\n" metadata["N_sim"]

            # write N_opt_bins
            @printf fout "Optimization bins = %d\n\n" metadata["N_opt_bins"]

            # write N_sim_bins
            @printf fout "Simulation bins = %d\n\n" metadata["N_sim_bins"]

            # write seed
            @printf fout "RNG seed = %d\n\n" metadata["seed"]

            # write total VMC time
            @printf fout "Total VMC time = %d\n\n" metadata["vmc_time"]

            # write opt_flags
            show(fout, MIME("text/plain"), metadata["opt_flags"])

            # # write model geometry out to file
            # show(fout, "text/plain", model_geometry)

            # # write tight binding model to file
            # show(fout, MIME("text/plain"), tight_binding_model)

        end
    end

    return nothing
end

