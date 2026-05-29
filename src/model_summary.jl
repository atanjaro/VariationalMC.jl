@doc raw"""

    model_summary(;
        simulation_info::SimulationInfo,
        model_geometry::ModelGeometry,
        particle_configuration::ParticleConfiguration,
        tight_binding_model::TightBindingModel,
        interactions::Union{Tuple,Nothing} = nothing,
        parameters::Union{Tuple, Nothing} = nothing
    ) 

Write model to summary to file. This includes information on the canonical ensemble configuration,
tight binding model, interactions, and variational parameters.

"""
function model_summary(;
    simulation_info::SimulationInfo,
    model_geometry::ModelGeometry,
    particle_configuration::ParticleConfiguration,
    tight_binding_model::TightBindingModel,
    interactions::Union{Tuple,Nothing} = nothing,
    parameters::Union{Tuple, Nothing} = nothing
) 
    # if process ID is 1
    if iszero(simulation_info.pID)

        # construct full filename, including filepath
        fn = joinpath(simulation_info.datafolder, "model_summary.toml")

        # open file to write to
        open(fn, "w") do fout
            # write model geometry out to file
            show(fout, "text/plain", model_geometry)

            # write particle configuration out to file
            show(fout, "text/plain", particle_configuration)

            # write tight binding model out to file
            show(fout, MIME("text/plain"), tight_binding_model)

            # write various interactions to file
            if !isnothing(interactions)
                for interaction in interactions
                    show(fout, "text/plain", interaction)
                end
            end

            # write various parameters to file
            if !isnothing(parameters)
                for parameter in parameters
                    show(fout, "text/plain", parameter)
                end
            end
        end
    end

    return nothing
end
