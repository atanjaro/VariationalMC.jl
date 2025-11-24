@doc raw"""

    model_summary( simulation_info::SimulationInfo,
                   nup::Int,
                   ndn::Int,
                   U::Float64,
                   dt::Float64, 
                   model_geometry::ModelGeometry,
                   tight_binding_model::TightBindingModel )::Nothing

Writes model summary to file.

"""
function model_summary(
    simulation_info::SimulationInfo,
    nup::Int,
    ndn::Int,
    U::Float64,
    dt::Float64, 
    N_equil::Int,
    N_opt::Int,
    N_opt_bins::Int,
    model_geometry::ModelGeometry,
    tight_binding_model::TightBindingModel
)::Nothing

    # if process ID is 1
    if iszero(simulation_info.pID)

        # construct full filename, including filepath
        fn = joinpath(simulation_info.datafolder, "model_summary.toml")

        # open file to write to
        open(fn, "w") do fout
            # write nup
            @printf fout "nup = %d\n\n" nup

            # write ndn
            @printf fout "ndn = %d\n\n" ndn

            # write dt
            @printf fout "dt = %.6f\n\n" dt

            # # write model geometry out to file
            # show(fout, "text/plain", model_geometry)

            # # write tight binding model to file
            # show(fout, MIME("text/plain"), tight_binding_model)

            # write U
            @printf fout "U = %.6f\n\n" U

            # write N_equil
            @printf fout "N_equil = %d\n\n" N_equil

            # write N_opt
            @printf fout "N_opt = %d\n\n" N_opt

            # write N_opt_bins
            @printf fout "N_opt_bins = %d\n\n" N_opt_bins

            # write seed
            @printf fout "seed = %d\n\n" seed
        end
    end

    return nothing
end


@doc raw"""

    parameter_summary( simulation_info::SimulationInfo,
                       determinantal_parameters::DeterminantalParameters)::Nothing
    
Writes parameter summary to file. 

"""
function parameter_summary(
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters
)::Nothing

    # if process ID is 1
    if iszero(simulation_info.pID)
        fn = joinpath(simulation_info.datafolder, "parameter_summary.toml")

        open(fn, "w") do fout
            @printf fout "[DeterminantalParameters]\n"
            @printf fout "num_det_pars = %d\n" determinantal_parameters.num_det_pars
            @printf fout "num_det_opts = %d\n\n" determinantal_parameters.num_det_opts

            @printf fout "[DeterminantalParameters.det_pars]\n"
            for (k, v) in determinantal_parameters.det_pars
                if isa(v, AbstractVector)
                    @printf fout "%s = %s\n" string(k), repr(v)
                else
                    @printf fout "%s = %.10g\n" string(k), v
                end
            end
            println(fout)
        end
    end

    return nothing
end


@doc raw"""

    parameter_summary( simulation_info::SimulationInfo,
                       determinantal_parameters::DeterminantalParameters,
                       jastrow_parametrs::JastrowParameters )::Nothing
    
Writes parameter summary to file. 

"""
function parameter_summary(
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parametrs::JastrowParameters
)::Nothing

    # if process ID is 1
    if iszero(simulation_info.pID)
        fn = joinpath(simulation_info.datafolder, "parameter_summary.toml")

        open(fn, "w") do fout
            @printf fout "[DeterminantalParameters]\n"
            @printf fout "num_det_pars = %d\n" determinantal_parameters.num_det_pars
            @printf fout "num_det_opts = %d\n\n" determinantal_parameters.num_det_opts

            @printf fout "[DeterminantalParameters.det_pars]\n"
            for (k, v) in determinantal_parameters.det_pars
                if isa(v, AbstractVector)
                    @printf fout "%s = %s\n" string(k), repr(v)
                else
                    @printf fout "%s = %.10g\n" string(k), v
                end
            end
            println(fout)

            @printf fout "[JastrowParameters]\n"
            @printf fout "jastrow_type = \"%s\"\n" jastrow_parametrs.jastrow_type
            @printf fout "num_jpars = %d\n" jastrow_parametrs.num_jpars
            @printf fout "num_jpar_opts = %d\n\n" jastrow_parametrs.num_jpar_opts

            @printf fout "[JastrowParameters.jpar_map]\n"
            for (k, v) in jastrow_parametrs.jpar_map
                key_str = string(k)
                val_str = isa(v, Tuple) && length(v) == 1 ? repr(v[1]) : repr(v)
                @printf fout "\"%s\" = %s\n" key_str, val_str
            end
        end
    end

    return nothing
end