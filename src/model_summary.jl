@doc raw"""

    model_summary(  simulation_info::SimulationInfo,
                    tight_binding_model::TightBindingModel,
                    hubbard_model::HubbardModel,
                    determinantal_parameters::DeterminantalParameters,
                    model_geometry::ModelGeometry,
                    pht::Bool ) where {T<:AbstractFloat}

Writes model summary to file.

- `simulation_info::SimulationInfo`: contains simulation information. 
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{T}`: Hubbard interaction parameters.
- `determinantal_parameters::DeterminantalParameters`: initial determinantal parameters values.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `pht::Bool`: whether model is particle-hole transformed. 
 
"""
function model_summary(
    simulation_info::SimulationInfo,
    tight_binding_model::TightBindingModel,
    hubbard_model::HubbardModel{T},
    determinantal_parameters::DeterminantalParameters,
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:AbstractFloat}

    # if process ID is 1
    if iszero(simulation_info.pID)

        # construct full filename, including filepath
        fn = joinpath(simulation_info.datafolder, "model_summary.toml")

        # open file to write to
        open(fn, "w") do fout
            # write model geometry out to file
            show(fout, "text/plain", model_geometry)

            # write tight binding model to file
            show(fout, MIME("text/plain"), tight_binding_model)

            # write Hubbard interaction
            show(fout, MIME("text/plain"), hubbard_model)

            # write pht
            @printf fout "[ParticleHole]\n\n"
            @printf(fout, "pht = %s\n\n", pht)

            # write initial values of determinantal parameters
            # Map parameter names → section labels
            name_map = Dict(
                :Δ_0   => "s-wave",
                :Δ_slo => "Larkin-Ovchinnikov",
                :Δ_sff => "Fulde-Ferrell",
                :Δ_d   => "d-wave",
                :Δ_dlo => "Larkin-Ovchinnikov",
                :Δ_dff => "Fulde-Ferrell",
                :Δ_sx  => "spin-x",
                :Δ_sz  => "spin-z",
                :Δ_ssd => "site-dependent spin",
                :μ     => "chemical potential",
                :Δ_cdw => "charge density wave",
                :Δ_csd => "site-dependent density"
            )

            @printf fout "[DeterminantalParameters]\n\n"

            for (key, val) in pairs(determinantal_parameters.det_pars)
                section = get(name_map, key, string(key))
                @printf fout "[[%s]]\n\n" section

                if isa(val, Number)
                    @printf fout "%s = %.10g\n\n" key val
                elseif isa(val, AbstractVector)
                    # flatten vector into string safely
                    @printf fout "%s = %s\n\n" key string(val)
                else
                    @printf fout "%s = %s\n\n" key string(val)
                end
            end
        end
    end

    return nothing
end


@doc raw"""

    model_summary(  simulation_info::SimulationInfo,
                    tight_binding_model::TightBindingModel,
                    hubbard_model::HubbardModel{T},
                    determinantal_parameters::DeterminantalParameters,
                    jastrow_parameters::JastrowParameters,
                    model_geometry::ModelGeometry,
                    pht::Bool ) where {T<:AbstractFloat}

Writes model summary to file.

- `simulation_info::SimulationInfo`: contains simulation information. 
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{T}`: Hubbard interaction parameters. 
- `determinantal_parameters::DeterminantalParameters`: initial determinantal parameters values.
- `jastrow_parameters::JastrowParameters`: initial set of Jastrow parameters. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `pht::Bool`: whether model is particle-hole transformed. 


"""
function model_summary(
    simulation_info::SimulationInfo,
    tight_binding_model::TightBindingModel,
    hubbard_model::HubbardModel{T},
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters::JastrowParameters,
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:AbstractFloat}

    # if process ID is 1
    if iszero(simulation_info.pID)

        # construct full filename, including filepath
        fn = joinpath(simulation_info.datafolder, "model_summary.toml")

        # open file to write to
        open(fn, "w") do fout
            # write model geometry out to file
            show(fout, "text/plain", model_geometry)

            # write tight binding model to file
            show(fout, MIME("text/plain"), tight_binding_model)

            # write Hubbard interaction
            show(fout, MIME("text/plain"), hubbard_model)

            # write pht
            @printf fout "[ParticleHole]\n\n"
            @printf(fout, "pht      = %s\n\n", pht)

            # write initial values of determinantal parameters
            # Map parameter names → section labels
            name_map = Dict(
                :Δ_0   => "s-wave",
                :Δ_slo => "Larkin-Ovchinnikov",
                :Δ_sff => "Fulde-Ferrell",
                :Δ_d   => "d-wave",
                :Δ_dlo => "Larkin-Ovchinnikov",
                :Δ_dff => "Fulde-Ferrell",
                :Δ_sx  => "spin-x",
                :Δ_sz  => "spin-z",
                :Δ_ssd => "site-dependent spin",
                :μ     => "chemical potential",
                :Δ_cdw => "charge density wave",
                :Δ_csd => "site-dependent density"
            )

            @printf fout "[DeterminantalParameters]\n\n"

            for (key, val) in pairs(determinantal_parameters.det_pars)
                section = get(name_map, key, string(key))
                @printf fout "[[%s]]\n\n" section

                if isa(val, Number)
                    @printf fout "%s    = %.10g\n\n" key val
                elseif isa(val, AbstractVector)
                    # flatten vector into string safely
                    @printf fout "%s    = %s\n\n" key string(val)
                else
                    @printf fout "%s    = %s\n\n" key string(val)
                end
            end

            # write Jastrow parameters
            @printf fout "[JastrowParameters]\n\n"
            @printf fout "type      = \"%s\"\n\n" jastrow_parameters.jastrow_type

            for (idx, (indices, value)) in jastrow_parameters.jpar_map
                @printf fout "[[Pseudopotential]]\n\n"
                @printf fout "index    = %d\n" idx
                @printf fout "v_%d      = %.17g\n\n" idx value
            end
        end
    end

    return nothing
end


@doc raw"""

    model_summary(  simulation_info::SimulationInfo,
                    tight_binding_model::TightBindingModel,
                    hubbard_model::HubbardModel{T},
                    determinantal_parameters::DeterminantalParameters,
                    jastrow_parameters_1::JastrowParameters,
                    jastrow_parameters_2::JastrowParameters,
                    model_geometry::ModelGeometry,
                    pht::Bool ) where {T<:AbstractFloat}

Writes model summary to file.

- `simulation_info::SimulationInfo`: contains simulation information.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{T}`: Hubbard interaction parameters. 
- `determinantal_parameters::DeterminantalParameters`: initial determinantal parameters values.
- `jastrow_parameters_1::JastrowParameters`: first set of Jastrow parameters.
- `jastrow_parameters_2::JastrowParameters`: second set of Jastrow parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function model_summary(
    simulation_info::SimulationInfo,
    tight_binding_model::TightBindingModel,
    hubbard_model::HubbardModel{T},
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters_1::JastrowParameters,
    jastrow_parameters_2::JastrowParameters,
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:AbstractFloat}

    # if process ID is 1
    if iszero(simulation_info.pID)

        # construct full filename, including filepath
        fn = joinpath(simulation_info.datafolder, "model_summary.toml")

        # open file to write to
        open(fn, "w") do fout
            # write model geometry out to file
            show(fout, "text/plain", model_geometry)

            # write tight binding model to file
            show(fout, MIME("text/plain"), tight_binding_model)

            # write Hubbard interaction
            show(fout, MIME("text/plain"), hubbard_model)

            # write pht
            @printf fout "[ParticleHole]\n\n"
            @printf(fout, "pht      = %s\n\n", pht)

            # write initial values of determinantal parameters
            # Map parameter names → section labels
            name_map = Dict(
                :Δ_0   => "s-wave",
                :Δ_slo => "Larkin-Ovchinnikov",
                :Δ_sff => "Fulde-Ferrell",
                :Δ_d   => "d-wave",
                :Δ_dlo => "Larkin-Ovchinnikov",
                :Δ_dff => "Fulde-Ferrell",
                :Δ_sx  => "spin-x",
                :Δ_sz  => "spin-z",
                :Δ_ssd => "site-dependent spin",
                :μ     => "chemical potential",
                :Δ_cdw => "charge density wave",
                :Δ_csd => "site-dependent density"
            )

            @printf fout "[DeterminantalParameters]\n\n"

            for (key, val) in pairs(determinantal_parameters.det_pars)
                section = get(name_map, key, string(key))
                @printf fout "[[%s]]\n\n" section

                if isa(val, Number)
                    @printf fout "%s    = %.10g\n\n" key val
                elseif isa(val, AbstractVector)
                    # flatten vector into string safely
                    @printf fout "%s    = %s\n\n" key string(val)
                else
                    @printf fout "%s    = %s\n\n" key string(val)
                end
            end

            # write first set of Jastrow parameters
            @printf fout "[JastrowParameters]\n\n"
            @printf fout "type      = \"%s\"\n\n" jastrow_parameters_1.jastrow_type

            for (idx, (indices, value)) in jastrow_parameters_1.jpar_map
                @printf fout "[[Pseudopotential]]\n\n"
                @printf fout "index    = %d\n" idx
                @printf fout "v_%d      = %.17g\n\n" idx value
            end

            # write second set of Jastrow parameters
            @printf fout "[JastrowParameters]\n\n"
            @printf fout "type      = \"%s\"\n\n" jastrow_parameters_2.jastrow_type

            for (idx, (indices, value)) in jastrow_parameters_2.jpar_map
                @printf fout "[[Pseudopotential]]\n\n"
                @printf fout "index    = %d\n" idx
                @printf fout "v_%d      = %.17g\n\n" idx value
            end
        end
    end

    return nothing
end

