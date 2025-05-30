@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction, 
                        tight_binding_model::TightBindingModel,
                        determinantal_parameters::DeterminantalParameters,
                        optimize::NamedTuple, 
                        model_geometry::ModelGeometry, 
                        Ne::Int64, 
                        pht::Bool )::Nothing

Measures global, optimization, and simulation observables.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current determinantal wavefunction.
- `tight_binding_model::TightBindingModel`: non-interacting tight-binding model. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: number of electrons.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    determinantal_parameters::DeterminantalParameters, 
    model_geometry::ModelGeometry,
    Ne::Int64, 
    pht::Bool
)::Nothing
    # # measure the variational parameters
    # measure_parameters!(
    #     measurement_container, 
    #     determinantal_parameters
    # )
    
    # measure the energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Ne, 
        pht
    )

    # measure the lograithmic derivatives
    measure_Δk!(
        measurement_container, 
        detwf, 
        determinantal_parameters,
        model_geometry, 
        Ne
    )
    measure_ΔkΔkp!(
        measurement_container, 
        detwf, 
        determinantal_parameters, 
        model_geometry, 
        Ne
    )
    measure_ΔkE!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        determinantal_parameters, 
        model_geometry, 
        Ne, 
        pht
    )

    # measure double occupancy
    measure_double_occ!(
        measurement_container, 
        detwf, 
        model_geometry, 
        pht
    )

    # measure average density
    measure_n!(
        measurement_container, 
        detwf, 
        model_geometry
    )

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction, 
                        tight_binding_model::TightBindingModel,
                        determinantal_parameters::DeterminantalParameters,
                        jastrow_parameters::JastrowParameters,
                        jastrow_factor::JastrowFactor, 
                        model_geometry::ModelGeometry, 
                        Ne::Int64, 
                        pht::Bool )::Nothing

Measures global, optimization, and simulation observables.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current determinantal wavefunction.
- `tight_binding_model::TightBindingModel`: non-interacting tight-binding model. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor`: current Jastrow factor. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: number of electrons.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    jastrow_parameters::JastrowParameters,
    jastrow_factor::JastrowFactor, 
    determinantal_parameters::DeterminantalParameters, 
    model_geometry::ModelGeometry,
    Ne::Int64, 
    pht::Bool
)::Nothing
    # # measure the variational parameters
    # measure_parameters!(
    #     measurement_container, 
    #     determinantal_parameters
    # )
    
    # measure the energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        jastrow_parameters,
        jastrow_factor, 
        model_geometry, 
        Ne, 
        pht
    )

    # measure the lograithmic derivatives
    measure_Δk!(
        measurement_container, 
        detwf,
        determinantal_parameters,
        jastrow_parameters,
        model_geometry,
        Ne,
        pht
    )
    measure_ΔkΔkp!(
        measurement_container, 
        detwf, 
        determinantal_parameters, 
        jastrow_parameters,
        model_geometry, 
        Ne,
        pht
    )
    measure_ΔkE!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        determinantal_parameters, 
        jastrow_parameters,
        model_geometry, 
        Ne, 
        pht
    )

    # measure double occupancy
    measure_double_occ!(
        measurement_container, 
        detwf, 
        model_geometry, 
        pht
    )

    # measure average density
    measure_n!(
        measurement_container, 
        detwf, 
        model_geometry
    )

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction, 
                        tight_binding_model::TightBindingModel, 
                        model_geometry::ModelGeometry, 
                        Ne::Int64, 
                        pht::Bool )

Measures simulation observables.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current determinantal wavefunction.
- `tight_binding_model::TightBindingModel`: non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: number of electrons.
- `pht::Bool`: whether or not model is particle-hoel transformed.

"""
function make_measurements!(measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing    
    # measure the local energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Ne, 
        pht
    )

    # measure double occupancy
    measure_double_occ!(
        measurement_container, 
        detwf, 
        model_geometry, 
        pht
    )

    # measure average density
    measure_n!(
        measurement_container, 
        detwf, 
        model_geometry
    )

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction, 
                        tight_binding_model::TightBindingModel, 
                        jastrow_parameters::JastrowParameters,
                        jastrow_factor::JastrowFactor, 
                        model_geometry::ModelGeometry, 
                        Ne::Int64, 
                        pht::Bool )

Measures simulation observables.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current determinantal wavefunction.
- `tight_binding_model::TightBindingModel`: non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor`: current Jastrow factor. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: number of electrons.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel,
    jastrow_parameters::JastrowParameters,
    jastrow_factor::JastrowFactor,  
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing    
    # measure the local energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        jastrow_parameters::JastrowParameters,
        jastrow_factor::JastrowFactor, 
        model_geometry, 
        Ne, 
        pht
    )

    # measure double occupancy
    measure_double_occ!(
        measurement_container, 
        detwf, 
        model_geometry, 
        pht
    )

    # measure average density
    measure_n!(
        measurement_container, 
        detwf, 
        model_geometry
    )

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


"""

    reset_measurements!( measurements::Dict{String, Any} )

Resets value of a dictionary (measurement container) to zero.

"""
function reset_measurements!(measurements::Dict{String, Any})
    for (key, value) in measurements
        if key == "pconfig"
            # Skip particle configurations
            continue
        elseif isa(value, Tuple)
            # Handle 2- or 3-element tuples
            new_vals = map(1:length(value)) do i
                if i == 1
                    value[1]  # Keep first element (accumulated sum)
                else
                    v = value[i]
                    if isa(v, AbstractVector)
                        empty!(v)
                    elseif isa(v, AbstractMatrix)
                        fill!(v, 0)
                    elseif isa(v, Number)
                        zero(v)
                    elseif isa(v, Any)  # Handle nested arrays like Any[ [a,b], [c,d] ]
                        if eltype(v) <: AbstractVector
                            foreach(empty!, v)
                            v
                        else
                            zero(v)
                        end
                    else
                        # Unknown type in tuple: fallback to zero
                        zero(v)
                    end
                end
            end
            measurements[key] = Tuple(new_vals)
        else
            @warn "Unhandled or unexpected entry at key `$key`: $(typeof(value))"
        end
    end

    return nothing
end