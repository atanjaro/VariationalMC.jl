@doc raw"""

    measure_parameters!( measurement_container::NamedTuple, 
                         determinantal_parameters::DeterminantalParameters )::Nothing

Measures all variational (determinantal) parameters and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.

"""
function measure_parameters!(
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters
)::Nothing
    # record current variational parameters
    parameters_current = collect(values(determinantal_parameters.det_pars))

    # add current parameter values to the container
    measurement_container.optimization_measurements["parameters"] += parameters_current

    # TODO: write final Hamiltonian to file?

    return nothing
end


@doc raw"""

    measure_parameters!( measurement_container::NamedTuple, 
                         determinantal_parameters::DeterminantalParameters, 
                         jastrow_parameters::JastrowParameters )::Nothing

Measures all initialized variational parameters and then writes them to the measurment container. 
The first 'p' are determinantal parameters and the rest are Jastrow parameters. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.

"""
function measure_parameters!(
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters
)::Nothing
    # record current variational parameters
    parameters_current = collect_parameters(determinantal_parameters, jastrow_parameters)

    # add current parameter values to the container
    measurement_container.optimization_measurements["parameters"] += parameters_current

    # TODO: write final Hamiltonian to file?

    return nothing
end


@doc raw"""

    measure_Δk!( measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction, 
                 determinantal_parameters::DeterminantalParameters, 
                 model_geometry::ModelGeometry, 
                 Ne::Int64 )::Nothing

Measures logarithmic derivatives for all variational parameters and then writes
them to the measurement container. The first 'p' are derivatives of determinantal parameters.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    determinantal_parameters::DeterminantalParameters,
    model_geometry::ModelGeometry, 
    Ne::Int64
)::Nothing
    # calculate variational parameter derivatives
    Δk_current = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Ne
    ) 

    # add current derivative to the accumulator
    measurement_container.optimization_measurements["Δk"] += Δk_current

    return nothing
end


@doc raw"""

    measure_Δk!( measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction, 
                 determinantal_parameters::DeterminantalParameters, 
                 jastrow_parameters::JastrowParameters, 
                 model_geometry::ModelGeometry, 
                 Ne::Int64, 
                 pht::Bool )::Nothing

Measures logarithmic derivatives for all variational parameters and then writes 
them to the measurement container. The first 'p' are derivatives of determinantal 
parameters and the rest are derivatives of Jastrow parameters. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
    # calculate variational parameter derivatives
    Δk_current = get_Δk(
        optimize::NamedTuple, 
        determinantal_parameters::DeterminantalParameters, 
        jastrow_parameters::JastrowParameters,
        detwf::DeterminantalWavefunction, 
        model_geometry::ModelGeometry, 
        Ne::Int,
        pht::Bool
    )

    # add current derivative to the accumulator
    measurement_container.optimization_measurements["Δk"] += Δk_current

    return nothing
end


@doc raw"""

    measure_ΔkΔkp!( measurement_container::NamedTuple, 
                    detwf::DeterminantalWavefunction, 
                    determinantal_parameters::DeterminantalParameters,
                    model_geometry::ModelGeometry,
                    Ne::Int64 )::Nothing

Measures the product of variational derivatives with other variational derivatives and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current variational determinantal parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    determinantal_parameters::DeterminantalParameters, 
    model_geometry::ModelGeometry, 
    Ne::Int64
)::Nothing
    # calculate variational parameter derivatives
    Δk = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Ne
    )

    # inner product of Δk and Δk′
    ΔkΔkp_current = Δk .* Δk'

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkΔkp"] += ΔkΔkp_current

    return nothing
end 


@doc raw"""

    measure_ΔkΔkp!( measurement_container::NamedTuple, 
                    detwf::DeterminantalWavefunction, 
                    determinantal_parameters::DeterminantalParameters, 
                    jastrow_parameters::JastrowParameters, 
                    model_geometry::ModelGeometry, 
                    Ne::Int64, 
                    pht::Bool )::Nothing

Measures the product of variational derivatives with other variational derivatives and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
    # calculate variational parameter derivatives
    Δk = get_Δk(
        optimize::NamedTuple, 
        determinantal_parameters::DeterminantalParameters, 
        jastrow_parameters::JastrowParameters,
        detwf::DeterminantalWavefunction, 
        model_geometry::ModelGeometry, 
        Ne::Int,
        pht::Bool
    )

    # inner product of Δk and Δk′
    ΔkΔkp_current = Δk .* Δk'

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkΔkp"] += ΔkΔkp_current

    return nothing
end 


@doc raw"""

    measure_ΔkE!( measurement_container::NamedTuple, 
                  detwf::DeterminantalWavefunction, 
                  tight_binding_model::TightBindingModel, 
                  determinantal_parameters::DeterminantalParameters, 
                  model_geometry::ModelGeometry, 
                  Ne::Int64, 
                  pht::Bool )::Nothing

Measures the product of variational derivatives with the local energy and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction. 
- `tight_binding_model::TightBindingModel`: parameter for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    determinantal_parameters::DeterminantalParameters, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
    # calculate variational parameter derivatives
    Δk = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Ne
    )

    # compute local energy
    E_loc = get_local_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Ne, 
        pht
    ) 

    # compute product of local derivatives with the local energy
    ΔkE_current = Δk * E_loc

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkE"] += ΔkE_current

    return nothing
end


@doc raw"""

    measure_ΔkE!( measurement_container::NamedTuple, 
                  detwf::DeterminantalWavefunction, 
                  tight_binding_model::TightBindingModel, 
                  determinantal_parameters::DeterminantalParameters, 
                  jastrow_parameters::JastrowParameters, 
                  model_geometry::ModelGeometry, 
                  Ne::Int64, 
                  pht::Bool )::Nothing

Measures the product of variational derivatives with the local energy and then 
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
    # calculate variational parameter derivatives
    Δk = get_Δk(
        optimize::NamedTuple, 
        determinantal_parameters::DeterminantalParameters, 
        jastrow_parameters::JastrowParameters,
        detwf::DeterminantalWavefunction, 
        model_geometry::ModelGeometry, 
        Ne::Int,
        pht::Bool
    )

    # compute local energy
    E_loc = get_local_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters,
        jastrow_factor, 
        model_geometry, 
        Ne, 
        pht
    ) 

    # compute product of local derivatives with the local energy
    ΔkE_current = Δk * E_loc

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkE"] += ΔkE_current

    return nothing
end