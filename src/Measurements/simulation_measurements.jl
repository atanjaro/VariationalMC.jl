@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction, 
                           tight_binding_model::TightBindingModel, 
                           model_geometry::ModelGeometry,
                           Ne::Int64,
                           pht::Bool )::Nothing

Measures the total local energy for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Ne, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction, 
                           tight_binding_model::TightBindingModel, 
                           jastrow_parameters::JastrowParameters,
                           jastrow_factor::JastrowFactor, 
                           model_geometry::ModelGeometry,
                           Ne::Int64, 
                           pht::Bool )::Nothing

Measures the total local energy for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor`: current Jastrow factor. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    jastrow_parameters::JastrowParameters,
    jastrow_factor::JastrowFactor, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters, 
        jastrow_factor, 
        model_geometry, 
        Ne, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_double_occ!( measurement_container::NamedTuple, 
                         detwf::DeterminantalWavefunction, 
                         model_geometry::ModelGeometry, 
                         pht::Bool )::Nothing

Measure the average double occupancy ⟨D⟩ = N⁻¹ ∑ᵢ ⟨nᵢ↑nᵢ↓⟩.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_double_occ!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Nothing
    # calculate the current double occupancy
    dblocc_current = get_double_occ(
        detwf, 
        model_geometry, 
        pht
    )

    # add the current measurement to the accumulator
    measurement_container.simulation_measurements["double_occ"] += dblocc_current

    return nothing
end


@doc raw"""

    measure_n!( measurement_container::NamedTuple, 
                detwf::DeterminantalWavefunction, 
                model_geometry::ModelGeometry )::Nothing

Measure the local particle density ⟨n⟩.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function measure_n!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry
)::Nothing
    # calculate current density
    density_current = get_n(
        detwf, 
        model_geometry
    )

    # add the current measurement to the accumulator
    measurement_container.simulation_measurements["density"] += density_current

    return nothing
end