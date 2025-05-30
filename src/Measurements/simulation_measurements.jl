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

    # get current values from the container
    energy_container = measurement_container.simulation_measurements["energy"]

    # update value for the current bin
    current_E_loc_bin = energy_container[2]
    current_E_loc_bin = E_loc_current

    # update accumuator for this bin
    thisbin_E_loc_sum = energy_container[1]
    thisbin_E_loc_sum += E_loc_current

    # combine the updated values 
    updated_values = (thisbin_E_loc_sum, current_E_loc_bin)

    # write the new values to the container
    measurement_container.simulation_measurements["energy"] = updated_values

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

    # get current values from the container
    energy_container = measurement_container.simulation_measurements["energy"]

    # update value for the current bin
    current_E_loc_bin = energy_container[2]
    current_E_loc_bin = E_loc_current

    # update accumuator for this bin
    thisbin_E_loc_sum = energy_container[1]
    thisbin_E_loc_sum += E_loc_current

    # combine the updated values 
    updated_values = (thisbin_E_loc_sum, current_E_loc_bin)

    # write the new values to the container
    measurement_container.simulation_measurements["energy"] = updated_values

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

    # get current values from the container
    dblocc_container = measurement_container.simulation_measurements["double_occ"]

    # update value for the current bin
    current_dblocc_bin = dblocc_container[2]
    current_dblocc_bin = dblocc_current

    # update accumuator for this bin
    thisbin_dblocc_sum = dblocc_container[1]
    thisbin_dblocc_sum += dblocc_current

    # combine the updated values 
    updated_values = (thisbin_dblocc_sum, current_dblocc_bin)

    # write the new values to the container
    measurement_container.simulation_measurements["double_occ"] = updated_values

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

    # get current values from the container
    density_container = measurement_container.simulation_measurements["density"]

    # update value for the current bin
    current_density_bin = density_container[2]
    current_density_bin = density_current

    # update accumuator for this bin
    thisbin_density_sum = density_container[1]
    thisbin_density_sum += density_current

    # combine the updated values 
    updated_values = (thisbin_density_sum, current_density_bin)

    # write the new values to the container
    measurement_container.simulation_measurements["density"] = updated_values

    return nothing
end