@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction{T, Q, E, I}, 
                           tight_binding_model::TightBindingModel{E}, 
                           model_geometry::ModelGeometry,
                           U::E,
                           Np::I,
                           pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the total local energy ``E_{\mathrm{loc}}`` for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `U::E`: Hubbard repulsion.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    model_geometry::ModelGeometry, 
    U::E,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        U,
        Np, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["local_energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction{T, Q, E, I}, 
                           tight_binding_model::TightBindingModel{E}, 
                           jastrow_parameters::JastrowParameters{S, K, V, I},
                           jastrow_factor::JastrowFactor{E}, 
                           model_geometry::ModelGeometry,
                           U::E,
                           Np::I, 
                           pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the total local energy ``E_{\mathrm{loc}}`` for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor{E}`: current Jastrow factor. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `U::E`: Hubbard repulsion.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    jastrow_factor::JastrowFactor{E}, 
    model_geometry::ModelGeometry, 
    U::E,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters, 
        jastrow_factor, 
        model_geometry, 
        U,
        Np, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["local_energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction{T, Q, E, I}, 
                           tight_binding_model::TightBindingModel{E}, 
                           jastrow_parameters_1::JastrowParameters{S, K, V, I},
                           jastrow_parameters_2::JastrowParameters{S, K, V, I},
                           jastrow_factor_1::JastrowFactor{E},
                           jastrow_factor_2::JastrowFactor{E}, 
                           model_geometry::ModelGeometry,
                           U::E,
                           Np::I, 
                           pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the total local energy ``E_{\mathrm{loc}}`` for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `jastrow_factor_1::JastrowFactor{E}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E}`: second Jastrow factor. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `U::E`: Hubbard repulsion.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    jastrow_factor_1::JastrowFactor{E},
    jastrow_factor_2::JastrowFactor{E}, 
    model_geometry::ModelGeometry, 
    U::E,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters_1,
        jastrow_parameters_2, 
        jastrow_factor_1,
        jastrow_factor_2, 
        model_geometry, 
        U,
        Np, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["local_energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_double_occ!( measurement_container::NamedTuple, 
                         detwf::DeterminantalWavefunction{T, Q, E, I}, 
                         model_geometry::ModelGeometry, 
                         pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the average double occupancy ``\langle D\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_double_occ!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
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
                detwf::DeterminantalWavefunction{T, Q, E, I}, 
                model_geometry::ModelGeometry ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the local particle density averaged over all sites ``\langle n\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function measure_n!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # calculate current density
    density_current = get_n(
        detwf, 
        model_geometry
    )

    # add the current measurement to the accumulator
    measurement_container.simulation_measurements["global_density"] += density_current

    return nothing
end


@doc raw"""

    measure_Sz!( measurement_container::NamedTuple, 
                detwf::DeterminantalWavefunction{T, Q, E, I}, 
                model_geometry::ModelGeometry ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the local spin z-component over averaged all sites ``\langle S_z\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function measure_Sz!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # calculate current Sz
    Sz_current = get_Sz(
        detwf::DeterminantalWavefunction{T, Q, E, I}, 
        model_geometry::ModelGeometry
    )

    # add the current measurement to the accumulator
    measurement_container.simulation_measurements["spin-z"] += Sz_current

    return nothing
end