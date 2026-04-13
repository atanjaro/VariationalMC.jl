@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                           tight_binding_model::TightBindingModel{E2}, 
                           hubbard_model::HubbardModel{E2},
                           model_geometry::ModelGeometry,
                           Np::I,
                           pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Measures the total local energy per site ``E_{\mathrm{loc}}`` for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        tight_binding_model, 
        hubbard_model,
        model_geometry, 
        Np, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["local_energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                           jastrow_factor::JastrowFactor{E2},
                           tight_binding_model::TightBindingModel{E2}, 
                           hubbard_model::HubbardModel{E2},
                           jastrow_parameters::JastrowParameters{S, K, V, I}, 
                           model_geometry::ModelGeometry,
                           Np::I, 
                           pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Measures the total local energy per site ``E_{\mathrm{loc}}`` for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor::JastrowFactor{E2}`: current Jastrow factor. 
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel`: Hubbard interaction parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    jastrow_parameters::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        jastrow_factor, 
        tight_binding_model, 
        hubbard_model,
        jastrow_parameters, 
        model_geometry, 
        Np, 
        pht
    )

    # add the current measured energy to the accumulator
    measurement_container.simulation_measurements["local_energy"] += E_loc_current

    return nothing
end


@doc raw"""

    measure_local_energy!( measurement_container::NamedTuple, 
                           detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                           jastrow_factor_1::JastrowFactor{E2},
                           jastrow_factor_2::JastrowFactor{E2},
                           tight_binding_model::TightBindingModel{E2}, 
                           hubbard_model::HubbardModel{E2},
                           jastrow_parameters_1::JastrowParameters{S, K, V, I},
                           jastrow_parameters_2::JastrowParameters{S, K, V, I},
                           model_geometry::ModelGeometry,
                           Np::I, 
                           pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Measures the total local energy per site ``E_{\mathrm{loc}}`` for a Hubbard model and writes to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor_1::JastrowFactor{E2}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E2}`: second Jastrow factor. 
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_local_energy!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor_1::JastrowFactor{E2},
    jastrow_factor_2::JastrowFactor{E2},
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
   # calculate the current local energy
    E_loc_current = get_local_energy(
        detwf, 
        jastrow_factor_1,
        jastrow_factor_2, 
        tight_binding_model, 
        hubbard_model,
        jastrow_parameters_1,
        jastrow_parameters_2, 
        model_geometry, 
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
) where {T<:Number, Q, E<:Number, I<:Integer}
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

    measure_n!( type::String,
                measurement_container::NamedTuple, 
                detwf::DeterminantalWavefunction{T, Q, E, I}, 
                model_geometry::ModelGeometry,
                pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the local or site-dependent particle density averaged over all sites ``\langle n\rangle``.

- `type::String`: type of measurement. Either `"local"` or `"site-dependent"`.
- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_n!(
    type::String,
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    if type == "local"
        # calculate current density
        density_current = get_n(
            detwf, 
            model_geometry,
            pht
        )

        # add the current measurement to the accumulator
        measurement_container.simulation_measurements["global_density"] += density_current
    elseif type == "site-dependent"
        # calculate current density
        density_current = get_site_dependent_n(
            detwf, 
            model_geometry,
            pht
        )

        # add the current measurement to the accumulator
        measurement_container.simulation_measurements["site-dependent_density"] += density_current
    end
    
    return nothing
end


@doc raw"""

    measure_Sz!( type::String,
                 measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction{T, Q, E, I}, 
                 model_geometry::ModelGeometry,
                 pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the local or site-dependent spin z-component over averaged all sites ``\langle S_z\rangle``.

- `type::String`: type of measurement. Either `"local"` or `"site-dependent"`.
- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_Sz!(
    type::String,
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    if type == "local"
        # calculate current Sz
        Sz_current = get_Sz(
            detwf, 
            model_geometry,
            pht
        )

        # add the current measurement to the accumulator
        measurement_container.simulation_measurements["local_spin-z"] += Sz_current
    elseif type == "site-dependent"
        # calculate current density
        Sz_current = get_site_dependent_s(
            detwf, 
            model_geometry,
            pht
        )

        # add the current measurement to the accumulator
        measurement_container.simulation_measurements["site-dependent_spin-z"] += Sz_current 
    end

    return nothing
end