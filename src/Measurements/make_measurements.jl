@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        tight_binding_model::TightBindingModel{E2}, 
                        hubbard_model::HubbardModel{E2},
                        determinantal_parameters::DeterminantalParameters{I}, 
                        model_geometry::ModelGeometry,
                        optimize::NamedTuple,
                        Np::I, 
                        pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Measures optimization and simulation observables including the local energy ``\langle E_{\mathrm{loc}}\rangle``, logarithmic
derivatives ``\langle\Delta_k\rangle``, ``\langle\Delta_{k}\Delta_{k}^\prime\rangle``, ``\langle\Delta_{k}E\rangle``, average double occupancy ``\langle D\rangle``, and average density ``\langle n\rangle``.
Also records the current particle configuration ``|x\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current determinantal wavefunction.
- `tight_binding_model::TightBindingModel{E}`: non-interacting tight-binding model. 
- `hubbard_model::HubbardModel`: Hubbard interaction parameters.
- `determinantal_parameters{I}::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    determinantal_parameters::DeterminantalParameters{I}, 
    model_geometry::ModelGeometry,
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
    # measure the energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        hubbard_model,
        model_geometry,
        Np, 
        pht
    )

    # measure the logarithmic derivatives
    measure_Δk!(
        measurement_container, 
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
        Np
    )
    measure_ΔkΔkp!(
        measurement_container, 
        detwf, 
        determinantal_parameters, 
        model_geometry,
        optimize, 
        Np
    )
    measure_ΔkE!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        hubbard_model,
        determinantal_parameters, 
        model_geometry, 
        optimize,
        Np, 
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
        "local",
        measurement_container, 
        detwf, 
        model_geometry,
        pht
    )

    # measure site_dependent density (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        measure_n!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure average Sz (if added)
    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        measure_Sz!(
            "local",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure site_dependent Sz (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        measure_Sz!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end
    
    # measure correlations (if added)
    if haskey(measurement_container.correlation_measurements, "density")
        measure_density_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end
    if haskey(measurement_container.correlation_measurements, "spin")
        measure_spin_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        jastrow_factor::JastrowFactor{E2, I},
                        tight_binding_model::TightBindingModel{E2},
                        hubbard_model::HubbardModel{E2},
                        determinantal_parameters::DeterminantalParameters{I},
                        jastrow_parameters::JastrowParameters{S, K, V, I},
                        model_geometry::ModelGeometry, 
                        optimize::NamedTuple,
                        Np::I, 
                        pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Measures optimization and simulation observables including the local energy ``\langle E_{\mathrm{loc}}\rangle``, logarithmic
derivatives ``\langle\Delta_k\rangle``, ``\langle\Delta_{k}\Delta_{k}^\prime\rangle``, ``\langle\Delta_{k}E\rangle``, average double occupancy ``\langle D\rangle``, and average density ``\langle n\rangle``.
Also records the current particle configuration ``|x\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current determinantal wavefunction.
- `jastrow_factor::JastrowFactor{E2, I}`: current Jastrow factor. 
- `tight_binding_model::TightBindingModel{E}`: non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    determinantal_parameters::DeterminantalParameters{I},
    jastrow_parameters::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry,
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    # measure the energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        jastrow_factor,
        tight_binding_model, 
        hubbard_model,
        jastrow_parameters,
        model_geometry, 
        Np, 
        pht
    )

    # measure the logarithmic derivatives
    measure_Δk!(
        measurement_container, 
        detwf,
        determinantal_parameters,
        jastrow_parameters,
        model_geometry,
        optimize,
        Np,
        pht
    )
    measure_ΔkΔkp!(
        measurement_container, 
        detwf, 
        determinantal_parameters, 
        jastrow_parameters,
        model_geometry, 
        optimize,
        Np,
        pht
    )
    measure_ΔkE!(
        measurement_container, 
        detwf, 
        jastrow_factor,
        tight_binding_model, 
        hubbard_model,
        determinantal_parameters, 
        jastrow_parameters,
        model_geometry, 
        optimize,
        Np, 
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
        "local",
        measurement_container, 
        detwf, 
        model_geometry,
        pht
    )

    # measure site_dependent density (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        measure_n!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure average Sz (if added)
    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        measure_Sz!(
            "local",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure site_dependent Sz (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        measure_Sz!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end
    
    # measure correlations (if added)
    if haskey(measurement_container.correlation_measurements, "density")
        measure_density_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end
    if haskey(measurement_container.correlation_measurements, "spin")
        measure_spin_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        jastrow_factor_1::JastrowFactor{E2, I},
                        jastrow_factor_2::JastrowFactor{E2, I},
                        tight_binding_model::TightBindingModel{E2},
                        hubbard_model::HubbardModel{E2},
                        determinantal_parameters::DeterminantalParameters{I},
                        jastrow_parameters_1::JastrowParameters{S, K, V, I},
                        jastrow_parameters_2::JastrowParameters{S, K, V, I},
                        model_geometry::ModelGeometry, 
                        optimize::NamedTuple,
                        Np::I, 
                        pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Measures optimization and simulation observables including the local energy ``\langle E_{\mathrm{loc}}\rangle``, logarithmic
derivatives ``\langle\Delta_k\rangle``, ``\langle\Delta_{k}\Delta_{k}^\prime\rangle``, ``\langle\Delta_{k}E\rangle``, average double occupancy ``\langle D\rangle``, and average density ``\langle n\rangle``.
Also records the current particle configuration ``|x\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current determinantal wavefunction.
- `jastrow_factor_1::JastrowFactor{E, I}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E, I}`: second Jastrow factor.
- `tight_binding_model::TightBindingModel{E}`: non-interacting tight-binding model. 
- `hubbard_model::HubbardModel`: Hubbard interaction parameters.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags. 
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor_1::JastrowFactor{E2, I},
    jastrow_factor_2::JastrowFactor{E2, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    determinantal_parameters::DeterminantalParameters{I},
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry,
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    # measure the energy
    measure_local_energy!(
        measurement_container, 
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

    # measure the lograithmic derivatives
    measure_Δk!(
        measurement_container, 
        detwf,
        determinantal_parameters,
        jastrow_parameters_1,
        jastrow_parameters_2,
        model_geometry,
        optimize,
        Np,
        pht
    )
    measure_ΔkΔkp!(
        measurement_container, 
        detwf, 
        determinantal_parameters, 
        jastrow_parameters_1,
        jastrow_parameters_2,
        model_geometry, 
        optimize,
        Np,
        pht
    )
    measure_ΔkE!(
        measurement_container, 
        detwf, 
        jastrow_factor_1,
        jastrow_factor_2,
        tight_binding_model, 
        hubbard_model,
        determinantal_parameters, 
        jastrow_parameters_1,
        jastrow_parameters_2,
        model_geometry, 
        optimize,
        Np, 
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
        "local",
        measurement_container, 
        detwf, 
        model_geometry,
        pht
    )

    # measure site_dependent density (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        measure_n!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure average Sz (if added)
    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        measure_Sz!(
            "local",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure site_dependent Sz (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        measure_Sz!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end
    
    # measure correlations (if added)
    if haskey(measurement_container.correlation_measurements, "density")
        measure_density_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end
    if haskey(measurement_container.correlation_measurements, "spin")
        measure_spin_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        tight_binding_model::TightBindingModel{E2}, 
                        hubbard_model::HubbardModel{E2},
                        model_geometry::ModelGeometry, 
                        Np::I, 
                        pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Measures simulation observables including the local energy ``\langle E_{\mathrm{loc}}\rangle``,
average double occupancy ``\langle D\rangle``, and average density ``\langle n\rangle``.
Also records the current particle configuration ``|x\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current determinantal wavefunction.
- `tight_binding_model::TightBindingModel{E}`: non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
    # measure the local energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        tight_binding_model, 
        hubbard_model,
        model_geometry, 
        Np, 
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
        "local",
        measurement_container, 
        detwf, 
        model_geometry,
        pht
    )

    # measure site_dependent density (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        measure_n!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure average Sz (if added)
    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        measure_Sz!(
            "local",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure site_dependent Sz (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        measure_Sz!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end
    
    # measure correlations (if added)
    if haskey(measurement_container.correlation_measurements, "density")
        measure_density_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end
    if haskey(measurement_container.correlation_measurements, "spin")
        measure_spin_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
                        detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        jastrow_factor::JastrowFactor{E2},
                        tight_binding_model::TightBindingModel{E2}, 
                        hubbard_model::HubbardModel{E2},
                        jastrow_parameters::JastrowParameters{S, K, V, I},
                        model_geometry::ModelGeometry, 
                        Np::I, 
                        pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Measures simulation observables including the local energy ``\langle E_{\mathrm{loc}}\rangle``,
average double occupancy ``\langle D\rangle``, and average density ``\langle n\rangle``.
Also records the current particle configuration ``|x\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current determinantal wavefunction.
- `jastrow_factor::JastrowFactor{E2}`: current Jastrow factor. 
- `tight_binding_model::TightBindingModel{E2}`: non-interacting tight-binding model. 
- `hubbard_model::HubbardMode{E2}`: Hubbard interaction parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2},
    tight_binding_model::TightBindingModel{E2},
    hubbard_model::HubbardModel{E2},
    jastrow_parameters::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    # measure the local energy
    measure_local_energy!(
        measurement_container, 
        detwf, 
        jastrow_factor, 
        tight_binding_model, 
        hubbard_model,
        jastrow_parameters,
        model_geometry, 
        Np, 
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
        "local",
        measurement_container, 
        detwf, 
        model_geometry,
        pht
    )

    # measure site_dependent density (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        measure_n!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure average Sz (if added)
    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        measure_Sz!(
            "local",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure site_dependent Sz (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        measure_Sz!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end
    
    # measure correlations (if added)
    if haskey(measurement_container.correlation_measurements, "density")
        measure_density_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end
    if haskey(measurement_container.correlation_measurements, "spin")
        measure_spin_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


@doc raw"""

    make_measurements!( measurement_container::NamedTuple, 
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

Measures simulation observables including the local energy ``\langle E_{\mathrm{loc}}\rangle``,
average double occupancy ``\langle D\rangle``, and average density ``\langle n\rangle``.
Also records the current particle configuration ``|x\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current determinantal wavefunction.
- `jastrow_factor_1::JastrowFactor{E2}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E2}`: second Jastrow factor. 
- `tight_binding_model::TightBindingModel{E2}`: non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether or not model is particle-hole transformed.

"""
function make_measurements!(measurement_container::NamedTuple, 
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
    # measure the local energy
    measure_local_energy!(
        measurement_container, 
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

    # measure double occupancy
    measure_double_occ!(
        measurement_container, 
        detwf, 
        model_geometry, 
        pht
    )

    # measure average density
    measure_n!(
        "local",
        measurement_container, 
        detwf, 
        model_geometry,
        pht
    )

    # measure site_dependent density (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        measure_n!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure average Sz (if added)
    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        measure_Sz!(
            "local",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end

    # measure site_dependent Sz (if added)
    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        measure_Sz!(
            "site-dependent",
            measurement_container, 
            detwf, 
            model_geometry,
            pht
        ) 
    end
    
    # measure correlations (if added)
    if haskey(measurement_container.correlation_measurements, "density")
        measure_density_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end
    if haskey(measurement_container.correlation_measurements, "spin")
        measure_spin_correlation!(
            measurement_container,
            detwf,
            model_geometry,
            pht
        )
    end

    # record the current particle configuration
    measurement_container.simulation_measurements["pconfig"] = detwf.pconfig

    return nothing
end


"""

    reset_measurements!( measurements::Dict{S, T}; 
                         in_place::Bool = false ) where {S<:AbstractString, T}

Resets values in the measurement container to zero.

- `measurements::Dict{S, T}`: measurement container.
- `in_place::Bool = false`: whether to perform in-place update of measurements. Best to set to `false` to prevent undesired mutation.

"""
function reset_measurements!(
    measurements::Dict{S, T}; 
    in_place::Bool = false
) where {S<:AbstractString, T}
    function reset_value(val)
        # Numbers & complex -> scalar zero
        if isa(val, Number) || isa(val, Complex)
            return zero(val)

        # 1D arrays / vectors
        elseif isa(val, AbstractVector)
            if in_place
                fill!(val, zero(eltype(val)))
                return val
            else
                return zeros(eltype(val), length(val))
            end

        # Matrices / N-d arrays
        elseif isa(val, AbstractMatrix)
            if in_place
                fill!(val, zero(eltype(val)))
                return val
            else
                return zeros(eltype(val), size(val)...)
            end

        # General <:AbstractArray (including Vector{Vector} etc.)
        elseif isa(val, AbstractArray)
            if in_place
                for i in eachindex(val)
                    val[i] = reset_value(val[i])
                end
                return val
            else
                return [ reset_value(x) for x in val ]
            end

        # Tuples: immutable -> build new tuple with reset elements
        elseif isa(val, Tuple)
            return tuple((reset_value(v) for v in val)...)

        # Dict: recurse and return a new dict of reset values
        elseif isa(val, Dict)
            newd = Dict{Any,Any}()
            for (k,v) in val
                newd[k] = reset_value(v)
            end
            return newd

        # fallback: unknown type -> return as-is
        else
            return val
        end
    end

    for (k, v) in measurements
        if k == "parameters"
            continue  # skip this key
        end
        measurements[k] = reset_value(v)
    end

    return nothing
end