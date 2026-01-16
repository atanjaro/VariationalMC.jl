@doc raw"""

    measure_Δk!( measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction{T, Q, E, I}, 
                 determinantal_parameters::DeterminantalParameters{I}, 
                 optimize::NamedTuple,
                 model_geometry::ModelGeometry, 
                 Np::I ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the logarithmic derivative ``\Delta_k`` for ``k`` variational parameters and then writes
them to the measurement container. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system. 

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    Np::I
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # calculate variational parameter derivatives
    Δk_current = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    ) 

    # add current derivative to the accumulator
    measurement_container.optimization_measurements["Δk"] += Δk_current

    return nothing
end


@doc raw"""

    measure_Δk!( measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction{T, Q, E, I}, 
                 determinantal_parameters::DeterminantalParameters{I}, 
                 jastrow_parameters::JastrowParameters{S, K, V, I}, 
                 optimize::NamedTuple,
                 model_geometry::ModelGeometry, 
                 Np::I, 
                 pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the logarithmic derivative ``\Delta_k`` for ``k` variational parameters and then writes
them to the measurement container. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    )  

    # calculate Jastrow parameter derivatives
    Δk_jastrow = get_Δk(
        optimize, 
        jastrow_parameters,
        detwf, 
        model_geometry,
        pht
    ) 

    # combine all derivatives
    Δk_current = vcat(Δk_determinantal, Δk_jastrow)

    # add current derivative to the accumulator
    measurement_container.optimization_measurements["Δk"] += Δk_current

    return nothing
end


@doc raw"""

    measure_Δk!( measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction{T, Q, E, I}, 
                 determinantal_parameters::DeterminantalParameters{I}, 
                 jastrow_parameters_1::JastrowParameters{S, K, V, I},
                 jastrow_parameters_2::JastrowParameters{S, K, V, I},
                 optimize::NamedTuple,
                 model_geometry::ModelGeometry, 
                 Np::I, 
                 pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the logarithmic derivative ``\Delta_k`` for ``k` variational parameters and then writes
them to the measurement container. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters`: second set of Jastrow variational parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    )  

    # calculate Jastrow parameter derivatives
    Δk_jastrow_1 = get_Δk(
        optimize, 
        jastrow_parameters_1,
        detwf, 
        model_geometry,
        pht
    ) 
    Δk_jastrow_2 = get_Δk(
        optimize, 
        jastrow_parameters_2,
        detwf, 
        model_geometry,
        pht
    ) 

    # combine all derivatives
    Δk_current = vcat(Δk_determinantal, Δk_jastrow_1, Δk_jastrow_2)

    # add current derivative to the accumulator
    measurement_container.optimization_measurements["Δk"] += Δk_current

    return nothing
end



@doc raw"""

    measure_ΔkΔkp!( measurement_container::NamedTuple, 
                    detwf::DeterminantalWavefunction{T, Q, E, I}, 
                    determinantal_parameters::DeterminantalParameters{I}, 
                    optimize::NamedTuple,
                    model_geometry::ModelGeometry, 
                    Np::I ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the product of logarithmic derivatives ``\Delta_{k}\Delta_{k}^\prime`` and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current variational determinantal parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::Int`: total number of particles in the system.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    Np::I
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # calculate variational parameter derivatives
    Δk = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    )

    # inner product of Δk and Δk′
    ΔkΔkp_current = Δk .* Δk'

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkΔkp"] += ΔkΔkp_current

    return nothing
end 


@doc raw"""

    measure_ΔkΔkp!( measurement_container::NamedTuple, 
                    detwf::DeterminantalWavefunction{T, Q, E, I}, 
                    determinantal_parameters::DeterminantalParameters{I}, 
                    jastrow_parameters::JastrowParameters{S, K, V, I},
                    optimize::NamedTuple,
                    model_geometry::ModelGeometry, 
                    Np::I, 
                    pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the product of logarithmic derivatives ``\Delta_{k}\Delta_{k}^\prime`` and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow = get_Δk(
        optimize, 
        jastrow_parameters,
        detwf, 
        model_geometry,
        pht
    ) 

    # combine all derivatives
    Δk_current = vcat(Δk_determinantal, Δk_jastrow)

    # inner product of Δk and Δk′
    ΔkΔkp_current = Δk_current .* Δk_current'

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkΔkp"] += ΔkΔkp_current

    return nothing
end 


@doc raw"""

    measure_ΔkΔkp!( measurement_container::NamedTuple, 
                    detwf::DeterminantalWavefunction{T, Q, E, I}, 
                    determinantal_parameters::DeterminantalParameters{I}, 
                    jastrow_parameters_1::JastrowParameters{S, K, V, I},
                    jastrow_parameters_2::JastrowParameters{S, K, V, I},
                    optimize::NamedTuple,
                    model_geometry::ModelGeometry, 
                    Np::I, 
                    pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the product of logarithmic derivatives ``\Delta_{k}\Delta_{k}^\prime`` and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters`: second set of Jastrow variational parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow_1 = get_Δk(
        optimize, 
        jastrow_parameters_1,
        detwf, 
        model_geometry,
        pht
    )
    Δk_jastrow_2 = get_Δk(
        optimize, 
        jastrow_parameters_2,
        detwf, 
        model_geometry,
        pht
    ) 

    # combine all derivatives
    Δk_current = vcat(Δk_determinantal, Δk_jastrow_1, Δk_jastrow_2)

    # inner product of Δk and Δk′
    ΔkΔkp_current = Δk_current .* Δk_current'

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkΔkp"] += ΔkΔkp_current

    return nothing
end 


@doc raw"""

    measure_ΔkE!( measurement_container::NamedTuple, 
                  detwf::DeterminantalWavefunction{T, Q, E, I}, 
                  tight_binding_model::TightBindingModel{E}, 
                  determinantal_parameters::DeterminantalParameters{I}, 
                  optimize::NamedTuple,
                  model_geometry::ModelGeometry, 
                  U::E,
                  Np::I, 
                  pht::Bool )where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the product of logarithmic derivatives with the local energy ``\Delta_{k}E`` and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction. 
- `tight_binding_model::TightBindingModel`: parameter for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `U::E`: Hubbard repulsion.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    U::E,
    Np::I, 
    pht::Bool
)where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # calculate variational parameter derivatives
    Δk = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    )

    # compute local energy
    E_loc = get_local_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        U,
        Np, 
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
                  detwf::DeterminantalWavefunction{T, Q, E, I}, 
                  tight_binding_model::TightBindingModel{E}, 
                  determinantal_parameters::DeterminantalParameters{I}, 
                  jastrow_parameters::JastrowParameters{S, K, V, I}, 
                  jastrow_factor::JastrowFactor{E, I},
                  optimize::NamedTuple,
                  model_geometry::ModelGeometry, 
                  U::E,
                  Np::I, 
                  pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the product of logarithmic derivatives with the local energy ``\Delta_{k}E`` and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor`: current Jastrow factor.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `U::E`: Hubbard repulsion.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I}, 
    jastrow_factor::JastrowFactor{E, I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    U::E,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow = get_Δk(
        optimize, 
        jastrow_parameters,
        detwf,
        model_geometry, 
        pht
    ) 

    # combine all derivatives
    Δk = vcat(Δk_determinantal, Δk_jastrow)

    # compute local energy
    E_loc = get_local_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters,
        jastrow_factor, 
        model_geometry, 
        U,
        Np, 
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
                  detwf::DeterminantalWavefunction{T, Q, E, I}, 
                  tight_binding_model::TightBindingModel{E}, 
                  determinantal_parameters::DeterminantalParameters{I}, 
                  jastrow_parameters_1::JastrowParameters{S, K, V, I},
                  jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
                  jastrow_factor_1::JastrowFactor{E, I},
                  jastrow_factor_2::JastrowFactor{E, I},
                  optimize::NamedTuple,
                  model_geometry::ModelGeometry, 
                  U::E,
                  Np::I, 
                  pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the product of logarithmic derivatives with the local energy ``\Delta_{k}E`` and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters`: second set of Jastrow variational parameters.
- `jastrow_factor_1::JastrowFactor`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor`: second Jastrow factor.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `U::E`: Hubbard repulsion.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
    jastrow_factor_1::JastrowFactor{E, I},
    jastrow_factor_2::JastrowFactor{E, I},
    optimize::NamedTuple,
    model_geometry::ModelGeometry, 
    U::E,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        optimize, 
        determinantal_parameters, 
        detwf, 
        model_geometry, 
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow_1 = get_Δk(
        optimize, 
        jastrow_parameters_1,
        detwf, 
        model_geometry,
        pht
    ) 
    Δk_jastrow_2 = get_Δk(
        optimize, 
        jastrow_parameters_2,
        detwf, 
        model_geometry,
        pht
    ) 

    # combine all derivatives
    Δk = vcat(Δk_determinantal, Δk_jastrow_1, Δk_jastrow_2)

    # compute local energy
    E_loc = get_local_energy(
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

    # compute product of local derivatives with the local energy
    ΔkE_current = Δk * E_loc

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkE"] += ΔkE_current

    return nothing
end