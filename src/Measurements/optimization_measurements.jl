@doc raw"""

    measure_Δk!( measurement_container::NamedTuple, 
                 detwf::DeterminantalWavefunction{T, Q, E, I}, 
                 determinantal_parameters::DeterminantalParameters{I}, 
                 model_geometry::ModelGeometry, 
                 optimize::NamedTuple,
                 Np::I ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the logarithmic derivative ``\Delta_k`` for ``k`` variational parameters and then writes
them to the measurement container. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::I`: total number of particles in the system. 

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I},
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I
) where {T<:Number, Q, E<:Number, I<:Integer}
    # calculate variational parameter derivatives
    Δk_current = get_Δk(
        detwf, 
        determinantal_parameters, 
        model_geometry, 
        optimize,
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
                 model_geometry::ModelGeometry, 
                 optimize::NamedTuple,
                 Np::I, 
                 pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the logarithmic derivative ``\Delta_k`` for ``k` variational parameters and then writes
them to the measurement container. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
        Np
    )  

    # calculate Jastrow parameter derivatives
    Δk_jastrow = get_Δk(
        detwf,  
        jastrow_parameters,
        model_geometry,
        optimize,
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
                 model_geometry::ModelGeometry, 
                 optimize::NamedTuple,
                 Np::I, 
                 pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the logarithmic derivative ``\Delta_k`` for ``k` variational parameters and then writes
them to the measurement container. 

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_Δk!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        detwf, 
        determinantal_parameters, 
        model_geometry, 
        optimize, 
        Np
    )  

    # calculate Jastrow parameter derivatives
    Δk_jastrow_1 = get_Δk(
        detwf, 
        jastrow_parameters_1,
        model_geometry,
        optimize,
        pht
    ) 
    Δk_jastrow_2 = get_Δk(
        detwf, 
        jastrow_parameters_2,
        model_geometry,
        optimize,
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
                    model_geometry::ModelGeometry,
                    optimize::NamedTuple, 
                    Np::I ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the product of logarithmic derivatives ``\Delta_{k}\Delta_{k}^\prime`` and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current variational determinantal parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I
) where {T<:Number, Q, E<:Number, I<:Integer}
    # calculate variational parameter derivatives
    Δk = get_Δk(
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
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
                    model_geometry::ModelGeometry, 
                    optimize::NamedTuple,
                    Np::I, 
                    pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the product of logarithmic derivatives ``\Delta_{k}\Delta_{k}^\prime`` and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        detwf, 
        determinantal_parameters, 
        model_geometry, 
        optimize,
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow = get_Δk(
        detwf, 
        jastrow_parameters,
        model_geometry,
        optimize,
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
                    model_geometry::ModelGeometry, 
                    optimize::NamedTuple,
                    Np::I, 
                    pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Measures the product of logarithmic derivatives ``\Delta_{k}\Delta_{k}^\prime`` and then writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `determinantal_parameters::DeterminantalParameters`: current set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkΔkp!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer, S<:AbstractString, K, V}
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow_1 = get_Δk(
        detwf,
        jastrow_parameters_1, 
        model_geometry,
        optimize,
        pht
    )
    Δk_jastrow_2 = get_Δk(
        detwf, 
        jastrow_parameters_2,
        model_geometry,
        optimize,
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
                  detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                  tight_binding_model::TightBindingModel{E2}, 
                  hubbard_model::HubbardModel{E2},
                  determinantal_parameters::DeterminantalParameters{I}, 
                  model_geometry::ModelGeometry, 
                  optimize::NamedTuple,
                  Np::I, 
                  pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Measures the product of logarithmic derivatives with the local energy ``\Delta_{k}E`` and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction. 
- `tight_binding_model::TightBindingModel{E2}`: parameter for a non-interacting tight-binding model.
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
    measurement_container::NamedTuple, 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    determinantal_parameters::DeterminantalParameters{I}, 
    model_geometry::ModelGeometry, 
    optimize::NamedTuple,
    Np::I, 
    pht::Bool
)where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
    # calculate variational parameter derivatives
    Δk = get_Δk(
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
        Np
    )

    # compute local energy
    E_loc = get_local_energy(
        detwf, 
        tight_binding_model, 
        hubbard_model,
        model_geometry, 
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

Measures the product of logarithmic derivatives with the local energy ``\Delta_{k}E`` and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor::JastrowFactor{E2, I}`: current Jastrow factor.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model.
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
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
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow = get_Δk(
        detwf,
        jastrow_parameters,
        model_geometry, 
        optimize,
        pht
    ) 

    # combine all derivatives
    Δk = vcat(Δk_determinantal, Δk_jastrow)

    # compute local energy
    E_loc = get_local_energy(
        detwf, 
        jastrow_factor,
        tight_binding_model,
        hubbard_model, 
        jastrow_parameters,
        model_geometry, 
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

Measures the product of logarithmic derivatives with the local energy ``\Delta_{k}E`` and then
writes them to the measurement container.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor_1::JastrowFactor{E2, I}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E2, I}`: second Jastrow factor.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model.
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`:: tuple of optimization flags.
- `Np::Int`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function measure_ΔkE!(
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
    # calculate determinantal parameter derivatives
    Δk_determinantal = get_Δk(
        detwf, 
        determinantal_parameters,
        model_geometry, 
        optimize,
        Np
    ) 

    # calculate Jastrow parameter derivatives
    Δk_jastrow_1 = get_Δk(
        detwf, 
        jastrow_parameters_1,
        model_geometry,
        optimize,
        pht
    ) 
    Δk_jastrow_2 = get_Δk(
        detwf, 
        jastrow_parameters_2,
        model_geometry,
        optimize,
        pht
    ) 

    # combine all derivatives
    Δk = vcat(Δk_determinantal, Δk_jastrow_1, Δk_jastrow_2)

    # compute local energy
    E_loc = get_local_energy(
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

    # compute product of local derivatives with the local energy
    ΔkE_current = Δk * E_loc

    # add current values to the accumulator
    measurement_container.optimization_measurements["ΔkE"] += ΔkE_current

    return nothing
end