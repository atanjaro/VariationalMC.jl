@doc raw"""

    get_local_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                      tight_binding_model::TightBindingModel{E2}, 
                      hubbard_model::HubbardModel{E2},
                      model_geometry::ModelGeometry, 
                      Np::I, 
                      pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Calculates the local variational energy per site ``E_{\mathrm{var}}/N`` for a Hubbard model.

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
    # number of lattice sites
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Np, 
        pht
    )

    # calculate on-site Hubbard energy
    E_hubb = get_local_hubbard_energy( 
        detwf, 
        hubbard_model.U₀,
        model_geometry, 
        pht
    )

    # calculate extended Hubbard energy
    if hubbard_model.U₁ != 0.0
        E_ext_hubb = get_extended_hubbard_energy(
            detwf,  
            hubbard_model.U₁,
            model_geometry,
            pht
        )
    else
        E_ext_hubb = 0.0
    end

    # calculate total local energy
    E_loc = E_k + E_hubb + E_ext_hubb
    
    return E_loc / N
end


@doc raw"""

    get_local_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                      jastrow_factor::JastrowFactor{E2},
                      tight_binding_model::TightBindingModel{E2},
                      hubbard_model::HubbardModel{E2},
                      jastrow_parameters::JastrowParameters{S, K, V, I}, 
                      model_geometry::ModelGeometry,
                      Np::I,
                      pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Calculates the local variational energy ``E_{\mathrm{var}}`` per site for a Hubbard model.

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor::JastrowFactor{E2}`: current Jastrow factor.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `hubbard_model::HubbardModel{E2}`: Hubbard interaction parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2},
    tight_binding_model::TightBindingModel{E2}, 
    hubbard_model::HubbardModel{E2},
    jastrow_parameters::JastrowParameters{S, K, V, I},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    # number of lattice sites
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        jastrow_factor, 
        tight_binding_model, 
        jastrow_parameters,
        model_geometry, 
        Np, 
        pht
    )

    # calculate on-site Hubbard energy
    E_hubb = get_local_hubbard_energy( 
        detwf, 
        hubbard_model.U₀,
        model_geometry, 
        pht
    )

    # calculate extended Hubbard energy
    if hubbard_model.U₁ != 0.0
        E_ext_hubb = get_extended_hubbard_energy(
            detwf,
            hubbard_model.U₁,
            model_geometry, 
            pht
        )
    else
        E_ext_hubb = 0.0
    end

    # calculate total local energy
    E_loc = E_k + E_hubb + E_ext_hubb
    
    return E_loc / N
end


@doc raw"""

    get_local_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                      jastrow_factor_1::JastrowFactor{E2},
                      jastrow_factor_2::JastrowFactor{E2},
                      tight_binding_model::TightBindingModel{E2},
                      hubbard_model::HubbardModel{E2},
                      jastrow_parameters_1::JastrowParameters{S, K, V, I},
                      jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
                      model_geometry::ModelGeometry,
                      Np::I,
                      pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Calculates the local variational energy ``E_{\mathrm{var}}`` per site for a Hubbard model.

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
function get_local_energy(
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
    # number of lattice sites
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        jastrow_factor_1, 
        jastrow_factor_2, 
        tight_binding_model, 
        jastrow_parameters_1,
        jastrow_parameters_2,
        model_geometry, 
        Np, 
        pht
    )

    # calculate on-site Hubbard energy
    E_hubb = get_local_hubbard_energy( 
        detwf, 
        hubbard_model.U₀,
        model_geometry, 
        pht
    )

    # calculate extended Hubbard energy
    if hubbard_model.U₁ != 0.0
        E_ext_hubb = get_extended_hubbard_energy( 
            detwf, 
            hubbard_model.U₁,
            model_geometry, 
            pht
        )
    else
        E_ext_hubb = 0.0
    end

    # calculate total local energy
    E_loc = E_k + E_hubb + E_ext_hubb
    
    return E_loc / N
end


@doc raw"""

    get_local_kinetic_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                              tight_binding_model::TightBindingModel{E2}, 
                              model_geometry::ModelGeometry, 
                              Np::I
                              pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Calculates the local kinetic energy ``E_{\mathrm{kin}}`` . 

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    tight_binding_model::TightBindingModel{E2}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
    dims = length(model_geometry.lattice.L)
    Lx   = model_geometry.lattice.L[1]
    if dims != 1
        Ly = model_geometry.lattice.L[2]
    else
        Ly = 1
    end

    # hopping amplitudes       
    t₀ = tight_binding_model.t₀
    t₁ = tight_binding_model.t₁

    # generate nearest neighbor table
    nbr_table0 = build_neighbor_table(
        model_geometry.bond[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # remove double counting if any lattice dimension equals 2
    if Lx == 2 || Ly == 2
        keep = trues(size(nbr_table0, 2))
        seen = Set{Tuple{Int,Int}}()

        for (j, col) in enumerate(eachcol(nbr_table0))
            key = Tuple(sort(col))
            if key in seen
                keep[j] = false
            else
                push!(seen, key)
            end
        end

        nbr_table0 = nbr_table0[:, keep]
    end

    # generate neighbor map
    nbr_map0 = map_neighbor_table(nbr_table0)

    E_loc_kinetic = 0.0

    for β in 1:Np
        # spindex of particle
        k = findfirst(x -> x == β, detwf.pconfig)

        # real position 
        ksite = get_index_from_spindex(k, model_geometry) 

        # check spin of particle  
        spin = get_spindex_type(k, model_geometry)
      
        # loop over nearest neighbors
        sum_nn = 0.0
        for lsite in nbr_map0[ksite][2]
            if spin == 1
                l = get_spindices_from_index(lsite, model_geometry)[1]
            else
                l = get_spindices_from_index(lsite, model_geometry)[2]
            end

            # check that neighboring site is unoccupied
            if detwf.pconfig[l] == 0
                sum_nn += detwf.W[l, β]
            end
        end

        sum_nnn = 0.0
        if t₁ != 0.0
            # generate next-nearest neighbor table
            nbr_table1 = build_neighbor_table(
                model_geometry.bond[2],
                model_geometry.unit_cell,
                model_geometry.lattice
            )

            # generate neighbor map
            nbr_map1 = map_neighbor_table(nbr_table1)

            # loop over next-nearest neighbors
            sum_nnn = 0.0
            for lsite in nbr_map1[ksite][2]
                if spin == 1
                    l = get_spindices_from_index(lsite, model_geometry)[1]
                else
                    l = get_spindices_from_index(lsite, model_geometry)[2]
                end

                # check that neighboring site is unoccupied
                if detwf.pconfig[l] == 0
                    sum_nnn += detwf.W[l, β]
                end
            end
        end

        if pht 
            if spin == 1
                E_loc_kinetic += (-t₀ * sum_nn) + (t₁ * sum_nnn)
            else
                E_loc_kinetic += (t₀ * sum_nn) - (t₁ * sum_nnn)
            end
        else
            E_loc_kinetic += (-t₀ * sum_nn) + (t₁ * sum_nnn)
        end
    end

    return real(E_loc_kinetic)
end


@doc raw"""

    get_local_kinetic_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                              jastrow_factor::JastrowFactor{E2},
                              tight_binding_model::TightBindingModel{E2}, 
                              jastrow_parameters::JastrowParameters{S, K, V, I},
                              model_geometry::ModelGeometry, 
                              Np::I,
                              pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Calculates the local kinetic energy ``E_{\mathrm{kin}}`` . 

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor::JastrowFactor{E2}`: current Jastrow factor.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2},
    tight_binding_model::TightBindingModel{E2}, 
    jastrow_parameters::JastrowParameters{S, K, V, I}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    dims = length(model_geometry.lattice.L)
    Lx   = model_geometry.lattice.L[1]
    if dims != 1
        Ly = model_geometry.lattice.L[2]
    else
        Ly = 1
    end

    # hopping amplitudes       
    t₀ = tight_binding_model.t₀
    t₁ = tight_binding_model.t₁

    # generate neighbor table
    nbr_table0 = build_neighbor_table(
        model_geometry.bond[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # remove double counting if any lattice dimension equals 2
    if Lx == 2 || Ly == 2
        keep = trues(size(nbr_table0, 2))
        seen = Set{Tuple{Int,Int}}()

        for (j, col) in enumerate(eachcol(nbr_table0))
            key = Tuple(sort(col))
            if key in seen
                keep[j] = false
            else
                push!(seen, key)
            end
        end

        nbr_table0 = nbr_table0[:, keep]
    end

    # generate neighbor maps
    nbr_map0 = map_neighbor_table(nbr_table0)

    E_loc_kinetic = 0.0

    for β in 1:Np
        # spindex of particle
        k = findfirst(x -> x == β, detwf.pconfig)

        # real position 
        ksite = get_index_from_spindex(k, model_geometry) 

        # check spin of particle  
        spin = get_spindex_type(k, model_geometry)
      
        # loop over nearest neighbors
        sum_nn = 0.0
        for lsite in nbr_map0[ksite][2]
            if spin == 1
                l = get_spindices_from_index(lsite, model_geometry)[1]
            else
                l = get_spindices_from_index(lsite, model_geometry)[2]
            end

            # check that neighboring site is unoccupied
            if detwf.pconfig[l] == 0
                @assert jastrow_parameters.jastrow_type == "e-den-den" || jastrow_parameters.jastrow_type == "e-spn-spn"

                Rⱼ = get_fermionic_jastrow_ratio(
                    k, 
                    l, 
                    jastrow_parameters, 
                    jastrow_factor, 
                    pht, 
                    spin, 
                    model_geometry
                )

                sum_nn += Rⱼ * detwf.W[l, β]
            end
        end

        sum_nnn = 0.0
        if t₁ != 0.0
            # generate next-nearest neighbor table
            nbr_table1 = build_neighbor_table(
                model_geometry.bond[2],
                model_geometry.unit_cell,
                model_geometry.lattice
            )

            # generate neighbor map
            nbr_map1 = map_neighbor_table(nbr_table1)

            # loop over next nearest neighbors
            sum_nnn = 0.0
            for lsite in nbr_map1[ksite][2]
                if spin == 1
                    l = get_spindices_from_index(lsite, model_geometry)[1]
                else
                    l = get_spindices_from_index(lsite, model_geometry)[2]
                end

                # check that neighboring site is unoccupied
                if detwf.pconfig[l] == 0
                    Rⱼ = get_fermionic_jastrow_ratio(
                        k, 
                        l, 
                        jastrow_parameters, 
                        jastrow_factor, 
                        pht, 
                        spin, 
                        model_geometry
                    )
                    sum_nnn += Rⱼ * detwf.W[l, β]
                end
            end
        end

        if pht 
            if spin == 1
                E_loc_kinetic += (-t₀ * sum_nn) + (t₁ * sum_nnn)
            else
                E_loc_kinetic += (t₀ * sum_nn) - (t₁ * sum_nnn)
            end
        else
            E_loc_kinetic += (-t₀ * sum_nn) + (t₁ * sum_nnn)
        end
    end

    return real(E_loc_kinetic)
end


@doc raw"""

    get_local_kinetic_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                              jastrow_factor_1::JastrowFactor{E2},
                              jastrow_factor_2::JastrowFactor{E2}, 
                              tight_binding_model::TightBindingModel{2E}, 
                              jastrow_parameters_1::JastrowParameters{S, K, V, I},
                              jastrow_parameters_2::JastrowParameters{S, K, V, I},
                              model_geometry::ModelGeometry, 
                              Np::I,
                              pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Calculates the local kinetic energy ``E_{\mathrm{kin}}`` . 

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `jastrow_factor_1::JastrowFactor{E2}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E2}`: second Jastrow factor.
- `tight_binding_model::TightBindingModel{E2}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor_1::JastrowFactor{E2},
    jastrow_factor_2::JastrowFactor{E2},
    tight_binding_model::TightBindingModel{E2}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    dims = length(model_geometry.lattice.L)
    Lx   = model_geometry.lattice.L[1]
    if dims != 1
        Ly = model_geometry.lattice.L[2]
    else
        Ly = 1
    end

    # hopping amplitudes       
    t₀ = tight_binding_model.t₀
    t₁ = tight_binding_model.t₁

    # generate neighbor table
    nbr_table0 = build_neighbor_table(
        model_geometry.bond[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # remove double counting if any lattice dimension equals 2
    if Lx == 2 || Ly == 2
        keep = trues(size(nbr_table0, 2))
        seen = Set{Tuple{Int,Int}}()

        for (j, col) in enumerate(eachcol(nbr_table0))
            key = Tuple(sort(col))
            if key in seen
                keep[j] = false
            else
                push!(seen, key)
            end
        end

        nbr_table0 = nbr_table0[:, keep]
    end

    # generate neighbor maps
    nbr_map0 = map_neighbor_table(nbr_table0)

    E_loc_kinetic = 0.0

    for β in 1:Np
        # spindex of particle
        k = findfirst(x -> x == β, detwf.pconfig)

        # real position 
        ksite = get_index_from_spindex(k, model_geometry) 

        # check spin of particle  
        spin = get_spindex_type(k, model_geometry)
      
        # loop over nearest neighbors
        sum_nn = 0.0
        for lsite in nbr_map0[ksite][2]
            if spin == 1
                l = get_spindices_from_index(lsite, model_geometry)[1]
            else
                l = get_spindices_from_index(lsite, model_geometry)[2]
            end

            # check that neighboring site is unoccupied
            if detwf.pconfig[l] == 0
                @assert jastrow_parameters_1.jastrow_type == "e-den-den" || jastrow_parameters_1.jastrow_type == "e-spn-spn" 
                @assert jastrow_parameters_2.jastrow_type == "e-den-den" || jastrow_parameters_2.jastrow_type == "e-spn-spn"

                Rⱼ₁ = get_fermionic_jastrow_ratio(
                    k, 
                    l, 
                    jastrow_parameters_1, 
                    jastrow_factor_1, 
                    pht, 
                    spin, 
                    model_geometry
                )

                Rⱼ₂ = get_fermionic_jastrow_ratio(
                    k, 
                    l, 
                    jastrow_parameters_2, 
                    jastrow_factor_2, 
                    pht, 
                    spin, 
                    model_geometry
                )
                sum_nn += Rⱼ₁ * Rⱼ₂ * detwf.W[l, β]
            end
        end

        sum_nnn = 0.0
        if t₁ != 0.0
            # generate next-nearest neighbor table
            nbr_table1 = build_neighbor_table(
                model_geometry.bond[2],
                model_geometry.unit_cell,
                model_geometry.lattice
            )

            # generate neighbor map
            nbr_map1 = map_neighbor_table(nbr_table1)

            # loop over next nearest neighbors
            sum_nnn = 0.0
            for lsite in nbr_map1[ksite][2]
                if spin == 1
                    l = get_spindices_from_index(lsite, model_geometry)[1]
                else
                    l = get_spindices_from_index(lsite, model_geometry)[2]
                end

                # check that neighboring site is unoccupied
                if detwf.pconfig[l] == 0
                    Rⱼ₁ = get_fermionic_jastrow_ratio(
                        k, 
                        l, 
                        jastrow_parameters_1, 
                        jastrow_factor_1, 
                        pht, 
                        spin, 
                        model_geometry
                    )

                    Rⱼ₂ = get_fermionic_jastrow_ratio(
                        k, 
                        l, 
                        jastrow_parameters_2, 
                        jastrow_factor_2, 
                        pht, 
                        spin, 
                        model_geometry
                    )

                    sum_nnn += Rⱼ₁ * Rⱼ₂ * detwf.W[l, β]
                end
            end
        end

        if pht 
            if spin == 1
                E_loc_kinetic += (-t₀ * sum_nn) + (t₁ * sum_nnn)
            else
                E_loc_kinetic += (t₀ * sum_nn) - (t₁ * sum_nnn)
            end
        else
            E_loc_kinetic += (-t₀ * sum_nn) + (t₁ * sum_nnn)
        end
    end

    return real(E_loc_kinetic)
end


@doc raw"""

    get_local_hubbard_energy( detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                              U₀::E2,
                              model_geometry::ModelGeometry, 
                              pht::Bool ) where {E1<:AbstractFloat, T<:Number, Q, E2<:Number, I<:Integer}

Calculates the energy due to on-site Hubbard interaction ``U``.  

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `U₀::E2`: on-site Hubbard interaction strength.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_hubbard_energy( 
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    U₀::E2,
    model_geometry::ModelGeometry, 
    pht::Bool
) where {E1<:AbstractFloat, T<:Number, Q, E2<:Number, I<:Integer}
    # number of sites
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells
        
    hubbard_sum = 0.0
    for i in 1:N
        occ_up, occ_dn, _ = get_onsite_fermion_occupation(i, detwf.pconfig, N)
        if pht
            hubbard_sum += occ_up .* (1 .- occ_dn)
        else
            hubbard_sum += occ_up .* occ_dn
        end
    end

    E_loc_hubbard = U₀ * hubbard_sum

    return E_loc_hubbard
end


@doc raw"""

    get_extended_hubbard_energy( detwf::DeterminantalWavefunction{T, Q, E1, I},
                                 U₁::E2,
                                 model_geometry::ModelGeometry,
                                 pht::Bool ) where {E1<:AbstractFloat, T<:Number, Q, E2<:Number, I<:Integer}

Calculates the energy due to the extended Hubbard interaction ``V``.  

- `detwf::DeterminantalWavefunction{T, Q, E1, I}`: current variational wavefunction.
- `U₁::E2`: extended Hubbard interaction strength.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_extended_hubbard_energy(
    detwf::DeterminantalWavefunction{T, Q, E1, I},
    U₁::E2,
    model_geometry::ModelGeometry,
    pht::Bool
) where {E1<:AbstractFloat, T<:Number, Q, E2<:Number, I<:Integer}
    # number of sites
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    # generate neighbor map
    nbr_map0 = map_neighbor_table(nbr_table0)

    ext_hubbard_sum = 0.0
    for i in 1:N
        # get occupations for site i
        occ_upᵢ, occ_dnᵢ, occ_eᵢ = get_onsite_fermion_occupation(i, detwf.pconfig, N)

        # get neighboring occupations for site j
        for j in nbr_map0[i][2]
            occ_upⱼ, occ_dnⱼ, occ_eⱼ = get_onsite_fermion_occupation(j, detwf.pconfig, N)
            
            # add to the energy sum
            if pht
                hubbard_sum += (occ_upᵢ - occ_dnᵢ) * (occ_upⱼ - occ_dnⱼ)
            else
                hubbard_sum += occ_eᵢ * occ_eⱼ
            end
        end
    end

    E_ext_hubbard = U₁ * ext_hubbard_sum

    return E_ext_hubbard
end


@doc raw"""

    get_double_occ( detwf::DeterminantalParameters{T, Q, E, I}, 
                    model_geometry::ModelGeometry, 
                    pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculates the average double occupancy``D = \sum_{i} n_{\boldsymbol{i},\uparrow}n_{\boldsymbol{i},\downarrow}``. 

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_double_occ(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry, 
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    nup_ndn = 0.0
    for site in 1:N
        occ_up, occ_dn, _ = get_onsite_fermion_occupation(site, detwf.pconfig, N)
        if pht
            nup_ndn += occ_up * (1 - occ_dn)
        else
            nup_ndn += occ_up * occ_dn
        end
    end
    
    return nup_ndn / N
end


@doc raw"""

    get_n( detwf::DeterminantalWavefunction{T, Q, E, I}, 
           model_geometry::ModelGeometry,
           pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculate the average particle density ``n = \sum_{i} (n_{\boldsymbol{i},\uparrow} + n_{\boldsymbol{i},\downarrow})``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_n(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    total_occ = 0.0
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    for i in 1:N
        occ_up, occ_dn, _ = get_onsite_fermion_occupation(i, detwf.pconfig, N)
        if pht
            local_occ = occ_up + (1 - occ_dn)
        else
            local_occ = occ_up + occ_dn
        end
        total_occ += local_occ
    end

    return total_occ / N
end


@doc raw"""

    get_Sz( detwf::DeterminantalWavefunction{T, Q, E, I}, 
            model_geometry::ModelGeometry,
            pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculate the average value of the operator 
``S_z = \frac{1}{2} \sum_{i} (n_{\boldsymbol{i},\uparrow} - n_{\boldsymbol{i},\downarrow})``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

""" 
function get_Sz(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    total_Sz = 0.0
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    for i in 1:N
        occ_up, occ_dn, _ = get_onsite_fermion_occupation(i, detwf.pconfig, N)
        if pht
            local_Sz = 0.5 * (occ_up + occ_dn - 1)
        else
            local_Sz = 0.5 * (occ_up - occ_dn)
        end
        total_Sz += local_Sz
    end

    return total_Sz / N
end