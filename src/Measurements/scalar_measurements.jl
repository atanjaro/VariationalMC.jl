@doc raw"""

    get_local_energy( detwf::DeterminantalWavefunction{T, Q, E, I}, 
                      tight_binding_model::TightBindingModel{E}, 
                      model_geometry::ModelGeometry, 
                      Np::I, 
                      pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculates the local variational energy ``E_{\mathrm{var}}`` per site for a Hubbard model.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # number of lattice sites
    N = model_geometry.lattice.N

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Np, 
        pht
    )

    # calculate Hubbard energy
    E_hubb = get_local_hubbard_energy(
        U, 
        detwf, 
        model_geometry, 
        pht
    )

    # calculate total local energy
    E_loc = E_k + E_hubb
    
    return E_loc / N
end


@doc raw"""

    get_local_energy( detwf::DeterminantalWavefunction{T, Q, E, I}, 
                      tight_binding_model::TightBindingModel{E},
                      jastrow_parameters::JastrowParameters{S, K, V, I}, 
                      jastrow_factor::JastrowFactor{E},
                      model_geometry::ModelGeometry
                      Np::I,
                      pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculates the local variational energy ``E_{\mathrm{var}}`` per site for a Hubbard model.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor{E}`: current Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    jastrow_factor::JastrowFactor{E}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # number of lattice sites
    N = model_geometry.lattice.N

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters,
        jastrow_factor, 
        model_geometry, 
        Np, 
        pht
    )

    # calculate Hubbard energy
    E_hubb = get_local_hubbard_energy(U, detwf, model_geometry, pht)

    # calculate total local energy
    E_loc = E_k + E_hubb
    
    return E_loc / N
end


@doc raw"""

    get_local_energy( detwf::DeterminantalWavefunction{T, Q, E, I}, 
                      tight_binding_model::TightBindingModel,
                      jastrow_parameters_1::JastrowParameters{S, K, V, I},
                      jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
                      jastrow_factor_1::JastrowFactor{E},
                      jastrow_factor_2::JastrowFactor{E},
                      model_geometry::ModelGeometry
                      Np::I,
                      pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Calculates the local variational energy ``E_{\mathrm{var}}`` per site for a Hubbard model.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `jastrow_factor_1::JastrowFactor{E}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E}`: second Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I},
    jastrow_factor_1::JastrowFactor{E},
    jastrow_factor_2::JastrowFactor{E}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
    # number of lattice sites
    N = model_geometry.lattice.N

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters_1,
        jastrow_parameters_2,
        jastrow_factor_1, 
        jastrow_factor_2, 
        model_geometry, 
        Np, 
        pht
    )

    # calculate Hubbard energy
    E_hubb = get_local_hubbard_energy(U, detwf, model_geometry, pht)

    # calculate total local energy
    E_loc = E_k + E_hubb
    
    return E_loc / N
end


@doc raw"""

    get_local_kinetic_energy( detwf::DeterminantalWavefunction{T, Q, E, I}, 
                              tight_binding_model::TightBindingModel{E}, 
                              model_geometry::ModelGeometry, 
                              Np::I
                              pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculates the local kinetic energy ``E_{\mathrm{kin}}`` . 

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
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

        # real position position 
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

    get_local_kinetic_energy( detwf::DeterminantalWavefunction{T, Q, E, I}, 
                              tight_binding_model::TightBindingModel{E}, 
                              jastrow_parameters::JastrowParameters{S, K, V, I},
                              jastrow_factor::JastrowFactor{E}, 
                              model_geometry::ModelGeometry, 
                              Np::I,
                              pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Calculates the local kinetic energy ``E_{\mathrm{kin}}`` . 

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor{E}`: current Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    jastrow_parameters::JastrowParameters{S, K, V, I}, 
    jastrow_factor::JastrowFactor{E},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
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

        # real position position 
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

    get_local_kinetic_energy( detwf::DeterminantalWavefunctionT, Q, E, I}, 
                              tight_binding_model::TightBindingModel{E}, 
                              jastrow_parameters_1::JastrowParameters{S, K, V, I},
                              jastrow_parameters_2::JastrowParameters{S, K, V, I},
                              jastrow_factor_1::JastrowFactor{E},
                              jastrow_factor_2::JastrowFactor{E}, 
                              model_geometry::ModelGeometry, 
                              Np::I,
                              pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}

Calculates the local kinetic energy ``E_{\mathrm{kin}}`` . 

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `jastrow_factor_1::JastrowFactor{E}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E}`: second Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    tight_binding_model::TightBindingModel{E}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
    jastrow_factor_1::JastrowFactor{E},
    jastrow_factor_2::JastrowFactor{E},
    model_geometry::ModelGeometry, 
    Np::I, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer, S<:AbstractString, K, V}
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

        # real position position 
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

    get_local_hubbard_energy( U::E, 
                              detwf::DeterminantalWavefunction{T, Q, E, I}, 
                              model_geometry::ModelGeometry, 
                              pht::Bool ) where {E<:AbstractFloat, T<:Number, Q, I<:Integer}

Calculates the energy due to on-site Hubbard repulsion ``U``.  

- `U::E`: Hubbard repulsion.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_hubbard_energy(
    U::E, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry, 
    pht::Bool
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # number of sites
    N = model_geometry.lattice.N
        
    hubbard_sum = 0.0
    for i in 1:N
        occ_up, occ_dn, occ_e = get_onsite_fermion_occupation(i, detwf.pconfig, N)
        if pht
            hubbard_sum += occ_up .* (1 .- occ_dn)
        else
            hubbard_sum += occ_up .* occ_dn
        end
    end

    E_loc_hubbard = U * hubbard_sum

    return E_loc_hubbard
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
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    N = model_geometry.lattice.N

    nup_ndn = 0.0
    for site in 1:N
        occ_up, occ_dn, _ = get_onsite_fermion_occupation(site, detwf.pconfig, N)
        if pht
            nup_ndn += occ_up .* (1 .- occ_dn)
        else
            nup_ndn += occ_up .* occ_dn
        end
    end
    
    return nup_ndn / N
end


@doc raw"""

    get_n( detwf::DeterminantalWavefunction{T, Q, E, I}, 
           model_geometry::ModelGeometry ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculate the average local density ``n = \sum_{i} (n_{\boldsymbol{i},\uparrow} + n_{\boldsymbol{i},\downarrow})``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_n(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    total_occ = 0.0
    N = model_geometry.lattice.N

    for i in 1:N
        local_occ = get_onsite_fermion_occupation(i, detwf.pconfig, N)[1] + 1 - get_onsite_fermion_occupation(i, detwf.pconfig, N)[2]
        total_occ += local_occ
    end

    return total_occ / N
end


@doc raw"""

    get_Sz( detwf::DeterminantalWavefunction{T, Q, E, I}, 
           model_geometry::ModelGeometry ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculate the average value of the operator 
``S_z = \frac{1}{2} \sum_{i} (n_{\boldsymbol{i},\uparrow} - n_{\boldsymbol{i},\downarrow})``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_Sz(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry
) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}
    total_Sz = 0.0
    N = model_geometry.lattice.N

    for i in 1:N
        local_Sz = 0.5 * (get_onsite_fermion_occupation(i, detwf.pconfig, N)[1] - get_onsite_fermion_occupation(i, detwf.pconfig, N)[2])
        total_Sz += local_Sz
    end

    return total_Sz / N
end