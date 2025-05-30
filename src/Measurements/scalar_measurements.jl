@doc raw"""
    get_local_energy( detwf::DeterminantalWavefunction, 
                      tight_binding_model::TightBindingModel, 
                      model_geometry::ModelGeometry, 
                      Ne::Int64, 
                      pht::Bool )::Float64

Calculates the local variational energy per site for a Hubbard model (without a Jastrow factor).

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Float64
    # number of lattice sites
    N = model_geometry.lattice.N

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        tight_binding_model, 
        model_geometry, 
        Ne, 
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
    
    return E_loc/N
end


@doc raw"""
    get_local_energy( detwf::DeterminantalWavefunction, 
                      tight_binding_model::TightBindingModel,
                      jastrow_parameters::JastrowParameters, 
                      jastrow_factor::JastrowParameters,
                      model_geometry::ModelGeometry
                      Ne::Int64,
                      pht::Bool )::Nothing

Calculates the local variational energy per site for a Hubbard model.

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor`: current Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_energy(
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    jastrow_parameters::JastrowParameters,
    jastrow_factor::JastrowFactor, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Nothing
    # number of lattice sites
    N = model_geometry.lattice.N

    # calculate kinetic energy
    E_k = get_local_kinetic_energy(
        detwf, 
        tight_binding_model, 
        jastrow_parameters,
        jastrow_factor, 
        model_geometry, 
        Ne, 
        pht
    )

    # calculate Hubbard energy
    E_hubb = get_local_hubbard_energy(U, detwf, model_geometry, pht)

    # calculate total local energy
    E_loc = E_k + E_hubb
    
    return E_loc/N
end


@doc raw"""
    get_local_kinetic_energy( detwf::DeterminantalWavefunction, 
                              tight_binding_model::TightBindingModel, 
                              model_geometry::ModelGeometry, 
                              Ne::Int64
                              pht::Bool )::Float64

Calculates the local electronic kinetic energy. 

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Float64
    # number of sites
    N = model_geometry.lattice.N

    # generate neighbor table
    nbr_table0 = build_neighbor_table(
        model_geometry.bond[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    nbr_table1 = build_neighbor_table(
        model_geometry.bond[2],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # generate neighbor maps
    nbr_map0 = map_neighbor_table(nbr_table0)
    nbr_map1 = map_neighbor_table(nbr_table1)

    E_loc_kinetic = 0.0

    for β in 1:Ne
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
                sum_nn += detwf.W[l, β]
            end
        end

        # hopping amplitudes       
        t₀ = tight_binding_model.t₀
        t₁ = tight_binding_model.t₁

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
    get_local_kinetic_energy( detwf::DeterminantalWavefunction, 
                              tight_binding_model::TightBindingModel, 
                              jastrow_parameters::JastrowParameters,
                              jastrow_factor::JastrowFactor, 
                              model_geometry::ModelGeometry, 
                              Ne::Int64,
                              pht::Bool )::Float64

Calculates the local electronic kinetic energy. 

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor`: current Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_kinetic_energy(
    detwf::DeterminantalWavefunction, 
    tight_binding_model::TightBindingModel, 
    jastrow_parameters::JastrowParameters, 
    jastrow_factor::JastrowFactor,
    model_geometry::ModelGeometry, 
    Ne::Int64, 
    pht::Bool
)::Float64
    # number of sites
    N = model_geometry.lattice.N

    # generate neighbor table
    nbr_table0 = build_neighbor_table(
        model_geometry.bond[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    nbr_table1 = build_neighbor_table(
        model_geometry.bond[2],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # generate neighbor maps
    nbr_map0 = map_neighbor_table(nbr_table0)
    nbr_map1 = map_neighbor_table(nbr_table1)

    E_loc_kinetic = 0.0

    for β in 1:Ne
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
                sum_nn += Rⱼ * detwf.W[l, β]
            end
        end

        # hopping amplitudes       
        t₀ = tight_binding_model.t₀
        t₁ = tight_binding_model.t₁

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
    get_local_hubbard_energy( U::Float64, 
                              detwf::DeterminantalWavefunction, 
                              model_geometry::ModelGeometry, 
                              pht::Bool )::Float64

Calculates the energy due to onsite Hubbard repulsion.  

- `U::Float64`: Hubbard repulsion.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_local_hubbard_energy(
    U::Float64, 
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Float64
    # number of sites
    N = model_geometry.lattice.N
        
    hubbard_sum = 0.0
    for i in 1:N
        occ_up, occ_dn, occ_e = get_onsite_fermion_occupation(i, detwf.pconfig)
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

    get_double_occ( detwf::DeterminantalParameters, 
                    model_geometry::ModelGeometry, 
                    pht::Bool )::Float64

Calculates the double occupancy. 

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_double_occ(
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Float64
    N = model_geometry.lattice.N

    nup_ndn = 0.0
    for site in 1:N
        occ_up, occ_dn, occ_e = get_onsite_fermion_occupation(site, detwf.pconfig)
        if pht
            nup_ndn += occ_up .* (1 .- occ_dn)
        else
            nup_ndn += occ_up .* occ_dn
        end
    end
    
    return nup_ndn / N
end


@doc raw"""

    get_n( detwf::DeterminantalWavefunction, 
           model_geometry::ModelGeometry )::Float64

Calculate the local density.

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_n(
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry
)::Float64
    total_occ = 0.0
    N = model_geometry.lattice.N

    for i in 1:N
        local_occ = get_onsite_fermion_occupation(i, detwf.pconfig)[1] + 1 - get_onsite_fermion_occupation(i, detwf.pconfig)[2]
        total_occ += local_occ
    end

    return total_occ / N
end