@doc raw"""

    abstract type AbstractHST{T<:Number, R<:AbstractFloat} end

Abstract type to represent a Jastrow factor.

"""
abstract type AbstractJastrowFactor{E<:AbstractFloat} end



@doc raw"""

    FermionJastrowFactor{E<:AbstractFloat} <: AbstractJastrowFactor{E}

    BosonJastrowFactor{E<:AbstractFloat} <: AbstractJastrowFactor{E}

    MixedJastrowFactor{E<:AbstractFloat} <: AbstractJastrowFactor{E}

A mutable struct containing all the parameters needed to characterize a Jastrow factor involving fermions and/or bosons (phonons).

# Fields

- `T::Vector{Vector{E}}`: Vector of Jastrow factor elements.
- `orbitals::Vector{Int}`: Orbital IDs subject to Jastrow pseudopotentials.
- `nq_updates_T::Int`: Tracker for the numbner of quick updates that have been performed to vector `T`.

"""
mutable struct FermionJastrowFactor{E<:AbstractFloat} <: AbstractJastrowFactor{E}
    # fermion T vector
    T::Vector{Vector{E}}

    # orbital IDS
    orbitals::Vector{Int}

    # number of T vector quick updates
    nq_updates_T::Int
end 


mutable struct BosonJastrowFactor{E<:AbstractFloat} <: AbstractJastrowFactor{E}
    # boson (phonon) T vector
    T::Vector{Vector{E}}

    # orbital IDS
    orbitals::Vector{Int}

    # number of T vector quick updates
    nq_updates_T::Int
end 


mutable struct MixedJastrowFactor{E<:AbstractFloat} <: AbstractJastrowFactor{E}
    # mixed (electron-phonon) T vector
    T::Union{Vector{Vector{E}}, Nothing}

    # orbital IDS
    orbitals::Vector{Int}

    # number of T vector quick updates
    nq_updates_T::Int
end 


@doc raw"""

    JastrowFactor(;
        # KEYWORD ARGUMENTS
        jastrow_parameters::JastrowParameters{E},
        particle_configuration::ParticleConfiguration,
        model_geometry::ModelGeometry
    ) where {E}

Initialize and return an instance of the `JastrowFactor` type depending upon the particle pair specified in `jastrow_parameters`.

"""
function JastrowFactor(;
    # KEYWORD ARGUMENTS
    jastrow_parameters::JastrowParameters{E},
    particle_configuration::ParticleConfiguration,
    model_geometry::ModelGeometry
) where {E}
    N_orbs = model_geometry.unit_cell.n
    N_cells = model_geometry.lattice.N
    N = N_orbs * N_cells

    (; particle_pair, order_pair, irr_indices, 
        irr_index_map, mean_v, orbitals) = jastrow_parameters
    (; pconfig, ph_transform) = particle_configuration

    # intialize tracker for quick updating of Jastrow T vector
    nq_updates_T = 0

    Ts = Vector{E}[]

    if particle_pair == "electron-electron"
        for o in orbitals
            fermion_T = zeros(E, N)
            
            for idx in eachindex(irr_indices[o])
                # get vₖ
                vₖ = mean_v[o][idx]
                # iterate over pairs for each irreducible index
                for (i,j) in irr_index_map[o][irr_indices[o][idx]]
                    build_fermion_T!(fermion_T, i+1, j+1, order_pair, vₖ, pconfig, N, ph_transform)
                end
            end

            push!(Ts, fermion_T)
        end

        return FermionJastrowFactor(Ts, orbitals, nq_updates_T)
    elseif particle_pair == "phonon-phonon"
        # TODO
        # return BosonJastrowFactor(Ts, unique(orbitals), nq_updates_T)
    elseif particle_pair == "electron-phonon"
        # TODO
        # return MixedJastrowFactor(Ts, unique(orbitals), nq_updates_T)
    end
end


@doc raw"""

    build_fermion_T!(
        T::Vector{E},
        i::Int,
        j::Int,
        order_pair::AbstractString,
        vₖ::E,
        pconfig::Vector{Int},
        N::Int,
        ph_transform::Bool
    ) where {E<:AbstractFloat}

Constructs the element of the fermion Jastrow vector `T` between lattice sites `i` and `j`.

"""
function build_fermion_T!(
    T::Vector{E},
    i::Int,
    j::Int,
    order_pair::AbstractString,
    vₖ::E,
    pconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
) where {E<:AbstractFloat}
    num_up, num_dn = get_fermion_occupations(j, pconfig, N)
    if order_pair == "density-density"
        T[i] += ph_transform ? vₖ * (num_up - num_dn) : vₖ * (num_up + num_dn)
    elseif order_pair == "spin-spin"
        T[i] += ph_transform ? 0.5 * vₖ * (num_up - num_dn) : 0.5 * vₖ * (num_up + num_dn)
    end 
    
    return nothing
end


@doc raw"""

    build_boson_T!(
        T::Vector{E},
        i::Int,
        j::Int,
        order_pair::AbstractString,
        vₖ::E,
        phconfig::Vector{Int},
        N::Int
    ) where {E<:AbstractFloat}

Constructs the element of the boson (phonon) Jastrow vector `T` between lattice sites `i` and `j`.

"""
function build_boson_T!(
    T::Vector{E},
    i::Int,
    j::Int,
    order_pair::AbstractString,
    vₖ::E,
    phconfig::Vector{Int},
    N::Int
) where {E<:AbstractFloat}
    if order_pair == "density-density"
        num_phs_j = get_boson_occupations(j, phconfig, N)
        T[i] += vₖ *  num_phs_j
    elseif order_pair == "displacement-displacement"
        # TODO 
    end 
    
    return nothing
end


@doc raw"""

    build_mixed_T!(
        T::Vector{E},
        i::Int,
        j::Int,
        order_pair::AbstractString,
        vₖ::E,
        pconfig::Vector{Int},
        phconfig::Vector{Int},
        N::Int
    ) where {E<:AbstractFloat}

Constructs the element of the electron-boson Jastrow vector `T` between lattice sites `i` and `j`.

"""
function build_mixed_T!(
    T::Vector{E},
    i::Int,
    j::Int,
    order_pair::AbstractString,
    vₖ::E,
    pconfig::Vector{Int},
    phconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
) where {E<:AbstractFloat}
    if order_pair == "density-density"
        num_up, num_dn = get_fermion_occupations(j, pconfig, N)
        num_phs_j = get_boson_occupations(j, phconfig, N)
        T[i] += ph_transform ? vₖ * (num_up - num_dn) * num_phs_j : vₖ * (num_up + num_dn) * num_phs_j
    elseif order_pair == "density-displacement"
        # TODO 
    end 
    
    return nothing
end


@doc raw"""

    calculate_fermion_jastrow_ratio(
        jfac::FermionJastrowFactor{E},
        jastrow_parameters::JastrowParameters{E},
        lattice::Lattice,
        o::Int,
        ksite::Int,
        lsite::Int,
        N::Int,
        ph_transform::Bool
    ) where {E}

Calculates ratio of fermionic Jastrow factors
```math
\frac{\mathcal{J}_1}{\mathcal{J}_2} = \exp[-s(T_{l} - T_{k}) + v_{ll} - v_{lk}],
```
after a particle `β` moves from site `k` to site `l`.

"""
function calculate_fermion_jastrow_ratio(
    jfac::FermionJastrowFactor{E},
    jastrow_parameters::JastrowParameters{E},
    lattice::Lattice,
    o::Int,
    ksite::Int,
    lsite::Int,
    N::Int,
    ph_transform::Bool
) where {E}
    (; T) = jfac
    (; irr_indices, mean_v) = jastrow_parameters

    # convert spindex to site index
    k = get_index_from_spindex(ksite, N)
    l = get_index_from_spindex(lsite, N)

    Tₖ = T[o][k]
    Tₗ = T[o][l]

    # convert site indices to irreducible index
    kl_idx = reduce_index(k-1, l-1, lattice)
    ll_idx = reduce_index(l-1, l-1, lattice)

    # lookup the correct index
    ll_pos = findfirst(==(ll_idx), irr_indices[o])
    kl_pos = findfirst(==(kl_idx), irr_indices[o])

    @assert kl_pos !== nothing "Invalid irreducible index $kl_idx for orbital $o"
    @assert ll_pos !== nothing "Invalid irreducible index $ll_idx for orbital $o"

    # get pseudopotential values for orbital o
    vₗₗ = mean_v[o][ll_pos]
    vₗₖ = mean_v[o][kl_pos]

    # check spin of β
    spin = ph_transform ? get_spindex_type(ksite, N) : 1.0

    # compute ratio
    jas_ratio = exp(spin * (Tₗ - Tₖ) + vₗₗ - vₗₖ)
    
    return jas_ratio
end


@doc raw"""

    update_fermion_jastrow_factor!(
        jfac::FermionJastrowFactor{E},
        jastrow_parameters::JastrowParameters{E},
        lattice::Lattice,
        o::Int,
        ksite::Int,
        lsite::Int,
        n_stab_T::Int,
        δT::AbstractFloat,
        pconfig::Vector{Int},
        N::Int,
        ph_transform::Bool
    ) where {E}

Updates the elements of the fermion Jastrow vector `T`.

"""
function update_fermion_jastrow_factor!(
    jfac::FermionJastrowFactor{E},
    jastrow_parameters::JastrowParameters{E},
    model_geometry::ModelGeometry,
    o::Int,
    ksite::Int,
    lsite::Int,
    n_stab_T::Int,
    δT::AbstractFloat,
    pconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
) where {E}
    (; particle_pair, order_pair, 
        irr_indices, irr_index_map, mean_v) = jastrow_parameters
    (; unit_cell, lattice) = model_geometry

    # convert spindex to site index
    k = get_index_from_spindex(ksite, N)
    l = get_index_from_spindex(lsite, N)

    # new elements for the T vector
    diff_fermion_T = zeros(E, N)

    # check spin of β once — loop invariant
    spin = ph_transform ? E(get_spindex_type(ksite, N)) : one(E)

    for i in 1:N
        # check that the site lives on the correct sublattice
        site_to_orbital(i, unit_cell) == o || continue
         
        # convert site indices to irreducible index
        ik_idx = reduce_index(i-1, k-1, lattice)
        il_idx = reduce_index(i-1, l-1, lattice)

        ik_pos = findfirst(==(ik_idx), irr_indices[o])
        il_pos = findfirst(==(il_idx), irr_indices[o])

        @assert ik_pos !== nothing "Invalid irreducible index $ik_idx for orbital $o"
        @assert il_pos !== nothing "Invalid irreducible index $il_idx for orbital $o"

        # get pseudopotential values for orbital o
        vᵢₖ = mean_v[o][ik_pos]
        vᵢₗ = mean_v[o][il_pos]

        # index into diff_fermion_T relative to start of orbital o's sites
        diff_fermion_T[i] += spin * (vᵢₗ - vᵢₖ)
    end

    # update the T vector in-place
    jfac.T[o] .+= diff_fermion_T

    if jfac.nq_updates_T >= n_stab_T
        jfac.nq_updates_T = 0

        # recalculate T vector — mirrors constructor algorithm exactly
        @assert particle_pair == "electron-electron"
        Tᵣ = zeros(E, N)
        for idx in eachindex(irr_indices[o])
            # get vₖ
            vₖ = mean_v[o][idx]
            # iterate over pairs for each irreducible index
            for (i,j) in irr_index_map[o][irr_indices[o][idx]]
                build_fermion_T!(Tᵣ, i+1, j+1, order_pair, vₖ, pconfig, N, ph_transform)
            end
        end

        # compute deviation and replace if necessary
        dev = check_deviation(jfac.T[o], Tᵣ)
        if dev > δT
            copyto!(jfac.T[o], Tᵣ)
        end
    else
        jfac.nq_updates_T += 1
    end

    return nothing
end


@doc raw"""

    check_deviation(
        T::Vector{E}, 
        Tᵣ::Vector{E}
    ) where {E<:AbstractFloat}

Checks floating point error accumulation in the fermionic T vector.

"""
function check_deviation(
    T::Vector{E}, 
    Tᵣ::Vector{E}
) where {E<:AbstractFloat}
    @assert size(T) == size(Tᵣ)

    exact_square_sum = sum(abs2, T)
    diff_square_sum  = sum(abs2, T .- Tᵣ)

    return exact_square_sum == 0.0 ?
           sqrt(diff_square_sum) :
           sqrt(diff_square_sum / exact_square_sum)
end