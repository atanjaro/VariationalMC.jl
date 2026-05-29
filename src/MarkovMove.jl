@doc raw"""

    MarkovMove

A type defining a Markov process where a particle is proposed to moved from spindex ``k``
to spindex ``l`` with a Boolean indicatinf whether said move is possible.

- `particle::Int`: Particle ID. 
- 'o::Int`: Orbital ID.
- `k::Int`: Initial spindex of the particle.
- `l::Int`: Final spindex of the particle. 
- `possible::B`: Whether a Markov move is possible.

"""
struct MarkovMove
    # particle
    β::Int

    # orbital
    o::Int

    # initial spindex site
    ksite::Int

    # neighboring spindex site
    lsite::Int

    # move possible
    possible::Bool
end


@doc raw"""

    function MarkovMove(
        pconfig::Vector{Int},
        Np::Int,
        N::Int,
        rng::AbstractRNG
    )

Returns an instance of `MarkovMove`.

"""
function MarkovMove(
    pconfig::Vector{Int},
    nbr_map::AbstractDict,
    unit_cell::UnitCell,
    Np::Int,
    N::Int,
    rng::AbstractRNG
)
    # select a random particle
    β = rand(rng, 1:Np)

    # spindex position of particle β
    ksite = findfirst(==(β), pconfig)

    # get spin of particle β
    spin = get_spindex_type(ksite, N)

    # real site of particle β
    k = get_index_from_spindex(ksite, N)

    # get the orbital ID of the particle
    orb = site_to_orbital(k, unit_cell)

    # choose a random neighbor
    nbr_site = rand(rng, nbr_map[k].neighbors)

    # get spindex of neighboring site with same spin
    lsite = spin == 1 ? get_spindices_from_index(nbr_site, N)[1] :
                    get_spindices_from_index(nbr_site, N)[2]

    @assert get_spindex_type(ksite, N) == get_spindex_type(lsite, N)

    # whether move is possible
    possible = iszero(pconfig[lsite])

    return MarkovMove(β, orb, ksite, lsite, possible)
end


@doc raw"""

    local_fermion_update!(;
        # KEYWORD ARGUMENTS
        detwf::DeterminantalWavefunction{T},
        jastrow_factors::Union{Tuple{<:AbstractJastrowFactor{T}},Nothing} = nothing,
        tight_binding_parameters::TightBindingParameters,
        jastrows_parameters::Union{Tuple{JastrowParameters{T}},Nothing} = nothing,
        particle_configuration::ParticleConfiguration,
        model_geometry::ModelGeometry,
        δW::AbstractFloat,
        n_stab_W::Int,
        rng::AbstractRNG
    ) where {T}

Performs a local update to the fermionic configuration by attempting to move a particle with ID ``\beta`` at a randomly 
selected lattice site ``k`` to another random site ``l``, with acceptance determined by the Metropolis-Hasting algorithm.

# KEYWORD ARGUMENTS

- `detwf::DeterminantalWavefunction{T}`: Current instance of `DeterminantalWavefunction`.
- `jastrow_factors::Union{Tuple{<:AbstractJastrowFactor{T}},Nothing} = nothing`: Tuple of `jastrow_factor`.
- `tight_binding_parameters::TightBindingParameters`: Instance of `TightBindingParameters`.
- `jastrows_parameters::Union{Tuple{JastrowParameters{T}},Nothing} = nothing`: Tuple of `jastrow_parameters`.
- `particle_configuration::ParticleConfiguration`: Instance of `ParticleConfiguration`.
- `model_geometry::ModelGeometry`: Instance of `ModelGeometry`.
- `δW::AbstractFloat`: Maximum allowed error in the equal-time Green's function matrix `W`.
- `δT::Union{AbstractFloat,Nothing} = nothing`: Maximum allowed error in the Jastrow vector `T`.
- `n_stab_W::Int`: Frequency of numerical stability checks performed on the matrix `W`.
- `n_stab_T::Union{Int,Nothing} = nothing`: Frequency of numerical stability checks performed on the vector `T`.
- `rng::AbstractRNG`: Random number generator.

"""
function local_fermion_update!(;
    # KEYWORD ARGUMENTS
    detwf::DeterminantalWavefunction{T},
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}},Nothing} = nothing,
    tight_binding_parameters::TightBindingParameters,
    jas_parameters::Union{Tuple{JastrowParameters{T}},Nothing} = nothing,
    particle_configuration::ParticleConfiguration,
    model_geometry::ModelGeometry,
    δW::AbstractFloat,
    δT::Union{AbstractFloat,Nothing} = nothing,
    n_stab_W::Int,
    n_stab_T::Union{Int,Nothing} = nothing,
    rng::AbstractRNG
) where {T}
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs*Ncells

    (; neighbor_map) = tight_binding_parameters
    (; Np, ph_transform) = particle_configuration

    # choose a particle to move
    markov_move = MarkovMove(particle_configuration.pconfig, neighbor_map, model_geometry.unit_cell, Np, N, rng)

    if markov_move.possible
        # get wavefunction ratio
        Rₛ = real(detwf.W[markov_move.lsite, markov_move.β])

        # calculate acceptance probablility
        acc_prob = Rₛ^2

        # get Jastrow ratios 
        if !isnothing(jas_factors) && !isnothing(jas_parameters)
            for (jasf, jasp) in zip(jas_factors, jas_parameters)
                Rⱼ =  calculate_fermion_jastrow_ratio(
                    jasf,
                    jasp,
                    model_geometry.lattice,
                    markov_move.o,
                    markov_move.ksite,
                    markov_move.lsite,
                    N,
                    ph_transform
                )

                acc_prob *= Rⱼ^2
            end
        end

        if acc_prob > rand(rng)
            # move accepted
            accepted = true

            # update the configuration to reflect the accepted move
            hopping_move!(particle_configuration.pconfig, markov_move.ksite, markov_move.lsite)

            # perform rank-1 update to the Green's function
            detwf.nq_updates_W = update_equal_time_greens_function!(
                detwf.W,
                detwf.M,
                markov_move.β,
                markov_move.lsite,
                Np,
                N,
                detwf.nq_updates_W,
                n_stab_W,
                δW,
                particle_configuration.pconfig
            )

            # update the Jastrow T vector
            if !isnothing(jas_factors) && !isnothing(jas_parameters)
                for (jasf, jasp) in zip(jas_factors, jas_parameters)
                    @assert !isnothing(n_stab_T) && !isnothing(δT)
                    update_fermion_jastrow_factor!(
                        jasf,
                        jasp,
                        model_geometry,
                        markov_move.o,
                        markov_move.ksite,
                        markov_move.lsite,
                        n_stab_T,
                        δT,
                        particle_configuration.pconfig,
                        N,
                        ph_transform
                    )
                end
            end
        else
            # move rejected
            accepted = false
        end 
    else
        # move rejected (not possible)
        accepted = false
    end
    
    return accepted
end


@doc raw"""

    hopping_move!(pconfig::Vector{Int}, ksite::Int, lsite::Int)

If a proposed hopping move is accepted, updates the particle positions in 
the current configuration.

# ARGUMENTS

- `pconfig::Vector{Int}`: Current particle configuration.
- `ksite::Int`: Initial site spindex of the hopping particle.
- `lsite::Int`: Final site spindex of the hopping particle.

"""
function hopping_move!(pconfig::Vector{Int}, ksite::Int, lsite::Int)
    @assert(pconfig[ksite] !== 0)
    @assert(pconfig[lsite] == 0)

    pconfig[lsite] = pconfig[ksite]
    pconfig[ksite] = 0

    return nothing
end























