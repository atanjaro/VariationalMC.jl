@doc raw"""

    MarkovMove{I<:Integer}

A type defining a Markov process where a particle is proposed to moved from spindex ``k``
to spindex ``l`` and a Boolean denoting whether said move is possible.

- `particle::I`: the particle undergoing a Markov process.
- `k::I`: initial spindex of the particle.
- `l::I`: final spindex of the particle. 
- `possible::B`: whether the Markov process is possible.

"""
struct MarkovMove{I<:Integer}
    # particle
    particle::I

    # initial spindex
    k::I

    # neighboring spindex
    l::I

    # move possible
    possible::Bool
end


@doc raw"""

    propose_random_move( Np::I, 
                         pconfig::AbstractVector{I}, 
                         model_geometry::ModelGeometry, 
                         rng::AbstractRNG ) where {I<:Integer}

Proposes randomly moving (via hopping or exchange) a particle from some intial spindex ``k`` 
to a neighboring spindex ``l`` and returns an instance of `MarkovMove`.

- `Np::I`: total number of particles in the system.  
- `pconfig::AbstractVector{I}`: current particle configuration. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `rng::AbstractRNG`: random number generator.

"""
function propose_random_move(
    Np::I, 
    pconfig::Vector{I}, 
    model_geometry::ModelGeometry, 
    rng::AbstractRNG
) where {I<:Integer}
    # create nearest neighbor table
    nbr_table = build_neighbor_table(
        model_geometry.bond[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # map nearest neighbor table to dictionary of bonds and neighbors                                
    nbr_map = map_neighbor_table(nbr_table)

    # # checks for next nearest neighbors
    # if length(bonds) == 2
    #     nbr_table_p = build_neighbor_table(model_geometry.bond[2],
    #                                       model_geometry.unit_cell,
    #                                       model_geometry.lattice)
    #     nbr_map_p = map_neighbor_table(nbr_table_p)
    # end

    # select a random particle
    β = rand(rng, 1:Np)

    # spindex position of particle β
    k = findfirst(x -> x == β, pconfig)

    # get spin of β
    spin = get_spindex_type(k, model_geometry)

    # real site of particle β
    ksite = get_index_from_spindex(k, model_geometry)

    # choose random neighboring site
    nbr_site = rand(rng, nbr_map[ksite][2]) 

    # get spindex of neighboring site
    if spin == 1
        l = get_spindices_from_index(nbr_site, model_geometry)[1]
    else
        l = get_spindices_from_index(nbr_site, model_geometry)[2]
    end

    @assert(get_spindex_type(k, model_geometry) == get_spindex_type(l, model_geometry))

    @debug """
    ParticleConfiguration::propose_random_move() :
    proposing random move =>
    particle: $(β)
    isite: $(k)
    jsite: $(l)
    """

    # whether move is possible
    if pconfig[l] == 0
        possible = true
    else
        possible = false
    end

    return MarkovMove(β, k, l, possible)
end


@doc raw"""

    hop!( markov_move::MarkovMove{I}, 
          pconfig::AbstractVector{I} ) where {I<:Integer}

If a proposed move (hopping) is accepted, updates the particle positions in 
the currenrt configuration.

- `markov_move::MarkovMove{I}`: quantities related to a Markov process.  
- `pconfig::AbstractVector{I}`: current particle configuration. 

"""
function hop!(
    markov_move::MarkovMove{I}, 
    pconfig::Vector{I}
) where {I<:Integer}
    @assert(markov_move.possible)

    # particle number
    β = markov_move.particle

    # initial site
    k = markov_move.k

    # final site
    l = markov_move.l

    @assert(pconfig[k] !== 0)
    @assert(pconfig[l] == 0)

    @debug """
    ParticleConfiguration::hop!() :
    preparing to hop =>
    particle $(β) from $(k) to $(l)
    """

    # update particle positions
    pconfig[l] = pconfig[k]
    pconfig[k] = 0

    @debug """
    ParticleConfiguration::hop!() :
    particle positions are =>
    $(pconfig)
    """

    return nothing
end


@doc raw"""

    exchange!( markov_move::MarkovMove{I}, 
               pconfig::AbstractVector{I},
               model_geometry::ModelGeometry ) where {I<:Integer}

If a proposed move (exchange) is accepted, updates the particle positions
in the current configuration.

- `markov_move::MarkovMove{I}`: quantities related to a Markov process.  
- `pconfig::AbstractVector{I}`: current particle configuration. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
# TODO: this method needs to be debugged!
function exchange!(
    markov_move::MarkovMove{I}, 
    pconfig::Vector{I},
    model_geometry::ModelGeometry
) where {I<:Integer}
    @assert(markov_move.possible)

    # particle number
    β = markov_move.particle

    # initial site
    k = markov_move.k

    # final site
    l = markov_move.l

    # get real site indices
    ksite = get_index_from_spindex(k, model_geometry)
    lsite = get_index_from_spindex(l, model_geometry)

    @debug """
    ParticleConfiguration::exchange!() :
    preparing to exchange =>
    particle: $(β₁)
    site1: $(k)
    particle2: $(β₂)
    site2: $(l)
    """

    # update particle positions
    pconfig[lsite] = pconfig[ksite]
    pconfig[lsite + N] = pconfig[ksite + N]
    pconfig[ksite] = 0
    pconfig[ksite + N] = 0

    @debug """
    ParticleConfiguration::exchange!() :
    particle positions are =>
    $(pconfig)
    """

    return nothing
end


@doc raw"""

    generate_initial_fermion_configuration!( pconfig::AbstractVector{I},
                                             nup::I, 
                                             ndn::I, 
                                             model_geometry::ModelGeometry, 
                                             rng::AbstractRNG ) where {I<:Integer}

Generates a random initial configuration of spin-up and spin-down fermions. The first `N` 
elements correspond to spin-up and the last `N` correspond to spin-down. Occupation is 
denoted by a positive integer corresponding to that particle's creation operator label. 

- `pconfig::Vector{I}`: vector to store configurations.
- `nup::I`: number of spin-up fermions.
- `ndn::I`: number of spin-down fermions.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `rng::AbstractRNG`: random number generator.

"""
function generate_initial_fermion_configuration!(
    pconfig::Vector{I},
    nup::I, 
    ndn::I, 
    model_geometry::ModelGeometry, 
    rng::AbstractRNG
) where {I<:Integer}
    # lattice sites
    N = model_geometry.lattice.N

    # reset configuration to zeros
    fill!(pconfig, 0)

    # assign spin-up electrons
    up_indices = shuffle(rng, 1:N)[1:nup]
    for (i, idx) in enumerate(up_indices)
        pconfig[idx] = i
    end

    # assign spin-down electrons
    down_indices = shuffle(rng, 1:N)[1:ndn]
    for (i, idx) in enumerate(down_indices)
        pconfig[idx + N] = i + nup
    end

    return pconfig
end


@doc raw"""

    get_onsite_fermion_occupation( site::I, 
                                   pconfig::AbstractVector{I},
                                   N::I ) where {I<:Integer}

Returns the number of spin-up and spin-down fermions occupying a real lattice site `i`.  

- `site::I`: lattice site. 
- `pconfig::AbstractVector{I}`: current particle configuration.
- `N::I`: total number of sites.

"""
function get_onsite_fermion_occupation(
    site::I, 
    pconfig::Vector{I},
    N::I
) where {I<:Integer}
    # count number of fermions
    num_up = pconfig[site] > 0 ? 1 : 0

    # count number of spin-down fermions
    num_dn = pconfig[site + N] > 0 ? 1 : 0

    # total number of fermions
    num_f = num_up + num_dn

    return num_up, num_dn, num_f
end


@doc raw"""

    get_particle_density( density::E, 
                          model_geometry::ModelGeometry,
                          pht::Bool ) where {E<:AbstractFloat}

Given a particle density, returns the total number of particles ``N_p``, the total number of 
electrons ``N_e`` number of spin-up fermions ``n_{\uparrow}``, number of spin-down fermions ``n_{\downarrow}``.
In cases where the given density is not a commensurate filling, a new density will be 
calculated to correspond with the new particle number. 

- `density::E`: desired electronic density.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_particle_density(
    density::E,
    model_geometry::ModelGeometry,
    pht::Bool
) where {E<:AbstractFloat}
    # total number of lattice sites
    N = model_geometry.lattice.N

    # total number of particles
    Np = 2 * round(Int, (density * N) / 2)
    
    # update to the appropriate density 
    if Np / N != density
        density = Np / N
        @debug """
        ParticleConfiguration::get_particle_numbers() :
        Particle density has been updated -> density = $density 
        """
    end

    # number of spin-up particles
    nup = Np / 2

    # number of spin-down particles
    if !pht
        ndn = nup
    else
        ndn = nup + N - Np
    end

    # total number of electrons
    Ne = nup + ndn

    return density, Int(Np), Int(Ne), Int(nup), Int(ndn)
end


@doc raw"""

    get_particle_density( Np::I,
                          model_geometry::ModelGeometry,
                          pht::Bool ) where {I<:Integer}

Given the total number of particles in the lattice ``N_p``, returns the particle density, 
the total number of electrons ``N_e``, the number of spin-up fermions ``n_{\uparrow}`` and 
the number of spin-down fermions ``n_{\downarrow}``.

- `Np::I`: total number of particles. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_particle_density(
    Np::I,
    model_geometry::ModelGeometry,
    pht::Bool
) where {I<:Integer}
    # check that the given Np is an integer
    @assert Np isa Integer && isinteger(Np) "Np must be an integer, but got $Np"

    # total number of lattice sites
    N = model_geometry.lattice.N

    # number of spin-up particles
    nup = Np / 2

    # number of spin-down particles
    if !pht
        ndn = nup
    else
        ndn = nup + N - Np
    end

    # total number of electrons
    Ne = nup + ndn

    # particle density
    density = Np / N

    return density, Np, Int(Ne), Int(nup), Int(ndn)
end


@doc raw"""

    get_particle_density( nup::I, 
                          ndn::I,
                          model_geometry::ModelGeometry, 
                          pht::Bool ) where {I<:Integer}

Given the number of spin-up fermions ``n_{\uparrow}``, and number of spin-down fermions ``n_{\downarrow}``, 
returns the particle density, total number of particles ``N_p``, and the total number of electrons ``N_e``.

- `nup::I`: number of spin-up fermions.
- `ndn::I`: number of spin-down fermions.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_particle_density(
    nup::I, 
    ndn::I,
    model_geometry::ModelGeometry,
    pht::Bool
) where {I<:Integer}
    # check that nup and ndn are integers
    @assert nup isa Integer && isinteger(nup) "nup must be an integer, but got $nup"
    @assert ndn isa Integer && isinteger(ndn) "ndn must be an integer, but got $ndn"

    # total number of lattice sites
    N = model_geometry.lattice.N

    # total number of particles and electrons
    if !pht
        Np = nup + ndn
        Ne = Np
    else
        Np = nup + N - ndn 
        Ne = nup + ndn
    end

     # particle density 
    density = Np / N

    return density, Int(Np), Int(Ne), nup, ndn
end


@doc raw"""

    get_spindex_type( spindex::I, 
                      model_geometry::ModelGeometry ) where {I<:Integer}

Returns the spin species at a given spindex.

- `spindex::I`: spin index.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_spindex_type(
    spindex::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
    @assert spindex < 2 * model_geometry.lattice.N + 1

    return spindex < model_geometry.lattice.N + 1 ? 1 : -1
end


@doc raw"""

    get_index_from_spindex( spindex::I, 
                            model_geometry::ModelGeometry ) where {I<:Integer}

Returns the lattice site ``i`` for a given spindex.

- `spindex::I`: spin index.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_index_from_spindex(
    spindex::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
    L = model_geometry.lattice.N
    @assert spindex < 2 * L + 1

    return spindex <= L ? spindex : spindex - L
end


@doc raw"""

    get_spindices_from_index( index::I, 
                              model_geometry::ModelGeometry ) where {I<:Integer}

Returns spin-up and spin-down indices from a given site index.

- `index::I`: lattice site index.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_spindices_from_index(
    index::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
    L = model_geometry.lattice.N
    @assert index <= L

    return index, index + L
end


@doc raw"""

    get_linked_spindex( i::I, 
                        N::I ) where {I<:Integer}

 Given an index ``i`` in the spin-up sector, returns an index in the spin-down sector.

- i::I: lattice index in the spin-up sector.
- N::I: total number of lattice sites. 

"""
function get_linked_spindex(
    i::I, 
    N::I
) where {I<:Integer}
    @assert i < 2 * N

    return i + (1 - 2 * (i ÷ N)) * N
end


