@doc raw"""

    MarkovMove( particle::Int, 
                k::Int, 
                l::Int, 
                possible::Bool )

A type defining a Markov move.

"""
struct MarkovMove
    # particle
    particle::Int

    # initial spindex
    k::Int

    # neighboring spindex
    l::Int

    # hop possible
    possible::Bool
end


@doc raw"""

    propose_random_move( Ne::Int64, 
                         pconfig::Vector{Int64}, 
                         model_geometry::ModelGeometry, 
                         rng::Xoshiro )::MarkovMove

Proposes randomly hopping or exchanging a particle from some intial site 'k' to a neighboring site 'l' 
and returns an instance of the MarkovMove type. 

- `Ne::Int64`: total number of electrons.  
- `pconfig::Vector{Int64}`: current particle configuration. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `rng::Xoshiro`: random number generator.

"""
function propose_random_move(
    Ne::Int64, 
    pconfig::Vector{Int64}, 
    model_geometry::ModelGeometry, 
    rng::Xoshiro
)::MarkovMove
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

    # select a particle, β
    β = rand(rng, 1:Ne)

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

    debug && println("ParticleConfiguration::propose_random_move() : proposing random move")
    debug && println("particle: ", β, ", isite: ", k, ", jsite: ", l)

    # whether move is possible
    if pconfig[l] == 0
        possible = true
    else
        possible = false
    end

    return MarkovMove(β, k, l, possible)
end


@doc raw"""
    hop!( markov_move::MarkovMove, 
          pconfig::Vector{Int} )::Nothing

If proposed hopping move is accepted, updates the particle positions.

- `markov_move::MarkovMove`: quantities related to a Markov move.  
- `pconfig::Vector{Int}`: current particle configuration. 

"""
function hop!(
    markov_move::MarkovMove, 
    pconfig::Vector{Int}
)::Nothing
    @assert(markov_move.possible)

    # particle number
    β = markov_move.particle

    # initial site
    k = markov_move.k

    # final site
    l = markov_move.l

    @assert(pconfig[k] !== 0)
    @assert(pconfig[l] == 0)

    debug && println("ParticleConfiguration::hop!() : preparing to hop")
    debug && println("particle ", β, " from ", k, " to ", l)

    # update particle positions
    pconfig[l] = pconfig[k]
    pconfig[k] = 0

    debug && println("ParticleConfiguration::hop!() : particle positions are")
    debug && println(pconfig)

    return nothing
end


@doc raw"""
    exchange!( markov_move::MarkovMove, 
               pconfig::Vector{Int},
               model_geometry::ModelGeometry )::Nothing

If proposed exchange move is accepted, updates the particle positions.

- `markov_move::MarkovMove`: quantities related to a Markov move.  
- `pconfig::Vector{Int}`: current particle configuration. 

"""
# TODO: need to debug!
function exchange!(
    markov_move::MarkovMove, 
    pconfig::Vector{Int},
    model_geometry::ModelGeometry
)::Nothing
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

    debug && println("ParticleConfiguration::exchange!() : preparing to exchange")
    debug && println("particle: ", β, ", isite: ", k, ", jsite ", l)

    # update particle positions
    pconfig[lsite] = pconfig[ksite]
    pconfig[lsite + N] = pconfig[ksite + N]
    pconfig[ksite] = 0
    pconfig[ksite + N] = 0

    debug && println("ParticleConfiguration::exchange!() : particle positions are")
    debug && println(pconfig)

    return nothing
end


@doc raw"""
    generate_initial_fermion_configuration( nup::Int64, 
                                            ndn::Int64, 
                                            model_geometry::ModelGeometry, 
                                            rng::Xoshiro )::Vector{Int64}

Generates a random initial configuration of spin-up and spin-down fermions. The first N elements correspond 
to spin-up and the last N correspond to spin-down. Occupation is denoted by a positive integer corresponding 
to that particle's creation operator label. 

- `nup::Int64`: number of spin-up electrons.
- `ndn::Int64`: number of spin-down electrons.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `rng::Xoshiro`: random number generator.

"""
function generate_initial_fermion_configuration(
    nup::Int64, 
    ndn::Int64, 
    model_geometry::ModelGeometry, 
    rng::Xoshiro
)::Vector{Int64}
    # lattice sites
    N = model_geometry.lattice.N

    # initialize configuration
    pconfig = fill(0, 2 * N)

    # Ensure unique random assignments for up-spin electrons
    up_indices = shuffle(rng, 1:N)[1:nup]
    for (i, idx) in enumerate(up_indices)
        pconfig[idx] = i
    end

    # Ensure unique random assignments for down-spin electrons
    down_indices = shuffle(rng, 1:N)[1:ndn]
    for (i, idx) in enumerate(down_indices)
        pconfig[idx + N] = i + nup
    end

    return pconfig
end


@doc raw"""

    get_onsite_fermion_occupation( site::Int, 
                                   pconfig::Vector{Int} )::Tuple{Int64, Int64, Int64}

Returns the number of spin-up and spin-down electrons occupying a real lattice site 'i'.  

- `site::Int`: lattice site. 
- `pconfig::Vector{Int}`: current particle configuration.

"""
function get_onsite_fermion_occupation(
    site::Int, 
    pconfig::Vector{Int}
)::Tuple{Int64, Int64, Int64}
    # count number of spin-up electrons
    num_up = pconfig[site] > 0 ? 1 : 0

    # count number of spin-down electrons
    num_dn = pconfig[site + model_geometry.lattice.N] > 0 ? 1 : 0

    # total number of electrons
    num_e = num_up + num_dn

    return num_up, num_dn, num_e
end


@doc raw"""

    get_particle_numbers( density::Float64 )::NTuple{4, Int64}

Given a particle density, returns the total number of particles Np, number of spin-up particles Npu, 
number of spin-down particles Npd, and total number of electrons Ne.

- `density::Float64`: desired electronic density.

"""
function get_particle_numbers(
    density::Float64
)::NTuple{4, Int64}
    N = model_geometry.lattice.N
    Np = density * N

    # Check if the total number of particles Np is an integer
    @assert isapprox(Np, round(Np), atol=1e-8) "Density does not correspond to a commensurate filling (Np must be an integer)"
    
    # Np must also be even
    @assert Np % 2 == 0 "Np must be even"
    
    if !pht
        # Number of spin-up and spin-down particles for particle-hole transformation off
        Npu = Np ÷ 2
        Npd = Np ÷ 2
        nup = Npu
        ndn = Npd
        Ne = Np
    else
        # Particle-hole transformation on
        nup = Np ÷ 2
        Npu = nup
        ndn = N + nup - Np
        Npd = N - ndn
        Ne = nup + ndn
    end
    
    return Int(Np), Int(Ne), Int(nup), Int(ndn)
end


@doc raw"""

    get_particle_density( nup::Int64, 
                          ndn::Int64 )::Tuple{Float64, Int64, Int64}

Given the number of spin-up electrons nup, and number of spin-down electrons ndn, returns 
the particle density, total number of particles Np, and the total number of electrons Ne.

- `nup::Int`: number of spin-up electrons.
- `ndn::Int`: number of spin-down electrons.

"""
function get_particle_density(
    nup::Int, 
    ndn::Int
)::Tuple{Float64, Int64, Int64}
    N = model_geometry.lattice.N

    if !pht
        # Number of spin-up and spin-down particles for particle-hole transformation off
        Npu = nup  
        Npd = ndn 
        Np = Npu + Npd
        Ne = Np
    else
        # Particle-hole transformation on
        Npu = nup
        Npd = N - ndn 
        Np = Npu + Npd
        Ne = nup + ndn
    end

    # calculate particle density
    density = Np / N

    return density, Int(Np), Int(Ne)
end


@doc raw"""
    get_spindex_type( spindex::Int, 
                      model_geometry::ModelGeometry )::Int

Returns the spin species at a given spindex.

- `spindex::Int`: spin index.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_spindex_type(
    spindex::Int, 
    model_geometry::ModelGeometry
)::Int
    @assert spindex < 2 * model_geometry.lattice.N + 1
    return spindex < model_geometry.lattice.N + 1 ? 1 : -1
end


@doc raw"""
    get_index_from_spindex( spindex::Int, 
                            model_geometry::ModelGeometry )::Int 

Returns the lattice site i for a given spindex.

- `spindex::Int`: spin index.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_index_from_spindex(
    spindex::Int, 
    model_geometry::ModelGeometry
)::Int
    L = model_geometry.lattice.N
    @assert spindex < 2 * L + 1
    return spindex <= L ? spindex : spindex - L
end


@doc raw"""
    get_spindices_from_index( index::Int, 
                              model_geometry::ModelGeometry )::Tuple{Int64, Int64}

Returns spin-up and spin-down indices from a given site index.

- `index::Int`: lattice site index.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_spindices_from_index(
    index::Int, 
    model_geometry::ModelGeometry
)::Tuple{Int64, Int64}
    L = model_geometry.lattice.N
    @assert index <= L
    return index, index + L
end


@doc raw"""
    get_linked_spindex( i::Int, 
                        N::Int )::Int

Returns an index in the spin-down sector, given an index in the spin-up sector.

- i::Int: lattice index in the spin-up sector.
- N::Int: total number of lattice sites. 

"""
function get_linked_spindex(
    i::Int, 
    N::Int
)::Int
    @assert i < 2 * N
    return i + (1 - 2 * (i ÷ N)) * N
end


##################################################### DEPRECATED FUNCTIONS #####################################################
# """
#     change_particle_number!( met_step, phconfig::Vector{Int}, model_geometry::ModelGeometry )

# If proposed particle addition or removal is accepted, add or remove a particle, and update the particle 
# configuration..

# """
# function change_particle_number!(met_step, phconfig, model_geometry)
#     # lattice site
#     N = model_geometry.lattice.N
#     i = met_step[3]

#     # addition or removal
#     cran = met_step[2]

#     if cran == 1
#         phconfig[i] += 1
#     elseif cran == -1
#         phconfig -= 1
#     end

#     return nothing
# end


# """
#     [DEPRECATED]
#     get_particle_positions( pconfig::Vector{Int}, model_geometry::ModelGeometry, Ne::Int )

# Given vector of spindex occupations, returns a vector κ of particle positions according to particle number.
# This vector is gives the order of creation operator which create the configuration. Initial configurations
# are automatically normal ordered.

# """
# function get_particle_positions(pconfig::Vector{Int}, model_geometry::ModelGeometry, Ne::Int)
#     N = model_geometry.lattice.N

#     # get current spindex occupations
#     config_indices = findall(x -> x == 1, pconfig)

#     κ = zeros(Int, 2 * N)

#     # fill κ according to position of particles
#     for i in 1:Ne
#         idx = config_indices[i]
#         κ[idx] = i 
#     end

#     return κ
# end

# """
#     [DEPRECATED]
#     get_onsite_fermion_occupation( site::Int, pconfig::Vector{Int} )

# Returns the number of spin-up and spin-down electrons occupying a real lattice site i.  

# """
# function get_onsite_fermion_occupation(site::Int, pconfig::Vector{Int})
#     nup = pconfig[site]
#     ndn = pconfig[site+model_geometry.lattice.N]
#     Ne = pconfig[site] + pconfig[site+model_geometry.lattice.N]
#     return nup, ndn, Ne
# end

# """

#     generate_initial_phonon_density_configuration()

# Returns initial vector to store phonon occupations.

# """
# function generate_initial_phonon_density_configuration(model_geometry::ModelGeometry)
#     N = model_geometry.lattice.N
#     init_phconfig = zeros(Int, N)
 
#     return init_phconfig
# end



# """

#     generate_initial_phonon_displacement_configuration( )

# Returns initial matrix to store phonon displacements. Each column represents displacements in each dimension. 

# """
# function generate_initial_phonon_displacement_configuration(loc::AbstractString, model_geometry::ModelGeometry)
#     # dimensions
#     dims = size(model_geometry.lattice.L)[1]
#     N = model_geometry.lattice.N

#     if loc == "bond"
#         # create neighbor table
#         nbr_table = build_neighbor_table(bonds[1],
#         model_geometry.unit_cell,
#         model_geometry.lattice)

#         # maps neighbor table to dictionary of bonds and neighbors                                
#         nbrs = map_neighbor_table(nbr_table)

#         # Collect all bonds into a Set to ensure uniqueness
#         unique_bonds = Set{Int}()

#         for (site, bond_info) in nbrs
#         # Add bonds to the set (automatically handles duplicates)
#         for bond in bond_info.bonds
#         push!(unique_bonds, bond)
#         end
#         end

#         # The total number of unique bonds is the length of the set
#         total_bonds = length(unique_bonds)

#         # initialize phonon configuration
#         init_phconfig = zeros(AbstractFloat, total_bonds, dims)        
#     elseif loc == "onsite"
#         # initialize phonon configuration
#         init_phconfig = zeros(AbstractFloat, N, dims)
#     else
#         println("ERROR: Not a valid location in the lattice!")
#     end

#     return init_phconfig
# end

# """
#     get_phonon_occupation( site::Int, phconfig::Vector{Int} )

# Returns the number of phonons occupying a real lattice site i or a lattice bond i.  

# """
# function get_phonon_occupation(i::Int, phconfig::Vector{Int})
#     nph = phconfig[i]

#     return nph
# end


# """
#     get_onsite_phonon_displacement( site::Int, phconfig::Vector{Int} )

# Returns the X and Y displacements of a phonon at site i.

# """
# function get_onsite_phonon_displacement(i::Int, phconfig::Matrix{AbstractFloat})
#     Xᵢ, Yᵢ = phconfig[i, :]

#     return (Xᵢ, Yᵢ)
# end


# """
#     get_bond_phonon_displacement( b::Int, phconfig::Vector{Int} )

# Given a bond in the lattice b, returns the displacement of phonon located on that bond. 

# """
# function get_bond_phonon_displacement(b::Int, phconfig::Matrix{AbstractFloat})
#     Xᵢ, Yᵢ = phconfig[i, :]

#     return [Xᵢ, Yᵢ]
# end












