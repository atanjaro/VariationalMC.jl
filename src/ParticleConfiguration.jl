@doc raw"""

    ParticleConfiguration

A type defining a particle configuration in the canonical ensemble.

# Comment

If the type field `ph_transform = true` then a (partial) particle-hole transformation is performed
on the spin-down sector. This is required when intorducing pair fields.

# Fields

- `density::E`: Density of particles on the lattice.
- `Np::I`: Total number of particles on the lattice.
- `nup::I`: Number of spin-up particles.
- `ndn::I`: Number of spin-down particles.
- `Ne::I`: Total number of electrons on the lattice.
- `pconfig::Vector{I}`: Vector of particle positions.
- `ph_transform::Bool`: Whether the particle configuration is particle-hole transformed.

"""

mutable struct ParticleConfiguration
    # particle density
    density::AbstractFloat
    
    # total number of particles
    Np::Int

    # number of spin-up particles
    nup::Int
    
    # number of spin-down particles
    ndn::Int
    
    # total number of electrons
    Ne::Int

    # vector of particle postions
    pconfig::Vector{Int}

    # whether configuration is particle-hole transformed
    ph_transform::Bool
end


@doc raw"""

    ParticleConfiguration( 
        density::AbstractFloat; 
        # KEYWORD ARGUMENTS
        model_geometry::ModelGeometry,
        particle_hole_transform::Bool 
    ) where {E<:AbstractFloat}

Given an initial particle density, initialize and return an instance of the type `ParticleConfiguration`.

- `density::AbstractFloat`: desired electronic density.

# KEYWORD ARGUMENTS

- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `ph_transform::Bool`: Whether configuration will be particle-hole transformed.

"""
function ParticleConfiguration(
    density::AbstractFloat;
    # KEYWORD ARGUMENTS
    model_geometry::ModelGeometry,
    ph_transform::Bool
)
    # total number of lattice sites
    N_orbs = model_geometry.unit_cell.n
    N_cells = model_geometry.lattice.N
    N = N_orbs * N_cells

    # total number of particles
    Np = 2 * round(Int, (density * N) / 2)
    
    # update to the appropriate density 
    if Np / N != density
        density = Np / N
    end

    # number of spin-up particles
    nup = Np / 2

    # number of spin-down particles
    if !ph_transform
        ndn = nup
    else
        ndn = nup + N - Np
    end

    # total number of electrons
    Ne = nup + ndn

    # allocated empty particle configuration
    pconfig = zeros(Int, 2*N)
    
    return ParticleConfiguration(density, Int(Np), Int(nup), Int(ndn), Int(Ne), pconfig, ph_transform)
end


@doc raw"""

    ParticleConfiguration( 
        Np::Int;
        # KEYWORD ARGUMENTS
        model_geometry::ModelGeometry,
        particle_hole_transform::Bool 
    ) where {I<:Integer}

Given an initial particle number, initialize and return an instance of the type `ParticleConfiguration`.

- `Np::Int`: total number of particles. 

# KEYWORD ARGUMENTS

- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `ph_transform::Bool`: Whether configuration is particle-hole transformed. 

"""
function ParticleConfiguration(
    Np::Int;
    # KEYWORD ARGUMENTS
    model_geometry::ModelGeometry,
    ph_transform::Bool
)
    # check that the given Np is an integer
    @assert Np isa Integer && isinteger(Np) "Np must be an integer, but got $Np"

    # total number of lattice sites
    N_orbs = model_geometry.unit_cell.n
    N_cells = model_geometry.lattice.N
    N = N_orbs * N_cells

    # number of spin-up particles
    nup = Np / 2

    # number of spin-down particles
    if !ph_transform
        ndn = nup
    else
        ndn = nup + N - Np
    end

    # total number of electrons
    Ne = nup + ndn

    # particle density
    density = Np / N

    # allocated empty particle configuration
    pconfig = zeros(Int, 2*N)

    return ParticleConfiguration(density, Np, Int(nup), Int(ndn), Int(Ne), pconfig, ph_transform)
end


@doc raw"""

    ParticleConfiguration(;
        # KEYWORD ARGUMENTS
        nup::Int, 
        ndn::Int,
        model_geometry::ModelGeometry, 
        particle_hole_transform::Bool 
    ) where {I<:Integer}

Given a certain number of spin-up and spin-down particles, initialize and return an instance of the type `ParticleConfiguration`.

# KEYWORD ARGUMENTS
- `nup::Int`: number of spin-up fermions.
- `ndn::Int`: number of spin-down fermions.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `ph_transform::Bool`: Whether configuration is particle-hole transformed. 

"""
function ParticleConfiguration(;
    # KEYWORD ARGUMENTS
    nup::Int, 
    ndn::Int,
    model_geometry::ModelGeometry,
    ph_transform::Bool
)
    # check that nup and ndn are integers
    @assert nup isa Integer && isinteger(nup) "nup must be an integer, but got $nup"
    @assert ndn isa Integer && isinteger(ndn) "ndn must be an integer, but got $ndn"

    # total number of lattice sites
    N_orbs = model_geometry.unit_cell.n
    N_cells = model_geometry.lattice.N
    N = N_orbs * N_cells

    # total number of particles and electrons
    if !ph_transform
        Np = nup + ndn
        Ne = Np
    else
        Np = nup + N - ndn 
        Ne = nup + ndn
    end

     # particle density 
    density = Np / N

    # allocated empty particle configuration
    pconfig = zeros(Int, 2*N)

    return ParticleConfiguration(density, Int(Np), nup, ndn, Int(Ne), pconfig, ph_transform)
end


# print struct info in TOML format
function Base.show(io::IO, ::MIME"text/plain", pconf::ParticleConfiguration)

    (; density, Np, nup, ndn, Ne, ph_transform) = pconf

    @printf io "[FermionConfiguration]\n\n"
    @printf io "DENSITY         = %.8f\n" density
    @printf io "PH_TRANSFORM    = %.d\n\n" ph_transform
    @printf io "[FermionConfiguration.ensemble]\n\n"
    @printf io "PARTICLE_NUM    = %.d\n" Np
    @printf io "SPIN-UP         = %.d\n" nup
    @printf io "SPIN-DOWN       = %.d\n" ndn
    @printf io "ELECTRON_NUM    = %.d\n" Ne

    return nothing
end


@doc raw"""

    generate_random_fermion_configuration!(
        pconfig::Vector{Int},
        nup::Int, 
        ndn::Int, 
        N::Int,
        rng::AbstractRNG
    ) 

Generates a random initial configuration of spin-up and spin-down fermions. The first `N` 
elements correspond to spin-up and the last `N` correspond to spin-down. Occupation is 
denoted by a positive integer corresponding to that particle's creation operator label. 

# Fields

- `pconfig::Vector{Int}`: Storage vector for configuration.
- `nup::Int`: Number of spin-up fermions.
- `ndn::Int`: Number of spin-down fermions.
- `N`: Total number of lattice sites.
- `rng::AbstractRNG`: Random number generator.

"""
function generate_random_fermion_configuration!(
    pconfig::Vector{Int},
    nup::Int, 
    ndn::Int, 
    N::Int,
    rng::AbstractRNG
) 
    # reset configuration to vacuum
    fill!(pconfig, 0)

    # assign spin-up fermions
    up_indices = shuffle(rng, 1:N)[1:nup]
    for (i, idx) in enumerate(up_indices)
        pconfig[idx] = i
    end

    # assign spin-down fermions
    down_indices = shuffle(rng, 1:N)[1:ndn]
    for (i, idx) in enumerate(down_indices)
        pconfig[idx + N] = i + nup
    end

    return pconfig
end


@doc raw"""

    get_fermion_occupations(
        site::I, 
        pconfig::Vector{I},
        N::I
    ) where {I<:Integer}

Returns the number of spin-up, spin-down, as well as the total fermion occupation at a lattice site.

"""
function get_fermion_occupations(
    site::I, 
    pconfig::Vector{I},
    N::I
) where {I<:Integer}
    # count number of spin-up fermions
    num_up = pconfig[site] > 0 ? 1 : 0

    # count number of spin-down fermions
    num_dn = pconfig[site + N] > 0 ? 1 : 0

    return num_up, num_dn
end


@doc raw"""

    calculate_chemical_potential(   
        N::I,
        n::I,
        Ne::I, 
        bond_slices::Vector{UnitRange{I}},
        neighbor_table::Matrix{I},
        t::Vector{E}
    ) where {I<:Integer, E<:AbstractFloat}

Calclates the exact chemical potential for a non-interacting tight binding model.

# Fields

- `N::I`: Number of unit cells in the lattice.
- `n::I`: Number of orbitals per unit cell.
- `Ne::I`: Number of electrons in the particle configuration.
- `bond_slices::Vector{UnitRange{I}}`: View into neighbor table for each bond ID.
- `neighbor_table::Matrix{I}`: Neighbor table for all pairs of orbitals connected by a bond in the lattice.
- `t::Vector{E}`: Hopping energies for all pairs of orbitals connected by a bond in the lattice.

"""
function calculate_chemical_potential(
    N::I,
    n::I,
    Ne::I, 
    bond_slices::Vector{UnitRange{I}},
    neighbor_table::Matrix{I},
    t::Vector{E}
) where {I<:Integer, E<:AbstractFloat}
    Nsites = n*N

    # allocate hopping matrix
    H_tb = zeros(Complex, 2*Nsites, 2*Nsites)

    for slice in bond_slices
        # get slice of neighbor table
        nbr_slice = neighbor_table[:,slice]

        # get slice of hopping parameters
        t_slice = t[slice]

        # spin-up sector
        for ((i,j),t_hop) in zip(eachcol(nbr_slice),t_slice)
            H_tb[i,j] += -t_hop
            H_tb[j,i] += -t_hop
        end

        # spin-down sector
        for ((i,j),t_hop) in zip(eachcol(nbr_slice),t_slice)
            H_tb[i+ Nsites,j+ Nsites] += -t_hop
            H_tb[j+ Nsites,i+ Nsites] += -t_hop
        end
    end

    # solve for eigenvalues
    ε_F, _ = diagonalize!(H_tb)

    # calculate the chemical potential
    μ = 0.5 * (ε_F[Ne + 1] + ε_F[Ne])

    return μ
end


# @doc raw"""

#     PhononConfiguration

# A type defining different phonon (boson) configurations.

# # Fields

# - `density::E`: Density of phonons on the lattice.
# - `Nph::Int`: Total number of phonons on the lattice.
# - `dconfig::Vector{Int}`: Density configuration of phonons for `N` lattice sites, where each element is the phonon number `nᵖʰᵢ` at site `i`.
# - `xconfig::Vector{AbstractFloat}`: Displacement configuration of phonons for `N` lattice sites, where each element is the phonon displacement `Xᵢ` at site `i`.

# """
# mutable struct PhononConfiguration
#     # phonon density
#     density::AbstractFloat

#     # total number of phonons
#     Nph::Int

#     # phonon density configuration
#     dconfig::Vector{Int}

#     # phonon displacement configuration
#     xconfig::Vector{AbstractFloat}
# end


# # print struct info in TOML format
# function Base.show(io::IO, ::MIME"text/plain", phconf::PhononConfiguration)

#     (; density, Nph, dconfig, xconfig) = phconf

#     @printf io "[BosonConfiguration]\n\n"
#     @printf io "DENSITY         = %.8f\n" density

#     @printf io "[BosonConfiguration.ensemble]\n\n"
    
#     if !isnothing(dconfig)
#         @printf io "TYPE    = DENSITY"
#         @printf io "PARTICLE_NUM    = %.d\n" Nph
#     elseif !isnothing(xconfig)
#         @printf io "TYPE    = DISPLACEMENT"
#     end

#     return nothing
# end































