@doc raw"""

    HubbardParameters{T<:AbstractFloat}

Hubbard parameters for finite lattice.

# Fields

- `U::Vector{T}`: On-site Hubbard interaction for each site with finite Hubbard interaction.
- `sites::Vector{Int}`: Site index associated with each finite Hubbard `U` interaction.
- `orbital_ids::Vector{Int}`: Orbital ID/species in unit cell with finite Hubbard interaction.

"""
struct HubbardParameters{T<:AbstractFloat}
    # Hubbard U for each orbital in the lattice
    U::Vector{T}

    # site index associated with each Hubbard U
    sites::Vector{Int}

    # orbital species in unit cell with finite hubbard interaction
    orbital_ids::Vector{Int}
end

@doc raw"""

    HubbardParameters(;
        hubbard_model::HubbardModel{T},
        model_geometry::ModelGeometry{D,T}
    ) where {D, T<:AbstractFloat}

Initialize an instance of `HubbardParameters`.

"""
function HubbardParameters(;
    hubbard_model::HubbardModel{T},
    model_geometry::ModelGeometry{D,T}
) where {D, T<:AbstractFloat}

    (; U_orbital_ids, U_mean) = hubbard_model
    (; lattice, unit_cell) = model_geometry

    # number of orbitals with finite hubbard interaction in unit cell
    n_hubbard = length(hubbard_model.U_orbital_ids)

    # number of unit cell in lattice
    N_unitcells = lattice.N

    # get the number of HS transformations per time-slice τ that is applied
    N_hubbard = N_unitcells * n_hubbard

    # allocate arrays
    U     = zeros(T, N_hubbard)
    sites = zeros(Int, N_hubbard)

    # reshape the allocated arrays
    U′     = reshape(U, (N_unitcells, n_hubbard))
    sites′ = reshape(sites, (N_unitcells, n_hubbard))

    # total number of orbitals in lattice
    N_orbitals = nsites(unit_cell, lattice)

    # iterate over orbitals in the unit cell with finite hubbard interaction
    for (n,o) in enumerate(U_orbital_ids)
        # iterate over unit cells in the lattice
        for u in 1:N_unitcells
            # calculate the site associated with the hubbard interaction
            sites′[u,n] = loc_to_site(u, o, unit_cell)
            # get the Hubbard U interaction on the site
            U′[u,n] = U_mean[n]
        end
    end

    return HubbardParameters(U, sites, U_orbital_ids)
end


@doc raw"""

    ExtendedHubbardParameters{T<:AbstractFloat}

Extended Hubbard interaction parameters for finite lattice.

# Fields

- `V::Vector{T}`: Extended Hubbard interaction strength for each pair neighbors in the lattice.
- `neighbor_table::Matrix{Int}`: Neighbor table for extended Hubbard interactions on lattice.
- `bond_ids::Vector{Int}`: Bond IDs used to define extended Hubbard interactions.

"""
struct ExtendedHubbardParameters{T<:AbstractFloat}
    # extended hubbard interaction strength between pair of orbitals
    V::Vector{T}

    # neighbor table for extended hubbard interactions
    neighbor_table::Matrix{Int}

    # bond IDs used to define extended hubbard interactions
    bond_ids::Vector{Int}
end


@doc raw"""

    ExtendedHubbardParameters(;
        # KEYWORD ARGUMENTS
        extended_hubbard_model::ExtendedHubbardModel{T},
        model_geometry::ModelGeometry{D,T}
    ) where {D, T<:AbstractFloat}

Initialize an instance of the `ExtendedHubbardParameters` type.

"""
function ExtendedHubbardParameters(;
    # KEYWORD ARGUMENTS
    extended_hubbard_model::ExtendedHubbardModel{T},
    model_geometry::ModelGeometry{D,T}
) where {D, T<:AbstractFloat}

    (; V_bond_ids, V_mean) = extended_hubbard_model
    (; bonds, unit_cell, lattice) = model_geometry

    # construct neighbor table to extended Hubbard interactions
    V_bonds = @view bonds[V_bond_ids]
    neighbor_table = build_neighbor_table(V_bonds, unit_cell, lattice)

    # number of types of extended hubbard interactions
    n = length(V_bond_ids)

    # number of unit cells with an extended hubbard interaction
    N = size(neighbor_table, 2) ÷ n

    # calculate the interaction stength between each pair of sites
    V  = zeros(T, N*n)
    V′ = reshape(V, (N, n))

    # iterate of types of extended Hubbard interactions
    for i in 1:n
        # iterate over unit cells
        for u in 1:N
            V′[u,i] = V_mean[i]
        end
    end

    return ExtendedHubbardParameters{T}(V, neighbor_table, V_bond_ids)
end