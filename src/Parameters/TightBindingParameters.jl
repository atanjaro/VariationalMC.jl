@doc raw"""

    TightBindingParameters{T<:Number, E<:AbstractFloat}

A mutable struct containing all the parameters needed to characterize a finite tight-binding Hamiltonian
for a single spin species ``\sigma`` on a finite lattice with periodic boundary conditions of the form
```math
\hat{H}_{0,\sigma}=-\sum_{\langle i,j\rangle}(t_{ij} \hat{c}_{\sigma,i}^{\dagger}\hat{c}_{\sigma,j}+\textrm{h.c.}),
```
where ``\hat{c}_{\sigma,i}^\dagger`` is the fermion creation operator for an electron with spin ``\sigma`` on orbital ``i,``
``t_{i,j}`` are the hopping energies, and ``\mu`` is the chemical potential.

# Fields

- `μ::E`: The chemical potential ``\mu.``
- `const t::Vector{T}`: The hopping energy ``t_{i,j}`` associated with each pair of neighboring orbitals connected by a bond in the lattice.
- `const neighbor_table::Matrix{Int}`: Neighbor table containing all pairs of orbitals in the lattices connected by a bond, with a non-zero hopping energy between them.
- `const neighbor_map::Dict`: Dictionary containing all sites in the lattice, mapping them to their valid neighboring sites and the bonds between them.
- `const bond_ids::Vector{Int}`: The bond ID definitions that define the types of hopping in the lattice.
- `const bond_slices::Vector{UnitRange{Int}}`: Slices of `neighbor_table` corresponding to given bond ID i.e. the neighbors `neighbor_table[:,bond_slices[i]]` corresponds the `bond_ids[i]` bond defintion.
- `const norbital::Int`: Number of orbitals per unit cell.

"""
mutable struct TightBindingParameters{T<:Number, E<:AbstractFloat}
    # chemical potential
    μ::E

    # hopping energies for all pairs of orbitals connected by a bond in the lattice
    const t::Vector{T}

    # neighbor table for all pairs of orbitals connected by a bond in the lattice
    const neighbor_table::Matrix{Int}

    # maps neighbor table to a dictionary of sites and their associated neighbors and bonds
    const neighbor_map::Dict

    # bond IDs that define hoppings
    const bond_ids::Vector{Int}

    # view into neighbor table for each bond ID
    const bond_slices::Vector{UnitRange{Int}}

    # number of orbitals per unit cell
    const norbital::Int
end


@doc raw"""

    TightBindingParameters(;
        # KEYWORD ARGUMENTS
        tight_binding_model::TightBindingModel{T,E,D},
        particle_configuration::ParticleConfiguration{E,I},
        model_geometry::ModelGeometry{D,E}
    ) where {T,E,D}

Initialize and return an instance of `TightBindingParameters`.

"""
function TightBindingParameters(;
    # KEYWORD ARGUMENTS
    tight_binding_model::TightBindingModel{T,E,D},
    particle_configuration::ParticleConfiguration,
    model_geometry::ModelGeometry{D,E}
) where {T,E,D}
    (; unit_cell, lattice) = model_geometry
    N = lattice.N # number of unit cells in lattice
    n = unit_cell.n # number of orbital per unit cell

    # construct neighbor table for hoppings
    t_bonds = tight_binding_model.t_bonds::Vector{Bond{D}}
    if length(t_bonds) > 0
        neighbor_table = build_neighbor_table(t_bonds, unit_cell, lattice)
    else
        neighbor_table = Matrix{Int}(undef,2,0)
    end

    # constuct the neighbor map for all sites
    neighbor_map = map_neighbor_table(neighbor_table)

    # get number of bond definitions in model
    nbonds = length(tight_binding_model.t_bonds)

    # get total number of bonds in lattice
    Nbonds = size(neighbor_table, 2)

    # get th slice of bonds in the neighbor_table associated with
    # hopping ID (and corresponding bond ID)
    bond_ids = copy(tight_binding_model.t_bond_ids)
    bond_slices = UnitRange{Int}[]
    offset = 0
    for bond in t_bonds
        nt = build_neighbor_table([bond], unit_cell, lattice)
        n_this = size(nt, 2)
        push!(bond_slices, offset+1:offset+n_this)
        offset += n_this
    end

    # set hopping energy for each bond in lattice
    t = zeros(T, Nbonds)
    if Nbonds > 0
        bonds_per_type = Nbonds ÷ nbonds
        t′ = reshape(t, (bonds_per_type, nbonds))
        for b in axes(t′,2)
            t_mean = tight_binding_model.t_mean[b]
            expniϕ = tight_binding_model.expniϕ[b]
            t_b = @view t′[:,b]
            @. t_b = expniϕ  * t_mean
        end
    end

    # set exact chemical potential
    μ = calculate_chemical_potential(    
        N,
        n,
        particle_configuration.Ne, 
        bond_slices,
        neighbor_table,
        t
    )

    return TightBindingParameters(μ, t, neighbor_table, neighbor_map, bond_ids, bond_slices, n)
end