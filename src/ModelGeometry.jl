@doc raw"""
    ModelGeometry{D, T<:AbstractFloat, N}

Contains all the information defining the lattice geometry for the model in `D` spatial dimensions.

# Comment

The bond ID associated with a `bond::Bond{D}` corresponds to the index associated with it into the `bonds` vector field.

# Fields

- `unit_cell::UnitCell{D,T,N}`: Defines unit cell.
- `lattice::Lattice{D}`: Defines finite lattice extent.
- `bonds::Vector{Bond{D}}`: All available bond definitions in simulation, with vector indices giving the bond ID.

"""
struct ModelGeometry{D, T<:AbstractFloat, N}
    # unit cell
    unit_cell::UnitCell{D,T,N}

    # lattice
    lattice::Lattice{D}

    # bonds
    bonds::Vector{Bond{D}}
end

@doc raw"""

    ModelGeometry(  unit_cell::UnitCell, 
                    lattice::Lattice ) where {D}

Initialize and return a [`ModelGeometry`](@ref) instance. Defines a "trivial" bond definition for each
orbital in the unit cell that connects an orbital to itself.
"""
function ModelGeometry(unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}
    # define trivial bond connecting each orbital in unit cell to itself
    n     = unit_cell.n
    bonds = Bond{D}[]
    for i in 1:n
        push!(bonds, Bond((i,i),zeros(Int,D)))
    end

    return ModelGeometry(unit_cell, lattice, bonds)
end


# print struct info in TOML format
function Base.show(io::IO, ::MIME"text/plain", model_geo::ModelGeometry{D,T}) where {D,T}

    (; unit_cell, lattice, bonds) = model_geo

    @printf io "[Geometry]\n\n"
    @printf io "dimensions = %d\n\n" D
    @printf io "[Geometry.UnitCell]\n\n"
    @printf io "orbitals = %d\n\n" unit_cell.n
    @printf io "[Geometry.UnitCell.LatticeVectors]\n\n"
    for d in 1:D
        a = @view unit_cell.lattice_vecs[:,d]
        @printf io "a_%d = %s\n" d string(round.(a, digits=6)) 
    end
    @printf io "\n"
    @printf io "[Geometry.UnitCell.ReciprocalVectors]\n\n"
    for d in 1:D
        b = @view unit_cell.reciprocal_vecs[:,d]
        @printf io "b_%d = %s\n" d string(round.(b, digits=6)) 
    end
    @printf io "\n"
    for i in 1:unit_cell.n
        r = unit_cell.basis_vecs[i]
        @printf io "[[Geometry.UnitCell.BasisVectors]]\n\n"
        @printf io "ORBITAL_ID = %d\n" i
        @printf io "r          = %s\n\n" string(round.(r, digits=6))
    end
    @printf io "\n"
    @printf io "[Geometry.Lattice]\n\n"
    @printf io "L        = %s\n" string(lattice.L)
    @printf io "periodic = [%s]\n\n" join(lattice.periodic, ", ")
    for i in eachindex(bonds)
        @printf io "[[Geometry.Bond]]\n\n"
        @printf io "BOND_ID      = %d\n" i
        @printf io "orbitals     = [%d, %d]\n" bonds[i].orbitals[1] bonds[i].orbitals[2]
        @printf io "displacement = %s\n\n" string(bonds[i].displacement)
    end

    return nothing
end


@doc raw"""

    add_bond!(  model_geometry::ModelGeometry{D,T}, 
                bond::Bond{D}) where {D, T}   

Add `bond` definition to `model_geometry`, returning the bond ID i.e. the index to `bond`
in the vector `model_geometry.bonds`.
This method first checks that `bond` is not already defined. If it is this method simply
returns the corresponding bond ID. If `bond` is not already defined, then it is appended
to the vector `model_geometry.bonds`.

"""
function add_bond!(model_geometry::ModelGeometry{D,T}, bond::Bond{D}) where {D, T}

    (; bonds) = model_geometry

    # get the bond ID
    bond_id = get_bond_id(model_geometry, bond)

    # if the bond is not already recorded, then record it and get its new bond ID
    if iszero(bond_id)
        # record the bond ID
        push!(bonds, bond)
        # get the ID of the new bond
        bond_id = length(bonds)
    end

    return bond_id
end


@doc raw"""

    get_bond_id(model_geometry::ModelGeometry{D,T}, bond::Bond{D}) where {D, T}

Return the bond ID associated with the bond defintion `bond`, returning `bond_id=0`
if the it is not a recorded bond.
    
"""
function get_bond_id(model_geometry::ModelGeometry{D,T}, bond::Bond{D}) where {D, T}

    (; bonds) = model_geometry

    if bond in bonds
        bond_id = findfirst(b -> b==bond, bonds)
    else
        bond_id = 0
    end
    
    return bond_id
end


@doc raw"""

    x(i::Int, lattice::Lattice) 

Convenience function for obtaining the ``x``-coordinate of a lattice site given a 
lattice spindex.

"""
function x(i::Int, lattice::Lattice) 
    L = lattice.L[1]
    return i % L
end


@doc raw"""

    y(i::Int, lattice::Lattice) 

Convenience function for obtaining the ``y``-coordinate of a lattice site given a 
lattice spindex.

"""
function y(i::Int, lattice::Lattice) 
    L = lattice.L[1]
    return div(i, L)
end


@doc raw"""

    z(i::Int, lattice::Lattice)

Convenience function for obtaining the z-coordinate of a lattice site given a 
lattice spindex.

"""
function z(i::Int, lattice::Lattice) 
    L = lattice.L[3]
    return div(i, L)
end

function coord(
    i::I,
    d::Int,
    lattice::Lattice
) where {I<:Integer}
    if d == 1
        return x(i, lattice)
    elseif d == 2
        return y(i, lattice)
    elseif d == 3
        return z(i, lattice)
    end
end


@doc raw"""

    reduce_index(i::Int, j::Int, lattice::Lattice) 

For two lattice sites ``i`` and ``j`` on a 1D, 2D, or 3D lattice, returns the 
irreducible index ``k`` between them, where ``k`` is an integer. 

# Note

The displacement components are sorted in descending order before index reduction to ensure that 
symmetry-equivalent pairs map to the same index.

"""
function reduce_index(i::Int, j::Int, lattice::Lattice)
    L   = lattice.L
    dim = length(L)

    # compute absolute displacements along each dimension
    disp = ntuple(dim) do n
        abs(d(coord(i, n, lattice), coord(j, n, lattice), lattice))
    end

    # sort components in descending order so symmetry-equivalent pairs share an index
    disp = Tuple(sort(collect(disp), rev=true))

    # reduce to scalar index using mixed-radix encoding
    strides = cumprod([1; collect(L[1:end-1])])
    return sum(disp[n] * strides[n] for n in 1:dim)
end



@doc raw"""

    d(p₁::Int, p₂::Int, lattice::Lattice) 

Given lattice points ``p\_1`` and ``p_2``, returns the distance between those two points, 
accounting for the lattice edges under different boundary conditions.

"""
function d(p₁::Int, p₂::Int, lattice::Lattice)
    L = lattice.L[1]
    dist = p₂ - p₁

    if dist >div(L, 2)
        dist -= L
    end
    if dist < -div(L,2)
        dist+= L
    end

    return dist
end


@doc raw"""

    max_dist( N::I, L::I ) where {I<:Integer}

For a lattice with ``N`` sites and extent ``L``, returns the maximum irreducible index ``k_{\mathrm{max}}``.

"""
function max_dist(N::Int, L::Int) 
    if L % 2 == 0
        return Int(N / 2 + L / 2)
    else
        return Int(N / 2)
    end
end

@doc raw"""

    map_spindex(model_geometry::ModelGeometry)

Creates a mapping between all spindices (spin-indices) to their real lattice site
coordinates. Works for 1D, 2D, and 3D geometries; dimensionality is inferred
automatically from the lattice.

"""
function map_spindex(model_geometry::ModelGeometry)
    (; unit_cell, lattice) = model_geometry
    Norbs = unit_cell.n
    Ncells = lattice.N
    N = Norbs * Ncells
    dims = length(lattice.L)

    locs = Matrix{Int}(undef, dims, 2*N)

    @inbounds for spin in 0:1
        for orb in 1:Norbs
            for cell in 1:Ncells
                # reconstruct spindex in desired orbital-major order
                s_out = spin * N + (orb - 1) * Ncells + cell
                # original spindex as laid out by get_index_from_spindex
                idx = (cell - 1) * Norbs + orb + spin * N
                loc = site_to_loc(idx, unit_cell, lattice)
                for d in 1:dims
                    locs[d, s_out] = loc[1][d]
                end
            end
        end
    end

    return locs
end


@doc raw"""

    get_spindex_type(spindex::Int, N::Int) 

Returns the spin species at a given spindex.

"""
function get_spindex_type(spindex::Int, N::Int) 
    @assert spindex < 2*N + 1

    return spindex < N + 1 ? 1 : -1
end


@doc raw"""

    get_index_from_spindex(spindex::Int, N::Int)

Returns the lattice site ``i`` for a given spindex.

"""
function get_index_from_spindex(spindex::Int, N::Int)
    @assert spindex < 2 * N + 1

    return spindex <= N ? spindex : spindex - N
end


@doc raw"""

    get_spindices_from_index(index::Int, N::Int) 

Returns spin-up and spin-down indices from a given site index.

"""
function get_spindices_from_index(index::Int, N::Int) 
    @assert index <= N

    return index, index + N
end


@doc raw"""

    get_linked_spindex(i::Int, N::Int) 

 Given an index ``i`` in the spin-up sector, returns an index in the spin-down sector.

"""
function get_linked_spindex(i::Int, N::Int) 
    @assert i < 2 * N

    return i + (1 - 2 * (i ÷ N)) * N
end
