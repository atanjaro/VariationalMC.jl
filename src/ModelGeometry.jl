@doc raw"""

    ModelGeometry( unit_cell::UnitCell, 
                   lattice::Lattice, 
                   bond::Vector{Vector{Any}} )

A type defining model geometry.

"""
struct ModelGeometry
    # unit cell
    unit_cell::UnitCell

    # extent of the lattice
    lattice::Lattice

    # lattice bonds
    bond::Vector{Vector{Any}}
end


@doc raw"""

    x( i::Int, model_geometry::ModelGeometry )

Convenience function for obtaining the x-coordinate of a lattice site given a 
lattice spindex.

"""
function x(i::Int, model_geometry::ModelGeometry)
    L = model_geometry.lattice.L[1]
    return i % L
end


@doc raw"""

    y( i::Int, model_geometry::ModelGeometry )

Convenience function for obtaining the y-coordinate of a lattice site given a 
lattice spindex.

"""
function y(i::Int, model_geometry::ModelGeometry)
    L = model_geometry.lattice.L[1]
    return div(i, L)
end


@doc raw"""

    d( p1::Int, p2::Int, model_geometry::ModelGeometry )

Given lattice indices i and j, returns the distances between those 2 points, accounting 
for the latticed edges with different boundary conditions.

"""
function d(p₁::Int, p₂::Int, model_geometry::ModelGeometry)
    L = model_geometry.lattice.L[1]
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

    reduce_index_2d( i::Int, j::Int, model_geometry::ModelGeometry )

Reduces the indices of 2 lattice sites (i,j) to irreducible indices (0,k), where k is an integer.

"""
function reduce_index_2d(i::Int, j::Int, model_geometry::ModelGeometry)
    L = model_geometry.lattice.L[1]

    dx = abs(d(x(i, model_geometry), x(j, model_geometry), model_geometry))
        dy = abs(d(y(i, model_geometry), y(j, model_geometry), model_geometry))

        if dy > dx
            dx, dy = dy, dx
        end
    
        return dx + L * dy
end


@doc raw"""

    reduce_index_1d( i::Int, j::Int, model_geometry::ModelGeometry )

Reduces the indices of 2 lattice sites (i,j) to irreducible indices (0,k), where k is an integer.

"""
function reduce_index_1d(i::Int, j::Int, model_geometry::ModelGeometry)
    L = model_geometry.lattice.L[1]

    dx = abs(d(x(i, model_geometry), x(j, model_geometry), model_geometry))

    return dx
end


@doc raw"""

    max_dist( N::Int, L::Int )

Obtains the maximum irreducible index given the total number of sites N and extent of the lattice L.

"""
function max_dist(N::Int, L::Int)
    if L % 2 == 0
        return Int(N / 2 + L / 2)
    else
        return Int(N / 2)
    end
end
