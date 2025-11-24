@doc raw"""

    ModelGeometry{T, B<:AbstractVector{<:AbstractVector{T}}}

A type defining the geometry of a lattice model, including the unit cell,
lattice, and bond definitions.

- `unit_cell::UnitCell`: contains unit cell definitions.
- `lattice::Lattice`: contains lattice definitions including boundary conditions.
- `bond::B`: vector of all bond definitions.

"""
struct ModelGeometry{T, B<:AbstractVector{<:AbstractVector{T}}}
    # unit cell
    unit_cell::UnitCell

    # extent of the lattice
    lattice::Lattice

    # lattice bonds
    bond::B
end


@doc raw"""

    apply_twist_angle!( H_t::Matrix{T},
                        θ_twist::E,
                        model_geometry::ModelGeometry ) where {T<:Number, E<:AbstractFloat}

Applies a twist angle to the hopping matrix by multiplying by an appropriate phase factor.

"""
function apply_twist_angle!(
    H_t::Matrix{T},
    θ_twist::E,
    model_geometry::ModelGeometry
) where {T<:Number, E<:AbstractFloat}
    θ_x, θ_y = θ_twist, θ_twist

    # dimensions
    Lx = model_geometry.lattice.L[1]
    Ly = model_geometry.lattice.L[2]

    # apply the twist in the x-direction 
    for y in 1:Lx
        idx1 = Lx * (y - 1) + 1  
        idx2 = Lx * y            

        # spin-up sector
        H_t[idx1, idx2] *= cis(θ_x) 
        H_t[idx2, idx1] *= cis(-θ_x) 

        # spin-down sector
        H_t[idx1 + N, idx2 + N] *= cis(-θ_x) 
        H_t[idx2 + N, idx1 + N] *= cis(θ_x) 
    end

    # apply the twist in the y-direction 
    for x in 1:Ly
        idx1 = x                  
        idx2 = Ly * (Ly - 1) + x  

        # spin-up sector
        H_t[idx1, idx2] *= cis(θ_y) 
        H_t[idx2, idx1] *= cis(-θ_y) 

        # spin-down sector
        H_t[idx1 + N, idx2 + N] *= cis(-θ_y) 
        H_t[idx2 + N, idx1 + N] *= cis(θ_y) 
    end

    return nothing
end


@doc raw"""

    x( i::I, 
       model_geometry::ModelGeometry ) where {I<:Integer}

Convenience function for obtaining the x-coordinate of a lattice site given a 
lattice spindex.

- `i::I`: spindex.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function x(
    i::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
    L = model_geometry.lattice.L[1]

    return i % L
end


@doc raw"""

    y( i::I, 
       model_geometry::ModelGeometry ) where {I<:Integer}

Convenience function for obtaining the y-coordinate of a lattice site given a 
lattice spindex.

- `i::I`: spindex.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function y(
    i::I, model_geometry::ModelGeometry
) where {I<:Integer}
    L = model_geometry.lattice.L[1]

    return div(i, L)
end


@doc raw"""

    d( p1::I, 
       p2::I, 
       model_geometry::ModelGeometry ) where {I<:Integer}

Given lattice points ``p\_1`` and ``p_2``, returns the distance between those two points, 
accounting for the lattice edges under different boundary conditions.

- `p1::I`: first lattice point.
- `p2::I`: second lattice point.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function d(
    p₁::I, 
    p₂::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
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

    reduce_index_2d( i::I, 
                     j::I, 
                     model_geometry::ModelGeometry ) where {I<:Integer}

For two lattice sites ``i`` and ``j`` on a 2D lattice, returns the irreducible index
``k`` between them, where ``k`` is an integer.

- `i::I`: first lattice site.
- `j::I`: second lattice site.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function reduce_index_2d(
    i::I, 
    j::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
    L = model_geometry.lattice.L[1]

    dx = abs(d(x(i, model_geometry), x(j, model_geometry), model_geometry))
    dy = abs(d(y(i, model_geometry), y(j, model_geometry), model_geometry))

    if dy > dx
        dx, dy = dy, dx
    end
    
    return dx + L * dy
end


@doc raw"""

    reduce_index_1d( i::I, 
                     j::I, 
                     model_geometry::ModelGeometry ) where {I<:Integer}

For two lattice sites ``i`` and ``j`` on a 1D lattice, returns the irreducible index
``k`` between them, where ``k`` is an integer.

- `i::I`: first lattice site.
- `j::I`: second lattice site.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function reduce_index_1d(
    i::I, 
    j::I, 
    model_geometry::ModelGeometry
) where {I<:Integer}
    L = model_geometry.lattice.L[1]

    dx = abs(d(x(i, model_geometry), x(j, model_geometry), model_geometry))

    return dx
end


@doc raw"""

    max_dist( N::I, 
              L::I ) where {I<:Integer}

For a lattice with ``N`` sites and extent ``L``, returns the maximum irreducible index ``k_{\mathrm{max}}``.

- `N::I`: total number of lattice sites.
- `L::I`: extent of the lattice.

"""
function max_dist(
    N::I, 
    L::I
) where {I<:Integer}
    if L % 2 == 0
        return Int(N / 2 + L / 2)
    else
        return Int(N / 2)
    end
end
