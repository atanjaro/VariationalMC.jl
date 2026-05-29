@doc raw"""

    HubbardModel{T<:AbstractFloat}

A type to represent a multi-orbital Hubbard model; specifically, the particle-hole asymmetric form for the Hubbard interaction
```math
\hat{H}_{U}=\sum_{\mathbf{i},\nu}U_{\nu,\mathbf{i}}\hat{n}_{\uparrow,\nu,\mathbf{i}}\hat{n}_{\downarrow,\nu,\mathbf{i}},
```
where ``\mathbf{i}`` specifies the unit cell, and ``\nu`` denotes the orbital in the unit cell.

# Fields

- `U_orbital_ids::Vector{Int}`: Orbital species/IDs in unit cell with finite Hubbard interaction.
- `U_mean::Vector{T}`: Average Hubbard interaction strength ``U_\nu`` for a given orbital species in the lattice.

"""
struct HubbardModel{T<:AbstractFloat}
    # orbital species
    U_orbital_ids::Vector{Int}

    # average Hubbard U
    U_mean::Vector{T}
end


@doc raw"""

    HubbardModel(;
        # KEYWORD ARGUMENTS
        U_orbital::AbstractVector{Int},
        U_mean::AbstractVector{T}
    ) where {T<:AbstractFloat}

Initialize and return an instance of the type `HubbardModel`.

# Keyword Arguments

- `U_orbital::Vector{Int}`: Orbital species/IDs in unit cell with finite Hubbard interaction.
- `U_mean::Vector{T}`: Average Hubbard interaction strength ``U_\nu`` for a given orbital species in the lattice.

"""
function HubbardModel(;
    # KEYWORD ARGUMENTS
    U_orbital::AbstractVector{Int},
    U_mean::AbstractVector{T}
) where {T<:AbstractFloat}
    
    return HubbardModel(U_orbital, U_mean)
end

# show struct info as TOML formatted string
function Base.show(io::IO, ::MIME"text/plain", hm::HubbardModel)

    (; U_orbital_ids, U_mean) = hm

    @printf io "[HubbardModel]\n\n"
    @printf io "HUBBARD_IDS = %s\n" string(collect(1:length(U_orbital_ids)))
    @printf io "ORBITAL_IDS = %s\n" string(U_orbital_ids)
    @printf io "U_mean      = %s\n" string(round.(U_mean, digits=6))

    return nothing
end


@doc raw"""

    ExtendedHubbardModel{T<:AbstractFloat}

A type to represent extended Hubbard interactions; specifically, the particle-hole asymmetric form of the extended Hubbard interaction
```math
\begin{align*}
\hat{H}_{V} = \sum_{\mathbf{j},\mathbf{r},\nu,\eta}V_{(\mathbf{j}+\mathbf{r},\nu),(\mathbf{j},\eta)} & \hat{n}_{\mathbf{j}+\mathbf{r},\nu}\hat{n}_{\mathbf{j},\eta} \\
    = \sum_{\mathbf{j},\mathbf{r},\nu,\eta}V_{(\mathbf{j}+\mathbf{r},\nu),(\mathbf{j},\eta)} & \bigg[\tfrac{1}{2}(\hat{n}_{\mathbf{j}+\mathbf{r},\nu}+\hat{n}_{\mathbf{j},\eta}-2)^{2}-1 \\
    & -\hat{n}_{\mathbf{j}+\mathbf{r},\nu,\uparrow}\hat{n}_{\mathbf{j}+\mathbf{r},\nu\downarrow}-\hat{n}_{\mathbf{j},\eta,\uparrow}\hat{n}_{\mathbf{j},\eta\downarrow}+\tfrac{3}{2}\hat{n}_{\mathbf{j}+\mathbf{r},\nu}+\tfrac{3}{2}\hat{n}_{\mathbf{j},\eta}\bigg],
\end{align*}
```
where ``\mathbf{j}`` specifies a unit cell in the lattice, ``\mathbf{r}`` is a displacement in units, and ``\nu`` and ``\eta`` specify the orbital in a given unit cell.
Here, ``\hat{n}_{\mathbf{j},\eta} = (\hat{n}_{\uparrow,\mathbf{j},\eta} + \hat{n}_{\downarrow,\mathbf{j},\eta})`` is the electron number operator for orbital
``\eta`` in unit cell ``\mathbf{j}`` in the lattice. ``V_{(\mathbf{j}+\mathbf{r},\nu),(\mathbf{j},\eta)}`` controls the strength of the
extended Hubbard interaction between orbital ``\eta`` in unit cell ``\mathbf{j}`` and orbital ``\nu`` in unit cell ``\mathbf{j}+\mathbf{r}``.

# Fields

- `V_bond_ids::Vector{Int}`: Bond IDs specifying bond definition that separates a pair of orbitals with an extended Hubbard interaction between them.
- `V_mean::Vector{T}`: Average extended Hubbard interaction strength ``V_{(\mathbf{j}+\mathbf{r},\nu),(\mathbf{j},\eta)}`` associated with bond definition.

"""
struct ExtendedHubbardModel{T<:AbstractFloat}
    # bond IDs
    V_bond_ids::Vector{Int}

    # average extend Hubbard interaction V
    V_mean::Vector{T}
end


@doc raw"""

    ExtendedHubbardModel(;
        # KEYWORD ARGUMENTS
        model_geometry::ModelGeometry{D,T},
        V_bond::Vector{Bond{D}},
        V_mean::Vector{T},
    ) where {T<:AbstractFloat, D}

Initialize and return an instance of the type `ExtendedHubbardModel`.

"""
function ExtendedHubbardModel(;
    # KEYWORD ARGUMENTS
    model_geometry::ModelGeometry{D,T},
    V_bond::Vector{Bond{D}},
    V_mean::Vector{T}
) where {T<:AbstractFloat, D}

    V_bond_ids = [add_bond!(model_geometry, bond) for bond in V_bond]

    return ExtendedHubbardModel{T}( V_bond_ids, V_mean)
end

# show struct info as TOML formatted string
function Base.show(io::IO, ::MIME"text/plain", ehm::ExtendedHubbardModel)

    (; V_bond_ids, V_mean) = ehm
    @printf io "[ExtendedHubbardModel]\n\n"
    @printf io "EXT_HUB_IDS = %s\n" string(collect(1:length(V_bond_ids)))
    @printf io "BOND_IDS    = %s\n" string(V_bond_ids)
    @printf io "V_mean      = %s\n" string(round.(V_mean, digits=6))

    return nothing
end
