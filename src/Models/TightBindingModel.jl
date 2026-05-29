@doc raw"""

    TightBindingModel{T<:Number, E<:AbstractFloat, D}

Defines a non-interacting tight binding model in `D` dimensions. 

# Fields

- `μ::E`: Initial value of the chemical potential.
- `t_bond_ids::Vector{Int}`: The bond ID for each bond/hopping definition.
- `t_bonds::Vector{Bond{D}}`: Bond definition for each type of hopping in the tight binding model.
- `t_mean::Vector{T}`: Mean hopping energy for each type of hopping.
- `η::Vector{E}`: Relative twist angle for each lattice vector direction ``d \in [1, D]`` such that ``\eta_d \in [0, 1).``
- `expniϕ::Vector{T}`: Twist angle phase ``\exp(i n_d \phi_d)`` for each hopping.

"""
struct TightBindingModel{T<:Number, E<:AbstractFloat,  D}
    # initial chemical potential
    μ::E

    # bond ID associated with each hopping in tight-binding model
    t_bond_ids::Vector{Int}

    # bond definition associated with each hopping
    t_bonds::Vector{Bond{D}}

    # mean hopping energy for each type of hopping in tight-binding model
    t_mean::Vector{T}

    # twist angle phase
    expniϕ::Vector{T}

    # relative twist angle
    η::Vector{E}
end


@doc raw"""

    TightBindingModel(;
        # KEYWORD ARGUMENTS
        model_geometry::ModelGeometry{D,E,N},
        μ::E,
        t_bonds::Vector{Bond{D}} = Bond{ndims(model_geometry.unit_cell)}[],
        t_mean::Vector{T} = T[],
        η::Union{Vector{E},Nothing} = nothing
    ) where {T<:Number, E<:AbstractFloat, D, N}

Initialize and return an instance of `TightBindingModel`, also adding/recording the bond defintions `t_bonds` to the
`ModelGeometry` instance `model_geometry`.

"""
function TightBindingModel(;
    # KEYWORD ARGUMENTS
    model_geometry::ModelGeometry{D,E,N},
    μ::E,
    t_bonds::Vector{Bond{D}} = Bond{ndims(model_geometry.unit_cell)}[],
    t_mean::Vector{T} = T[],
    η::Union{Vector{E},Nothing} = nothing
) where {T<:Number, E<:AbstractFloat, D, N}

    # get the number of orbitals per unit cell
    (; bonds, unit_cell, lattice) = model_geometry
    (; n, reciprocal_vecs) = unit_cell
    (; L) = lattice

    # record bonds associated with hoppings
    t_bond_ids = zeros(Int, length(t_bonds))
    for i in eachindex(t_bonds)
        t_bond_ids[i] = add_bond!(model_geometry, t_bonds[i])
    end

    # set default twist angle to zero
    η = isnothing(η) ? zeros(E, D) : η

    # get the type of the hopping energy
    H = all(i->iszero(i), η) ? T : Complex{E}

    # set twist angle phase to unity if twist angles are all zero
    if H<:Real
        expniϕ = ones(H, length(t_bonds))
    # calculate the twist angle phase for each hopping
    else
        expniϕ = ones(H, length(t_bonds))
        # iterate over hoppings
        for h in eachindex(t_bonds)
            # calculate the displacement vector associated the bond
            R = bond_to_vec(t_bonds[h], unit_cell)
            # iterate of reciprocal lattice vectors
            for d in 1:D
                # get the current reciprocal lattice vector
                Kd = @view reciprocal_vecs[:,d]
                # update the twist angle phase
                expniϕ[h] = exp(-1im * η[d] * dot(Kd,R) / L[d]) * expniϕ[h]
            end
        end
    end

    return TightBindingModel{H,E,D}(μ, t_bond_ids, t_bonds, Vector{H}(t_mean), expniϕ, η)
end


# show struct info as TOML formatted string for real hopping energies
function Base.show(io::IO, ::MIME"text/plain", tbm::TightBindingModel{T,E,D}; spin::Int=0) where {T<:AbstractFloat,E,D}

    if iszero(spin)
        @printf io "[TightBindingModel]\n\n"
    elseif isone(spin)
        @printf io "[TightBindingModelUp]\n\n"
    else
        @printf io "[TightBindingModelDown]\n\n"
    end
    @printf io "relative_twist_angle = %s\n"  string(round.(tbm.η, digits=6))
    @printf io "initial chemical_potential = %.8f\n\n" tbm.μ
    for i in eachindex(tbm.t_bonds)
        if iszero(spin)
            @printf io "[[TightBindingModel.hopping]]\n\n"
        elseif isone(spin)
            @printf io "[[TightBindingModelUp.hopping]]\n\n"
        else
            @printf io "[[TightBindingModelDown.hopping]]\n\n"
        end
        @printf io "HOPPING_ID   = %d\n" i
        @printf io "BOND_ID      = %d\n" tbm.t_bond_ids[i]
        @printf io "orbitals     = [%d, %d]\n" tbm.t_bonds[i].orbitals[1] tbm.t_bonds[i].orbitals[2]
        @printf io "displacement = %s\n" string(tbm.t_bonds[i].displacement)
        @printf io "t_mean       = %.8f\n" tbm.t_mean[i]
    end

    return nothing
end