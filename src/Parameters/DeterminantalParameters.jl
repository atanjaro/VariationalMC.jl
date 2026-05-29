@doc raw"""

    DeterminantalParameters{S<:AbstractString, T<:Number}

Determinantal parameters allowing for different types of ordering. This include charge order in the form of `"μ"` 
and ``"density"``, spin order in the form of `"spin-x"` and `"spin-x"`, and pairing in the form of `"s-wave"` 
and `"d-wave"`.

# Fields

- `order_type::Vector{S}`: The class of each order parameter, being `"charge"`, `"spin"`, or `"pair"`.
- `param_name::Vector{S}`: Specific order parameters names.
- `symmetry::Vector{S}`: Symmetry of each order parameter, being `"uniform"` or `"site-dependent"`.
- `p::Vector{Vector{T}}`: Vector of order parameters.
- `qₚ::Vector{AbstractFloat}`: Center-of-mass pairing momentum when introducing `site-dependent` pairing. 
- `optimize::Vector{Bool}`: Determines whether an order parameter will be optimized.
- `orbital::Vector{Int}`: All orbital species where each order parameter applies.

"""
mutable struct DeterminantalParameters{S<:AbstractString, T<:Number}
    # ordering type
    order_type::Vector{S}

    # parameter name
    param_name::Vector{S}

    # parameter symmetry
    symmetry::Vector{S}

    # determinantal parameter values
    p::Vector{Vector{T}}

    # pairing momentum
    qₚ::Union{Vector{AbstractFloat}, Nothing}

    # optimization flags
    optimize::Vector{Bool}

    # spindex table
    spdx_table::Matrix{Int}

    # orbital IDs
    orbital::Vector{Vector{Int}}
end


@doc raw"""

    DeterminantalParameters(;
        tight_binding_parameters::TightBindingParameters{T},
        model_geometry::ModelGeometry{D,T},
        rng::AbstractRNG
    ) where {D, T<:AbstractFloat}

Initialize and return an instance of `DeterminantalParameters'.

"""
function DeterminantalParameters(;
    tight_binding_parameters::TightBindingParameters{T, E},
    model_geometry::ModelGeometry{D,E},
    rng::AbstractRNG
) where {D, T<:Number, E<:AbstractFloat}
    # these represent the minimal order parameters for any model in arbitrary dimensions
    param_name  = ["μ", "density", "spin-x", "spin-z"]
    order_type  = ["charge", "charge", "spin", "spin"]
    symmetry    = ["uniform", "uniform", "uniform", "uniform"]
    p           = [[tight_binding_parameters.μ], [1e-4], [1e-4], [1e-4]]

    # orbitals for which the symmetries will apply
    Norbitals = model_geometry.unit_cell.n
    orbital  = Vector{Int}[]

    for par in eachindex(param_name)
        par_orb = zeros(Int,Norbitals)
        for o in 1:Norbitals
            par_orb[o] = o
        end
        push!(orbital, par_orb)
    end

    # create spindex table for building parameter matrices
    spdx_table = map_spindex(model_geometry)

    # initialize as static parameters
    optimize = Bool[false, false, false, false]

    # no pairing present in the system yet
    qₚ = nothing

    return DeterminantalParameters(order_type, param_name, symmetry, p, qₚ, optimize, spdx_table, orbital)
end


# print struct info in TOML format
function Base.show(io::IO, ::MIME"text/plain", det_par::DeterminantalParameters{T,E}) where {T,E}

    (; order_type, param_name, symmetry, p, orbital) = det_par

    println(io, "[DeterminantalParameters]\n")

    # preserve order instead of sorting
    tracking = Set{String}()
    ordered_types = String[]

    for ot in order_type
        if !(ot in tracking)
            push!(tracking, ot)
            push!(ordered_types, ot)
        end
    end

    for ord in ordered_types
        println(io, "[[DeterminantalParameters.$ord]]")

        inds = findall(==(ord), order_type)

        for i in inds
            name = param_name[i]
            sym  = symmetry[i]
            orb  = orbital[i]
            mean_val = p[i]

            println(io, "ORDER       = \"", name, "\"")
            println(io, "SYMMETRY    = \"", sym, "\"")
            println(io, "ORBITAL_IDS = ", repr(orb))
            println(io, "mean        = ", repr(mean_val))
            println(io)
        end
    end
end


@doc raw"""

    add_parameter!(
        determinantal_parameters::DeterminantalParameters;
        param_name::Union{S, Nothing} = nothing,
        order_type::Union{S, Nothing} = nothing,
        symmetry::Union{S, Nothing} = nothing,
        optimize::Union{Bool, Nothing} = nothing,
        orbital::Union{Vector{Int}, Nothing} = nothing,
        p::Union{Vector{AbstractFloat}, Nothing} = nothing,
        model_geometry::Union{ModelGeometry, Nothing} = nothing,
        rng::Union{AbstractRNG, Nothing} = nothing
    ) where {S<:AbstractString}  

Adds a parameters definition to `determinantal_parameters`.

This method first checks that the given parameter is not already defined. If it is, this method simply
flags the parameter with that name for optimization. If the given parameter is not already defined, then 
it is added to the definitons in `determinantal_parameters`.

"""
function add_parameter!(
    determinantal_parameters::DeterminantalParameters;
    param_name::Union{S, Nothing} = nothing,
    order_type::Union{S, Nothing} = nothing,
    symmetry::Union{S, Nothing} = nothing,
    optimize::Union{Bool, Nothing} = nothing,
    orbital::Union{Vector{Int}, Nothing} = nothing,
    p::Union{Vector{AbstractFloat}, Nothing} = nothing,
    qₚ::Union{Vector{AbstractFloat}, Nothing} = nothing,
    model_geometry::Union{ModelGeometry, Nothing} = nothing,
    rng::Union{AbstractRNG, Nothing} = nothing
) where {S<:AbstractString}

    SYMMETRY_TYPES = ("uniform", "site-dependent")

    # valid param_names are all values across all keys
    valid_params = reduce(vcat, values(DETERMINANTAL_PARAMETERS))

    @assert param_name in valid_params """
    Invalid param_name: "$param_name". Must be one of: $(join(valid_params, ", "))
    """

    @assert order_type in keys(DETERMINANTAL_PARAMETERS) """
    Invalid order_type: "$order_type". Must be one of: $(join(keys(DETERMINANTAL_PARAMETERS), ", "))
    """

    @assert symmetry in SYMMETRY_TYPES """
    Invalid symmetry: "$symmetry". Must be one of: $(join(SYMMETRY_TYPES, ", "))
    """

    @assert param_name in DETERMINANTAL_PARAMETERS[order_type] """
    Invalid combination: param_name "$param_name" is not compatible with order_type "$order_type".
    Expected one of: $(join(DETERMINANTAL_PARAMETERS[order_type], ", "))
    """

    name_idx = findfirst(==(param_name), determinantal_parameters.param_name)
    if !isnothing(name_idx) && determinantal_parameters.symmetry[name_idx] == symmetry
        # parameter already exists, just enable optimization
        determinantal_parameters.optimize[name_idx] = true
        # initial parameter value override
        if !isnothing(p)
            determinantal_parameters.p[name_idx] = p
        end
    else
        # add the order parameter identifiers
        push!(determinantal_parameters.order_type, order_type)
        push!(determinantal_parameters.param_name, param_name)
        push!(determinantal_parameters.symmetry, symmetry)

        if optimize
            # flag parameter for optimization
            push!(determinantal_parameters.optimize, true)
        else
            # initialize as a static parameter
            push!(determinantal_parameters.optimize, false)
        end

        # add orbitals
        push!(determinantal_parameters.orbital, orbital)

        if !isnothing(p)
            # manual override of default initial parameter values
            push!(determinantal_parameters.p, p)
        else
            # default value
            if symmetry == "uniform"
                push!(determinantal_parameters.p, [1e-4])
            elseif symmetry == "site-dependent" && order_type == "pair"
                N = model_geometry.lattice.N
                # Fulde-Ferrell parameters
                push!(determinantal_parameters.p, zeros(AbstractFloat, N))

                # Larkin-Ovchinnikov parameters
                push!(determinantal_parameters.p, zeros(AbstractFloat, N))
                # set pairing momenutum
                if isnothing(qₚ)
                    determinantal_parameters.qₚ = [0.0, 0.0]    # default pairing momentum
                else
                    determinantal_parameters.qₚ = qₚ
                end
            else
                Lx = model_geometry.lattice.L[1]
                push!(determinantal_parameters.p, zeros(AbstractFloat, Lx))
            end
        end
    end
    
    return nothing
end