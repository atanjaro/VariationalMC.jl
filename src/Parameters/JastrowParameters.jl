@doc raw"""

    JastrowParameters{E<:AbstractFloat}

Jastrow pseudopotentials for different types of Jastrow factor. These include density-density, spin-spin, phonon-phonon,
and electron-phonon.

# Fields

- `particle_pair::AbstractString`: Pair of particles being subject to Jastrow pseudopotentials, being `"electron-electron"`, `"phonon-phonon"`, and `"electron-phonon"`.
- `order_pair::AbstractString`: Pair of particle ordering types that are being reduced with Jastrow pseupotentials, being `"density-density"`, `"spin-spin"`, `"density-displacement"`, and `"displacement-displacement"`.
- `irr_indices::Vector{Int}`: Set of (irreducible) single indices identifying equivalence classes of lattice site index pairs. Note that irreducible indices start at 0.
- `mean_v::Vector{E}`: Strength of each pseudopotential for each irreducible index.
- `optimize::Bool`: Optimization flag for whether the pseudopotentials will be optimized. If `true`, all parameters will be optimized excpet for ones associated with the largest irreducible index.
- `orbitals::Vector{Int}`: Orbital IDs for each irreducible index.

"""
mutable struct JastrowParameters{E<:AbstractFloat}
    # name of particle pairs
    particle_pair::AbstractString

    # name of order pairs
    order_pair::AbstractString

    # irreducible site indices
    irr_indices::Vector{Vector{Int}}

    # maps each irreducible index to every site index pair which generates it
    irr_index_map::Vector{Dict{Int, Vector{Tuple{Int,Int}}}}

    # mean value of Jastrow pseudopotentials
    mean_v::Vector{Vector{E}}

    # optimization flag
    optimize::Bool

    # orbital IDs 
    orbitals::Vector{Int}
end


@doc raw"""

    JastrowParameters(;
        particle_pair::AbstractString,
        order_pair::AbstractString,
        orbitals::Vector{Int},
        optimize::Bool,
        model_geometry::ModelGeometry,
        rng::AbstractRNG,
        v::Union{Vector{E}, Nothing} = nothing
    ) where {E<:AbstractFloat}

Initialize and return an instance of `JastrowParameters'.

"""
function JastrowParameters(;
    particle_pair::AbstractString,
    order_pair::AbstractString,
    orbitals::Vector{Int},
    optimize::Bool,
    model_geometry::ModelGeometry,
    rng::AbstractRNG,
    v::Union{Vector{Vector{E}}, Nothing} = nothing
) where {E<:AbstractFloat}
    (; unit_cell, lattice) = model_geometry
    # Norbs = unit_cell.n
    Norbs = length(orbitals)
    Ncells = lattice.N

    valid_orders = reduce(vcat, values(JASTROW_PARAMETERS))

    @assert particle_pair in keys(JASTROW_PARAMETERS) """
    Invalid particle_pair: "$particle_pair". Must be one of: $(join(keys(JASTROW_PARAMETERS), ", "))
    """

    @assert order_pair in valid_orders """
    Invalid order_pair: "$order_pair". Must be one of: $(join(valid_orders, ", "))
    """

    @assert order_pair in JASTROW_PARAMETERS[particle_pair] """
    Invalid combination: order_pair "$order_pair" is not compatible with particle_pair "$particle_pair".
    Expected one of: $(join(JASTROW_PARAMETERS[particle_pair], ", "))
    """
    
    # pre-allocate with known size hint where possible
    irr_indices = [Vector{Int}() for _ in 1:Norbs]
    mean_v = [Vector{Float64}() for _ in 1:Norbs]
    irr_index_map = [Dict{Int, Vector{Tuple{Int,Int}}}() for _ in 1:Norbs]

    # OLD METHOD
    # for o in 1:Norbs
    #     d = irr_index_map[o]
    #     ii = irr_indices[o]
    #     mv = mean_v[o]

    #     for u in 1:Ncells
    #         i = loc_to_site(u, o, unit_cell) - 1
    #         red_idx = reduce_index(0, i, lattice)
    #         pairs = get!(d, red_idx) do
    #             push!(ii, red_idx)
    #             push!(mv, 0.01 * rand(rng))
    #             Vector{Tuple{Int,Int}}()
    #         end
    #         push!(pairs, (0, i))
    #     end

    #     for u in 1:Ncells
    #         i = loc_to_site(u, o, unit_cell) - 1
    #         for u2 in 1:Ncells
    #             j = loc_to_site(u2, o, unit_cell) - 1
    #             i == 0 && continue
    #             j == 0 && continue
    #             red_idx = reduce_index(i, j, lattice)
    #             pairs = get!(d, red_idx) do
    #                 push!(ii, red_idx)
    #                 push!(mv, 0.01 * rand(rng))
    #                 Vector{Tuple{Int,Int}}()
    #             end
    #             push!(pairs, (i, j))
    #         end
    #     end
    # end

    for o in 1:Norbs
        d = irr_index_map[o]   
        ii = irr_indices[o]   
        mv = mean_v[o]

        for u in 1:Ncells
            i = loc_to_site(u, o, unit_cell) - 1
            for u2 in 1:Ncells
                j = loc_to_site(u2, o, unit_cell) - 1
                # j < i && continue
                red_idx = reduce_index(i, j, lattice)
                pairs = get!(d, red_idx) do
                    push!(ii, red_idx)
                    push!(mv, 0.01 * rand(rng))
                    Vector{Tuple{Int,Int}}()
                end
                push!(pairs, (i, j))
            end
        end
    end

    # override mean_v with provided values if given
    if !isnothing(v)
        @assert length(v) == Norbs "Length of v ($(length(v))) must match number of orbitals ($Norbs)"
        for o in 1:Norbs
            @assert length(v[o]) == length(irr_indices[o]) "Length of v[$o] ($(length(v[o]))) must match number of irreducible indices for orbital $o ($(length(irr_indices[o])))"
            mean_v[o] = collect(Float64, v[o])
        end
    end

    # set maximum distance parameter to 0 and sort, per orbital
    for o in 1:Norbs
        ii = irr_indices[o]
        isempty(ii) && continue

        mean_v[o][argmax(ii)] = 0.0

        sort_perm = sortperm(ii)
        irr_indices[o] = ii[sort_perm]
        mean_v[o] = mean_v[o][sort_perm]
    end

    return JastrowParameters(particle_pair, order_pair, irr_indices, irr_index_map, mean_v, optimize, orbitals)
end


# print struct info in TOML format
function Base.show(io::IO, ::MIME"text/plain", jas_par::JastrowParameters{E}) where {E}

    (; particle_pair, order_pair, irr_indices, mean_v, orbitals) = jas_par

    println(io, "[JastrowParameters]")
    println(io, "PARTICLES = $particle_pair\n")
    println(io, "ORDERS = $order_pair\n")

    unique_orbitals = unique(orbitals)
    for o in unique_orbitals
        println(io, "[[JastrowParameters.pseudopotentials]]")
        println(io, "ORBITAL_ID = $o\n")

        for (idx,v) in zip(irr_indices[o], mean_v[o])
            println(io, "k       = $(idx)")
            println(io, "mean    = $(v)\n")
        end
    end

    return nothing
end


