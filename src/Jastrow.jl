@doc raw"""

    JastrowFactor{E<:AbstractFloat, I<:Integer}

A type defining quantities related to a Jastrow factor.  

- `Tvec_f::Vector{E}`: fermionic T vector.
- `Tvec_b::Vector{E}`: bosonic or phononic T vector. 
- `nq_updates_T::I`: tracker for the number of quick updates to the T vector.

"""
mutable struct JastrowFactor{E<:AbstractFloat, I<:Integer}
    # fermionic T vector
    Tvec_f::Vector{E}

    # bosonic (phononic) T vector
    Tvec_b::Vector{E}

    # number of T vector quick updates
    nq_updates_T::I
end 


@doc raw"""

    get_jastrow_factor( jastrow_parameters::JastrowParameters{S, K, V, I}, 
                        detwf::DeterminantalWavefunction{T, Q, E, I}, 
                        model_geometry::ModelGeometry, 
                        pht::Bool )::JastrowFactor

Given a set of Jastrow parameters, constructs the relevant Jastrow factor ``\mathcal{J}_n`` and 
returns a instance of the `JastrowFactor` type. 

- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow parameters. 
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function get_jastrow_factor(
    jastrow_parameters::JastrowParameters{S, K, V, I}, 
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry, 
    pht::Bool
) where {S<:AbstractString, K, V, I<:Integer, T<:Number, Q, E<:AbstractFloat}
    # extent of the lattice
    N = model_geometry.lattice.N

    # Jastrow type
    jastrow_type = jastrow_parameters.jastrow_type

    # map of Jastrow parameters
    jpar_map = jastrow_parameters.jpar_map

    if jastrow_type == "e-den-den" || jastrow_type == "e-spn-spn"
        # generate fermionic T vector
        init_Tvec_f = get_fermionic_Tvec(
            jastrow_type,
            jpar_map, 
            detwf,
            N, 
            pht
        )

        # create null bosonic T vector
        init_Tvec_b = zeros(Float64, length(init_Tvec_f))
    end

    # intialize quick updating tracker
    nq_updates_T = 0

    return JastrowFactor(init_Tvec_f, init_Tvec_b, nq_updates_T) 
end


@doc raw"""

    get_fermionic_Tvec( jastrow_type::S,  
                        jpar_map::OrderedDict{Any, Any},
                        detwf::DeterminantalWavefunction{T, Q, E, I}, 
                        N::I,
                        pht::Bool )

Returns a fermionic T vector with elements ``T_{i} = \sum_{j} v_{ij} n_{i}`` where ``v\_{ij}`` are the   
Jastrow peseudopotentials and ``n_{i}`` are the total fermion occupations.

- `jastrow_type::S`: either "e-den-den" or "e-spn-spn.
- `jpar_map::OrderedDict{Any, Any}`: current map of Jastrow parameters.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction. 
- `N::I`: number of lattice sites.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_fermionic_Tvec(
    jastrow_type::S,
    jpar_map::OrderedDict{Any, Any},
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    N::I,
    pht::Bool
) where {S<:AbstractString, I<:Integer, T<:Number, Q, E<:AbstractFloat}
    # initialize T vector
    Tvec_f = zeros(N) 

    for (irr_index, _) in jpar_map
        # get Jastrow parameter
        v_ij = jpar_map[irr_index][2]

        # loop over 
        for (i,j) in jpar_map[irr_index][1]
            # get site occupations
            num_up = get_onsite_fermion_occupation(j+1, detwf.pconfig)[1]
            num_dn = get_onsite_fermion_occupation(j+1, detwf.pconfig)[2]

            if jastrow_type == "e-den-den"
                if pht
                    Tvec_f[i+1] += v_ij * (num_up - num_dn)
                else
                    Tvec_f[i+1] += v_ij * (num_up + num_dn)
                end
            elseif jastrow_type == "e-spn-spn"
                if pht
                    Tvec_f[i+1] += 0.5 * v_ij * (num_up - num_dn)
                else
                    Tvec_f[i+1] += 0.5 * v_ij * (num_up + num_dn)
                end
            end
        end
    end
    
    return Tvec_f
end


@doc raw"""

    get_fermionic_jastrow_ratio( markov_move::MarkovMove{I}, 
                                 jastrow_parameters::JastrowParameters{S, K, V, I},
                                 jastrow_factor::JastrowFactor{E}, 
                                 pht::Bool, 
                                 spin::I, 
                                 model_geometry::ModelGeometry ) where {I<:Integer, S<:AbstractString, K, V, E<:AbstractFloat}

Calculates ratio of fermionic Jastrow factors
```math
\frac{\mathcal{J}_1}{\mathcal{J}_2} = \exp[s(T_{l} - T_{k}) + v_{ll} - v_{lk}],
```
after a single particle completes a move from site ``k`` to site ``l`` using the corresponding T vectors ``T_k`` and ``T_l``.

- `markov_move::MarkovMove{I}`: quantities related to a Markov process.  
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `jastrow_factor::JastrowFactor{E}`: current Jastrow factor.
- `pht::Bool`: whether model is particle-hole transformed.
- `spin::I`: spin of the current particle.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_fermionic_jastrow_ratio(
    markov_move::MarkovMove{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    jastrow_factor::JastrowFactor{E}, 
    pht::Bool, 
    spin::I, 
    model_geometry::ModelGeometry
) where {I<:Integer, S<:AbstractString, K, V, E<:AbstractFloat}
    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # initial site of particle
    k = markov_move.k

    # final site of particle
    l = markov_move.l

    # convert spindices to real sites
    ksite = get_index_from_spindex(k, model_geometry)
    lsite = get_index_from_spindex(l, model_geometry)

    # current fermionic T vector
    Tvec_f = jastrow_factor.Tvec_f

    # get T vector elements 
    Tₖ = Tvec_f[ksite]
    Tₗ = Tvec_f[lsite]

    # map of Jastrow parameters
    jpar_map = jastrow_parameters.jpar_map

    # obtain reduced indices
    if dims == 1
        red_ll = reduce_index_1d(lsite-1, lsite-1, model_geometry)
        red_lk = reduce_index_1d(lsite-1, ksite-1, model_geometry)
    elseif dims == 2
        red_ll = reduce_index_2d(lsite-1, lsite-1, model_geometry)
        red_lk = reduce_index_2d(lsite-1, ksite-1, model_geometry)
    end

    @assert haskey(jpar_map, red_ll)
    @assert haskey(jpar_map, red_lk)

    # get Jastrow parameters
    vₗₗ = jpar_map[red_ll][2]
    vₗₖ = jpar_map[red_lk][2]

    # compute Jastrow ratio
    if pht  
        jas_ratio_f = exp(spin * (Tₗ - Tₖ) + vₗₗ - vₗₖ)      
    else
        jas_ratio_f = exp(Tₗ - Tₖ + vₗₗ - vₗₖ)              
    end

    return jas_ratio_f 
end


@doc raw"""

    get_fermionic_jastrow_ratio( k::I, 
                                 l::I, 
                                 jastrow_parameters::JastrowParameters{S, K, V, I},
                                 jastrow_factor::JastrowFactor{E}, 
                                 pht::Bool, 
                                 spin::I, 
                                 model_geometry::ModelGeometry ) where {I<:Integer, S<:AbstractString, K, V, E<:AbstractFloat}

Calculates ratio of fermionic Jastrow factors
```math
\frac{\mathcal{J}_1}{\mathcal{J}_2} = \exp[-s(T_{l} - T_{k}) + v_{ll} - v_{lk}],
```
after a single particle completes a move from site ``k`` to site ``l`` using the corresponding T vectors ``T_k`` and ``T_l``.

- `k::Int`: initial site of the current particle. 
- `l::Int`: final site of the current particle. 
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters.
- `jastrow_factor`::JastrowFactor: current Jastrow factor.
- `pht`::Bool: whether model is particle-hole transformed.
- `spin`::Int: spin of the current particle.
- `model_geometry`::ModelGeometry: contains unit cell and lattice quantities.

"""
function get_fermionic_jastrow_ratio(
    k::I, 
    l::I, 
    jastrow_parameters::JastrowParameters{S, K ,V, I},
    jastrow_factor::JastrowFactor{E}, 
    pht::Bool, 
    spin::I,
    model_geometry::ModelGeometry
) where {I<:Integer, S<:AbstractString, K, V, E<:AbstractFloat}
    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # convert spindices to real sites
    ksite = get_index_from_spindex(k, model_geometry)
    lsite = get_index_from_spindex(l, model_geometry)

    # current fermionic T vector
    Tvec_f = jastrow_factor.Tvec_f

    # get T vector elements 
    Tₖ = Tvec_f[ksite]
    Tₗ = Tvec_f[lsite]

    # map of Jastrow parameters
    jpar_map = jastrow_parameters.jpar_map

    # obtain reduced indices
    if dims == 1
        red_ll = reduce_index_1d(lsite-1, lsite-1, model_geometry)
        red_lk = reduce_index_1d(lsite-1, ksite-1, model_geometry)
    elseif dims == 2
        red_ll = reduce_index_2d(lsite-1, lsite-1, model_geometry)
        red_lk = reduce_index_2d(lsite-1, ksite-1, model_geometry)
    end

    @assert haskey(jpar_map, red_ll)
    @assert haskey(jpar_map, red_lk)

    # get Jastrow parameters
    vₗₗ = jpar_map[red_ll][2]
    vₗₖ = jpar_map[red_lk][2]

    # compute Jastrow ratio
    if pht  
        jas_ratio_f = exp(spin * (Tₗ - Tₖ) + vₗₗ - vₗₖ)      
    else
        jas_ratio_f = exp(Tₗ - Tₖ + vₗₗ - vₗₖ)               
    end

    return jas_ratio_f 
end



@doc raw"""

    update_fermionic_Tvec!( markov_move::MarkovMove{I}, 
                            spin::I, 
                            jastrow_parameters::JastrowParameters{S, K, V, I},
                            jastrow_factor::JastrowFactor{E}, 
                            detwf::DeterminantalWavefunction,
                            model_geometry::ModelGeometry, 
                            n_stab_T::I, 
                            δT::E, 
                            pht::Bool ) where {S<:AbstractString, K, V, I<:Integer, E<:AbstractFloat}

Updates elements ``T_{i}`` of the fermion-type T vector after an accepted Metropolis step.

- `markov_move::MarkovMove{I}`: quantities related to a Markov process. 
- `spin::I`: spin of the current particle. 
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `jastrow::JastrowFactor{E}`: current Jastrow factor.
- `detwf::DeterminantalWavefunction`: current determinantal wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `n_stab_T::I`: frequency of T vector stabilization setps.
- `δT::E`: deviation threshold for the T vector.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function update_fermionic_Tvec!(
    markov_move::MarkovMove{I}, 
    spin::I, 
    jastrow_parameters::JastrowParameters{S, K, V, I},
    jastrow_factor::JastrowFactor{E}, 
    detwf::DeterminantalWavefunction,
    model_geometry::ModelGeometry, 
    n_stab_T::I, 
    δT::E, 
    pht::Bool
) where {S<:AbstractString, K, V, I<:Integer, E<:AbstractFloat}
    # number of lattice sites
    N = model_geometry.lattice.N

    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # initial site of particle
    k = markov_move.k

    # final site of particle
    l = markov_move.l

    # convert spindices to real sites
    ksite = get_index_from_spindex(k, model_geometry)
    lsite = get_index_from_spindex(l, model_geometry)

    # different T vectors
    Tvec_f_diff = zeros(N)

    # map of Jastrow parameters
    jpar_map = jastrow_parameters.jpar_map

    for i in 1:N
        if dims == 1
            red_ik = reduce_index_1d(i-1, ksite-1, model_geometry)
            red_il = reduce_index_1d(i-1, lsite-1, model_geometry)
        elseif dims == 2
            red_ik = reduce_index_2d(i-1, ksite-1, model_geometry)
            red_il = reduce_index_2d(i-1, lsite-1, model_geometry)
        end
        
        @assert haskey(jpar_map, red_ik)
        @assert haskey(jpar_map, red_il)

        # get Jastrow parameters
        vᵢₖ = jpar_map[red_ik][2]
        vᵢₗ = jpar_map[red_il][2]

        if pht
            Tvec_f_diff[i] += spin * (vᵢₗ - vᵢₖ)
        else
            Tvec_f_diff[i] += (vᵢₗ - vᵢₖ)
        end
    end

    if pht
        jastrow_factor.Tvec_f += spin * Tvec_f_diff
    else
        jastrow_factor.Tvec_f += Tvec_f_diff
    end

    if jastrow_factor.nq_updates_T >= n_stab_T
        @debug """
        Jastrow::update_fermionic_Tvec!() :
        recalculating T!
        """

        # reset counter 
        jastrow_factor.nq_updates_T = 0

        # recalculate T vector from scratch
        # Jastrow type
        jastrow_type = jastrow_parameters.jastrow_type

        # re-calculate the fermionic T vector 
        Tvec_f_r = get_fermionic_Tvec(jastrow_type, jpar_map, detwf, N, pht)    

        # compute deviation of the original T vector and the recomputed T vector
        dev = check_deviation(jastrow_factor.Tvec_f, Tvec_f_r)

        @debug """
        Jastrow::update_fermionic_Tvec!() :
        recalculated T with deviation = $(dev)
        """

        if dev > δT
            @debug """
            Jastrow::update_fermionic_Tvec!() :
            Deviation goal for vector T (fermionic) not met!
            """

            # replace original T vector with new one
            jastrow_factor.Tvec_f = Tvec_f_r

        else
            @debug """
            Jastrow::update_fermionic_Tvec!() :
            Deviation goal for T met! Jastrow T vector is stable
            """

            @assert dev < δT
        end

        return nothing
    else
        @debug """
        Jastrow::update_fermionic_Tvec!() :
        performing quick update of T!
        """

        jastrow_factor.nq_updates_T += 1 

        return nothing
    end
end


@doc raw"""

    map_jastrow_parameters( model_geometry::ModelGeometry, 
                            rng::AbstractRNG )

Generates a dictionary of irreducible index keys ``k`` with values being a tuple of all 
lattice index pairs ``(i,j)``, which generate ``k``, and their associated Jastrow parameter 
``v_{ij}``. The parameter corresponding to the largest irreducible index is initialized to zero.

- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.  
- `rng::AbstractRNG`: random number generator.

"""
function map_jastrow_parameters(
    model_geometry::ModelGeometry, 
    rng::AbstractRNG
)
    # number of lattice sites
    N = model_geometry.lattice.N
    
    # vector to store reduced indices
    reduced_indices = []

    # initialize map of Jastrow parameters
    jpar_map = OrderedDict()

    for i in 0:N-1
        # get irreducible lattice index 
        red_idx = reduce_index_2d(0, i, model_geometry)
        push!(reduced_indices, red_idx)
        if haskey(jpar_map, red_idx)
            indices, init_val = jpar_map[red_idx]
            push!(indices, (0, i))
            jpar_map[red_idx] = (indices, init_val)
        else
            jpar_map[red_idx] = ([(0, i)], 0.01*rand(rng)) #0.01*rand(rng)
        end
    end

    for i in 1:N-1
        for j in 1:N-1
            # get irreducible lattice index 
            red_idx = reduce_index_2d(i, j, model_geometry)
            push!(reduced_indices, red_idx)
            if haskey(jpar_map, red_idx)
                indices, init_val = jpar_map[red_idx]
                push!(indices, (i, j))
                jpar_map[red_idx] = (indices, init_val)
            else
                jpar_map[red_idx] = ([(i, j)], 0.01*rand(rng))#0.01*rand(rng)
            end
        end
    end

    # set the parameter corresponding to the maximum distance to 0
    max_idx = maximum(keys(jpar_map))
    if haskey(jpar_map, max_idx)
        indices, _ = jpar_map[max_idx]
        jpar_map[max_idx] = (indices, 0.0)
    end

    # sort the dictionary by irreducible indices
    sorted_jpar_map = OrderedDict(sort(collect(jpar_map)))

    return sorted_jpar_map
end


@doc raw"""

    map_jastrow_parameters( model_geometry::ModelGeometry, 
                            init_jpars::Vector{E} ) where {E<:AbstractFloat}

Given an initial set of Jastrow parameters from file, generates a dictionary of irreducible 
index keys ``k`` with values being a tuple of all lattice index pairs ``(i,j)``, which generate 
``k``, and their associated Jastrow parameter ``v_{ij}``. The parameter corresponding to the 
largest irreducible index is initialized to zero.

- `model_geometry::ModelGeometry`: contains unit cell and lattice qunatities. 
- `init_jpars::Vector{E}`: initial Jastrow parameters from file.

"""
function map_jastrow_parameters(
    model_geometry::ModelGeometry, 
    init_jpars::Vector{E}
) where {E<:AbstractFloat}
    # number of lattice sites
    N = model_geometry.lattice.N
    
    # vector to store reduced indices
    reduced_indices = []

    # initialize map of Jastrow parameters
    jpar_map = OrderedDict()

    param_index = 1
    for i in 0:N-1
        red_idx = reduce_index_2d(0, i, model_geometry)
        push!(reduced_indices, red_idx)
        if haskey(jpar_map, red_idx)
            indices, mean_vij = jpar_map[red_idx]   
            push!(indices, (0, i))
            jpar_map[red_idx] = (indices, mean_vij)
        else
            mean_vij = init_jpars[param_index]
            jpar_map[red_idx] = ([(0, i)], mean_vij)
            param_index += 1
        end
    end

    for i in 1:N-1
        for j in 1:N-1
            red_idx = reduce_index_2d(i, j, model_geometry)
            push!(reduced_indices, red_idx)
            if haskey(jpar_map, red_idx)
                indices, mean_vij = jpar_map[red_idx]
                push!(indices, (i, j))
                jpar_map[red_idx] = (indices, mean_vij)
            else
                mean_vij = init_jpars[param_index]
                jpar_map[red_idx] = ([(i, j)], mean_vij)
                param_index += 1
            end
        end
    end

    # set the parameter corresponding to the maximum distance to 0
    max_idx = maximum(keys(jpar_map))
    if haskey(jpar_map, max_idx)
        indices, _ = jpar_map[max_idx]
        jpar_map[max_idx] = (indices, 0.0)
    end

    # Sort the dictionary by key
    sorted_jpar_map = OrderedDict(sort(collect(jpar_map)))

    return sorted_jpar_map
end


@doc raw"""

    check_deviation( jastrow_Tvec::Vector{E}, 
                     Tvec_r::Vector{E} ) where {E<:AbstractFloat}

Checks floating point error accumulation in the fermionic T vector.

- `jastrow_Tvec::Vector{E}`: current T vector. 
- `Tvec_r::Vector{E}`: recalculated T vector. 

"""
function check_deviation(
    jastrow_Tvec::Vector{E}, 
    Tvec_r::Vector{E}
) where {E<:AbstractFloat}
    @assert size(jastrow_Tvec) == size(Tvec_r)

    exact_square_sum = sum(abs2, jastrow_Tvec)
    diff_square_sum  = sum(abs2, jastrow_Tvec .- Tvec_r)

    return exact_square_sum == 0.0 ?
           sqrt(diff_square_sum) :
           sqrt(diff_square_sum / exact_square_sum)
end




