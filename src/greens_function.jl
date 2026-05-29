@doc raw"""

    calculate_equal_time_greens_function!(
        W::AbstractMatrix{<:Number},
        D::AbstractMatrix{<:Number},
        M::AbstractMatrix{<:Number},
        Dt::AbstractMatrix{<:Number},
        tmp::AbstractMatrix{<:Number},
        pconfig::Vector{Int},
        config_indices::Vector{Int},
        Np::Int
    )
    
Computes the equal-time Green's function by solving the matrix equation  ``D^{T}W^{T} = M^{T}`` using LU decomposition, 
where ``W`` is the equal-time Green's function, ``D`` is the Slater matrix, and ``M`` is the reduced ``U`` matrix.

# ARGUMENTS

- `W::AbstractMatrix{T}`: Equal-time Green's function matrix.
- `D::AbstractMatrix{T}`: Slater matrix.
- `M::AbstractMatrix{T}`: Reduced unitary matrix.
- `Dt::AbstractMatrix{T}`: Temporary storage matrix.
- `tmp::AbstractMatrix{T}`: Temporary storage matrix.
- `pconfig::Vector{Int}`: Particle configuration.
- `config_indices::Vector{Int}`: Temporary storage vector.
- `Np::Int`: Total number of particles in the configuration. 

"""
function calculate_equal_time_greens_function!(
    W::AbstractMatrix{<:Number},
    D::AbstractMatrix{<:Number},
    M::AbstractMatrix{<:Number},
    Dt::AbstractMatrix{<:Number},
    tmp::AbstractMatrix{<:Number},
    pconfig::Vector{Int},
    config_indices::Vector{Int},
    Np::Int
)
    # get indices for the current particle configuration
    @inbounds for i in 1:Np
        config_indices[i] = findfirst(==(i), pconfig)
    end

    # build Slater matrix D in-place using views
    @inbounds for (j, row_idx) in enumerate(config_indices)
        @views D[j, :] .= M[row_idx, :]
    end

    if !is_invertible(D)
        return false
    end

    # LU decompose D' in-place
    transpose!(Dt, D)
    lu_decomp = lu!(Dt)

    # solve and compute W = M * D^{-1} in-place
    transpose!(tmp, M)
    ldiv!(lu_decomp, tmp)
    transpose!(W, tmp)

    return true
end


@doc raw"""

    update_equal_time_greens_function!(
        W::Matrix{T},
        M::Matrix{T},
        β::Int,
        lsite::Int,
        Np::Int, 
        N::Int,
        nq_updates_W::Int,
        n_stab_W::Int, 
        δW::AbstractFloat
    ) where {T}

Updates the equal-time Green's function matrix `W` while performing a numerical stabilzation check after
`n_stab_W` steps. If the calculated deviation exceeds the set threshold ``\delta W``, then the current Green's 
function is replaced by one calculated from scratch.

# ARGUMENTS

- `W::Matrix{T}`: Equal-time Green's function matrix.
- `M::Matrix{T}`: Reduced unitary matrix.
- `β::Int`: ID of the hopping particle.
- `lsite::Int`: Final site spindex of the hopping particle with ID `β`.
- `Np::Int`: Total number of particles in the configuration.
- `N::Int`: Total number of lattice sites.
- `nq_updates_W::Int`: Quick update tracker for the eqaul-time Green's function.
- `n_stab_W::Int`: Frequency of numerical stability checks done on matrix `W`.
- `δW::AbstractFloat`: Maximum allowed error in the equal-time Green's function.

"""
function update_equal_time_greens_function!(
    W::Matrix{T},
    M::Matrix{T},
    β::Int,
    lsite::Int,
    Np::Int, 
    N::Int,
    nq_updates_W::Int,
    n_stab_W::Int, 
    δW::AbstractFloat,
    pconfig::Vector{Int}
) where {T}
    # perform rank-1 update
    rank1_update!(W, β, lsite)

    if nq_updates_W >= n_stab_W
        # preallocate matrices for stabilization
        Dᵣ  = Matrix{T}(undef, Np, Np)
        Wᵣ  = Matrix{T}(undef, 2*N, Np)
        Dt  = Matrix{T}(undef, Np, Np)
        tmp = Matrix{T}(undef, Np, 2*N)
        config_indices = Vector{Int}(undef, Np)

        # recalculate the Green's function
        calculate_equal_time_greens_function!(Wᵣ, Dᵣ, M, Dt, tmp, pconfig, config_indices, Np)

        # compute deviation
        dev = check_deviation(W, Wᵣ)

        if dev > δW
            copyto!(W, Wᵣ)  # in-place copy instead of rebinding
        end

        return 0  # reset counter
    else
        return nq_updates_W + 1  # increment counter
    end
end


@doc raw"""

    rank1_update!(W::Matrix{E}, β::Int, lsite::Int) where {E}

Performs a rank-1 update on matrix `W`. 

# ARGUMENTS

- `W::Matrix{T}`: Equal-time Green's function matrix.
- `β::Int`: ID of the hopping particle.
- `lsite::Int`: Final site spindex of the hopping particle `β`.

"""
function rank1_update!(W::Matrix{T}, β::Int, lsite::Int) where {T}
    # get lth row of the Green's function
    rₗ  = W[lsite, :]

    # subtract 1 from the βth element
    rₗ[β] -= 1.0

    # get the βth column of the Green's function
    cᵦ = W[:, β]

    # # Perform rank-1 update (slower)
    # W -= (cᵦ / W[lsite, β]) * rₗ' 

    # Perform rank-1 update
    BLAS.geru!(
        -1.0 / W[lsite, β], 
        cᵦ, 
        rₗ, 
        W
    )

    return nothing
end


@doc raw"""

    check_deviation(W::Matrix{T}, Wᵣ::Matrix{T}) where {T} 
    
Checks floating point error accumulation in the equal-time Green's function.

# ARGUMENTS

-`W::Matrix{T}`: Updated Green's function matrix.
-`Wᵣ::Matrix{T}`: Recalculated Green's function matrix.


"""
function check_deviation(W::Matrix{T}, Wᵣ::Matrix{T}) where {T} 
    @assert size(W) == size(Wᵣ)

    exact_square_sum = sum(abs2, W)
    diff_square_sum  = sum(abs2, W .- Wᵣ)

    return exact_square_sum == 0.0 ?
           sqrt(diff_square_sum) :
           sqrt(diff_square_sum / exact_square_sum)
end





