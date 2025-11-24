@doc raw"""

    initialize_equal_time_greens( W::Matrix{T}, 
                                  D::Matrix{T}, 
                                  M::AbstractMatrix{T}, 
                                  pconfig::Vector{I}, 
                                  Np::I,
                                  config_indices::Vector{I},
                                  Dt::Matrix{T},
                                  tmp::Matrix{T} ) where {T<:Number, I<:Integer}
    
Computes the equal-time Green's function by solving the matrix equation  ``D^{T}W^{T} = M^{T}`` 
using LU decomposition, where ``W`` is the equal-time Green's function, ``D`` is the Slater
matrix, and ``M`` is the reduced ``U`` matrix.

- `W::Matrix{T}`: equal-time Green's function matrix.
- `D::Matrix{T}`: Slater determinant matrix.
- `M::AsbtractMatrix{T}`: reduced U_int matrix.
- `pconfig::Vector{Int}`: current particle configuration.
- `Np::I`: total number of particles in the system.
- `config_indices::Vector{I}`: vector to store particle indices.
- `Dt::Matrix{T}`: helper matrix for the calculation of transpose(D)
- `tmp::Matrix{T}`: helper matrix for performing the LU decomposition

"""
function initialize_equal_time_greens!(
    W::Matrix{T}, 
    D::Matrix{T}, 
    M::AbstractMatrix{T}, 
    pconfig::Vector{I}, 
    Np::I,
    config_indices::Vector{I},
    Dt::Matrix{T},
    tmp::Matrix{T}
) where {T<:Number, I<:Integer}
    # get indices for the current particle configuration
    @inbounds for i in 1:Np
        config_indices[i] = findfirst(==(i), pconfig)
    end

    # build Slater matrix D
    @inbounds for (j, row_idx) in enumerate(config_indices)
        @views D[j, :] .= M[row_idx, :]
    end

    if is_invertible(D)   
        # LU decomposition of D'
        @views Dt .= D'
        lu_decomp = lu!(Dt)

        @views tmp .= transpose(M)
        ldiv!(lu_decomp, tmp)
        
        # calculate the equal-time Green's function
        permutedims!(W, tmp, (2, 1))    

        return true
    else        
        # former invertibility test
        # abs(det(D)) < 1e-12 * size(D, 1) 

        @debug """
        Greens::initialize_equal_time_greens() : 
        state has no overlap with the determinantal wavefunction => 
        determinant of D = $(det(D))
        """

        return false
    end            
end


@doc raw"""

    update_equal_time_greens!( markov_move::MarkovMove{I}, 
                               detwf::DeterminantalWavefunction{T, V, E, I}, 
                               model_geometry::ModelGeometry,
                               Np::I, 
                               n_stab_W::I, 
                               δW::E ) where {I<:Integer, T<:Number, V, E<:AbstractFloat}

Updates the equal-time Green's function while performing a numerical stabilzation check after
`n_stab_W` steps. If the calculated deviation exceeds the set threshold ``\delta W``, then the current 
Green's function is replaced by one calculated from scratch.

- `markov_move::MarkovMove{I}`: quantities related to a Markov process.
- `detwf::DeterminantalWavefunction{T, V, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `n_stab_W::I`: frequency of Green's function stability steps.
- `δW::E`: Green's function stability threshold.

"""
function update_equal_time_greens!(
    markov_move::MarkovMove{I}, 
    detwf::DeterminantalWavefunction{T, V, E, I}, 
    model_geometry::ModelGeometry,
    Np::I, 
    n_stab_W::I, 
    δW::E
) where {I<:Integer, T<:Number, V, E<:AbstractFloat}
    if detwf.nq_updates_W >= n_stab_W
        @debug """
        Greens::update_equal_greens!() : 
        recalculating W!
        """

        # reset counter 
        detwf.nq_updates_W = 0

        # perform rank-1 update
        rank1_update!(
            markov_move, 
            detwf
        )

        # number of lattice sites
        N = model_geometry.lattice.N

        # re-initialize W matrix
        Wᵣ = zeros(ComplexF64, 2*N, Np)

        # re-initialize Slater matrix
        Dᵣ = zeros(ComplexF64, Np, Np)

        # vector to store configuration indices
        config_indices = Vector{Int}(undef, Np)

        # helper matrices
        Dt = zeros(ComplexF64, Np, Np)            
        tmp = zeros(ComplexF64, Np, 2*N) 

        # recalclate Green's function from scratch
        recalculate_equal_time_greens!(
            Wᵣ, 
            Dᵣ, 
            detwf.M, 
            detwf.pconfig, 
            Np, 
            config_indices,
            Dt,
            tmp
        )

        # compute deviation between current Green's function and the recalculated Green's function
        dev = check_deviation(
            detwf, 
            Wᵣ
        )

        @debug """
        Greens:update_equal_greens!() : 
        recalculated W with deviation = $(dev)
        """

        if dev > δW
            @debug """
            Greens::update_equal_greens!() : 
            deviation goal for matrix W not met!
            """

            # replace original W matrix with new one
            detwf.W = Wᵣ
        else
            @debug """
            Greens::update_equal_greens!() :
            deviation goal for matrix W met! Green's function is stable.
            """
            @assert dev < δW
        end  

        return nothing
    else
        @debug """
        Greens::update_equal_greens!() : 
        performing rank-1 update of W!
        """ 

        # perform rank-1 update
        rank1_update!(
            markov_move, 
            detwf
        )

        detwf.nq_updates_W += 1

        return nothing
    end
end


@doc raw"""

    rank1_update!( markov_move::MarkovMove{I}, 
                   detwf::DeterminantalWavefunction{T, V, E, I} ) where {I<:Integer, T<:Number, V, E<:AbstractFloat}
    
Performs an in-place rank-1 update of the equal-time Green's function. 

- `markov_move::MarkovMove{I}`: quantities related to a Markov process. 
- `detwf::DeterminantalWavefunction{T, V, E, I}`: current variational wavefunction. 

"""
function rank1_update!(
    markov_move::MarkovMove{I}, 
    detwf::DeterminantalWavefunction{T, V, E, I}
) where {I<:Integer, T<:Number, V, E<:AbstractFloat}
    # particle 
    β = markov_move.particle

    # final site of the hopping particle
    l = markov_move.l

    # get lth row of the Green's function
    rₗ  = detwf.W[l, :]

    # subtract 1 from the βth element
    rₗ[β] -= 1.0

    # get the βth column of the Green's function
    cᵦ = detwf.W[:, β]

    # # Perform rank-1 update (by hand)
    # detwf.W -= (cᵦ / detwf.W[l, β]) * rₗ' 

    # Perform rank-1 update
    BLAS.geru!(
        -1.0 / detwf.W[l, β], 
        cᵦ, 
        rₗ, 
        detwf.W
    )

    return nothing
end


@doc raw"""

    recalculate_equal_time_greens( Wᵣ::Matrix{T}, 
                                   Dᵣ::Matrix{T}, 
                                   M::AbstractMatrix{T}, 
                                   pconfig::Vector{I}, 
                                   Np::I,
                                   config_indices::Vector{I},
                                   Dt::Matrix{T},
                                   tmp::Matrix{T} ) where {T<:Number, I<:Integer}
    
Recomputes the equal-time Green's function.

- `Wᵣ::Matrix{T}`: equal-time Green's function matrix.
- `Dᵣ::Matrix{T}`: Slater determinant matrix.
- `M::AbstractMatrix{T}`: reduced U_int matrix.
- `pconfig::Vector{I}`: current particle configuration.
- `Np::I`: total number of particles in the system. 
- `config_indices::Vector{I}`: vector to store particle indices.
- `Dt::Matrix{T}`: helper matrix for the calculation of transpose(D).
- `tmp::Matrix{T}`: helper matrix for performing the LU decomposition.

"""
function recalculate_equal_time_greens!(
    Wᵣ::Matrix{T}, 
    Dᵣ::Matrix{T}, 
    M::AbstractMatrix{T}, 
    pconfig::Vector{I}, 
    Np::I,
    config_indices::Vector{I},
    Dt::Matrix{T},
    tmp::Matrix{T}
) where {T<:Number, I<:Integer}
    # get indices for the current particle configuration
    @inbounds for i in 1:Np
        config_indices[i] = findfirst(==(i), pconfig)
    end

    # build Slater matrix
    @inbounds for (j, row_idx) in enumerate(config_indices)
        @views Dᵣ[j, :] .= M[row_idx, :]
    end

    # LU decomposition of D'
    @views Dt .= Dᵣ'
    lu_decomp = lu!(Dt)

    @views tmp .= transpose(M)
    ldiv!(lu_decomp, tmp)
        
    # calculate the equal-time Green's function
    permutedims!(Wᵣ, tmp, (2, 1))    

    # former method for calculating W
    # lu_decomp = lu(Dᵣ')
    # Wᵣ .= transpose(lu_decomp \ transpose(M))

    return true      
end


@doc raw"""

    check_deviation!( detwf::DeterminantalWavefunction{T, V, E, I}, 
                      Wᵣ::Matrix{T} ) where {T<:Number, V, E<:AbstractFloat, I<:Integer}
    
Checks floating point error accumulation in the equal-time Green's function.

- `detwf::DeterminantalWavefunction{T, V, E, I}`: current variational wavefunction. 
- `Wᵣ::Matrix{T}`: reclculated Green's function matrix

"""
## UPDATED
function check_deviation(
    detwf::DeterminantalWavefunction{T, V, E, I}, 
    Wᵣ::Matrix{T}
) where {T<:Number, V, E<:AbstractFloat, I<:Integer} 
    @assert size(detwf.W) == size(Wᵣ)

    exact_square_sum = sum(abs2, detwf.W)
    diff_square_sum  = sum(abs2, detwf.W .- Wᵣ)

    return exact_square_sum == 0.0 ?
           sqrt(diff_square_sum) :
           sqrt(diff_square_sum / exact_square_sum)
end


# @doc raw"""

#     is_invertible( D::Matrix{T} ) where {T<:Number}
    
# Checks the invertibility of matrix `D`.

# - `D::Matrix{T}`: Slater matrix. 

# """
# function is_invertible(
#     D::Matrix{T}
# ) where {T<:Number}
#     try
#         lu_fact = lu(D)
#         return true   
#     catch e
#         return false  
#     end
# end


"""

    is_invertible_( D::AbstractMatrix{T}; 
                    rcond_tol=1e-12 ) where {T<:Number}

Checks the invertibility of matrix `D`. Returns `true` if D is factorable 
and well conditioned.

- `D::Matrix{T}`: Slater matrix. 
- `rcond_tol::Float64=1e-12`: condition number.
- `det_tol::Float64=0.0`: determinant tolerance. 

"""
## NEW 
function is_invertible(
    D::AbstractMatrix{T}; 
    rcond_tol::Float64=1e-12,
    det_tol::Float64=0.0 
) where {T<:Number}
    # check the determinant
    info = (rcond = nothing, det = nothing)

    detD = try
        det(D)
    catch e
        @debug """
        Greens::is_invertible(): determinant computation failed!
        """
        return false
    end

    info = merge(info, (det = detD,))

    # check magnitiude of the determinant
    if detD == 0.0 || abs(detD) <= det_tol
        @debug """
        Greens::is_invertible(): |det(D)| = 0 or is below det_tol = $det_tol 
        """
        @debug info
        return false
    end

    # # check LU factorization
    # info = (rcond = nothing, cond_est = nothing, svmin = nothing)
    # try
    #     F = lu(D) 
    # catch e
    #     @debug rcond = 0.0
    #     @debug cond_est = Inf
    #     @debug svmin = 0.0
    #     @debug """
    #     Greens::is_invertible() : 
    #     LU factorization of D failed! ->
    #     rcond = 0.0
    #     cond_est = Inf
    #     svmin = 0.0
    #     """
    #     return false
    # end

    # fast reciprocal condition estimate
    rc = rcond(D)     
    info = merge(info, (rcond = rc, cond_est = nothing, svmin = nothing))

    if rc < rcond_tol
        # optionally attach smallest singular value for debugging (expensive)
        @debug sv = svdvals(D)
        @debug info = (rcond = rc, cond_est = sv[1]/sv[end], svmin = sv[end])
        @debug """
        Greens::is_invertible() : 
        Slater matrix is ill-conditioned! ->
        $info
        """

        return false
    else
        return true
    end
end


"""
    rcond( D::AbstractMatrix{T} ) where {T<:Number}

Return an SVD-based reciprocal condition estimate: σ_min / σ_max.
May be slower than some LAPACK implementations but is robust.

- `D::Matrix{T}`: Slater matrix. 

"""
#TODO: faster/ more efficeint method?
function rcond(D::AbstractMatrix{T}) where {T<:Number}
    s = svdvals(D)        # returns singular values sorted descending
    if isempty(s) || s[1] == 0.0
        return 0.0
    else
        return s[end] / s[1]
    end
end


@doc raw"""

    inverse( D::Matrix{T} ) where {T<:Number}

Computes the full matrix inverse of matrix `D`.

- `D::Matrix{T}`: Slater matrix. 

"""
function inverse(
    D::Matrix{T}
) where {T<:Number}
    D_copy = copy(D)
    
    # perform LU factorization
    lu_fact = lu(D_copy)
    Dinv = inv(lu_fact)
        
    return Dinv
end