@doc raw"""

    diagonalize!(H::Matrix{T}) where {T<:Number}
    
Returns all eigenenergies ``\varepsilon`` and all eigenstates ``\phi`` of any Hamiltonian matrix. 

# Note 

The columns of the returned unitary matrix are the eigenstates of the Hamiltonian matrix.


"""
function diagonalize!(H::Matrix{T}) where {T<:Number}
    if T <: Real
        # in-place path exists
        F = eigen!(H)
    else
        # Complex Hermitian must use non-mutating eigen
        F = eigen(H)
    end

    return F.values, F.vectors  
end


"""

    is_invertible(
        D::AbstractMatrix{T}; 
        # KEYWORD ARGUMENTS
        rcond_tol::Float64=1e-12,
        det_tol::Float64=0.0 
    ) where {T<:Number}

Checks the invertibility of matrix `D`. Returns `true` if `D` is factorable and well-conditioned and will return `false` if
matrix inversion fails due to the failing to compute the determinant, the determinant value is outside of the tolerance, or
if the Slater matrix is ill-conditioned.

"""
function is_invertible(
    D::AbstractMatrix{T}; 
    # KEYWORD ARGUMENTS
    rcond_tol::Float64=1e-12,
    det_tol::Float64=0.0 
) where {T<:Number}
    info = (rcond = nothing, det = nothing)

    detD = try
        det(D)
    catch e
        #@warn "Matrix inversion failure: cannot compute determinant."

        return false
    end

    info = merge(info, (det = detD,))

    if detD == 0.0 || abs(detD) <= det_tol
        # @warn """
        # Matrix inversion failure: |det(D)| = 0 or |det(D)| < $det_tol.
        # """ info

        return false
    end

    rc = rcond(D)
    info = merge(info, (rcond = rc, cond_est = nothing, svmin = nothing))

    if rc < rcond_tol
        sv = svdvals(D)
        info = (rcond = rc, cond_est = sv[1]/sv[end], svmin = sv[end])
        # @warn """
        # Matrix inversion failure: Slater matrix is ill-conditioned.
        # """ info

        return false
    else
        return true
    end
end


"""
    
    rcond(D::AbstractMatrix{T}) where {T<:Number}

Return an SVD-based reciprocal condition estimate.

# Note

This may be slower than some LAPACK implementations but is robust.

"""
function rcond(D::AbstractMatrix{T}) where {T<:Number}
    s = svdvals(D)       
    if isempty(s) || s[1] == 0.0
        return 0.0
    else
        return s[end] / s[1]
    end
end


@doc raw"""

    inverse(D::AbstractMatrix{T}) where {T<:Number}

Computes the full matrix inverse of matrix `D`.

"""
function inverse(D::AbstractMatrix{T}) where {T<:Number}
    D_copy = copy(D)
    
    # perform LU factorization
    lu_fact = lu(D_copy)
    Dinv = inv(lu_fact)
        
    return Dinv
end


# rebin a data array along a specified array dimension
function rebin(
    data::AbstractArray{T},
    N_bins::Int,
    dim::Int = 1
) where {T<:Number}

    # Get the length of the dimension to be rebinned
    N_data = size(data,dim)
    @assert iszero(mod(N_data, N_bins))
    # if no rebinning is required
    if N_data == N_bins
        # relabel the original data array the rebinned array
        rebinned_data =  data
    # perform rebinning if necessary
    else
        # Calculate the bin size
        N_binsize = N_data ÷ N_bins
        # get size of data array
        data_dims = size(data)
        # calculated reshaped dims
        reshaped_dims = (data_dims[1:dim-1]..., N_binsize, N_bins, data_dims[dim+1:end]...)
        # reshape the data array
        reshaped_data = reshape(data, reshaped_dims)
        # calculate the average of each new bin and write the rebinned data to a new array
        rebinned_data = dropdims(mean(reshaped_data, dims=dim), dims=dim)
    end

    return rebinned_data
end