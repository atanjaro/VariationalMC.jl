@doc raw"""

    Optimizer{T<:Number}

A mutable struct containing all the parameters needed to characterize a Stochastic Reconfiguration optimizer.

# Fields

- `var_params::Vector{T}`: All variational parameters in the model.
- `η::AbstractFloat`: Parameter which ensures the stability of the optimization.
- `dt::AbstractFloat`: Parameter which controls the rate of optimization.
- `S::Matrix{T}`: Covariance matrix.
- `f::Vector{T}`: Forces.

"""
mutable struct Optimizer{T<:Number}
    # variational parameters
    var_params::Vector{T}

    # optimization stability factor
    η::AbstractFloat

    # optimization rate
    dt::AbstractFloat

    # covariance matrix
    S::Matrix{T}

    # force vector
    f::Vector{T}
end


@doc raw"""

    Optimizer(;
        determinantal_parameters::DeterminantalParameters{S, T},
        η::AbstractFloat, dt::AbstractFloat,
        jas_parameters::Union{Tuple{JastrowParameters{T}}, Nothing} = nothing
    ) where {S, T}

Initialize and return an instance of `Optimizer`.

"""
function Optimizer(;
    determinantal_parameters::DeterminantalParameters{S, T},
    η::AbstractFloat, dt::AbstractFloat,
    jas_parameters::Union{Tuple{JastrowParameters{T}}, Nothing} = nothing
) where {S, T}
    # get the total number of variational parameters
    n_params = sum(length, determinantal_parameters.p)

    # collect all of the parameters
    var_params = reduce(vcat, determinantal_parameters.p)

    if !isnothing(jas_parameters)
        for jasp in jas_parameters
            for o in jasp.orbitals
                n_params += length(jasp.mean_v[o])
                var_params = vcat(var_params, jasp.mean_v[o])
            end
        end
    end

    @assert length(var_params) == n_params

    # initialize covariance matrix
    Sk = zeros(T, n_params, n_params)

    # initialize force vector
    fk = zeros(T, n_params)

    return Optimizer(var_params, η, dt, Sk, fk)
end


@doc raw"""

    update_optimizer!(;
        optimizer::Optimizer{T},
        bin_size::Int,
        determinantal_parameters::DeterminantalParameters{S,T},
        jas_parameters::Union{Tuple{JastrowParameters{T}}, Nothing} = nothing
    ) where {S, T<:Number}

Updates the optimizer by performing Stochastic Reconfiguration. The updated parameters are then pushed back to their respective containers.

# NOTE

During the optimization step, this function normalizes all measurements so that the optimizer
can construct `S` and `f`. Ensure that the `opt` field in the `write_measurements!` method is 
set to `true`.

"""
function update_optimizer!(;
    optimizer::Optimizer{T},
    measurement_container::NamedTuple,
    bin_size::Int,
    determinantal_parameters::DeterminantalParameters{S,T},
    jas_parameters::Union{Tuple{JastrowParameters{T}}, Nothing} = nothing
) where {S, T<:Number}
    # (; global_measurements, optimization_measurements) = measurement_container

    # normalize all measurements by the bin size
    normalize_measurements!(measurement_container, bin_size)

    # calculate S
    calculate_covariance_matrix!(optimizer, measurement_container.optimization_measurements)

    # calculate f
    calculate_forces!(optimizer, measurement_container.optimization_measurements, measurement_container.global_measurements)

    # perform gradient descent
    gradient_descent!(optimizer)

    # push back parameters to their respective containers
    idx = 1
    for pv in determinantal_parameters.p
        n = length(pv)
        pv .= @view optimizer.var_params[idx:idx+n-1]
        idx += n
    end

    if !isnothing(jas_parameters)
        for jasp in jas_parameters
            for orb in jasp.orbitals
                n = length(jasp.mean_v[orb])
                jasp.mean_v[orb] .= @view optimizer.var_params[idx:idx+n-1]
                idx += n
            end
        end
    end

    # add parameters to the measurement container
    measurement_container.optimization_measurements["parameters"] += optimizer.var_params

    return nothing
end


@doc raw"""

    calculate_covariance_matrix!(
        optimizer::Optimizer,
        optimization_measurements::Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}}
    ) where {T<:Number}

Calculates the covariance matrix ``S_{kk}' = \langle \Delta_{k}\Delta_{k^\prime}\rangle - \langle \Delta_{k} \rangle\langle \Delta_{k^\prime} \rangle``.

"""
function calculate_covariance_matrix!(
    optimizer::Optimizer,
    optimization_measurements::Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}}
) where {T<:Number}
    logDk       = real(optimization_measurements["logDk"])
    logDklogDkp = real(optimization_measurements["logDklogDkp"])

    # reset Sk before accumulating
    fill!(optimizer.S, zero(T))

    # Sk = logDklogDkp - logDk ⊗ logDk, written in-place
    copyto!(optimizer.S, logDklogDkp)
    BLAS.ger!(-one(T), logDk, logDk, optimizer.S)

    return nothing
end


@doc raw"""

    calculate_forces!(
        optimizer::Optimizer,
        optimization_measurements::Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}},
        global_measurements::Dict{String, Complex{T}}
    ) where {T<:Number}

Calculates the forces ``f_k = \langle \Delta_{k} \rangle\langle H\rangle - \langle \Delta_{k}H\rangle``.

"""
function calculate_forces!(
    optimizer::Optimizer,
    optimization_measurements::Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}},
    global_measurements::Dict{String, Complex{T}}
) where {T<:Number}
    logDk  = real(optimization_measurements["logDk"])
    logDkE = real(optimization_measurements["logDkE"])
    E_global = real(global_measurements["energy_per_site"])

    # reset f before accumulating
    fill!(optimizer.f, zero(T))

    @inbounds @simd for i in eachindex(logDk)
        optimizer.f[i] = logDk[i]*E_global - logDkE[i]
    end

    return nothing
end


@doc raw"""

    gradient_descent!(optimizer::Optimizer{T}) where {T}

Performs gradient descent on the set of variational parameters.

"""
function gradient_descent!(optimizer::Optimizer{T}) where {T}
    (; η, dt, S, f) = optimizer

    # stabilize S
    S_stable = copy(S)
    @inbounds for i in axes(S_stable, 1)
        S_stable[i, i] += η
    end

    # solve for parameter variation
    F = cholesky!(Symmetric(S_stable))  # cholesky!(Symmetric(S_stable, :L))
    δp = F \ f
    optimizer.var_params += dt * δp

    return nothing
end


