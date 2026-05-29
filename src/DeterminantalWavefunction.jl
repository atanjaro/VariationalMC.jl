@doc raw"""

    DeterminantalWavefunction{T<:Number, E<:AbstractMatrix{T}}

A mutable struct containing all the parameters needed to characterize the determinantal part of the trial wavefunction.

# Fields

- `W::E`: Equal-time Green's function matrix.
- `D::E`: Slater matrix.
- `M::E`: Reduced form of the unitary matrix `U`.
- `A::Vector{E}`: Set of logarithmic derivative operator matrices for all parameters that will be optimized. 
- `nq_updates_W::Int`: Tracker for the number of quick updates that have been performed on the equal-time Green's function.

"""
mutable struct DeterminantalWavefunction{T<:Number}
    # equal-time Green's function
    W::AbstractMatrix{T}

    # Slater matrix
    D::AbstractMatrix{T}
    
    # M matrix
    M::AbstractMatrix{T}

    # logarithmic derivative matrices
    A::Vector{Matrix{T}}
    
    # number of W matrix quick updates
    nq_updates_W::Int
end


@doc raw"""

    DeterminantalWavefunction(;
        # KEYWORD ARGUMENTS
        determinantal_parameters::DeterminantalParameters,
        hamiltonian::Hamiltonian{T,E},
        particle_configuration::ParticleConfiguration,
        model_geometry::ModelGeometry
    ) where {T, E}

Initialize and return an instance of `DeterminantalWavefunction`.

"""
function DeterminantalWavefunction(;
    # KEYWORD ARGUMENTS
    determinantal_parameters::DeterminantalParameters{S,T},
    hamiltonian::Hamiltonian{T,E},
    particle_configuration::ParticleConfiguration,
    model_geometry::ModelGeometry,
    rng::AbstractRNG
) where {S,T,E}
    Norbs = model_geometry.unit_cell.n
    Ncells = model_geometry.lattice.N
    N = Norbs * Ncells

    (; p, optimize) = determinantal_parameters
    (; V_p, U, ε₀) = hamiltonian
    (; Np, pconfig) = particle_configuration

    # previously, only certain `V` matrices were calculated along with their associated `A` matrix.
    # However, since the calculation of `V` is relatively cheap and is used to contruct `H`, now
    # all `V` matrices are calculated and saved, while only optimized parameters receive an `A` matrix.

    # allocate temporary matrices
    Nt = size(U, 1)
    T_buf   = Matrix{T}(undef, Nt, Np)
    Q_buf   = Matrix{T}(undef, Nt-Np, Np)
    Tmp_buf = Matrix{T}(undef, Nt, Np)
    A_buf   = Matrix{T}(undef, Nt, Nt)

    As = Matrix{T}[]

    V_p_idx = 1
    for (opt, params) in zip(optimize, p)
        if opt
            for _ in eachindex(params)
                A = calculate_derivative_operator_matrix(V_p[V_p_idx], U, T_buf, Q_buf, Tmp_buf, A_buf, ε₀, Np)
                push!(As, copy(A))
                V_p_idx += 1
            end
        else
            # skip over the matrices belonging to this parameter set
            V_p_idx += length(params)
        end
    end

    # generate a random particle configuration is one wasn't provided
    if iszero(pconfig)
        generate_random_fermion_configuration!(pconfig, particle_configuration.nup, particle_configuration.ndn, N, rng)
    end 

    # get reduced unitary matrix
    M = U[:, 1:Np]

    # preallocate matrices for the Slater and Green's function matrices
    D   = Matrix{T}(undef, Np, Np)
    W   = Matrix{T}(undef, 2*N, Np)

    # temporary storage
    Dt  = Matrix{T}(undef, Np, Np)
    tmp = Matrix{T}(undef, Np, 2*N)
    config_indices = Vector{Int}(undef, Np)

    overlap = calculate_equal_time_greens_function!(W, D, M, Dt, tmp, pconfig, config_indices, Np)

    while !overlap
        #@info "Proposed initial state has no overlap with the wavefunction. Retrying..."

        # get a new random configuration
        generate_random_fermion_configuration!(pconfig, particle_configuration.nup, particle_configuration.ndn, N, rng)

        # re-try getting Green's function
        overlap = calculate_equal_time_greens_function!(W, D, M, Dt, tmp, pconfig, config_indices, Np)
    end
    
    # intialize tracker for quick updating of the Green's function
    nq_updates_W = 0

    return DeterminantalWavefunction(W, D, M, As, nq_updates_W)
end