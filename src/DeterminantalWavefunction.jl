@doc raw"""

    DeterminantalWavefunction{T<:Number, Q, E<:AbstractFloat, I<:Integer} 

A type defining quantities related to a determinantal wavefunction.

- `W::Matrix{T}`: equal-time Green's function matrix.
- `D::Matrix{T}`: Slater matrix.
- `M::AbstractMatrix{T}`: reduced U_aux matrix.
- `U_aux::Matrix{T}`: unitary matrix that diagonalizes the auxiliary Hamiltonian.
- `A::Vector{Q}`: variational parameter matrices.
- `ε::Vector{E}`: vector of mean-field energies.
- `pconfig::Vector{I}`: particle configuration.
- `nq_updates_W::I`: tracker for the number for quick updates to the W matrix.

"""
mutable struct DeterminantalWavefunction{T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # equal-time Green's function
    W::Matrix{T}

    # Slater matrix
    D::Matrix{T}
    
    # M matrix
    M::AbstractMatrix{T}
    
    # U matrix 
    U_aux::Matrix{T}
    
    # variational parameter matrices
    A::Vector{Q}
    
    # initial energies
    ε::Vector{E}

    # particle configuration
    pconfig::Vector{I}

    # number of W matrix quick updates
    nq_updates_W::I
end


@doc raw"""

    DeterminantalWavefunctionTABC{T<:Number, Q, E<:AbstractFloat, I<:Integer} 

A type defining quantities related to a determinantal wavefunction over a range 
of twist angles.

- `W_θ::Vector{Matrix{T}}`: set of equal-time Green's function matrices.
- `D_θ::Vector{Matrix{T}}`: set of Slater matrices.
- `M_θ::Vector{AbstractMatrix{T}}`: set of reduced U_aux matrices.
- `U_θ::Vector{Matrix{T}}`: set of unitary matrices that diagonalizes the auxiliary Hamiltonian.
- `A_θ::Vector{Vector{Q}}`: set of variational parameter matrices.
- `ε_θ::Vector{Vector{E}}`: vector of mean-field energies.
- `pconfig::Vector{I}`: particle configuration.
- `N_θ::I`: number of twist angles to average over.
- `twist_angles::AbstractRange{E}`: set of twist angles.
- `nq_updates_W::I`: tracker for the number for quick updates to the W matrix.

"""
mutable struct DeterminantalWavefunctionTABC{T<:Number, Q, E<:AbstractFloat, I<:Integer}
    # equal-time Green's functions
    W_θ::Vector{Matrix{T}}

    # Slater matrices
    D_θ::Vector{Matrix{T}}
    
    # M matrices
    M_θ::Vector{AbstractMatrix{T}}
    
    # U matrices
    U_θ::Vector{Matrix{T}}
    
    # sets variational parameter matrices
    A_θ::Vector{Vector{Q}}
    
    # initial energies
    ε_θ::Vector{Vector{E}}

    # particle configuration
    pconfig::Vector{I}

    # number of twiste angles
    N_θ::I
    
    # set of twist angles
    twist_angles::AbstractRange{E}

    # number of W matrix quick updates
    nq_updates_W::I
end


@doc raw"""

    get_determinantal_wavefunction( tight_binding_model::TightBindingModel{E}, 
                                    determinantal_parameters::DeterminantalParameters{I}, 
                                    optimize::NamedTuple, 
                                    Np::I, 
                                    nup::I, 
                                    ndn::I, 
                                    model_geometry::ModelGeometry, 
                                    rng::AbstractRNG, 
                                    pht::Bool,
                                    pconfig::Vector{I} = Int[] ) where {E<:AbstractFloat, I<:Integer}

Constructs a variational wavefunction ``|\Phi_0\rangle`` based on parameters given by the tight-binding model 
and determinantal parameter and returns an instance of the `DeterminantalWavefunction` type. If no initial 
particle configuration is specified, a random configuration will be generated.                            

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `Np::I`: total number of particles in the system.
- `nup::I`: number of spin-up electrons.
- `ndn::I`: number of spin-down electrons.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `rng::AbstractRNG`: random number generator. 
- `pht::Bool`: whether model is particle-hole transformed.
- `pconfig::Vector{I} = Int[]`: optional initial particle configuration.

"""
function get_determinantal_wavefunction(
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    Np::I, 
    nup::I, 
    ndn::I, 
    model_geometry::ModelGeometry, 
    rng::AbstractRNG, 
    pht::Bool,
    pconfig::Vector{I} = Int[]
) where {E<:AbstractFloat, I<:Integer}
    # number of lattice sites
    N = model_geometry.lattice.N 

    # build auxiliary (mean-field) Hamiltonian and variational operators
    (H, V) = build_auxiliary_hamiltonian(
        tight_binding_model, 
        determinantal_parameters, 
        optimize, 
        model_geometry, 
        pht
    )

    # diagonalize Hamiltonian
    (ε, U_aux) = diagonalize!(H)

    if is_openshell(ε, Np, pht) 
        # @debug """
        # DeterminantalWavefunction::build_determinantal_wavefunction() : 
        # WARNING! Open shell detected!
        # """
        println("DeterminantalWavefunction::build_determinantal_wavefunction() :")
        println("WARNING: Open shell is detected!")
    else
        @debug """
        DeterminantalWavefunction::build_determinantal_wavefunction() : 
        Closed shell is detected! Forming shell..."""
    end

    # initialize variational parameter matrices
    A = get_variational_matrices(
        V, 
        U_aux, 
        ε, 
        Np,
        model_geometry
    )

    # get M matrix
    @views M = U_aux[:, 1:Np]
    # M = Matrix{ComplexF64}(view(U_aux, 1:size(U_aux,1), 1:Np))

    # initialize Slater matrix
    D = zeros(ComplexF64, Np, Np)

    # initialize W matrix
    W = zeros(ComplexF64, 2*N, Np)

    # vector to store configuration indices
    config_indices = Vector{Int}(undef, Np)     

    # helper matrices
    Dt = zeros(ComplexF64, Np, Np)            
    tmp = zeros(ComplexF64, Np, 2*N)     
    
    # generate a new configuration at random if none was provided
    if isempty(pconfig)
        resize!(pconfig, 2 * N)
        generate_initial_fermion_configuration!(pconfig, nup, ndn, model_geometry, rng)
    end 

    # initialize equal-time Green's function and Slater matrix
    overlap = initialize_equal_time_greens!(
        W, 
        D, 
        M, 
        pconfig, 
        Np,
        config_indices,
        Dt,
        tmp
    )

    # if starting configuration was used, check whether accepted
    if !isempty(pconfig) && overlap
        @debug """DeterminantalWavefunction::build_determinantal_wavefunction() : 
        Initial configuration accepted!
        """
    elseif !isempty(pconfig) && !overlap
        @debug """DeterminantalWavefunction::build_determinantal_wavefunction() : 
        Initial configuration rejected!
        """
    end

    while !overlap
        @debug """
        DeterminantalWavefunction::build_determinantal_wavefunction() : 
        configuration does not have a finite overlap with the current determinantal wavefunction => 
        generating a new configuration"
        """

        # re-generate a random particle configuration
        generate_initial_fermion_configuration!(
            pconfig,
            nup, 
            ndn, 
            model_geometry, 
            rng
        )
 
        # re-initialize equal-time Green's function
        overlap = initialize_equal_time_greens!(
            W, 
            D, 
            M, 
            pconfig, 
            Np,
            config_indices,
            Dt,
            tmp
        )
    end

    # intialize quick updating tracker
    nq_updates_W = 0

    return DeterminantalWavefunction(W, D, M, U_aux, A, ε, pconfig, nq_updates_W)
end

# TABC prototypes
N_θ = 10
twist_angles = range(-π, π; length=N_θ+1)[1:end-1] .+ π/N_θ

@doc raw"""

    get_determinantal_wavefunction( tight_binding_model::TightBindingModel{E}, 
                                    determinantal_parameters::DeterminantalParameters{I}, 
                                    optimize::NamedTuple, 
                                    Np::I, 
                                    nup::I, 
                                    ndn::I, 
                                    model_geometry::ModelGeometry, 
                                    N_θ::I,
                                    twist_angles::AbstractRange{E},
                                    rng::AbstractRNG, 
                                    pht::Bool,
                                    pconfig::Vector{I} = Int[] ) where {E<:AbstractFloat, I<:Integer}

Constructs a variational wavefunction ``|\Phi_{0}^{\theta}\rangle`` for `N_\theta` twist angles based on parameters 
given by the tight-binding model and determinantal parameter and returns an instance of the `DeterminantalWavefunctionTABC` 
type. If no initial particle configuration is specified, a random configuration will be generated.                            

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `Np::I`: total number of particles in the system.
- `nup::I`: number of spin-up electrons.
- `ndn::I`: number of spin-down electrons.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `N_θ::I`: number of twist angles.
- `twist_angles::AbstractRange{E}`: list of twist angles.
- `rng::AbstractRNG`: random number. 
- `pht::Bool`: whether model is particle-hole transformed.
- `pconfig::Vector{I} = Int[]`: optional initial particle configuration.

"""
function get_determinantal_wavefunction(
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    Np::I, 
    nup::I, 
    ndn::I, 
    model_geometry::ModelGeometry, 
    N_θ::I,
    twist_angles::AbstractRange{E},
    rng::AbstractRNG, 
    pht::Bool,
    pconfig::Vector{I} = Int[]
) where {E<:AbstractFloat, I<:Integer}
    N = model_geometry.lattice.N 
    ε_θ = Vector{Vector{AbstractFloat}}()
    U_θ = Vector{Matrix{<:Number}}()
    A_θ = Vector{Vector{Matrix{<:Number}}}()
    D_θ = Vector{Matrix{<:Number}}()
    W_θ = Vector{Matrix{<:Number}}()
    M_θ = Vector{Matrix{<:Number}}()

    # build auxiliary Hamiltonian matrices for each twist angle
    H_θ, V = build_auxiliary_hamiltonian(
        tight_binding_model, 
        determinantal_parameters, 
        optimize, 
        model_geometry, 
        twist_angles, 
        pht
    )

    # select starting configuration 
    config_indices = Vector{Int}(undef, Np)
    if isempty(pconfig)
        resize!(pconfig, 2 * N)
        generate_initial_fermion_configuration!(pconfig, nup, ndn, model_geometry, rng)
    end 

    # get wavefunction quantities for each twist angle
    for n in 1:N_θ
        # get mean-field energy and unitary matrix
        ε, U_aux = diagonalize!(H_θ[n])
        push!(ε_θ, ε)
        push!(U_θ, U_aux)

        if is_openshell(ε, Np, pht)
            # @debug """
            # DeterminantalWavefunction::build_determinantal_wavefunction() : 
            # WARNING! Open shell detected!
            # """
            println("DeterminantalWavefunction::build_determinantal_wavefunction() :")
            println("WARNING! Open shell is detected in state $(n)!")
        else
            @debug """
            DeterminantalWavefunction::build_determinantal_wavefunction() : 
            Closed shell is detected! Forming shell..."""
        end

        # get variaitonal matrices
        A = get_variational_matrices(
            V, 
            U_aux, 
            ε, 
            Np,
            model_geometry
        )
        push!(A_θ, A)

        # get M matrix
        @views M = U_aux[:, 1:Np]
        push!(M_θ, M)

        D = zeros(ComplexF64, Np, Np)
        W = zeros(ComplexF64, 2*N, Np)
        Dt = zeros(ComplexF64, Np, Np)            
        tmp = zeros(ComplexF64, Np, 2*N)     
        
        # get equal-time Green's function
        overlap = initialize_equal_time_greens!(
            W, 
            D, 
            M, 
            pconfig, 
            Np,
            config_indices,
            Dt,
            tmp
        )
        push!(D_θ, D)
        push!(W_θ, W)

        # warn if there is no overlap
        if !overlap
            println("WARNING! Configuration has no overlap with wavefunction number $n.")
        end
        
        # if starting configuration was used, check whether accepted
        if !isempty(pconfig) && overlap
            @debug """DeterminantalWavefunction::build_determinantal_wavefunction() : 
            Initial configuration accepted!
            """
        elseif !isempty(pconfig) && !overlap
            @debug """DeterminantalWavefunction::build_determinantal_wavefunction() : 
            Initial configuration rejected!
            """
        end
    end

    # intialize quick updating tracker
    nq_updates_W = 0
    
    return DeterminantalWavefunctionTABC(W_θ, D_θ, M_θ, U_θ, A_θ, ε_θ, pconfig, N_θ, twist_angles, nq_updates_W)

end



