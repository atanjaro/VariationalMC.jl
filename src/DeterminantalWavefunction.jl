@doc raw"""

    DeterminantalWavefunction( W::Matrix{ComplexF64}, 
                               D::Matrix{ComplexF64}, 
                               M::Matrix{ComplexF64}
                               U_int::Matrix{ComplexF64}, 
                               A::Vector{Any}, 
                               ε::Vector{Float64}, 
                               pconfig::Vector{Int64} )

A type defining quantities related to a determinantal wavefunction.

"""
mutable struct DeterminantalWavefunction
    # equal-time Green's function
    W::Matrix{ComplexF64}

    # Slater matrix
    D::Matrix{ComplexF64}
    
    # M matrix
    M::Matrix{ComplexF64}
    
    # U matrix (that diagonalizes H)
    U_int::Matrix{ComplexF64}
    
    # variational parameter matrices
    A::Vector{Any}
    
    # initial energies
    ε::Vector{Float64}

    # particle configuration
    pconfig::Vector{Int64}

    # number of W matrix quick updates
    nq_updates_W::Int64

    # number of T vector quick updates
    nq_updates_T::Int64
end


@doc raw"""

    get_determinantal_wavefunction( tight_binding_model::TightBindingModel, 
                                    determinantal_parameters::DeterminantalParameters, 
                                    Ne::Int64, 
                                    nup::Int64, 
                                    ndn::Int64,
                                    model_geometry::ModelGeometry, 
                                    rng::Xoshiro)::DeterminantalWavefunction

Constructs a variational wavefunction based on parameters given by the tight binding model and determinantal parameters. 
Returns an instances of the DeterminantalWavefunction type. If no particle configuration is specified, a random
configuration will be generated.                            

- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `Ne::Int`: total number of electrons.
- `nup::Int`: number of spin-up electrons.
- `ndn::Int`: number of spin-down electrons.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `rng::Xoshiro`: random number.
- `pconfig::Union{Nothing, Vector{Int}}=nothing`: current particle configuration. 

"""
function get_determinantal_wavefunction(
    tight_binding_model::TightBindingModel, 
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    Ne::Int, 
    nup::Int, 
    ndn::Int, 
    model_geometry::ModelGeometry, 
    rng::Xoshiro, 
    pconfig::Union{Nothing, Vector{Int}}=nothing
)::DeterminantalWavefunction
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
    (ε, U_int) = diagonalize(H)

    if is_openshell(ε, Ne)
        if isnothing(pconfig)
            error("   ERROR! Open shell detected! Exiting...")           
        else
            debug && println("   WARNING! Open shell detected!")
        end
    else
        debug && println("   Forming shell...")
    end

    # initialize variational parameter matrices
    A = get_variational_matrices(
        V, 
        U_int, 
        ε, 
        model_geometry
    )

    # get M matrix
    M = Matrix{ComplexF64}(view(U_int, 1:size(U_int,1), 1:Ne))

    # initialize Slater matrix
    D = zeros(ComplexF64, Ne, Ne)

    # initialize W matrix
    W = zeros(ComplexF64, 2*N, Ne)

     # use provided pconfig or generate a new one
     if isnothing(pconfig)
        pconfig = generate_initial_fermion_configuration(
            nup, 
            ndn, 
            model_geometry, 
            rng
        )
    end

    # # generate a random particle configuration
    # pconfig = generate_initial_fermion_configuration(nup, ndn, model_geometry, rng)

    # initialize equal-time Green's function and Slater matrix
    overlap = initialize_equal_time_greens!(
        W, 
        D, 
        M, 
        pconfig, 
        Ne
    )

    while overlap == false
        debug && println("Wavefunction::build_determinantal_wavefunction() : ")
        debug && println("configuration does not have ")
        debug && println("an overlap with the determinantal wavefunction => ")
        debug && println("generating a new configuration")

        # re-generate a random particle configuration
        pconfig = generate_initial_fermion_configuration(
            nup, 
            ndn, 
            model_geometry, 
            rng
        )

        # initialize W matrix
        W = zeros(ComplexF64, 2*N, Ne)

        # initialize Slater matrix
        D = zeros(ComplexF64, Ne, Ne)

        # re-initialize equal-time Green's function and Slater matrix
        overlap = initialize_equal_time_greens!(
            W, 
            D, 
            M, 
            pconfig, 
            Ne
        )
    end

    # intialize quick updating tracker
    nq_updates_W = 0
    nq_updates_T = 0

    return DeterminantalWavefunction(W, D, M, U_int, A, ε, pconfig, nq_updates_W, nq_updates_T)
end