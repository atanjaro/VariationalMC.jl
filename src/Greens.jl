@doc raw"""

    initialize_equal_time_greens( W::Matrix{ComplexF64}, 
                                  D::Matrix{ComplexF64}, 
                                  M::Matrix{ComplexF64}, 
                                  pconfig::Vector{Int64}, 
                                  N::Int64, 
                                  Ne::Int64 )::Bool
    
Computes the equal-time Green's function by solving the matrix equation DᵀWᵀ = Mᵀ through
LU decomposition.

- `W::Matrix{ComplexF64}`: equal-time Green's function matrix.
- `D::Matrix{ComplexF64}`: Slater determinant matrix.
- `M::Matrix{ComplexF64}`: reduced U_int matrix.
- `pconfig::Vector{Int64}`: current particle configuration.
- `Ne::Int64`: total number of electrons.

"""
function initialize_equal_time_greens!(
    W::Matrix{ComplexF64}, 
    D::Matrix{ComplexF64}, 
    M::Matrix{ComplexF64}, 
    pconfig::Vector{Int64}, 
    Ne::Int64
)::Bool
    # get indices from the particle configuration
    config_indices = [findfirst(==(i), pconfig) for i in 1:Ne]

    # get Slater matrix
    D .= M[config_indices, :] 

    if abs(det(D)) < 1e-12 * size(D, 1) 
        debug && println("Wavefunction::initialize_equal_time_greens() : state has no")
        debug && println("overlap with the determinantal wavefunction, ")
        debug && println("D = ")
        debug && display(D)
        debug && println("determinant of D = ", det(D))

        return false
    else        
        # LU decomposition of D'
        lu_decomp = lu(D')
        
        # calculate the equal-time Green's function
        W .= transpose(lu_decomp \ transpose(M))

        return true
    end            
end


@doc raw"""

    function update_equal_time_greens!( markov_move::MarkovMove, 
                                        detwf::DeterminantalWavefunction, 
                                        model_geometry::ModelGeometry,
                                        Ne::Int64, 
                                        n_stab_W::Int64, 
                                        δW::Float64 )::Nothing

Updates the equal-time Green's function while performing a numerical stabilzation check. If the calculated 
deviation exceeds the set threshold, then the current Green's function is replaced by one calculated from scratch.

- `markov_move::MarkovMove`: quantities related to a Markov move.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int64`: total number of electrons.
- `n_stab_W::Int64`: frequency of Green's function stability steps.
- `δW::Float64`: Green's function stability threshold.

"""
function update_equal_time_greens!(
    markov_move::MarkovMove, 
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry,
    Ne::Int64, 
    n_stab_W::Int64, 
    δW::Float64
)::Nothing
    if detwf.nq_updates_W >= n_stab_W
        debug && println("Wavefunction::update_equal_greens!() : recalculating W!")

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
        Wᵣ = zeros(ComplexF64, 2*N, Ne)

        # re-initialize Slater matrix
        Dᵣ = zeros(ComplexF64, Ne, Ne)

        # recalclate Green's function from scratch
        recalculate_equal_time_greens!(
            Wᵣ, 
            Dᵣ, 
            detwf.M, 
            detwf.pconfig, 
            Ne
        )

        # compute deviation between current Green's function and the recalculated Green's function
        dev = check_deviation(
            detwf, 
            Wᵣ
        )

        debug && println("Wavefunction:update_equal_greens!() : recalculated W with deviation = ", dev)

        debug && println("Wavefunction::update_equal_greens!() : deviation goal for matrix")

        if dev > δW
            debug && println("W not met!")
            debug && println("Wavefunction::update_equal_greens!() : updated W = ")
            debug && display(detwf.W)
            debug && println("Wavefunction::update_equal_greens!() : exact W = ")
            debug && display(Wᵣ)

            # replace original W matrix with new one
            detwf.W = Wᵣ
        else
            debug && println("W met! Green's function is stable")
            @assert dev < δW
        end  

        return nothing
    else
        debug && println("Wavefunction::update_equal_greens!() : performing rank-1 update of W!") 

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

    rank1_update!( markov_move::MarkovMove, 
                   detwf::DeterminantalWavefunction )::Nothing
    
Performs and in-place rank-1 update of the equal-time Green's function. The default method 
uses BLAS.geru!(). Also available is a DEBUG method which performs the update by hand.

- `markov_move::MarkovMove`: quantities related to a Markov move. 
- `detwf::DeterminantalWavefunction`: current variational wavefunction. 

"""
function rank1_update!(
    markov_move::MarkovMove, 
    detwf::DeterminantalWavefunction
)::Nothing
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

    # # Perform rank-1 update by hand
    # detwf.W -= (cᵦ / detwf.W[l, β]) * rₗ' 

    # Perform rank-1 update using BLAS.geru!
    BLAS.geru!(
        -1.0 / detwf.W[l, β], 
        cᵦ, 
        rₗ, 
        detwf.W
    )

    return nothing
end


@doc raw"""

    recalculate_equal_time_greens( Wᵣ::Matrix{ComplexF64}, 
                                   Dᵣ::Matrix{ComplexF64}, 
                                   M::Matrix{ComplexF64}, 
                                   pconfig::Vector{Int64}, 
                                   Ne::Int64 )::Bool
    
Recomputes the equal-time Green's function.

- `Wᵣ::Matrix{ComplexF64}`: equal-time Green's function matrix.
- `Dᵣ::Matrix{ComplexF64}`: Slater determinant matrix.
- `M::Matrix{ComplexF64}`: reduced U_int matrix.
- `pconfig::Vector{Int64}`: current particle configuration.
- `Ne::Int64`: total number of electrons. 

"""
function recalculate_equal_time_greens!(
    Wᵣ::Matrix{ComplexF64}, 
    Dᵣ::Matrix{ComplexF64}, 
    M::Matrix{ComplexF64}, 
    pconfig::Vector{Int64}, 
    Ne::Int64
)::Bool
    # get indices from the particle configuration
    config_indices = [findfirst(==(i), pconfig) for i in 1:Ne] 

    # get Slater matrix
    Dᵣ .= M[config_indices, :]

    # LU decomposition of D'
    lu_decomp = lu(Dᵣ')

    # calculate the equal-time Green's function
    Wᵣ .= transpose(lu_decomp \ transpose(M))

    return true      
end


@doc raw"""

    check_deviation!( detwf::DeterminantalWavefunction, 
                      Wᵣ::Matrix{ComplexF64} )::Float64 
    
Checks floating point error accumulation in the equal-time Green's function.

- `detwf::DeterminantalWavefunction`: current variational wavefunction. 
- `Wᵣ::Matrix{ComplexF64}`: reclculated Green's function matrix

"""
function check_deviation(
    detwf::DeterminantalWavefunction, 
    Wᵣ::Matrix{ComplexF64}
)::Float64    
    # Difference in updated Green's function and recalculated Green's function
    difference = detwf.W .- Wᵣ

    # Sum the absolute differences and the recalculated Green's function elements
    diff_sum = sum(abs.(difference))
    W_sum = sum(abs.(Wᵣ))

    # condition for recalculation
    ΔW = sqrt(diff_sum / W_sum)

    return ΔW
end