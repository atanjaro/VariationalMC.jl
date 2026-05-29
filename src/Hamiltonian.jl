@doc raw"""

    Hamiltonian{T<:Number}

A mutable struct containing all the parameters needed to characterize an auxiliary Hamiltonian used to 
construct a determinantal wavefunction. 

# Fields

- `V_t::AbstractMatrix{T}`: Hopping matrix constructed from the tight binding model.
- `V_p::Vector{AbstractMatrix{T}}`: Contains all matrices constructed from the determinantal order parameters.
- `U::AbstractMatrix{T}`: Matrix that diagonalizes the auxiliary Hamiltonian.
- `ε₀::Vector{T}`: Initial energies from the diagnolization of the auxiliary Hamiltonian.

"""
mutable struct Hamiltonian{T<:Number, M<:AbstractMatrix{T}}
    # hopping matrix
    V_t::M

    # vector of variational matrices
    V_p::Vector{M}

    # unitary matrix that diagonalizes H
    U::M

    # auxiliary energies
    ε₀::Vector{T}
end


@doc raw"""

    Hamiltonian(
        tight_binding_parameters::TightBindingParameters{T},
        determinantal_parameters::DeterminantalParameters{S,T},
        particle_configuration::ParticleConfiguration,
        model_geometry::ModelGeometry
    ) where {S,T}

Initialize and return an instance of `Hamiltonian`.

"""
function Hamiltonian(;
    # KEYWORD ARGUMENTS
    tight_binding_parameters::TightBindingParameters{T},
    determinantal_parameters::DeterminantalParameters{S, T},
    particle_configuration::ParticleConfiguration,
    model_geometry::ModelGeometry
) where {S, T}
    (; unit_cell, lattice) = model_geometry
    (; neighbor_table, bond_slices, t) = tight_binding_parameters
    (; order_type, param_name, symmetry, p, spdx_table) = determinantal_parameters
    (; ph_transform) = particle_configuration

    Norbs = unit_cell.n
    Ncells = lattice.N
    N =  Norbs * Ncells

    V_t = zeros(T, 2*N, 2*N)
    V_p = Matrix{T}[]

    # construct the hopping matrix
    build_hopping_matrix!(V_t, neighbor_table, bond_slices, t, N, ph_transform)

    # construct variational matrices, while preserving the order of the parameters as located in `determinantal_parameters`
    for (otype, pname, sym) in zip(order_type, param_name, symmetry)
        if otype == "charge"
            if sym == "site-dependent"
                nshifts = lattice.L[1]
                for shift in 0:(nshifts - 1)
                    V = zeros(T, 2*N, 2*N)
                    build_charge_ordering_matrix!(V, pname, sym, spdx_table, lattice, N, ph_transform, shift=shift)
                    push!(V_p, V)
                end
            else
                V = zeros(T, 2*N, 2*N)
                build_charge_ordering_matrix!(V, pname, sym, spdx_table, lattice, N, ph_transform)
                push!(V_p, V)
            end
        elseif otype == "spin"
            if sym == "site-dependent"
                nshifts = lattice.L[1]
                for shift in 0:(nshifts - 1)
                    V = zeros(T, 2*N, 2*N)
                    build_spin_ordering_matrix!(V, pname, sym, spdx_table, lattice, N, ph_transform, shift=shift)
                    push!(V_p, V)
                end
            else
                V = zeros(T, 2*N, 2*N)
                build_spin_ordering_matrix!(V, pname, sym, spdx_table, lattice, N, ph_transform)
                push!(V_p, V)
            end
        elseif otype == "pair"
            @assert ph_transform == true """
            ph_transform = $ph_transform. This must be enabled to add pairing matrices.
            """
            if sym == "site-dependent"
                for s in 1:N
                    # Fulde-Ferrell matrix for site s
                    V = zeros(T, 2*N, 2*N)
                    build_pairing_matrix!(V, pname, sym, neighbor_table, bond_slices, spdx_table, lattice, N,
                        site        = s,
                        qₚ          = determinantal_parameters.qₚ,
                        sym_subtype = "FF")
                    push!(V_p, V)
                end
                for s in 1:N
                    # Larkin-Ovchinnikov matrix for site s
                    V = zeros(T, 2*N, 2*N)
                    build_pairing_matrix!(V, pname, sym, neighbor_table, bond_slices, spdx_table, lattice, N,
                        site        = s,
                        qₚ          = determinantal_parameters.qₚ,
                        sym_subtype = "LO")
                    push!(V_p, V)
                end
            else
                V = zeros(T, 2*N, 2*N)
                build_pairing_matrix!(V, pname, sym, neighbor_table, bond_slices, spdx_table, lattice, N)
                push!(V_p, V)
            end
        else
            @error "Unknown order_type: \"$otype\". Expected \"charge\", \"spin\", or \"pair\"."
        end
    end

    # get full Hamiltonian
    H = copy(V_t)
    V_p_idx = 1
    for params in p
        for param in params
            H .+= param .* V_p[V_p_idx]
            V_p_idx += 1
        end
    end

    # solve for initial energies
    ε₀, U = diagonalize!(H)

    # check for open shell configuration
    is_openshell!(ε₀, particle_configuration.Np, ph_transform)

    return Hamiltonian(V_t, V_p, U, ε₀)
end 


@doc raw"""

    build_hopping_matrix!(
        V_t::AbstractMatrix{T},
        neighbor_table::Matrix{Int},
        bond_slices::Vector{UnitRange{Int}},
        t::Vector{E},
        N::Int,
        ph_transform::Bool
    ) where {T<:Number, E<:AbstractFloat}

Constructs the hopping matrix using parameters defined in `tight_binding_parameters`.

# ARGUMENTS

- `V_t::AbstractMatrix{T}`: Empty hopping matrix.
- `neighbor_table::Matrix{Int}`: Neighbor table containing all pairs of orbitals in the lattices connected by a bond, with a non-zero hopping energy between them.
- `bond_slices::Vector{UnitRange{Int}}`: Slices of `neighbor_table` corresponding to given bond ID.
- `t::Vector{E}`: The hopping energy ``t_{i,j}`` associated with each pair of neighboring orbitals connected by a bond in the lattice.
- `N::Int`: Total number of lattice sites.
- `ph_transform::Bool`: Whether the model is particle-hole transformed.

"""
function build_hopping_matrix!(
    V_t::AbstractMatrix{T},
    neighbor_table::Matrix{Int},
    bond_slices::Vector{UnitRange{Int}},
    t::Vector{E},
    N::Int,
    ph_transform::Bool
) where {T<:Number, E<:AbstractFloat}
    for slice in bond_slices
        # get slice of neighbor table
        nbr_slice = neighbor_table[:,slice]

        # get slice of hopping parameters
        t_slice = t[slice]

        # spin-up sector
        for ((i,j),t_hop) in zip(eachcol(nbr_slice),t_slice)
            V_t[i,j] += -t_hop
            V_t[j,i] += -t_hop
        end

        # spin-down sector
        for ((i,j),t_hop) in zip(eachcol(nbr_slice),t_slice)
            V_t[i+ N,j+ N] += ph_transform ? t_hop : -t_hop
            V_t[j+ N,i+ N] += ph_transform ? t_hop : -t_hop
        end
    end

    return nothing
end


@doc raw"""

    build_charge_ordering_matrix!(
        V::AbstractMatrix{T},
        pname::AbstractString,
        sym::AbstractString,
        spdx_table::Matrix{Int},
        lattice::Lattice,
        N::Int,
        ph_transform::Bool
    ) where {T}

Constructs a matrix representing a charge ordering operator in an auxiliary Hamiltonian.

# ARGUMENTS

- `V::AbstractMatrix{T}`: Temporary storage matrix.
- `pname::AbstractString`: Name of the charge order parameter.
- `sym::AbstractString`: Symmetry of the charge order.
- `spdx_table::Matrix{Int}`: Table of spindex pairs.
- `lattice::Lattice`: Instance of the `Lattice` type.
- `N::Int`: Total number of sites on the lattice.
- `ph_transform::Bool`: Whether the model is particle-hole transformed.

# KEYWORD ARGUMENTS

- `shift::Int = 0`: Accounts for the necessary shift required to `site-dependent` parameters such that each parameter will receive its own `V` matrix.

"""
function build_charge_ordering_matrix!(
    V::AbstractMatrix{T},
    pname::AbstractString,
    sym::AbstractString,
    spdx_table::Matrix{Int},
    lattice::Lattice,
    N::Int,
    ph_transform::Bool;
    shift::Int = 0
) where {T}
    if (pname, sym) == ("density", "uniform")
        dims = length(lattice.L)
        @inbounds for s in 1:N
            parity = sum(spdx_table[d, s] for d in 1:dims)
            v = T(isodd(parity) ? -1 : 1)
            V[s, s]     +=  v
            V[s+N, s+N] += ph_transform ? v : -v
        end
    elseif (pname, sym) == ("density", "site-dependent")
        @inbounds for idx in (1 + shift):lattice.L[1]:N
            V[idx, idx]     += -one(T)
            V[idx+N, idx+N] += ph_transform ? one(T) : -one(T)
        end
    elseif (pname, sym) == ("μ", "uniform")
        @inbounds for s in 1:2*N
            V[s, s] += ph_transform && s > N ? one(T) : -one(T)
        end
    else
        @error "Parameter = \"$pname\" and symmetry = \"$sym\" combination not recognized."
    end

    return nothing
end


@doc raw"""

    build_spin_ordering_matrix!(
        V::AbstractMatrix{T},
        pname::AbstractString,
        sym::AbstractString,
        spdx_table::Matrix{Int},
        lattice::Lattice,
        N::Int,
        ph_transform::Bool
    ) where {T}

Constructs a matrix representing a charge ordering operator in an auxiliary Hamiltonian.

# ARGUMENTS

- `V::AbstractMatrix{T}`: Temporary storage matrix.
- `pname::AbstractString`: Name of the spin order parameter.
- `sym::AbstractString`: Symmetry of the spin order.
- `spdx_table::Matrix{Int}`: Table of spindex pairs.
- `lattice::Lattice`: Instance of the `Lattice` type.
- `N::Int`: Total number of sites on the lattice.
- `ph_transform::Bool`: Whether the model is particle-hole transformed.

# KEYWORD ARGUMENTS

- `shift::Int = 0`: Accounts for the necessary shift required to `site-dependent` parameters such that each parameter will receive its own `V` matrix.

"""
function build_spin_ordering_matrix!(
    V::AbstractMatrix{T},
    pname::AbstractString,
    sym::AbstractString,
    spdx_table::Matrix{Int},
    lattice::Lattice,
    N::Int,
    ph_transform::Bool;
    shift::Int = 0
) where {T}
    if (pname, sym) == ("spin-x", "uniform")
        @inbounds for s in 1:N
            V[s, s + N] += one(T)
            V[s + N, s] += one(T)
        end
    elseif (pname, sym) == ("spin-z", "uniform")
        dims = length(lattice.L)
        @inbounds for s in 1:N
            parity = sum(spdx_table[d, s] for d in 1:dims)
            v = T(isodd(parity) ? -1 : 1)
            V[s, s]     +=  v
            V[s+N, s+N] += ph_transform ? -v : v
        end
    elseif (pname, sym) == ("spin-z", "site-dependent")
        @inbounds for idx in (1 + shift):lattice.L[1]:N
            V[idx, idx]     +=  one(T)
            V[idx+N, idx+N] += ph_transform ? -one(T) : one(T)
        end
    else
        @error "Combination of parameter \"$pname\" and symmetry \"$sym\" not recognized."
    end

    return nothing
end


@doc raw"""

    function build_pairing_matrix!(
        V::AbstractMatrix{T},
        pname::AbstractString,
        sym::AbstractString,
        neighbor_table::Matrix{Int},
        bond_slices::Vector{UnitRange{Int}},
        spdx_table::Matrix{Int},
        lattice::Lattice,
        N::Int;
        qₚ::Vector{<:AbstractFloat} = [0.0, 0.0],
        sym_subtype::Union{AbstractString, Nothing} = nothing
    ) where {T}

Constructs a matrix representing a charge ordering operator in an auxiliary Hamiltonian.

# ARGUMENTS

- `V::AbstractMatrix{T}`: Temporary storage matrix.
- `pname::AbstractString`: Name of the charge order parameter.
- `sym::AbstractString`: Symmetry of the charge order.
- `neighbor_table::Matrix{Int}`: Table of neighbor pairs.
- `bond_slices::Vector{UnitRange{Int}}`: Slices of `neighbor_table`.
- `spdx_table::Matrix{Int}`: Table of spindex pairs.
- `lattice::Lattice`: Instance of the `Lattice` type.
- `N::Int`: Total number of sites on the lattice.

# KEYWORD ARGUMENTS

- `site::Union{Int, Nothing} = nothing`: Site where the pair-density-wave parameter applies.
- `qₚ::Vector{<:AbstractFloat} = [0.0, 0.0]`: Center-of-mass pairing momentum.
- `sym_subtype::Union{AbstractString, Nothing} = nothing`: Symmetry type corresponding to pair-density-wave order, being either `FF` for Fulde-Furrell and `LO` for Larkin-Ovchinnikov.

"""
function build_pairing_matrix!(
    V::AbstractMatrix{T},
    pname::AbstractString,
    sym::AbstractString,
    neighbor_table::Matrix{Int},
    bond_slices::Vector{UnitRange{Int}},
    spdx_table::Matrix{Int},
    lattice::Lattice,
    N::Int;
    site::Union{Int, Nothing} = nothing,
    qₚ::Vector{<:AbstractFloat} = [0.0, 0.0],
    sym_subtype::Union{AbstractString, Nothing} = nothing
) where {T}
    dims = length(lattice.L)

    if (pname, sym) == ("s-wave", "uniform")
        @inbounds for r in 1:2*N
            c = get_linked_spindex(r - 1, N) + 1
            V[r, c] = one(T)
        end

    elseif (pname, sym) == ("s-wave", "site-dependent")
        @assert site !== nothing "site must be provided for site-dependent s-wave"
        θ = sum(qₚ[d] * spdx_table[d, site] for d in 1:dims)

        phase = if sym_subtype == "LO"
            cis(θ) + cis(-θ)
        elseif sym_subtype == "FF"
            cis(θ)
        else
            error("Symmetry sub-type \"$sym_subtype\" not recognized. Expected \"LO\" or \"FF\".")
        end

        si_up, si_dn = get_spindices_from_index(site, N)
        V[si_up, si_dn] += phase
        V[si_dn, si_up] += phase

    elseif (pname, sym) == ("d-wave", "uniform")
        for (slice, dsgn) in zip(bond_slices[1:2], (one(T), -one(T)))
            @inbounds for col in slice
                i = neighbor_table[1, col]
                j = neighbor_table[2, col]

                si_up, si_dn = get_spindices_from_index(i, N)
                sj_up, sj_dn = get_spindices_from_index(j, N)

                V[si_up, sj_dn] += dsgn
                V[sj_up, si_dn] += dsgn
                V[sj_dn, si_up] += dsgn
                V[si_dn, sj_up] += dsgn
            end
        end

    elseif (pname, sym) == ("d-wave", "site-dependent")
        @assert site !== nothing "site must be provided for site-dependent d-wave"
        θ = sum(qₚ[d] * spdx_table[d, site] for d in 1:dims)

        phase = if sym_subtype == "LO"
            cis(θ) + cis(-θ)
        elseif sym_subtype == "FF"
            cis(θ)
        else
            error("Symmetry sub-type \"$sym_subtype\" not recognized. Expected \"LO\" or \"FF\".")
        end

        for (slice, dsgn) in zip(bond_slices[1:2], (one(T), -one(T)))
            @inbounds for col in slice
                i = neighbor_table[1, col]
                j = neighbor_table[2, col]
                (i == site || j == site) || continue

                si_up, si_dn = get_spindices_from_index(i, N)
                sj_up, sj_dn = get_spindices_from_index(j, N)

                V[si_up, sj_dn] += dsgn * phase
                V[sj_dn, si_up] += dsgn * phase
                V[sj_up, si_dn] += dsgn * phase
                V[si_dn, sj_up] += dsgn * phase
            end
        end

    else
        error("Combination of parameter \"$pname\" and symmetry \"$sym\" not recognized.")
    end

    return nothing
end


@doc raw"""

    is_openshell!(ε₀::Vector{E}, Np::Int, ph_transform::Bool) where {E<:AbstractFloat}

Checks whether auxiliary energy configuration is open shell. If detected, will warn the user.

"""
function is_openshell!(ε₀::Vector{E}, Np::Int, ph_transform::Bool) where {E<:AbstractFloat}
    if ph_transform
        if ε₀[Np] - ε₀[Np - 1] < 1e-4
            @warn "Open shell detected!"
        end  
    else
        if ε₀[Np + 1] - ε₀[Np] < 1e-4 
            @warn "Open shell detected!"
        end 
    end 
    return nothing
end


@doc raw"""

    calculate_derivative_operator_matrix(
        v::M,
        U::M,
        T::M,      
        Q::M,
        Tmp::M,
        A::M,
        ε::Vector{E},
        Np::Int
    ) where {E<:Number, M<:AbstractMatrix{E}}

Computes the relevant lograithmic derivative operator matrix `A` from a variational parameter matrix `v`.

"""
function calculate_derivative_operator_matrix(
    v::AbstractMatrix{E},
    U::AbstractMatrix{E},
    T::AbstractMatrix{E},      
    Q::AbstractMatrix{E},
    Tmp::AbstractMatrix{E},
    A::AbstractMatrix{E},
    ε::Vector{E},
    Np::Int
) where {E<:Number}
    Nt  = size(U, 1)
    Np1 = Np + 1

    Uocc  = @view U[:, 1:Np]
    Uvirt = @view U[:, Np1:Nt]

    mul!(T, v, Uocc)                  # T   = v * Uocc        [Nt × Np]
    mul!(Q, Uvirt', T)                # Q   = Uvirt' * T      [Nt-Np × Np]

    # apply ptmask: zero out degenerate pairs, divide by energy difference otherwise
    @inbounds for j in axes(Q,2), i in axes(Q,1)
        εν = ε[j]
        εη = ε[i+Np]
        Q[i,j] = εν == εη ? zero(E) : Q[i,j] / (εν - εη)
    end

    mul!(Tmp, Uvirt, Q)               # Tmp = Uvirt * Q       [Nt × Np]
    mul!(A, Tmp, Uocc')               # A   = Tmp * Uocc'     [Nt × Nt]

    @inbounds @simd for i in eachindex(A)
        A[i] += conj(A[i])
    end

    return A
end

























































