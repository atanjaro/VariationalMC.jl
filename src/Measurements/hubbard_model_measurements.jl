@doc raw"""

    measure_hubbard_energy(
        hubbard_parameters::HubbardParameters{E},
        hubbard_id::Int,
        pconfig::Vector{Int},
        ph_transform::Bool
    ) where {E<:AbstractFloat}

Calculates the average Hubbard energy ``U \langle \hat{n}_\uparrow \hat{n}_\downarrow \rangle``.

# ARGUMENTS

- `hubbard_parameters::HubbardParameters{E}`: Instance of `HubbardParameters`.
- `hubbard_id::Int`: ID of the Hubbard interaction.
- `pconfig::Vector{Int}`: Current particle configuration.
- `ph_transform::Bool`: Whether the model is particle-hole transformed.

"""
function measure_hubbard_energy(
    hubbard_parameters::HubbardParameters{E},
    hubbard_id::Int,
    pconfig::Vector{Int}, ph_transform::Bool
) where {E<:AbstractFloat}
    (; U, sites, orbital_ids) = hubbard_parameters

    # initialize the Hubbard energy
    ε_hubb = zero(E)

    # number of orbitals in the unit cell with finite hubbard interaction
    n_hubbard = length(orbital_ids)

    # number of Hubbard interactions in the lattice
    N_hubbard = length(U)

    # get number of unit cells in the lattice
    N_unitcells = N_hubbard ÷ n_hubbard

    # reshape so each column corresponds to a given orbital with finite hubbard interaction
    U′ = reshape(U, (N_unitcells, n_hubbard))
    sites′ = reshape(sites, (N_unitcells, n_hubbard))

    @fastmath @inbounds for i in axes(U′,1)
        site = sites′[i,hubbard_id]
        nup, ndn = get_fermion_occupations(site, pconfig, Int(length(pconfig)/2))
        ε_hubb += ph_transform ? U′[i, hubbard_id] * nup * (1 - ndn) : U′[i, hubbard_id] * nup * ndn
    end
    
    # normalize the measurement
    ε_hubb /= N_unitcells

    return ε_hubb 
end


@doc raw"""

    measure_ext_hubbard_energy(
        extended_hubbard_parameters::ExtendedHubbardParameters{E},
        ext_hubbard_id::Int,
        pconfig::Vector{Int},
        ph_transform::Bool
    ) where {E<:AbstractFloat}

Calculates the average extended Hubbard interaction energy ``V \hat{n}_i \hat{n}_j``.

# ARGUMENTS

- `extended_hubbard_parameters::ExtendedHubbardParameters{E}`: Instance of `ExtendedHubbardParameters`.
- `ext_hubbard_id::Int`: ID of the extended Hubbard interaction.
- `pconfig::Vector{Int}`: Current particle configuration.
- `ph_transform::Bool`: Whether model is particle-hole transformed.

"""
function measure_ext_hubbard_energy(
    extended_hubbard_parameters::ExtendedHubbardParameters{E},
    ext_hubbard_id::Int,
    pconfig::Vector{Int}, ph_transform::Bool
) where {E<:AbstractFloat}
    (; V, neighbor_table, bond_ids) = extended_hubbard_parameters

    # number of types of extended hubbard interactions
    n = length(bond_ids)

    # number of unit cells
    N = size(neighbor_table, 2) ÷ n

    # get appropriate view in interaction
    V′ = reshape(V, (N, n))
    V′ = @view V′[:,ext_hubbard_id]

    # get appropriate view into neighbor table
    nt = reshape(neighbor_table, (2, N, n))
    nt = @view nt[:,:,ext_hubbard_id]

    # initialize measurement to zero
    ε_ext_hubb = zero(E)

    # iterate over coupled sites
    for n in 1:N
        # get the pair of sites that are coupled
        i = nt[1,n]
        j = nt[2,n]

        # get the relevant densities
        n_i_up, n_i_dn = get_fermion_occupations(i, pconfig, Int(length(pconfig)/2))
        n_j_up, n_j_dn = get_fermion_occupations(j, pconfig, Int(length(pconfig)/2))

        # get interaction strength
        Vij = V′[n]

        # calculate interaction energy
        ε_ext_hubb += ph_transform ? Vij * (n_i_up - n_i_dn) * (n_j_up - n_j_dn) : Vij * (n_i_up + n_i_dn) * (n_j_up + n_j_dn)
    end

    # normalize the energy
    ε_ext_hubb /= N

    return ε_ext_hubb
end
