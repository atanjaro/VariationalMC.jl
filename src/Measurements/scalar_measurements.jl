@doc raw"""

    measure_double_occ(
        pconfig::Vector{Int},
        N::Int,
        ph_transform::Bool
    )

Measures the average double occupancy ``\langle \hat{n}_\uparrow \hat{n}_\downarrow \rangle`` given the current particle configuration. 

# ARGUMENTS

-`pconfig::Vector{Int}`: Current particle configuration.
-`N::Int`: Total number of sites on the lattice.
-`ph_transform::Bool`: Whether the model is particle-hole transformed.

"""
function measure_double_occ(
    pconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
)
    nup_ndn = 0.0
    for site in 1:N
        nup, ndn = get_fermion_occupations(site, pconfig, N)
        nup_ndn += ph_transform ? nup * (1 - ndn) : nup * ndn
    end

    nup_ndn /= N

    return nup_ndn
end


@doc raw"""

    measure_double_occ(
        pconfig::Vector{Int},
        o::Int,
        Ncells::Int,
        unit_cell::UnitCell,
        ph_transform::Bool
    )

Measures the average double occupancy ``\langle \hat{n}_\uparrow \hat{n}_\downarrow \rangle`` for a 
given orbital species `o` in the current particle configuration. 

# ARGUMENTS

-`pconfig::Vector{Int}`: Current particle configuration.
-`unit_cell::UnitCell`: Instance of `UnitCell`.
-`o::Int`: Orbital ID.
-`Ncells::Int`: Total number of unit cells on the lattice.
-`N::Int`: Total number of sites on the lattice.
-`ph_transform::Bool`: Whether the model is particle-hole transformed.

"""
function measure_double_occ(
    pconfig::Vector{Int},
    unit_cell::UnitCell,
    o::Int,
    Ncells::Int,
    N::Int,
    ph_transform::Bool
)
    nup_ndn = 0.0
    for uc in 1:Ncells
        site = loc_to_site(uc, o, unit_cell)
        nup, ndn = get_fermion_occupations(site, pconfig, N)
        nup_ndn += ph_transform ? nup * (1 - ndn) : nup * ndn
    end

    nup_ndn /= N

    return nup_ndn
end


@doc raw"""

    measure_n(
        type::AbstractString,
        pconfig::Vector{Int},
        N::Int,
        ph_transform::Bool
    )

Measures the average (or site-dependent) density ``\langle \hat{n}_\sigma \rangle`` given the current particle configuration. 

# ARGUMENTS

-`type::AbstractString`: Type of density observable being calculated.
-`pconfig::Vector{Int}`: Current particle configuration.
-`N::Int`: Total number of sites on the lattice.
-`ph_transform::Bool`: Whether the model is particle-hole transformed.

# NOTE 

The valid options for `type` are `"mean"` and `"site-dependent"`.

"""
function measure_n(
    type::AbstractString,
    pconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
)
    if type == "mean"
        n̄ = 0.0
        for site in 1:N
            nup, ndn = get_fermion_occupations(site, pconfig, N)
            n̄ += ph_transform ? nup + (1 - ndn) : nup + ndn
        end
        n̄ /= N
    elseif type == "site-dependent"
        up_occs = Int.(pconfig[1:N] .!= 0)
        dn_occs = Int.(pconfig[N+1:2*N] .!= 0)

        n̄ = ph_transform ? up_occs .- dn_occs .+ 1 : up_occs .+ dn_occs
    end

    return n̄
end


@doc raw"""

    measure_n(
        type::AbstractString,
        pconfig::Vector{Int},
        o::Int,
        Ncells::Int,
        unit_cell::UnitCell,
        ph_transform::Bool
    )

Measures the average density ``\langle \hat{n}_\sigma \rangle`` for a given orbital species `o` in the current particle configuration. 

# ARGUMENTS

-`pconfig::Vector{Int}`: Current particle configuration.
-`unit_cell::UnitCell`: Instance of `UnitCell`.
-`o::Int`: Orbital ID.
-`Ncells::Int`: Total number of unit cells on the lattice.
-`N::Int`: Total number of sites on the lattice.
-`ph_transform::Bool`: Whether the model is particle-hole transformed.

# NOTE 

The valid options for `type` are `"mean"` and `"site-dependent"`.

"""
function measure_n(
    type::AbstractString,
    pconfig::Vector{Int},
    unit_cell::UnitCell,
    o::Int,
    Ncells::Int,
    N::Int,
    ph_transform::Bool
)
    if type == "mean"
        n̄ = 0.0
        for uc in 1:Ncells
            site = loc_to_site(uc, o, unit_cell)
            nup, ndn = get_fermion_occupations(site, pconfig, N)
            n̄ += ph_transform ? nup + (1 - ndn) : nup + ndn
        end
        n̄ /= N
    elseif type == "site-dependent"
        # get sites corresponding to orbital ID
        orb_sites = [(uc - 1) * unit_cell.n + o for uc in 1:Ncells]

        up_occs = Int.(pconfig[orb_sites] .!= 0)
        dn_occs = Int.(pconfig[orb_sites .+ N] .!= 0)

        n̄ = ph_transform ? up_occs .- dn_occs .+ 1 : up_occs .+ dn_occs
    end

    return n̄
end


@doc raw"""

    measure_spinz(
        type::AbstractString,
        pconfig::Vector{Int},
        N::Int,
        ph_transform::Bool
    )   

Measures the average (or site-dependent) z-component of the spin ``\langle \hat{n}_{\uparrow} - \hat{n}_{\downarrow} \rangle/2`` given the current particle configuration. 

# ARGUMENTS

-`type::AbstractString`: Type of spin-z observable being calculated.
-`pconfig::Vector{Int}`: Current particle configuration.
-`N::Int`: Total number of sites on the lattice.
-`ph_transform::Bool`: Whether the model is particle-hole transformed.

# NOTE 

The valid options for `type` are `"mean"` and `"site-dependent"`.

"""
function measure_spinz(
    type::AbstractString,
    pconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
)   
    if type =="mean"
        Sz = 0.0
        for site in 1:N
            nup, ndn = get_fermion_occupations(site, pconfig, N)
            Sz += ph_transform ? 0.5 * (nup + ndn - 1) : 0.5 * (nup - ndn)
        end
        Sz /= N
    elseif type == "site-dependent"
        up_occs = Int.(pconfig[1:N] .!= 0)
        dn_occs = Int.(pconfig[N+1:2*N] .!= 0)

        Sz = ph_transform ? 0.5 * (up_occs .+ dn_occs .- 1) : 0.5 * (up_occs .- dn_occs)
    end

    return Sz
end


@doc raw"""

    measure_spinz(
        pconfig::Vector{Int},
        unit_cell::UnitCell,
        o::Int,
        Ncells::Int,
        N::Int,
        ph_transform::Bool
    )   

Measures the average z-component of the spin ``\langle \hat{n}_{\uparrow} - \hat{n}_{\downarrow} \rangle`` for a given orbital species `o` in the current particle configuration. 

# ARGUMENTS

-`pconfig::Vector{Int}`: Current particle configuration.
-`unit_cell::UnitCell`: Instance of `UnitCell`.
-`o::Int`: Orbital ID.
-`Ncells::Int`: Total number of unit cells on the lattice.
-`N::Int`: Total number of sites on the lattice.
-`ph_transform::Bool`: Whether the model is particle-hole transformed.

# NOTE 

The valid options for `type` are `"mean"` and `"site-dependent"`.

"""
function measure_spinz(
    type::AbstractString,
    pconfig::Vector{Int},
    unit_cell::UnitCell,
    o::Int,
    Ncells::Int,
    N::Int,
    ph_transform::Bool
)   
    if type == "mean"
        Sz = 0.0
        for uc in 1:Ncells
            site = loc_to_site(uc, o, unit_cell)
            nup, ndn = get_fermion_occupations(site, pconfig, N)
            Sz += ph_transform ? 0.5 * (nup + ndn - 1) : 0.5 * (nup - ndn)
        end

        Sz /= N
    elseif type == "site-dependent"
        # get sites corresponding to orbital ID
        orb_sites = [(uc - 1) * unit_cell.n + o for uc in 1:Ncells]

        up_occs = Int.(pconfig[orb_sites] .!= 0)
        dn_occs = Int.(pconfig[orb_sites .+ N] .!= 0)

        Sz = ph_transform ? 0.5 * (up_occs .+ dn_occs .- 1) : 0.5 * (up_occs .- dn_occs)
    end

    return Sz
end


@doc raw"""

    measure_variational_energy(
        detwf::DeterminantalWavefunction{T},
        tight_binding_parameters::TightBindingParameters{T, E},
        coupling_parameters::Tuple,
        pconfig::Vector{Int},
        ph_transform::Bool;
        model_geometry::ModelGeometry,
        jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
        jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing
    ) where {T<:Number, E<:AbstractFloat}

Measures the average variational energy per site.

"""
function measure_variational_energy(
    detwf::DeterminantalWavefunction{T},
    tight_binding_parameters::TightBindingParameters{T, E},
    coupling_parameters::Tuple,
    pconfig::Vector{Int}, N::Int, ph_transform::Bool;
    model_geometry::ModelGeometry,
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
    jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing
) where {T<:Number, E<:AbstractFloat}
    # number of orbitals per unit cell
    norbital = tight_binding_parameters.norbital

    # number of types of hopping
    bond_ids = tight_binding_parameters.bond_ids
    nhopping = length(tight_binding_parameters.bond_ids)

    E_var = zero(Complex{T})

    # iterate over the hoppings in the lattice
    if nhopping > 0
        for hopping_id in 1:nhopping
            # measure the hopping energy
            e_hop = measure_hopping_energy(
                detwf,
                tight_binding_parameters,
                pconfig,
                hopping_id,
                N,
                ph_transform,
                model_geometry = model_geometry,
                jas_factors = jas_factors,
                jas_parameters = jas_parameters
            )

            E_var += e_hop
        end
    end

    # iterate over the couplings
    for coupling in coupling_parameters
        e_coupling = measure_coupling_energy(
            coupling,
            pconfig, ph_transform
        )
        E_var += e_coupling
    end

    return E_var
end


@doc raw"""

    measure_coupling_energy(
        hubbard_parameters::HubbardParameters{E},
        pconfig::Vector{Int}, ph_transform::Bool
    ) where {E<:AbstractFloat}

Measures the average energy per site due to coupling in the Hubbard model.

"""
function measure_coupling_energy(
    hubbard_parameters::HubbardParameters{E},
    pconfig::Vector{Int}, ph_transform::Bool
) where {E<:AbstractFloat}
    E_hub = zero(Complex{E})
    for hubbard_id in eachindex(hubbard_parameters.orbital_ids)
        e_hub = measure_hubbard_energy(hubbard_parameters, hubbard_id, pconfig, ph_transform)
        E_hub += e_hub
    end

    return E_hub
end


@doc raw"""

    measure_coupling_energy(
        extended_hubbard_parameters::ExtendedHubbardParameters{E},
        pconfig::Vector{Int}, ph_transform::Bool
    ) where {E<:AbstractFloat}

Measures the average energy per site due to coupling in the extended Hubbard model.

"""
function measure_coupling_energy(
    extended_hubbard_parameters::ExtendedHubbardParameters{E},
    pconfig::Vector{Int}, ph_transform::Bool
) where {E<:AbstractFloat}
    E_ext_hub = zero(Complex{E})
    for ext_hub_id in eachindex(extended_hubbard_parameters.bond_ids)
        e_ext_hub = measure_ext_hubbard_energy(extended_hubbard_parameters, ext_hub_id, pconfig, ph_transform)
        E_ext_hub += e_ext_hub
    end

    return E_ext_hub
end