@doc raw"""

    measure_hopping_energy(
        detwf::DeterminantalWavefunction{T},
        tight_binding_parameters::TightBindingParameters{T,E},
        pconfig::Vector{Int},
        hopping_id::Int,
        N::Int,
        ph_transform::Bool;
        # KEYWORD ARGUMENTS
        model_geometry::ModelGeometry = model_geometry,
        jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}},Nothing} = nothing,
        jas_parameters::Union{Tuple{JastrowParameters{T}},Nothing} = nothing
    ) where {T<:Number, E<:AbstractFloat}

Calculate the average hopping energy ``-\langle t_{\langle i,j \rangle} \hat{c}^\dagger_{\sigma,i} \hat{c}_{\sigma,j} + {\rm h.c.} \rangle``.


"""
function measure_hopping_energy(
    detwf::DeterminantalWavefunction{T},
    tight_binding_parameters::TightBindingParameters{T,E},
    pconfig::Vector{Int},
    hopping_id::Int,
    N::Int,
    ph_transform::Bool;
    # KEYWORD ARGUMENTS
    model_geometry::ModelGeometry = model_geometry,
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}},Nothing} = nothing,
    jas_parameters::Union{Tuple{JastrowParameters{T}},Nothing} = nothing
) where {T<:Number, E<:AbstractFloat}

    (; t, neighbor_table, bond_slices) = tight_binding_parameters
    (; W) = detwf

    # initialize hopping energy to zero
    ε_hop = zero(E)

    # get the neighbor table associated with the bond/hopping in question
    nt = @view neighbor_table[:, bond_slices[hopping_id]]

    # get the hopping associated with the bond/hopping in question
    t′ = @view t[bond_slices[hopping_id]]

    Np = size(W, 2)

    # sweep over particles (guarantees β is always a valid column index into W)
    for β in 1:Np
        # find the spindex of this particle
        ksite = findfirst(x -> x == β, pconfig)

        # determine spin (1 = up, 2 = down)
        spin = ksite > N ? 2 : 1

        # get the real-space site index
        k = spin == 1 ? ksite : ksite - N

        # loop over bonds involving site k
        for (col, t_hop) in zip(axes(nt, 2), t′)
            local_k, local_l = nt[1, col], nt[2, col]

            # check if this bond involves site k (in either direction)
            if local_k == k
                l = local_l
                lsite = get_spindices_from_index(l, N)[spin]

                if iszero(pconfig[lsite])
                    t_eff = ph_transform ? spin == 2 ? t_hop : -t_hop : -t_hop
                    if isnothing(jas_factors) && isnothing(jas_parameters)
                        ε_hop += t_eff * W[lsite, β]
                    else
                        o = site_to_orbital(l, model_geometry.unit_cell)
                        for (jasf, jasp) in zip(jas_factors, jas_parameters)
                            Rⱼ = calculate_fermion_jastrow_ratio(
                                jasf, jasp, model_geometry.lattice,
                                o, ksite, lsite, N, ph_transform
                            )
                            ε_hop += t_eff * Rⱼ * W[lsite, β]
                        end
                    end
                end

            elseif local_l == k
                l = local_k
                lsite = get_spindices_from_index(l, N)[spin]

                if iszero(pconfig[lsite])
                    t_eff = ph_transform ? spin == 2 ? t_hop : -t_hop : -t_hop
                    if isnothing(jas_factors) && isnothing(jas_parameters)
                        ε_hop += t_eff * W[lsite, β]
                    else
                        o = site_to_orbital(l, model_geometry.unit_cell)
                        for (jasf, jasp) in zip(jas_factors, jas_parameters)
                            Rⱼ = calculate_fermion_jastrow_ratio(
                                jasf, jasp, model_geometry.lattice,
                                o, ksite, lsite, N, ph_transform
                            )
                            ε_hop += t_eff * Rⱼ * W[lsite, β]
                        end
                    end
                end
            end
        end
    end

    # normalize the measurement
    ε_hop /= length(t′)

    return ε_hop
end
