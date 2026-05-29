@doc raw"""

    measure_Dk(
        detwf::DeterminantalWavefunction{T},
        pconfig::Vector{Int},
        determinantal_parameters::DeterminantalParameters,
        N::Int,
        Np::Int,
    ) where {T}

Calculates the logarithmic derivative ``\Delta_k = \frac{\partial\ln\Psi_{T}}{\partial\alpha_k}`` for each determinantal parameter ``\alpha_k``.

"""
function measure_Dk(
    detwf::DeterminantalWavefunction{T},
    determinantal_parameters::DeterminantalParameters,
    pconfig::Vector{Int},
    N::Int,
    Np::Int,
) where {T}
    (; W, A) = detwf
    (; p, optimize) = determinantal_parameters

    Dk = zeros(T, sum(length, p))
    G = zeros(Complex{T}, 2*N, 2*N)

    for β in 1:Np
        k = findfirst(x -> x == β, pconfig)
        G[k, :] .= W[:, β]
    end

    A_idx   = 1
    res_idx = 1

    for (opt, params) in zip(optimize, p)
        if opt
            for _ in eachindex(params)
                Dk[res_idx] = sum(A[A_idx] .* G)
                A_idx   += 1
                res_idx += 1
            end
        else
            # skip over non-optimized parameters
            res_idx += length(params)
        end
    end

    return Dk
end


@doc raw"""

    measure_Dk(
        jastrow_parameters::JastrowParameters{T},
        o::Int,
        pconfig::Vector{Int},
        N::Int,
        ph_transform::Bool
    ) where {T}

Calculates the logarithmic derivative ``\Delta_k = \frac{\partial\ln\Psi_{T}}{\partial\alpha_k}`` for each Jastrow parameter ``\alpha_k``.

"""
function measure_Dk(
    jastrow_parameters::JastrowParameters{T},
    o::Int,
    pconfig::Vector{Int},
    N::Int,
    ph_transform::Bool
) where {T}
    (; particle_pair, order_pair, irr_indices, irr_index_map, mean_v, optimize) = jastrow_parameters

    Dk = zeros(T, length(mean_v[o]))

    max_idx = maximum(irr_indices[o])

    if optimize
        for idx in eachindex(irr_indices[o])
            irr_indices[o][idx] == max_idx && continue
            for (i,j) in irr_index_map[o][irr_indices[o][idx]]
                if particle_pair == "electron-electron"
                    ni_up, ni_dn = get_fermion_occupations(i+1, pconfig, N)
                    nj_up, nj_dn = get_fermion_occupations(j+1, pconfig, N)
                    if order_pair == "density-density"
                        pfc = 0.5
                    elseif order_pair == "spin-spin"
                        pfc = 0.25 
                    end
                    Dk[idx] += ph_transform ?
                        pfc * (ni_up - ni_dn) * (nj_up - nj_dn) :
                        pfc * (ni_up + ni_dn) * (nj_up + nj_dn)
                elseif particle_pair == "phonon-phonon"
                    # TODO
                elseif particle_pair == "electron-phonon"
                    # TODO
                end
            end
        end
    end

    return Dk
end


