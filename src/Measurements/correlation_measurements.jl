@doc raw"""

    measure_density_correlation!( measurement_container::NameTuple,
                                  detwf::DeterminantalWavefunction{T, Q, E, I},
                                  model_geometry::ModelGeometry,
                                  pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the equal-time density-density correlation function ``\langle \hat{n}_{i}\hat{n}_{j}\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stroed.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transform.

"""
function measure_density_correlation!(
    measurement_container::NamedTuple,
    detwf::DeterminantalWavefunction{T, Q, E, I},
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    # get site-dependent density
    n = get_site_dependent_n(detwf, model_geometry, pht)

    # calculate density-density correlation
    n_corr_current = n * n'

    # add current correlation to the accumulator
    measurement_container.correlation_measurements["density"] += n_corr_current

    return nothing
end


@doc raw"""

    measure_spin_correlation!( measurement_container::NamedTuple,
                               detwf::DeterminantalWavefunction{T, Q, E, I},
                               model_geometry::ModelGeometry,
                               pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Measures the equal-time spin-spin correlation function ``\langle \hat{S}_{i}\hat{S}_{j}\rangle``. 

- `measurement_container::NamedTuple`: container where measurements are stroed.
- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transform.

"""
function measure_spin_correlation!(
    measurement_container::NamedTuple,
    detwf::DeterminantalWavefunction{T, Q, E, I},
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    # get site-dependent spin
    s = get_site_dependent_s(detwf, model_geometry, pht)

    # calculate spin-spin correlation
    s_corr_current = s * s'

    # add current correlation to the accumulator
    measurement_container.correlation_measurements["spin"] += s_corr_current

    return nothing
end


@doc raw"""

    get_site_dependent_n( detwf::DeterminantalWavefunction{T, Q, E, I}, 
                          model_geometry::ModelGeometry,
                          pht::Bool ) where {T<:Number, Q, E<:AbstractFloat, I<:Integer}

Calculates the local density ``n_{\boldsymbol{i}}`` on each lattice site.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_site_dependent_n(
    detwf::DeterminantalWavefunction{T, Q, E, I}, 
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    # number of lattice sites
    N = model_geometry.lattice.N

    # current particle configuration
    pconfig = detwf.pconfig

    # get spin-up and spin-down occupations
    up_occs = Int.(pconfig[1:N] .!= 0)
    dn_occs = Int.(pconfig[N+1:2*N] .!= 0)

    if pht
        # n_occs = up_occs .- dn_occs .+ 1
        n_occs = up_occs .+ dn_occs .- 1
    else
        n_occs = up_occs .+ dn_occs 
    end

    return n_occs
end


@doc raw"""

    get_site_dependent_s( detwf::DeterminantalParameters{T, Q, E, I}
                          model_geometry::ModelGeometry,
                          pht::Bool )

Calculates the local spin ``S_{\boldsymbol{i}}`` on each lattice site.

- `detwf::DeterminantalParameters{T, Q, E, I}`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_site_dependent_s(
    detwf::DeterminantalWavefunction{T, Q, E, I},
    model_geometry::ModelGeometry,
    pht::Bool
) where {T<:Number, Q, E<:Number, I<:Integer}
    # number of lattice sites
    N = model_geometry.lattice.N

    # current particle configuration
    pconfig = detwf.pconfig

    # get spin-up and spin-down occupations
    up_occs = Int.(pconfig[1:N] .!= 0)
    dn_occs = Int.(pconfig[N+1:2*N] .!= 0)

    if pht
        # s_occs = up_occs .+ dn_occs .- 1
        s_occs = up_occs .- dn_occs .+ 1
    else
        s_occs = up_occs .- dn_occs
    end

    return s_occs
end