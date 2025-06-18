@doc raw"""

    measure_density_correlation!( measurement_container::NameTuple,
                                  detwf::DeterminantalWavefunction,
                                  model_geometry::ModelGeometry )::Nothing

Measures the equal-time density-density correlation function.

- `measurement_container::NamedTuple`: container where measurements are stroed.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function measure_density_correlation!(
    measurement_container::NamedTuple,
    detwf::DeterminantalWavefunction,
    model_geometry::ModelGeometry
)::Nothing
    # get site-dependent density
    n = get_site_dependent_n(detwf, model_geometry)

    # calculate density-density correlation
    n_corr_current = n * n'

    # get current values from the container
    n_corr_container = measurement_container.correlation_measurements["density-density"]

    # update values for the current bin
    current_n_corr_bin = n_corr_container[2]
    current_n_corr_bin = n_corr_current

    # update accumulator for this bin
    thisbin_n_corr_sum = n_corr_container[1]
    thisbin_n_corr_sum += n_corr_current

    # combine the updated values 
    updated_values = (thisbin_n_corr_sum, current_n_corr_bin)

    # write the new values to the container
    measurement_container.correlation_measurements["density-density"] = updated_values

    return nothing
end


@doc raw"""

    measure_spin_correlation!( measurement_container::NamedTuple,
                               detwf::DeterminantalWavefunction,
                               model_geometry::ModelGeometry )::Nothing

Measures the equal-time spin-spin correlation function. 

- `measurement_container::NamedTuple`: container where measurements are stroed.
- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function measure_spin_correlation!(
    measurement_container::NamedTuple,
    detwf::DeterminantalWavefunction,
    model_geometry::ModelGeometry
)
    # get site-dependent spin
    s = get_site_dependent_s(detwf, model_geometry)

    # calculate spin-spin correlation
    s_corr_current = s * s'

    # get current values from the container
    s_corr_container = measurement_container.correlation_measurements["spin-spin"]

    # update values for the current bin
    current_s_corr_bin = s_corr_container[2]
    current_s_corr_bin = s_corr_current

    # update accumulator for this bin
    thisbin_s_corr_sum = s_corr_container[1]
    thisbin_s_corr_sum += s_corr_current

    # combine the updated values 
    updated_values = (thisbin_s_corr_sum, current_s_corr_bin)

    # write the new values to the container
    measurement_container.correlation_measurements["spin-spin"] = updated_values

    return nothing
end


@doc raw"""

    get_site_dependent_n( detwf::DeterminantalWavefunction, 
                          model_geometry::ModelGeometry )::Vector{Int}

Calculates the density on each lattice site.

- `detwf::DeterminantalWavefunction`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_site_dependent_n(
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry
)::Vector{Int}
    # number of lattice sites
    N = model_geometry.lattice.N

    # current particle configuration
    pconfig = detwf.pconfig

    # get spin-up and spin-down occupations
    up_occs = Int.(pconfig[1:N] .!= 0)
    dn_occs = Int.(pconfig[N+1:2*N] .!= 0)

    n_occs = up_occs .+ 1 .- dn_occs 

    return n_occs
end


@doc raw"""

    get_site_dependent_s( detwf::DeterminantalParameters
                          model_geometry::ModelGeometry )::Vector{Int}

Calculates the spin on each lattice site.

- `detwf::DeterminantalParameters`: current variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_site_dependent_s(
    detwf::DeterminantalWavefunction,
    model_geometry::ModelGeometry
)::Vector{Int}
    # number of lattice sites
    N = model_geometry.lattice.N

    # current particle configuration
    pconfig = detwf.pconfig

    # get spin-up and spin-down occupations
    up_occs = Int.(pconfig[1:N] .!= 0)
    dn_occs = Int.(pconfig[N+1:2*N] .!= 0)

    s_occs = up_occs .- 1 + dn_occs

    return s_occs
end


# @doc raw"""

#     calculate_structure_factor()

# Computes the static structure factor from either density or spin correlation data.

# """
# function calculate_structure_factor()

# end