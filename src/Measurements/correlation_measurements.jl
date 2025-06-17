# @doc raw"""

#     measure_density_correlation!( measurement_container::NameTuple,
#                                   detwf::DeterminantalWavefunction,
#                                   model_geometry::ModelGeometry )::Nothing

# Measures the density-density correlation function.

# """
# function measure_density_correlation!(
#     measurement_container::NamedTuple,
#     detwf::DeterminantalWavefunction,
#     model_geometry::ModelGeometry
# )::Nothing


# end


# @doc raw"""

#     measure_spin_correlation!( measurement_container::NamedTuple,
#                                detwf::DeterminantalWavefunction,
#                                model_geometry::ModelGeometry )::Nothing

# Measures the spin-spin correlation function. 

# """
# function measure_spin_correlation!(
#     measurement_container::NamedTuple,
#     detwf::DeterminantalWavefunction,
#     model_geometry::ModelGeometry
# )

# end


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