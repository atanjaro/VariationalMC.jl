@doc raw"""

    initialize_measurement_container(
        determinantal_parameters::DeterminantalParameters,
        model_geometry::ModelGeometry{D,T,N},
        N_opt::Int,
        N_sim::Int,
        opt_bin_size::Int,
        sim_bin_size::Int;
        jas_parameters::Union{Tuple{JastrowParameters},Nothing}=nothing
    ) where {T<:AbstractFloat, D, N}

Initialize and return a measurement container of type `NamedTuple`.

"""
function initialize_measurement_container(
    determinantal_parameters::DeterminantalParameters,
    model_geometry::ModelGeometry{D,T,N},
    jas_parameters::Union{Tuple{JastrowParameters},Nothing} = nothing
) where {T<:AbstractFloat, D, N}
    lattice   = model_geometry.lattice::Lattice{D}
    unit_cell = model_geometry.unit_cell::UnitCell{D,T,N}
    # bonds     = model_geometry.bonds::Vector{Bond{D}}

    # number of orbitals per unit cell
    norbitals = unit_cell.n

    # size of lattice in unit cells in direction of each lattice vector
    L = lattice.L

    # initialize global measurements
    global_measurements = Dict{String, Complex{T}}(k => zero(Complex{T}) for k in GLOBAL_MEASUREMENTS)

    # initialize local measurements
    local_measurements = Dict{String, Union{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}}(
        "density"                   => zeros(Complex{T}, norbitals),    # average density for each orbital species
        "double_occ"                => zeros(Complex{T}, norbitals),    # average double occupancy for each orbital
        "spin-z"                    => zeros(Complex{T}, norbitals),    # average z-component of spin for each orbital
    )

    # get the total number of variational parameters
    n_params = sum(length, determinantal_parameters.p)

    # count optimized determinantal parameters
    n_opt_params = sum(length(p) for (opt, p) in zip(determinantal_parameters.optimize, determinantal_parameters.p) if opt)

    if !isnothing(jas_parameters)
        for jasp in jas_parameters
            for o in jasp.orbitals
                n_params += length(jasp.mean_v[o]) 
                n_opt_params += length(jasp.mean_v[o]) - 1
            end
        end
    end

    # initialize optimization measurements
    optimization_measurements = Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}}(
        "logDk"         => zeros(Complex{T}, n_params),
        "logDklogDkp"   => zeros(Complex{T}, n_params, n_params),
        "logDkE"        => zeros(Complex{T}, n_params),
        "parameters"    => zeros(T, n_params)
    )

    # initialize equal-time correlation measurement dictionary
    equaltime_correlations = Dict{String, CorrelationContainer{D,T}}()

    # initialize fft plan
    pfft! = plan_fft!(zeros(Complex{T}, prod(L), prod(L)); flags=FFTW.PATIENT)

    # initialize measurement container
    measurement_container = (
        global_measurements = global_measurements,
        local_measurements = local_measurements,
        optimization_measurements = optimization_measurements,
        equaltime_correlations = equaltime_correlations,
        hopping_to_bond_id = Int[],
        L = L,
        n_params = n_params,
        n_opt_params = n_opt_params,
        pfft! = pfft!
    )

    return measurement_container
end


@doc raw"""

    initialize_measurements!(
        measurement_container::NamedTuple,
        tight_binding_model::TightBindingModel{T,E}
    ) where {T<:Number, E<:AbstractFloat}

Initialize tight-binding model related measurements.

# Initialized Measurements

- `hopping_energy`: Refer to [`measure_hopping_energy`](@ref).

"""

function initialize_measurements!(
    measurement_container::NamedTuple,
    tight_binding_model::TightBindingModel{T,E}
) where {T<:Number, E<:AbstractFloat}

    (; local_measurements, hopping_to_bond_id) = measurement_container

    # number of types of hoppings
    nhopping = length(tight_binding_model.t_bond_ids)

    # initialize hopping related measurements
    if nhopping > 0
        local_measurements["hopping_energy"] = zeros(Complex{E}, nhopping)
    end

    # record bond ID associated with each hopping ID
    for id in tight_binding_model.t_bond_ids
        push!(hopping_to_bond_id, id)
    end

    return nothing
end


@doc raw"""

    initialize_measurements!(
        measurement_container::NamedTuple,
        hubbard_model::HubbardModel{T}
    ) where {T<:AbstractFloat}

Initialize Hubbard model related measurements.

# Initialized Measurements

- `hubbard_energy`: Refer to [`measure_hopping_energy`](@ref).

"""
function initialize_measurements!(
    measurement_container::NamedTuple,
    hubbard_model::HubbardModel{T}
) where {T<:AbstractFloat}

    (; local_measurements) = measurement_container
    (; U_orbital_ids) = hubbard_model

    # number of orbitals in unit cell
    n_hubbard = length(U_orbital_ids)

    # initialize hubbard energy measurement U⋅nup⋅ndn
    local_measurements["hubbard_energy"] = zeros(Complex{T}, n_hubbard)

    return nothing
end


@doc raw"""

    initialize_measurements!(
        measurement_container::NamedTuple,
        extended_hubbard_model::ExtendedHubbardModel{T}
    ) where {T<:AbstractFloat}

Initialize Extended Hubbard model related measurements.

# Initialized Measurements

- `ext_hubbard_energy`: Refer to [`measure_ext_hubbard_energy`](@ref).

"""
function initialize_measurements!(
    measurement_container::NamedTuple,
    extended_hubbard_model::ExtendedHubbardModel{T}
) where {T<:AbstractFloat}

    (; local_measurements) = measurement_container
    (; V_bond_ids) = extended_hubbard_model

    # number of extended hubbard interactions
    n_ehi = length(V_bond_ids)

    # initialize extended hubbard energy measurement
    local_measurements["ext_hub_energy"] = zeros(Complex{T}, n_ehi)

    return nothing
end


@doc raw"""

    initialize_site_dependent_measurements!(;
        measurement_container::NamedTuple,
        model_geometry::ModelGeometery,
        observable::String
    )

Initialize site-dependent observable measurements.

# Initialized Measurements

- `site-dependent density`
- `site-dependent spin-z`

"""
function initialize_site_dependent_measurements!(;
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry{D,T,N},
    observable::String
) where {T<:AbstractFloat, D, N}
    (; local_measurements, L) = measurement_container
    (; unit_cell) = model_geometry

    @assert observable == "density" || observable == "spin-z"

    # number of orbitals per unit cell
    norbitals = unit_cell.n

    # initialize site-dependent measurement
    local_measurements["site-dependent_" * observable] = [zeros(Complex{T}, prod(L)) for _ in 1:norbitals]

    return nothing
end


@doc raw"""

    initialize_correlation_measurements!(;
        measurement_container::NamedTuple,
        correlation::String,
        pairs::AbstractVector{NTuple{2,Int}},
        time_displaced::Bool,
        integrated::Bool = false
    )  where {T<:AbstractFloat, D, N}

Initialize measurements of `correlation` for all ID pairs; refer to [`CORRELATION_FUNCTIONS`](@ref) for ID type associated
with each correlation measurement.
The name `correlation` must therefore also appear in [`CORRELATION_FUNCTIONS`]@ref.

"""
function initialize_correlation_measurements!(;
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry{D,T,N},
    correlation::String,
    pairs::AbstractVector{NTuple{2,Int}},
    integrated::Bool = false
)  where {T<:AbstractFloat, D, N}

    # iterate over all bond/orbial ID pairs
    for pair in pairs
        initialize_correlation_measurement!(
            measurement_container = measurement_container,
            model_geometry = model_geometry,
            correlation = correlation,
            pair = pair,
            integrated = integrated
        )
    end

    return nothing
end

# initialize a single correlation measurement
function initialize_correlation_measurement!(;
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry{D,T,N},
    correlation::String,
    pair::NTuple{2,Int},
    integrated::Bool = false
)  where {T<:AbstractFloat, D, N}

    (; equaltime_correlations) = measurement_container

    # check to make sure valid correlation measurement
    @assert correlation in keys(CORRELATION_FUNCTIONS)

    # extent of lattice in unit cells
    L = measurement_container.L

    # add equal-time correlation key, if not present
    if !haskey(equaltime_correlations, correlation)
        equaltime_correlations[correlation] = CorrelationContainer(D, T)
    end

    # add equal-time correlation measurement
    push!(equaltime_correlations[correlation].id_pairs, pair)
    push!(equaltime_correlations[correlation].correlations, zeros(Complex{T}, prod(L), prod(L)))

    return nothing
end
