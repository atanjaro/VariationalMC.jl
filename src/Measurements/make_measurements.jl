@doc raw"""

    make_measurements!(
        measurement_container::NamedTuple;
        detwf::DeterminantalWavefunction{T},
        jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
        tight_binding_parameters::TightBindingParameters{T, E},
        jas_parameters::Union{Tuple{JastrowParameters{E}}, Nothing} = nothing,
        coupling_parameters::Tuple,
        particle_configuration::ParticleConfiguration,
        model_geometry::ModelGeometry;
        opt_step::Bool = false
    ) where {T<:Number, E<:AbstractFloat}

Make measurements. If used during the optimization phase of the simulation, set `opt_step = true`, which will activate the 
measurements of the logarithmic derivatives.

"""
function make_measurements!(
    measurement_container::NamedTuple;
    detwf::DeterminantalWavefunction{T},
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
    tight_binding_parameters::TightBindingParameters{T, E},
    determinantal_parameters::DeterminantalParameters,
    jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing,
    coupling_parameters::Tuple,
    particle_configuration::ParticleConfiguration,
    model_geometry::ModelGeometry,
    opt_step::Bool = false
) where {T<:Number, E<:AbstractFloat}
    (; equaltime_correlations) = measurement_container
    (; unit_cell, lattice) = model_geometry
    (; pconfig, Np, ph_transform) = particle_configuration

    # total number of lattice sites
    Norbs = unit_cell.n 
    Ncells = lattice.N
    N = Norbs * Ncells 

    # make local measurements
    local_measurements = measurement_container.local_measurements
    make_local_measurements!(
        local_measurements,
        detwf,
        tight_binding_parameters,
        coupling_parameters,
        model_geometry,
        pconfig, N, ph_transform,
        jas_factors = jas_factors,
        jas_parameters = jas_parameters
    )

    # make global_measurements
    global_measurements = measurement_container.global_measurements
    make_global_measurements!(
        global_measurements,
        detwf,
        tight_binding_parameters,
        coupling_parameters,
        pconfig, N, ph_transform;
        model_geometry = model_geometry,
        jas_factors = jas_factors,
        jas_parameters = jas_parameters
    )

    if opt_step
        # make optimization measurements
        optimization_measurements = measurement_container.optimization_measurements
        make_optimization_measurements!(
            optimization_measurements,
            detwf,
            tight_binding_parameters,
            determinantal_parameters,
            coupling_parameters,
            pconfig, N, Np,
            model_geometry = model_geometry,
            jas_factors = jas_factors,
            jas_parameters = jas_parameters,
            ph_transform = ph_transform
        )
    end

    # make equal-time correlation measurements
    if !isempty(equaltime_correlations)
        make_equaltime_measurements!(
            equaltime_correlations,
            unit_cell,
            pconfig, N, Ncells, ph_transform
        )
    end

    return nothing
end


# make local measurements
function make_local_measurements!(
    local_measurements::Dict{String, Union{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}},
    detwf::DeterminantalWavefunction{T},
    tight_binding_parameters::TightBindingParameters{T, E},
    coupling_parameters::Tuple,
    model_geometry::ModelGeometry,
    pconfig::Vector{Int}, N::Int, ph_transform::Bool;
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
    jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing
) where {T<:Number, E<:AbstractFloat}
    # number of orbitals per unit cell
    unit_cell = model_geometry.unit_cell
    lattice = model_geometry.lattice
    norbital = unit_cell.n
    Ncells = lattice.N

    for n in 1:norbital
        # measure density
        local_measurements["density"][n] += measure_n("mean", pconfig, unit_cell, n, Ncells, N, ph_transform)
        if haskey(local_measurements, "site-dependent_density")
            local_measurements["site-dependent_density"][n] += measure_n("site-dependent", pconfig, unit_cell, n, Ncells, N, ph_transform)
        end
        # measure double occupancy
        local_measurements["double_occ"][n] += measure_double_occ(pconfig, unit_cell, n, Ncells, N, ph_transform)
        # measure spin-z
        local_measurements["spin-z"][n] += measure_spinz("mean", pconfig, unit_cell, n, Ncells, N, ph_transform)
        if haskey(local_measurements, "site-dependent_spin-z")
            local_measurements["site-dependent_spin-z"][n] += measure_spinz("site-dependent", pconfig, unit_cell, n, Ncells, N, ph_transform)
        end
    end

    # make tight binding measurements
    make_local_measurements!(
        local_measurements,
        detwf,
        tight_binding_parameters,
        pconfig, N, ph_transform,
        model_geometry = model_geometry,
        jas_factors = jas_factors,
        jas_parameters = jas_parameters
    )

    # make local measurements associated with couplings
    for coupling_parameter in coupling_parameters
        make_local_measurements!(
            local_measurements,
            coupling_parameter,
            pconfig,
            ph_transform
        )
    end

    return nothing
end


# make local measurements associated with the tight binding model
function make_local_measurements!(
    local_measurements::Dict{String, Union{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}},
    detwf::DeterminantalWavefunction{T},
    tight_binding_parameters::TightBindingParameters{T, E},
    pconfig::Vector{Int}, N::Int, ph_transform::Bool;
    model_geometry::ModelGeometry = model_geometry,
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
    jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing
) where {T<:Number, E<:AbstractFloat}
    # number of orbitals per unit cell
    norbital = tight_binding_parameters.norbital

    # number of types of hopping
    bond_ids = tight_binding_parameters.bond_ids
    nhopping = length(tight_binding_parameters.bond_ids)

    if nhopping > 0
        for hopping_id in 1:nhopping
            # measure the hopping hopping
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

            local_measurements["hopping_energy"][hopping_id] += e_hop
        end
    end

    return nothing
end


# make local measurements associated with the Hubbard model
function make_local_measurements!(
    local_measurements::Dict{String, Union{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}},
    hubbard_parameters::HubbardParameters{E},
    pconfig::Vector{Int}, ph_transform::Bool
) where {T<:Number, E<:AbstractFloat}
    # measure hubbard energy for each orbital in unit cell
    hubbard_energies = local_measurements["hubbard_energy"]
    for hubbard_id in eachindex(hubbard_energies)
        hubbard_energies[hubbard_id] += measure_hubbard_energy(hubbard_parameters, hubbard_id, pconfig, ph_transform)
    end

    return nothing
end


# make local measurements associated with the extended Hubbard model
function make_local_measurements!(
    local_measurements::Dict{String, Union{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}},
    extended_hubbard_parameters::ExtendedHubbardParameters{E},
    pconfig::Vector{Int}, ph_transform::Bool
) where {T<:Number, E<:AbstractFloat}

    # measure hubbard energy for each orbital in unit cell
    ext_hub_energies = local_measurements["ext_hub_energy"]
    for ext_hub_id in eachindex(ext_hub_energies)
        ext_hub_energies[ext_hub_id] += measure_ext_hubbard_energy(extended_hubbard_parameters, ext_hub_id, pconfig, ph_transform)
    end

    return nothing
end


# make global measurements
function make_global_measurements!(
    global_measurements::Dict{String, Complex{T}},
    detwf::DeterminantalWavefunction{T},
    tight_binding_parameters::TightBindingParameters{T, E},
    coupling_parameters::Tuple,
    pconfig::Vector{Int}, N::Int, ph_transform::Bool;
    model_geometry::ModelGeometry,
    jas_factors::Union{Tuple{<:AbstractJastrowFactor{T}}, Nothing} = nothing,
    jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing
) where {T<:Number, E<:AbstractFloat}
    # measure the average density
    global_measurements["density"] += measure_n("mean", pconfig, N, ph_transform)

    # measure the double occupancy
    global_measurements["double_occ"] += measure_double_occ(pconfig, N, ph_transform)

    # measure the average spin-z
    global_measurements["spin-z"] += measure_spinz("mean", pconfig, N, ph_transform)
    
    # measure the total energy per site 
    global_measurements["energy_per_site"] += measure_variational_energy(
                                                    detwf, 
                                                    tight_binding_parameters, 
                                                    coupling_parameters,
                                                    pconfig, N, ph_transform,
                                                    model_geometry = model_geometry, 
                                                    jas_factors = jas_factors,
                                                    jas_parameters = jas_parameters
                                                )

    return nothing
end


# make optimization measurements
function make_optimization_measurements!(
    optimization_measurements::Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}},
    detwf::DeterminantalWavefunction{T},
    tight_binding_parameters::TightBindingParameters{T, E},
    determinantal_parameters::DeterminantalParameters{S, T},
    coupling_parameters::Tuple,
    pconfig::Vector{Int}, N::Int, Np::Int;
    model_geometry::ModelGeometry,
    jas_factors::Union{Tuple{<:AbstractJastrowFactor}, Nothing} = nothing,
    jas_parameters::Union{Tuple{JastrowParameters}, Nothing} = nothing,
    ph_transform::Bool = false
) where {T<:Number, E<:AbstractFloat, S}
    logD = measure_Dk(detwf, determinantal_parameters, pconfig, N, Np)

    if !isnothing(jas_parameters)
        for jasp in jas_parameters
            for o in jasp.orbitals
                logDj = measure_Dk(jasp, o, pconfig, N, ph_transform)
                logD = vcat(logD, logDj)
            end
        end
    end

    # measure the global energy
    E_glob = measure_variational_energy(
        detwf, 
        tight_binding_parameters, 
        coupling_parameters,
        pconfig, N, ph_transform,
        model_geometry = model_geometry, 
        jas_factors = jas_factors,
        jas_parameters = jas_parameters
    )

    optimization_measurements["logDk"]          += logD
    optimization_measurements["logDklogDkp"]    += logD .* logD'
    optimization_measurements["logDkE"]         += logD .* E_glob

    return nothing
end


# make equal-time measurements
function make_equaltime_measurements!(
    equaltime_correlations::Dict{String, CorrelationContainer{D,E}},
    unit_cell::UnitCell,
    pconfig::Vector{Int}, N::Int, Ncells::Int, ph_transform::Bool
) where {D, E}
    # iterate over equal-time correlation function getting measured
    for correlation in keys(equaltime_correlations)

        correlation_container = equaltime_correlations[correlation]::CorrelationContainer{D,E}
        id_pairs = correlation_container.id_pairs::Vector{NTuple{2,Int}}
        correlations = correlation_container.correlations::Vector{Array{Complex{E}, D}}

        if correlation == "density"
            for k in eachindex(id_pairs)
                (i,j) = id_pairs[k]
                n_i = measure_n("site-dependent", pconfig, unit_cell, i, Ncells, N, ph_transform) 
                n_j = measure_n("site-dependent", pconfig, unit_cell, j, Ncells, N, ph_transform) 
                
                correlations[k] += n_i * n_j'
            end
        elseif correlation == "spin-z"
            for k in eachindex(id_pairs)
                (i,j) = id_pairs[k]
                s_i = measure_spinz("site-dependent", pconfig, unit_cell, i, Ncells, N, ph_transform) 
                s_j = measure_spinz("site-dependent", pconfig, unit_cell, j, Ncells, N, ph_transform) 
  
                correlations[k] += s_i * s_j'
            end
        elseif correlation == "pair"
            @assert ph_transform == true
            # TODO
        end
    end

    return nothing
end