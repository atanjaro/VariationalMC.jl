@doc raw"""

    write_measurements!(;
        # KEYWORD ARGUMENTS
        measurement_container::NamedTuple,
        simulation_info::SimulationInfo,
        model_geometry::ModelGeometry,
        bin::Int,
        bin_size::Int,
        opt_step::Bool = false
    )

Write the measurements contained in `measurement_container` to file.
Measurements are written to file in HDF5 format.

"""
function write_measurements!(;
    # KEYWORD ARGUMENTS
    measurement_container::NamedTuple,
    simulation_info::SimulationInfo,
    model_geometry::ModelGeometry,
    bin::Int,
    bin_size::Int,
    opt_step::Bool = false
)
    (; datafolder, pID, write_bins_concurrent, opt_bin_files, sim_bin_files) = simulation_info
    lattice   = model_geometry.lattice
    unit_cell = model_geometry.unit_cell
    bonds     = model_geometry.bonds

    (; global_measurements, local_measurements, optimization_measurements, 
        equaltime_correlations, n_params, n_opt_params, pfft!) = measurement_container

    bin_files = sim_bin_files
    if opt_step
        bin_files = opt_bin_files
    end

    # normalize all measurements by the bin size (if this is not during optimization)
    if !opt_step
        normalize_measurements!(measurement_container, bin_size)
    end

    # construct filename
    filename = opt_step ? joinpath(datafolder, "opt_bins", @sprintf("pID-%d", pID), @sprintf("bin-%d.h5", bin)) : 
        joinpath(datafolder, "sim_bins", @sprintf("pID-%d", pID), @sprintf("bin-%d.h5", bin))


     # open HDF5 file to write binned data to
    if write_bins_concurrent
        file = h5open(filename, "w")
    else
        # Use in-memory file
        file = h5open(filename, "w"; driver=Drivers.Core(; backing_store=false))
    end

    # if first bin record system info
    if isone(bin)
        # record total number of variational parameters
        attributes(file)["NUM_PARAMS"] = n_params
        # records total number of optimized parameters
        attributes(file)["NUM_OPT_PARAMS"] = n_opt_params
        # record total number of orbitals in lattice
        attributes(file)["N_ORBITALS"] = nsites(unit_cell, lattice)
    end

    # write global measurements to group
    Global = create_group(file, "GLOBAL")
    for (measurement, value) in global_measurements
        Global[measurement] = value
    end

    # write local measurements to group
    Local = create_group(file, "LOCAL")
    for (measurement, value) in local_measurements
        if value isa Vector{<:Vector}
            grp = create_group(Local, measurement)
            for (i, v) in enumerate(value)
                grp[string(i)] = v
            end
        else
            Local[measurement] = value
        end
    end
    
    # write optimization measurements to group
    if opt_step
        Optimization = create_group(file, "OPTIMIZATION")
        Optimization["parameters"] = optimization_measurements["parameters"]
    end

    # create group to contain correlation measurements
    Correlations = create_group(file, "CORRELATIONS")

    # create group for standard equal-time correlation measurements
    EqualTime = create_group(Correlations, "EQUAL-TIME")

     # iterate over equal-time correlation measurements
    for correlation in keys(equaltime_correlations)
        # get the correlation container for current standard equal-time correlation measurement
        correlation_container = equaltime_correlations[correlation]
        id_pairs = correlation_container.id_pairs::Vector{NTuple{2,Int}}
        id_type = CORRELATION_FUNCTIONS[correlation]
        correlations = correlation_container.correlations::Vector{Matrix{Complex{T}}} where T<:AbstractFloat

         # create a group for correlation measurement
        EqualTimeCorrelation = create_group(EqualTime, correlation)

        # record ID pairs that were measured
        attributes(EqualTimeCorrelation)["ID_PAIRS"] = id_pairs

        # record the ID type corresponding to correlation measurement
        attributes(EqualTimeCorrelation)["ID_TYPE"] = id_type

        # record the position space correlations
        EqualTimeCorrelation["POSITION"] = copy(stack(correlations))

        # Fourier transform correlations to momentum space
        for i in eachindex(correlations)
            @assert (id_type == "ORBITAL_ID") || (id_type == "BOND_ID")         # this works for now

            bond_b_id, bond_a_id = id_pairs[i]
            a = bonds[bond_a_id].orbitals[1]
            b = bonds[bond_b_id].orbitals[1]

            # perform Fourier transform
            fourier_transform!(correlations[i], a, b, unit_cell, lattice, pfft!)
        end

        # record the momentum space correlations
        EqualTimeCorrelation["MOMENTUM"] = copy(stack(correlations))
    end

    # If writing binned data to file.
    if write_bins_concurrent
        close(file)
        file_bytes = read(filename)  
        push!(bin_files, file_bytes)
        # If using in-memory files.
    else
        flush(file)
        file_bytes = Vector{UInt8}(file)  
        close(file)
        push!(bin_files, file_bytes)
    end

    # ZERO THE MEASUREMENT CONTAINER

    # reset global measurements to zero
    reset_global_measurements!(global_measurements)

    # reset local measurements to zero
    reset_local_measurements!(local_measurements)

    # reset optimization measurements to zero
    reset_optimization_measurements!(optimization_measurements)

    # reset equal-time correlation measurements to zero
    for correlation in keys(equaltime_correlations)
        correlation_container = equaltime_correlations[correlation]
        reset!(correlation_container)
    end

    return nothing
end


# normalize measurements that will be written to file by the bin size
function normalize_measurements!(measurement_container::NamedTuple, bin_size::Int)
    # normalize global measurements by bin size
    global_measurements = measurement_container.global_measurements
    normalize_global_measurements!(global_measurements, bin_size)

    # normalize local measurements by bin size
    local_measurements = measurement_container.local_measurements
    normalize_local_measurements!(local_measurements, bin_size)

    # normalize optimization measurments by bin size
    optimization_measurements = measurement_container.optimization_measurements
    normalize_optimization_measurements!(optimization_measurements, bin_size)

    # normalize equal-time correlation function measurement
    equaltime_correlations = measurement_container.equaltime_correlations
    normalize_correlation_measurements!(equaltime_correlations, bin_size)

    return nothing
end


# normalize global measurements by bin size
function normalize_global_measurements!(
    global_measurements::Dict{String, Complex{T}}, bin_size::Int
) where {T<:AbstractFloat}

    for global_measurement in keys(global_measurements)
        global_measurements[global_measurement] /= bin_size
    end

    return nothing
end


# normalize local measurements by bin size
function normalize_local_measurements!(
    local_measurements::Dict{String, Union{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}}, bin_size::Int
) where {T<:AbstractFloat}

    for local_measurement in keys(local_measurements)
        @. local_measurements[local_measurement] /= bin_size
    end

    return nothing
end


# normalize optimization measurments by bin size
function normalize_optimization_measurements!(
    optimization_measurements::Dict{String, Union{Vector{T}, Vector{Complex{T}}, Matrix{Complex{T}}}}, bin_size::Int
) where {T<:AbstractFloat}

    for optimization_measurement in keys(optimization_measurements)
        @. optimization_measurements[optimization_measurement] /= bin_size
    end

    return nothing
end


# normalize correlation measurement
function normalize_correlation_measurements!(
    correlation_measurements::Dict{String, CorrelationContainer{D, T}}, bin_size::Int
) where {D, T<:AbstractFloat}

    for measurement in keys(correlation_measurements)
        correlation_container = correlation_measurements[measurement]
        pairs = correlation_container.id_pairs::Vector{NTuple{2,Int}}
        correlations = correlation_container.correlations::Vector{Array{Complex{T}, D}}
        for i in eachindex(pairs)
            @. correlations[i] /= bin_size
        end
    end

    return nothing
end


# zero the global measurements
function reset_global_measurements!(global_measurements)
    for key in keys(global_measurements)
        global_measurements[key] = zero(typeof(global_measurements[key]))
    end
    return nothing
end


# zero the local measurements
function reset_local_measurements!(local_measurements)
    for value in values(local_measurements)
        _reset_container!(value)
    end
    return nothing
end


# zero the optimization_measurements
function reset_optimization_measurements!(optimization_measurements)
    for value in values(optimization_measurements)
        _reset_container!(value)
    end
    return nothing
end

# for numeric arrays
function _reset_container!(A::AbstractArray{<:Number})
    fill!(A, zero(eltype(A)))
    return nothing
end

# for arrays of arrays
function _reset_container!(A::AbstractArray)
    for x in A
        _reset_container!(x)
    end
    return nothing
end

