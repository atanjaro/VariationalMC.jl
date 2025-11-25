@doc raw"""

    write_measurements!( step::AbstractString
                         bin::I, 
                         bin_size::I,
                         measurement_container::NamedTuple, 
                         simulation_info::SimulationInfo;
                         write_parameters=false ) where {I<:Integer}

Writes simulation and correlation (if measured) measurements in the current bin to file. 
Files created are in HDF5 format. 

- `step::AbstractString`: VMC step: `opt` or `sim`, `optimization` or `simulation`.
- `bin::I`: current bin.
- `bin_size::I`: size of the current bins.
- `measurement_container::NamedTuple`: container where measurements are stored.
- `simulation_info::SimulationInfo`: contains datafolder names.
- `write_parameters=false`: whether parameters are written to file. Set to `true` during optimization.

Data in each file can be accessed by doing:

h5open(sim_file, "r") do f
    # list all top-level measurements
    println(keys(f))  # e.g., ["global_density", "double_occ", ...]

    # access a specific dataset
    data = read(f["/local_energy/bin-1/"])  # returns a Julia array
    println(data)
end

for num in 1:N_opt_bins
    h5open(opt_file, "r") do f
        data = read(f["/parameters/bin-$(num)/"]) 
        println(data)
    end
end

""" 
function write_measurements!(
    step::AbstractString,
    bin::I,
    bin_size::I,
    measurement_container::NamedTuple,
    simulation_info::SimulationInfo;
    write_parameters=false
) where {I<:Integer}
    (; datafolder, pID) = simulation_info
    (; simulation_measurements, optimization_measurements, correlation_measurements) = measurement_container

    # ensure directories exist
    mkpath(datafolder)

    sim_dir  = joinpath(datafolder, "simulation")
    mkpath(sim_dir)

    opt_dir = joinpath(datafolder, "optimization")
    mkpath(opt_dir)

    corr_dir = joinpath(datafolder, "correlation")
    mkpath(corr_dir)

    step_filename = step * "_bin_measurements.h5"
    sim_file  = joinpath(sim_dir,  step_filename)
    opt_file = joinpath(opt_dir, step_filename)
    corr_file = joinpath(corr_dir, step_filename)

    bin_group = @sprintf("bin-%d", bin)

    # write simulation data
    sim_mode = isfile(sim_file) ? "r+" : "w"
    h5open(sim_file, sim_mode) do f
        if !haskey(attrs(f), "created_by")
            attrs(f)["created_by"] = "VariationalMC"
            attrs(f)["proc_id"]    = pID
        end

        expected = ["global_density", "double_occ", "local_energy", "pconfig"]
        for key in expected
            if haskey(simulation_measurements, key)
                data = simulation_measurements[key]
                data_to_write = normalize_measurements(data, bin_size, key)
                dset_path = "/" * key * "/" * bin_group

                if haskey(f, dset_path)
                    delete!(f, dset_path)
                end
                f[dset_path] = data_to_write
            end
        end
    end

    # write variational parameters to file (if optimizing)
    if write_parameters
        opt_mode = isfile(opt_file) ? "r+" : "w"
        h5open(opt_file, opt_mode) do f
            if !haskey(attrs(f), "created_by")
                attrs(f)["created_by"] = "VariationalMC"
                attrs(f)["proc_id"]    = pID
            end

            data_to_write = optimization_measurements["parameters"]
            dset_path = "/" * "parameters" * "/" * bin_group
            if haskey(f, dset_path)
                delete!(f, dset_path)
            end
            f[dset_path] = data_to_write
        end
    end

    # write correlation data (if measuring)
    if !isempty(correlation_measurements)
        corr_mode = isfile(corr_file) ? "r+" : "w"
        h5open(corr_file, corr_mode) do f
            if !haskey(attrs(f), "created_by")
                attrs(f)["created_by"] = "VariationalMC"
                attrs(f)["proc_id"]    = pID
            end

            for key in keys(correlation_measurements)
                data = correlation_measurements[key]
                if !isempty(data)
                    data_to_write = normalize_measurements(data, bin_size, string(key))
                    dset_path = "/" * string(key) * "/" * bin_group
                    if haskey(f, dset_path)
                        delete!(f, dset_path)
                    end
                    f[dset_path] = data_to_write
                end
            end
        end
    end

    @debug """
    Measurements::write_measurements!() :
    Simulation data written to HDF5: $(abspath(sim_file))
    """

    if isfile(corr_file)
        @debug """
        Measurements::write_measurements!() :
        Correlation data written to HDF5: $(abspath(corr_file))
        """
    else
        @debug """
        Measurements::write_measurements!() :
        No correlation file written (no correlation measurements).
        """
    end

    # reset measurement containers
    reset_measurements!(simulation_measurements)
    if write_parameters
        reset_measurements!(optimization_measurements)
    end
    if !isempty(correlation_measurements)
        reset_measurements!(correlation_measurements)
    end

    return nothing
end


@doc raw"""

    normalize_measurements( data, 
                            bin_size::I, 
                            key::AbstractString )

Normalizes different types of data by the length of a bin.

- `data`: data of abstract type.
- `bin_size::I`: length of the bin.
- `key::AbstractString`: name of measurement.

"""
function normalize_measurements(
    data, 
    bin_size::I, 
    key::AbstractString
) where {I<:Integer}
    if key == "pconfig"
        return data
    end

    if isa(data, Number)
        return data / bin_size
    elseif isa(data, AbstractArray)
        return data ./ bin_size
    elseif isa(data, Dict)
        return Dict(k => normalize_measurements(v, bin_size, string(k)) for (k, v) in data)
    elseif isa(data, NamedTuple)
        return NamedTuple{keys(data)}(Tuple(normalize_measurements(v, bin_size, string(k)) for (k, v) in pairs(data)))
    elseif isa(data, Tuple)
        return tuple((normalize_measurements(x, bin_size, key) for x in data)...)
    else
        try
            return data ./ bin_size
        catch
            return data
        end
    end
end








