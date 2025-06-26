@doc raw"""

    process_scalar_measurements( datafolder::String,
                                 measurement::String,
                                 output_csv.String )::Nothing

For energy, density, or double occupancy measurements, averages over all measurements
in a certain bin. 

"""
function process_scalar_measurements(
    datafolder::String,
    measurement::String,        # "energy", "density", or "double_occ"
    output_csv::String
)::Nothing
    # All now live in 'simulation' directory
    measurement_dir = joinpath(datafolder, "simulation", measurement)

    # Find and sort JLD2 files
    files = sort(glob("bin-*_pID-*.jld2", measurement_dir))

    # Prepare DataFrame for results
    results = DataFrame(index=Int[], local_mean=Float64[], local_std=Float64[])

    for file in files
        bin_match = match(r"bin-(\d+)_pID-\d+\.jld2", basename(file))
        bin_idx = parse(Int, bin_match.captures[1])

        # Gather scalar measurements from all steps
        values = Float64[]
        jldopen(file, "r") do jld
            for key in keys(jld)
                val = jld[key]
                if isa(val, Number)
                    push!(values, val)
                else
                    @warn "Skipping non-numeric value in $file at $key"
                end
            end
        end

        stats = jackknife_stats(values)
        push!(results, (bin_idx, stats.mean, stats.std))
    end

    sort!(results, :index)
    rename!(results, [:index, :local_mean, :local_std])
    CSV.write(output_csv, results)

    return nothing
end


@doc raw"""

    process_optimization_measurements()

For logarithmic derivative measurements, calculates the mean value.


"""
function process_optimization_measurements()

end


@doc raw"""

    process_correlation_measurements()

For either density-density or spin-spin correlation data, calculates the static structure factor.

"""
function process_correlation_measurements(correlation_type::String, model_geometry::ModelGeometry)
    # calculate all momentum points that will used in the FT
    q_points = calc_k_points(model_geometry.unit_cell, model_geometry.lattice)

    if correlation_type == "density"
        # read in the site-dependent density data 

    elseif correlation_type == "spin"
        # read in the site-dependent spin data

    end

    # calculate_structure_factor()

    # write structure factor info to file
end


function process_measurements()

    # set the walkers to iterate over
    # this will be used for MPI
    if isempty(pIDs)

        # get the number of MPI walkers
        N_walkers = get_num_walkers(folder)

        # get the pIDs
        pIDs = collect(0:(N_walkers-1))
    end

    process_scalar_measurements(datafolder, "density", output_csv)
    process_scalar_measurements(datafolder, "double_occ", output_csv)
    process_scalar_measurements(datafolder, "energy", output_csv)

    # # only call if using
    # process_optimization_measurements()
    # process_correlation_measurements()
end