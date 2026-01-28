@doc raw"""

    process_measurements( measurement_container::NamedTuple,
                          simulation_info::SimulationInfo,
                          determinantal_parameters::DeterminantalParameters,
                          model_geometry::ModelGeometry )

Processes all simulation and optimization measurements by organinzing them into CSV files.

- `measurement_container::NamedTuple`: contains measurement quantities.
- `simulation_info::SimulationInfo`: contains all simulation info.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 

"""
function process_measurements(
    measurement_container, 
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters,
    model_geometry::ModelGeometry
)
    (; datafolder, pID) = simulation_info
    (; N_opt, opt_bin_size) = measurement_container

    bins_to_avg = div(N_opt, opt_bin_size)

    # merge `opt` and `sim` bins
    merge_bin_measurements!(datafolder * "/simulation/", pID)
    if haskey(measurement_container.correlation_measurements, "density") || haskey(measurement_container.correlation_measurements, "spin")
        merge_bin_measurements!(datafolder * "/correlation/", pID)
    end


    # process all scalar measurements
    process_scalar_measurements(
        datafolder, 
        pID,
        "local_energy";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder, 
        pID,
        "double_occ";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder, 
        pID,
        "global_density";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder, 
        pID,
        "pconfig"
    )

    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        process_scalar_measurements(
            datafolder, 
            pID,
            "local_spin-z";
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        process_scalar_measurements(
            datafolder, 
            pID,
            "site-dependent_spin-z";
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        process_scalar_measurements(
            datafolder, 
            pID,
            "site-dependent_density";
            bins_to_avg = bins_to_avg
        )
    end

   # process optimization measurements
    process_optimization_measurements(
        datafolder, 
        pID,
        determinantal_parameters
    )

    if haskey(measurement_container.correlation_measurements, "density")
        process_correlation_measurements(
            datafolder, 
            pID,
            "density", 
            model_geometry;
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.correlation_measurements, "spin")
        process_correlation_measurements(
            datafolder, 
            pID,
            "spin", 
            model_geometry;
            bins_to_avg = bins_to_avg
        )
    end

   return nothing
end


@doc raw"""

    process_measurements( measurement_container::NamedTuple,
                          simulation_info::SimulationInfo,
                          determinantal_parameters::DeterminantalParameters,
                          jastrow_parameters::JastrowParameters,
                          model_geometry::ModelGeometry )

Processes all simulation and optimization measurements by organinzing them into CSV files.

- `measurement_container::NamedTuple`: contains measurement quantities.
- `simulation_info::SimulationInfo`: contains all simulation info.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: set of Jastrow parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 

"""
function process_measurements(
    measurement_container, 
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters::JastrowParameters,
    model_geometry::ModelGeometry
)
    (; datafolder, pID) = simulation_info
    (; N_opt, opt_bin_size) = measurement_container

    bins_to_avg = div(N_opt, opt_bin_size)

    # merge `opt` and `sim` bins
    merge_bin_measurements!(datafolder * "/simulation/", pID)
    if haskey(measurement_container.correlation_measurements, "density") || haskey(measurement_container.correlation_measurements, "spin")
        merge_bin_measurements!(datafolder * "/correlation/", pID)
    end

    # process all scalar measurements
    process_scalar_measurements(
        datafolder, 
        pID,
        "local_energy";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder,
        pID, 
        "double_occ";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder,
        pID, 
        "global_density";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder,
        pID, 
        "pconfig"
    )

    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        process_scalar_measurements(
            datafolder,
            pID, 
            "local_spin-z";
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        process_scalar_measurements(
            datafolder, 
            pID,
            "site-dependent_spin-z";
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        process_scalar_measurements(
            datafolder, 
            pID,
            "site-dependent_density";
            bins_to_avg = bins_to_avg
        )
    end

   # process optimization measurements
    process_optimization_measurements(
        datafolder, 
        pID,
        determinantal_parameters,
        jastrow_parameters
    )

    if haskey(measurement_container.correlation_measurements, "density")
        process_correlation_measurements(
            datafolder, 
            pID,
            "density", 
            model_geometry;
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.correlation_measurements, "spin")
        process_correlation_measurements(
            datafolder, 
            pID,
            "spin", 
            model_geometry;
            bins_to_avg = bins_to_avg
        )
    end

    return nothing
end


@doc raw"""

    process_measurements( measurement_container::NamedTuple,
                          simulation_info::SimulationInfo,
                          determinantal_parameters::DeterminantalParameters,
                          jastrow_parameters_1::JastrowParameters,
                          jastrow_parameters_2::JastrowParameters,
                          model_geometry::ModelGeometry )

Processes all simulation and optimization measurements by organinzing them into CSV files.

- `measurement_container::NamedTuple`: contains measurement quantities.
- `simulation_info::SimulationInfo`: contains all simulation info.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters`: first set of Jastrow parameters.
- `jastrow_parameters_2::JastrowParameters`: second set of Jastrow parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 

"""
function process_measurements(
    measurement_container, 
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters_1::JastrowParameters,
    jastrow_parameters_2::JastrowParameters,
    model_geometry::ModelGeometry
)
    (; datafolder, pID) = simulation_info
    (; N_opt, opt_bin_size) = measurement_container

    bins_to_avg = div(N_opt, opt_bin_size)

    # merge `opt` and `sim` bins
    merge_bin_measurements!(datafolder * "/simulation/", pID)
    if haskey(measurement_container.correlation_measurements, "density") || haskey(measurement_container.correlation_measurements, "spin")
        merge_bin_measurements!(datafolder * "/correlation/", pID)
    end

    # process all scalar measurements
    process_scalar_measurements(
        datafolder, 
        pID,
        "local_energy";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder,
        pID, 
        "double_occ";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder,
        pID, 
        "global_density";
        bins_to_avg = bins_to_avg
    )
    process_scalar_measurements(
        datafolder,
        pID, 
        "pconfig"
    )

    if haskey(measurement_container.simulation_measurements, "local_spin-z")
        process_scalar_measurements(
            datafolder,
            pID, 
            "local_spin-z";
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.simulation_measurements, "site-dependent_spin-z")
        process_scalar_measurements(
            datafolder, 
            pID,
            "site-dependent_spin-z";
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.simulation_measurements, "site-dependent_density")
        process_scalar_measurements(
            datafolder, 
            pID,
            "site-dependent_density";
            bins_to_avg = bins_to_avg
        )
    end

   # process optimization measurements
    process_optimization_measurements(
        datafolder, 
        pID,
        determinantal_parameters,
        jastrow_parameters_1,
        jastrow_parameters_2
    )

    if haskey(measurement_container.correlation_measurements, "density")
        process_correlation_measurements(
            datafolder, 
            pID,
            "density", 
            model_geometry;
            bins_to_avg = bins_to_avg
        )
    end

    if haskey(measurement_container.correlation_measurements, "spin")
        process_correlation_measurements(
            datafolder, 
            pID,
            "spin", 
            model_geometry;
            bins_to_avg = bins_to_avg
        )
    end

    return nothing
end


@doc raw"""

    process_scalar_measurements( datafolder::T,
                                 pID::I,
                                 measurement::T;
                                 N_bins::Union{I, Nothing}=nothing;
                                 bins_to_avg::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

Write binned simulation measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `pID::I`: processor ID/MPI rank
- `measurement::T`: "local_energy", "double_occ", "global_density",
  "site-dependent_density", "pconfig", "local_spin-z", or "site-dependent_spin-z".
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.
- `bins_to_avg::Union{I, Nothing}=nothing`: number of bins to average over when performing jackknife resampling.

"""
function process_scalar_measurements(
    datafolder::T,
    pID::I,
    measurement::T,
    N_bins::Union{I, Nothing}=nothing;
    bins_to_avg::Union{I, Nothing}=nothing
) where {T<:AbstractString, I<:Integer}

    sim_file = joinpath(datafolder, "simulation", "bin_measurements_rank-$(pID).h5")
    @assert isfile(sim_file) "HDF5 file not found: $sim_file"

    allowed = Set([
        "local_energy",
        "double_occ",
        "global_density",
        "site-dependent_density",
        "pconfig",
        "local_spin-z",
        "site-dependent_spin-z"
    ])
    @assert measurement in allowed "Unknown measurement: $measurement"

    results = Vector{Any}()

    h5open(sim_file, "r") do f
        path_for_bin(bin_num) = "/$measurement/bin-$bin_num"

        if N_bins === nothing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    bins = collect(1:length(results))

    if measurement == "local_energy"
        # Complex scalar per bin
        mean_r = real.(results)
        mean_i = imag.(results)
        df = DataFrame(BIN = bins, MEAN_R = mean_r, MEAN_I = mean_i)

    elseif measurement in ("site-dependent_spin-z", "site-dependent_density")
        # Vector-valued per bin (one value per site)
        @assert all(!ismissing(r) for r in results) """
        Missing bins are not supported for site-dependent measurements
        """

        n_sites = length(results[1])
        @assert all(length(r) == n_sites for r in results) """
        Inconsistent site vector lengths across bins
        """

        # bins × sites matrix
        data = reduce(hcat, results)'

        df = DataFrame(BIN = bins)
        for s in 1:n_sites
            df[!, Symbol("SITE_$s")] = data[:, s]
        end

    else
        # Real scalar per bin
        df = DataFrame(BIN = bins, MEAN_R = results)
    end

    output_csv = joinpath(
        datafolder,
        "simulation",
        measurement * "_stats_rank-$(pID).csv"
    )

    CSV.write(output_csv, df)

    ## JACKKNIFE RESAMPLING OF SELECT OBSERVABLES ##
    # allowable observables
    jk_measurements = Set([
        "local_energy",
        "double_occ",
        "global_density",
        "local_spin-z"
    ])

    if measurement in jk_measurements
        if measurement == "local_energy"
            samples = real.(results)
        else
            samples = Float64.(results)
        end

        # collect samples
        samples = @view samples[bins_to_avg+1:end]

        # perform jackknife resampling
        mean, err = jackknife(identity, samples)

        # write to file
        summary_file = joinpath(datafolder, "vmc_summary_rank-$(pID).csv")

        summary_df = DataFrame(
            MEASUREMENT = [measurement],
            MEAN_R      = [mean],
            STD         = [err]
        )

        if isfile(summary_file)
            CSV.write(summary_file, summary_df; append=true)
        else
            CSV.write(summary_file, summary_df)
        end
    end

    ## JACKKNIFE RESAMPLING FOR SITE-DEPENDENT OBSERVABLES ##
    site_jk_measurements = Set([
        "site-dependent_density",
        "site-dependent_spin-z"
    ])

    if measurement in site_jk_measurements
        n_sites = length(results[1])

        # bins × sites matrix
        data = reduce(hcat, results)'  # size: (nbins, n_sites)

        # discard thermalization
        data = @view data[bins_to_avg+1:end, :]

        # output file
        site_summary_file = joinpath(
            datafolder,
            "site-dependent_summary_rank-$(pID).csv"
        )

        # load existing file if present
        if isfile(site_summary_file)
            site_df = CSV.read(site_summary_file, DataFrame)
        else
            site_df = DataFrame(
                SITE   = 1:n_sites,
                MEAN_D = fill(NaN, n_sites),
                STD_D  = fill(NaN, n_sites),
                MEAN_S = fill(NaN, n_sites),
                STD_S  = fill(NaN, n_sites)
            )
        end

        for s in 1:n_sites
            samples = data[:, s]
            mean, err = jackknife(identity, samples)

            if measurement == "site-dependent_density"
                site_df.MEAN_D[s] = mean
                site_df.STD_D[s]  = err
            elseif measurement == "site-dependent_spin-z"
                site_df.MEAN_S[s] = mean
                site_df.STD_S[s]  = err
            end
        end

        CSV.write(site_summary_file, site_df)
    end

    return nothing
end


@doc raw"""

    process_optimization_measurements( datafolder::T,
                                       pID::I,
                                       determinantal_parameters::DeterminantalParameters;
                                       N_bins::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

Writes binned optimization measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `pID::I`: processor ID/MPI rank
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.

"""
function process_optimization_measurements(
    datafolder::T,
    pID::I,
    determinantal_parameters::DeterminantalParameters;
    N_bins::Union{I, Nothing}=nothing
    
) where {T<:AbstractString, I<:Integer}
    opt_file = joinpath(datafolder, "optimization", "opt_bin_measurements_rank-$(pID).h5")
    @assert isfile(opt_file) "HDF5 file not found: $opt_file"

    # helper function to convert filenames to ASCII
    convert_name(name::AbstractString) = replace(name, "Δ" => "delta", "μ" => "mu", " " => "_")

    # extract parameter information 
    current_det_pars = determinantal_parameters.det_pars
    param_names = Tuple(keys(current_det_pars))  

    # helper function to reconstruct a parameter NamedTuple
    function reconstruct_from_vpars(vpars::AbstractVector)
        reconstructed = Dict{Symbol, Any}()
        pos = 1
        for pname in param_names
            old_val = current_det_pars[pname]
            if old_val isa AbstractVector
                n = length(old_val)
                if pos + n - 1 > length(vpars)
                    throw(ArgumentError("vpars too short for parameter $(pname): need $n elements starting at $pos"))
                end
                slice = vpars[pos:pos + n - 1]
                reconstructed[pname] = collect(slice) 
                pos += n
            else
                if pos > length(vpars)
                    throw(ArgumentError("vpars too short for scalar parameter $(pname) at position $pos"))
                end
                reconstructed[pname] = vpars[pos]
                pos += 1
            end
        end
        if pos - 1 != length(vpars)
            throw(ArgumentError("Length mismatch: consumed $(pos - 1) elements but vpars has length $(length(vpars))."))
        end
        return reconstructed
    end

    results = Vector{Any}()

    h5open(opt_file, "r") do f
        path_for_bin(bin_num) = "/parameters/bin-$bin_num"

        if N_bins === nothing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    bins = collect(1:length(results))

    # reconstruct parameter dictionary
    per_bin_params = Vector{Union{Dict{Symbol,Any}, Missing}}(undef, length(results))
    for (i, entry) in enumerate(results)
        if entry === missing
            per_bin_params[i] = missing
        else
            vpars = collect(entry)
            per_bin_params[i] = reconstruct_from_vpars(vpars)
        end
    end

    outdir = joinpath(datafolder, "optimization")
    isdir(outdir) || mkpath(outdir)

    for pname in param_names
        converted = convert_name(String(pname))
        outcsv = joinpath(outdir, "parameter_$(converted)_stats_rank-$(pID).csv")

        template_val = current_det_pars[pname]
        if template_val isa AbstractVector
            # for vector parameters
            n = length(template_val)
            el_is_complex = Base.eltype(template_val) <: Complex || any(x->x isa Complex, template_val)

            df_cols = Dict{Symbol, Any}()
            df_cols[:BIN] = bins
            for j in 1:n
                df_cols[Symbol("MEAN_R_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
            end
            if el_is_complex
                for j in 1:n
                    df_cols[Symbol("MEAN_I_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
                end
            end

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    for j in 1:n
                        df_cols[Symbol("MEAN_R_$(j)")][i] = missing
                    end
                    if el_is_complex
                        for j in 1:n
                            df_cols[Symbol("MEAN_I_$(j)")][i] = missing
                        end
                    end
                else
                    val = entry[pname]
                    if !(val isa AbstractVector)
                        throw(ArgumentError("Expected vector for parameter $(pname) but got scalar in bin $i"))
                    end
                    for j in 1:n
                        comp = val[j]
                        if comp isa Complex
                            df_cols[Symbol("MEAN_R_$(j)")][i] = real(comp)
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = imag(comp)
                            end
                        else
                            df_cols[Symbol("MEAN_R_$(j)")][i] = comp
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = 0.0
                            end
                        end
                    end
                end
            end

            df = DataFrame(df_cols)
            CSV.write(outcsv, df)

        else
            # for scalar parameters
            el_is_complex = (template_val isa Complex) ||
                            any(b -> (b !== missing && (b[pname] isa Complex)), per_bin_params)

            MEAN_R = Vector{Union{Float64, Missing}}(undef, length(bins))
            MEAN_I = el_is_complex ? Vector{Union{Float64, Missing}}(undef, length(bins)) : nothing

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    MEAN_R[i] = missing
                    if el_is_complex
                        MEAN_I[i] = missing
                    end
                else
                    val = entry[pname]
                    if val isa Complex
                        MEAN_R[i] = real(val)
                        if el_is_complex
                            MEAN_I[i] = imag(val)
                        end
                    else
                        MEAN_R[i] = val
                        if el_is_complex
                            MEAN_I[i] = 0.0
                        end
                    end
                end
            end

            if el_is_complex
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R, MEAN_I = MEAN_I)
            else
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R)
            end
            CSV.write(outcsv, df)
        end
    end

    return nothing
end


@doc raw"""

    process_optimization_measurements( datafolder::T,
                                       pID::I,
                                       determinantal_parameters::DeterminantalParameters,
                                       jastrow_parameters::JastrowParameters;
                                       N_bins::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

Writes binned optimization measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `pID::I`: processor ID/MPI rank
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: set of Jastrow parameters.
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.

"""
function process_optimization_measurements(
    datafolder::T,
    pID::I,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters::JastrowParameters;
    N_bins::Union{I, Nothing}=nothing
) where {T<:AbstractString, I<:Integer}
    opt_file = joinpath(datafolder, "optimization", "opt_bin_measurements_rank-$(pID).h5")
    @assert isfile(opt_file) "HDF5 file not found: $opt_file"

    # helper function to convert filenames to ASCII
    convert_name(name::AbstractString) = replace(name, "Δ" => "delta", "μ" => "mu", " " => "_")

    # extract parameter information 
    current_det_pars = determinantal_parameters.det_pars
    param_names = Tuple(keys(current_det_pars))

    # helper function to reconstruct a parameter NamedTuple and Jastrow parameters
    function reconstruct_from_vpars(vpars::AbstractVector)
        reconstructed = Dict{Symbol, Any}()
        pos = 1
        for pname in param_names
            old_val = current_det_pars[pname]
            if old_val isa AbstractVector
                n = length(old_val)
                if pos + n - 1 > length(vpars)
                    throw(ArgumentError("vpars too short for parameter $(pname): need $n elements starting at $pos"))
                end
                slice = vpars[pos:pos + n - 1]
                reconstructed[pname] = collect(slice)
                pos += n
            else
                if pos > length(vpars)
                    throw(ArgumentError("vpars too short for scalar parameter $(pname) at position $pos"))
                end
                reconstructed[pname] = vpars[pos]
                pos += 1
            end
        end
        # If there are extra elements after the determinantal parameters, treat them as jastrow tail
        if pos <= length(vpars)
            jtail = collect(vpars[pos:end])
        else
            jtail = missing
        end
        return reconstructed, jtail, (pos - 1)
    end

    results = Vector{Any}()

    h5open(opt_file, "r") do f
        path_for_bin(bin_num) = "/parameters/bin-$bin_num"

        if N_bins === nothing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    bins = collect(1:length(results))

    # reconstruct parameter dictionary and collect Jastrow parameters
    per_bin_params = Vector{Union{Dict{Symbol,Any}, Missing}}(undef, length(results))
    per_bin_jpars  = Vector{Union{Vector{Float64}, Vector{Any}, Missing}}(undef, length(results))
    consumed_lengths = Vector{Int}(undef, length(results))

    for (i, entry) in enumerate(results)
        if entry === missing
            per_bin_params[i] = missing
            per_bin_jpars[i]  = missing
            consumed_lengths[i] = 0
        else
            vpars = collect(entry)
            det_dict, jtail, consumed = reconstruct_from_vpars(vpars)
            per_bin_params[i] = det_dict
            per_bin_jpars[i] = jtail
            consumed_lengths[i] = consumed
        end
    end

    # get number of Jastrow parameters
    n_jpars = nothing
    try
        if hasproperty(jastrow_parameters, :num_jpar_opts)
            n_jpars = Int(getproperty(jastrow_parameters, :num_jpar_opts))
        end
    catch
        n_jpars = nothing
    end

    if n_jpars === nothing
        try
            if hasproperty(jastrow_parameters, :jpars)
                n_jpars = length(getproperty(jastrow_parameters, :jpars))
            elseif hasproperty(jastrow_parameters, :vpars)
                n_jpars = length(getproperty(jastrow_parameters, :vpars))
            elseif hasproperty(jastrow_parameters, :parameters)
                n_jpars = length(getproperty(jastrow_parameters, :parameters))
            end
        catch
            n_jpars = nothing
        end
    end

    if n_jpars === nothing
        for jtail in per_bin_jpars
            if jtail !== missing
                n_jpars = length(jtail)
                break
            end
        end
    end
    n_jpars = n_jpars === nothing ? 0 : n_jpars

    outdir = joinpath(datafolder, "optimization")
    isdir(outdir) || mkpath(outdir)

    # write determinantal parameters
    for pname in param_names
        converted = convert_name(String(pname))
        outcsv = joinpath(outdir, "parameter_$(converted)_stats_rank-$(pID).csv")

        template_val = current_det_pars[pname]
        if template_val isa AbstractVector
            # for vector parameters
            n = length(template_val)
            el_is_complex = Base.eltype(template_val) <: Complex || any(x->x isa Complex, template_val)

            df_cols = Dict{Symbol, Any}()
            df_cols[:BIN] = bins
            for j in 1:n
                df_cols[Symbol("MEAN_R_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
            end
            if el_is_complex
                for j in 1:n
                    df_cols[Symbol("MEAN_I_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
                end
            end

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    for j in 1:n
                        df_cols[Symbol("MEAN_R_$(j)")][i] = missing
                    end
                    if el_is_complex
                        for j in 1:n
                            df_cols[Symbol("MEAN_I_$(j)")][i] = missing
                        end
                    end
                else
                    val = entry[pname]
                    if !(val isa AbstractVector)
                        throw(ArgumentError("Expected vector for parameter $(pname) but got scalar in bin $i"))
                    end
                    for j in 1:n
                        comp = val[j]
                        if comp isa Complex
                            df_cols[Symbol("MEAN_R_$(j)")][i] = real(comp)
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = imag(comp)
                            end
                        else
                            df_cols[Symbol("MEAN_R_$(j)")][i] = comp
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = 0.0
                            end
                        end
                    end
                end
            end

            df = DataFrame(df_cols)
            CSV.write(outcsv, df)

        else
            # for scalar parameters
            el_is_complex = (template_val isa Complex) ||
                            any(b -> (b !== missing && (b[pname] isa Complex)), per_bin_params)

            MEAN_R = Vector{Union{Float64, Missing}}(undef, length(bins))
            MEAN_I = el_is_complex ? Vector{Union{Float64, Missing}}(undef, length(bins)) : nothing

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    MEAN_R[i] = missing
                    if el_is_complex
                        MEAN_I[i] = missing
                    end
                else
                    val = entry[pname]
                    if val isa Complex
                        MEAN_R[i] = real(val)
                        if el_is_complex
                            MEAN_I[i] = imag(val)
                        end
                    else
                        MEAN_R[i] = val
                        if el_is_complex
                            MEAN_I[i] = 0.0
                        end
                    end
                end
            end

            if el_is_complex
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R, MEAN_I = MEAN_I)
            else
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R)
            end
            CSV.write(outcsv, df)
        end
    end

    # write Jastrow parameters
    jtype = try
        string(getproperty(jastrow_parameters, :jastrow_type))
    catch
        "jastrow"
    end
    jfilename = joinpath(outdir, "$(jtype)_jastrow_parameters_rank-$(pID).csv")

    if n_jpars > 0
        jcolnames = [Symbol("MEAN_V_$(k)") for k in 1:n_jpars]
        jcols = Dict{Symbol, Any}()
        jcols[:BIN] = bins
        for sym in jcolnames
            jcols[sym] = Vector{Union{Float64, Missing}}(undef, length(bins))
        end

        for i in 1:length(bins)
            jtail = per_bin_jpars[i]
            if jtail === missing
                for sym in jcolnames
                    jcols[sym][i] = missing
                end
            else
                for k in 1:n_jpars
                    if k <= length(jtail)
                        val = jtail[k]
                        if val isa Complex
                            jcols[jcolnames[k]][i] = real(val)
                        else
                            jcols[jcolnames[k]][i] = Float64(val)
                        end
                    else
                        jcols[jcolnames[k]][i] = missing
                    end
                end
            end
        end

        jdf = DataFrame(jcols)
        CSV.write(jfilename, jdf)
    else
        CSV.write(jfilename, DataFrame(BIN = bins))
    end

    return nothing
end


@doc raw"""

    process_optimization_measurements( datafolder::T,
                                       pID::I,
                                       determinantal_parameters::DeterminantalParameters,
                                       jastrow_parameters_1::JastrowParameters,
                                       jastrow_parameters_2::JastrowParameters;
                                       N_bins::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

Writes binned optimization measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `pID::I`: processor ID/MPI rank
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters`: first set of Jastrow parameters.
- `jastrow_parameters_2::JastrowParameters`: second set of Jastrow parameters.
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.

"""
function process_optimization_measurements(
    datafolder::T,
    pID::I,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters_1::JastrowParameters,
    jastrow_parameters_2::JastrowParameters;
    N_bins::Union{I, Nothing}=nothing
) where {T<:AbstractString, I<:Integer}

    opt_file = joinpath(datafolder, "optimization", "opt_bin_measurements_rank-$(pID).h5")
    @assert isfile(opt_file) "HDF5 file not found: $opt_file"

    convert_name(name::AbstractString) = replace(name, "Δ" => "delta", "μ" => "mu", " " => "_")

    current_det_pars = determinantal_parameters.det_pars
    param_names = collect(keys(current_det_pars))

    # determine Jastrow lengths
    function infer_n_jpars(jp)
        try
            if hasproperty(jp, :num_jpar_opts)
                return Int(getproperty(jp, :num_jpar_opts))
            elseif hasproperty(jp, :jpars)
                return length(getproperty(jp, :jpars))
            elseif hasproperty(jp, :vpars)
                return length(getproperty(jp, :vpars))
            elseif hasproperty(jp, :parameters)
                return length(getproperty(jp, :parameters))
            end
        catch
        end
        return 0
    end

    n_j1 = infer_n_jpars(jastrow_parameters_1)
    n_j2 = infer_n_jpars(jastrow_parameters_2)

    function reconstruct_from_vpars(vpars::AbstractVector,
                                    current_det_pars,
                                    param_names::Vector{Symbol},
                                    n_j1::Int,
                                    n_j2::Int)

        reconstructed = Dict{Symbol, Any}()
        pos = 1

        for pname in param_names
            old_val = current_det_pars[pname]
            if old_val isa AbstractVector
                n = length(old_val)
                reconstructed[pname] = collect(vpars[pos:pos+n-1])
                pos += n
            else
                reconstructed[pname] = vpars[pos]
                pos += 1
            end
        end

        j1 = pos <= length(vpars) ?
            collect(vpars[pos : min(pos + n_j1 - 1, length(vpars))]) :
            missing
        pos += n_j1

        j2 = pos <= length(vpars) ?
            collect(vpars[pos : min(pos + n_j2 - 1, length(vpars))]) :
            missing

        return reconstructed, j1, j2
    end

    results = Vector{Any}()

    h5open(opt_file, "r") do f
        path_for_bin(bin_num) = "/parameters/bin-$bin_num"

        if N_bins === nothing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    bins = collect(1:length(results))

    per_bin_params = Vector{Union{Dict{Symbol,Any}, Missing}}(undef, length(results))
    per_bin_j1 = Vector{Union{Vector{Any}, Missing}}(undef, length(results))
    per_bin_j2 = Vector{Union{Vector{Any}, Missing}}(undef, length(results))

    for (i, entry) in enumerate(results)
        if entry === missing
            per_bin_params[i] = missing
            per_bin_j1[i] = missing
            per_bin_j2[i] = missing
        else
            vpars = collect(entry)
            det_dict, j1, j2 = reconstruct_from_vpars(
                vpars, current_det_pars, param_names, n_j1, n_j2
            )
            per_bin_params[i] = det_dict
            per_bin_j1[i] = j1
            per_bin_j2[i] = j2
        end
    end

    outdir = joinpath(datafolder, "optimization")
    isdir(outdir) || mkpath(outdir)

    # --- determinantal CSVs ---
    for pname in param_names
        converted = convert_name(String(pname))
        outcsv = joinpath(outdir, "parameter_$(converted)_stats_rank-$(pID).csv")

        template_val = current_det_pars[pname]
        if template_val isa AbstractVector
            n = length(template_val)
            el_is_complex = Base.eltype(template_val) <: Complex ||
                any(b -> b !== missing && any(x -> x isa Complex, b[pname]),
                    per_bin_params)

            df_cols = Dict{Symbol, Any}()
            df_cols[:BIN] = bins
            for j in 1:n
                df_cols[Symbol("MEAN_R_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
            end
            if el_is_complex
                for j in 1:n
                    df_cols[Symbol("MEAN_I_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
                end
            end

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    for j in 1:n
                        df_cols[Symbol("MEAN_R_$(j)")][i] = missing
                        if el_is_complex
                            df_cols[Symbol("MEAN_I_$(j)")][i] = missing
                        end
                    end
                else
                    val = entry[pname]
                    for j in 1:n
                        comp = val[j]
                        df_cols[Symbol("MEAN_R_$(j)")][i] = comp isa Complex ? real(comp) : comp
                        if el_is_complex
                            df_cols[Symbol("MEAN_I_$(j)")][i] = comp isa Complex ? imag(comp) : 0.0
                        end
                    end
                end
            end

            CSV.write(outcsv, DataFrame(df_cols))
        else
            el_is_complex = template_val isa Complex ||
                any(b -> b !== missing && (b[pname] isa Complex), per_bin_params)

            MEAN_R = Vector{Union{Float64, Missing}}(undef, length(bins))
            MEAN_I = el_is_complex ? Vector{Union{Float64, Missing}}(undef, length(bins)) : nothing

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    MEAN_R[i] = missing
                    if el_is_complex; MEAN_I[i] = missing; end
                else
                    val = entry[pname]
                    MEAN_R[i] = val isa Complex ? real(val) : val
                    if el_is_complex
                        MEAN_I[i] = val isa Complex ? imag(val) : 0.0
                    end
                end
            end

            df = el_is_complex ?
                DataFrame(BIN=bins, MEAN_R=MEAN_R, MEAN_I=MEAN_I) :
                DataFrame(BIN=bins, MEAN_R=MEAN_R)

            CSV.write(outcsv, df)
        end
    end

    # --- write Jastrow CSVs ---
    function write_jastrow(per_bin, n_j, jp)
        jtype = try string(getproperty(jp, :jastrow_type)) catch; "jastrow" end
        jfile = joinpath(outdir, "$(jtype)_jastrow_parameters_rank-$(pID).csv")

        if n_j > 0
            jcols = Dict{Symbol,Any}()
            jcols[:BIN] = bins
            for k in 1:n_j
                jcols[Symbol("MEAN_V_$(k)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
            end

            for i in 1:length(bins)
                row = per_bin[i]
                for k in 1:n_j
                    if row === missing || k > length(row)
                        jcols[Symbol("MEAN_V_$(k)")][i] = missing
                    else
                        v = row[k]
                        jcols[Symbol("MEAN_V_$(k)")][i] = v isa Complex ? real(v) : Float64(v)
                    end
                end
            end

            CSV.write(jfile, DataFrame(jcols))
        else
            CSV.write(jfile, DataFrame(BIN=bins))
        end
    end

    write_jastrow(per_bin_j1, n_j1, jastrow_parameters_1)
    write_jastrow(per_bin_j2, n_j2, jastrow_parameters_2)

    return nothing
end


@doc raw"""

    process_correlation_measurements( datafolder::T,
                                      pID::I
                                      correlation_type::T,
                                      model_geometry::ModelGeometry;
                                      bins_to_avg::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

For either density-density or spin-spin correlation data, calculates the static structure factor 
``N(\mathbf{q}) = \langle \hat{n}_{-\mathbf{q}\hat{n}_{\mathbf{q}}\rangle`` or 
``S(\mathbf{q}) = \langle \hat{S}_{-\mathbf{q}\hat{S}_{\mathbf{q}}\rangle``, respectively.

- `datafolder::T`: path to folder where simulation files are written.
- `pID::I`: processor ID/MPI rank
- `correlation_type::T`: either "density" or "spin".
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `bins_to_avg::Union{I, Nothing}=nothing`: number of bins to average over when performing jackknife resampling.

"""
function process_correlation_measurements(
    datafolder::T,
    pID::I,
    correlation_type::T,
    model_geometry::ModelGeometry;
    bins_to_avg::Union{I, Nothing}=nothing
) where {T<:AbstractString, I<:Integer}
    unit_cell = model_geometry.unit_cell
    lattice = model_geometry.lattice
    N = model_geometry.lattice.N

    # compute q-points used in FT
    q_points = calc_k_points(unit_cell, lattice)
    nq = length(q_points)

    corr_file = joinpath(datafolder, "correlation", "bin_measurements_rank-$(pID).h5")
    @assert isfile(corr_file) "HDF5 file not found: $corr_file"

    function parse_bin_number(name::AbstractString)::Int
        if (m = match(r"bin-(\d+)", name)) !== nothing
            return parse(Int, m.captures[1])
        end
        try
            return parse(Int, split(name, r"[-_]") |> last)
        catch err
            return 0
        end
    end

    rows = Vector{Vector{Float64}}()
    binnums = Int[]

    h5open(corr_file, "r") do fh
        grp_path = "/" * correlation_type
        if !haskey(fh, grp_path)
            @warn "Group $grp_path not found in $corr_file"
            return nothing
        end
        grp = fh[grp_path]

        # iterate bins (keys(grp) yields group names like "bin-0", "bin-1", ...)
        for bin_group in keys(grp)
            dset = grp[bin_group]
            # read dataset — expect a matrix; try to coerce to Matrix{Float64}
            raw = read(dset)
            corr = try
                Array{Float64}(raw)
            catch err
                # try reshaping if raw is a vector representing a square matrix
                nd = ndims(raw)
                if nd == 1
                    len = length(raw)
                    s = Int(round(sqrt(len)))
                    if s * s == len
                        corr_mat = reshape(collect(raw), s, s)
                        Array{Float64}(corr_mat)
                    else
                        rethrow(err)
                    end
                else
                    rethrow(err)
                end
            end

            # ensure corr is a 2D matrix
            if ndims(corr) != 2
                @warn "Dataset for $bin_group does not produce a 2D matrix; skipping"
                continue
            end

            # calculate structure factor vector
            Sq = calculate_structure_factor(corr, q_points, N, unit_cell, lattice)

            # parse bin number robustly
            binnum = parse_bin_number(String(bin_group))

            push!(binnums, binnum)
            push!(rows, vcat(Float64(binnum), Sq))  # first column is bin number
        end
    end

    if isempty(rows)
        @warn "No bins processed under /$correlation_type in $corr_file"
        return nothing
    end

    # Build dataframe: first column "bin", then q columns
    colnames = ["bin"]
    for i in 1:nq
        push!(colnames, "q$(i)")
    end

    # transpose rows -> columns for DataFrame constructor
    mat = reduce(hcat, rows)'   # rows is Vector{Vector}, produce matrix rows×cols then transpose
    df = DataFrame(mat, Symbol.(colnames))

    # sort by bin number
    sort!(df, :bin)

    # write CSV
    out_csv = joinpath(datafolder, "correlation", "$(correlation_type)_static_structure_factor_stats_rank-$(pID).csv")
    CSV.write(out_csv, df)

    ## JACKKNIFE RESAMPLING OF STRUCTURE FACTORS ##
    # Determine which columns to write
    is_density = correlation_type == "density"
    is_spin    = correlation_type == "spin"

    # Extract q-columns as a matrix: bins × nq
    qcols = Symbol.("q" .* string.(1:nq))
    data = Matrix(df[:, qcols])

    # discard thermalization bins
    data = @view data[bins_to_avg+1:end, :]

    # output file
    summary_file = joinpath(
        datafolder,
        "structure_factor_summary_rank-$(pID).csv"
    )

    # load or initialize summary DataFrame
    if isfile(summary_file)
        sum_df = CSV.read(summary_file, DataFrame)
    else
        sum_df = DataFrame(
            Q      = 1:nq,
            MEAN_D = fill(NaN, nq),
            STD_D  = fill(NaN, nq),
            MEAN_S = fill(NaN, nq),
            STD_S  = fill(NaN, nq)
        )
    end

    # sanity check
    @assert nrow(sum_df) == nq "Mismatch in number of q-points in structure factor summary"

    # jackknife per q
    for iq in 1:nq
        samples = data[:, iq]
        mean, err = jackknife(identity, samples)

        if is_density
            sum_df.MEAN_D[iq] = mean
            sum_df.STD_D[iq]  = err
        elseif is_spin
            sum_df.MEAN_S[iq] = mean
            sum_df.STD_S[iq]  = err
        end
    end

    # write back
    CSV.write(summary_file, sum_df)

    return nothing
end


@doc raw"""

    calculate_structure_factor( corr::Matrix{T},
                                q_points::Matrix{SVector{2, Float64}},
                                N::I,
                                unit_cell::UnitCell,
                                lattice::Lattice ) where {T<:AbstractFloat, I<:Integer}

Calculates the static structure factor from either density-density or spin-spin correlations.

- `corr::Matrix{T}`: correlation data.
- `q_points::Matrix{SVector{2, Float64}}`: set of momentum points.
- `N::I`: number of lattice sites.
- `unit_cell::UnitCell`: unit cell.
- `lattice::Lattice`: lattice.

"""
function calculate_structure_factor(
    corr::Matrix{T},
    q_points::Matrix{SVector{2, Float64}},
    N::I,
    unit_cell::UnitCell,
    lattice::Lattice
) where {T<:AbstractFloat, I<:Integer}
    nq = length(q_points)
    Sq = zeros(Float64, nq)

    # onsite contribution
    onsite = tr(corr) / N

    # Precompute all Δr_{ij} for i < j
    # Store them in a flat vector aligned with (i,j) pairs
    Δrs = Vector{SVector{2,Float64}}()
    Cs  = Vector{T}()
    sizehint!(Δrs, N*(N-1) ÷ 2)
    sizehint!(Cs,  N*(N-1) ÷ 2)

    @inbounds for i in 1:N
        for j in i+1:N
            Δl = sites_to_displacement(i, j, unit_cell, lattice)
            Δr = displacement_to_vec(Δl, unit_cell.n, unit_cell.n, unit_cell)
            push!(Δrs, Δr)
            push!(Cs, corr[i, j])
        end
    end

    # Now loop over q only
    @inbounds for (k, q) in enumerate(q_points)
        s = 0.0
        for n in eachindex(Δrs)
            s += cos(dot(q, Δrs[n])) * Cs[n]
        end
        Sq[k] = (2.0 / N) * s + onsite
    end

    return Sq
end


"""

    merge_bin_measurements!( datafolder::AbstractString,
                             pID::Int;
                             clear_src::Bool=false )

Merge the HDF5 files `opt_name` and `sim_name` (located inside `datafolder`)
into a single HDF5 file named `out_name` (also inside `datafolder`).

Bins from the sim file are offset by the maximum bin index found in the opt file,
so the result has contiguous, non-overlapping `bin-%d` subgroup names.

- `datafolder::AbstractString`: path to data files.
- `pID::Int`: processor ID/MPI rank
- `clear_src::Bool=false`: optionallyy clear source bin files after merger.

"""
function merge_bin_measurements!(
    datafolder::AbstractString,
    pID::Int;
    clear_src::Bool=false
)
    opt_name = "opt_bin_measurements_rank-$(pID).h5"
    sim_name = "sim_bin_measurements_rank-$(pID).h5"
    out_name = "bin_measurements_rank-$(pID).h5"

    opt_file = joinpath(datafolder, opt_name)
    sim_file = joinpath(datafolder, sim_name)
    out_file = joinpath(datafolder, out_name)

    bin_regex = r"^bin-(\d+)$"
    parse_bin(s::AbstractString) = begin
        m = match(bin_regex, s)
        m === nothing ? nothing : parse(Int, m.captures[1])
    end

    function top_level_groups(fname)
        if !isfile(fname)
            @warn "File not found: $fname"
            return String[]
        end
        groups = String[]
        h5open(fname, "r") do f
            append!(groups, collect(keys(f)))   # use keys() instead of names()
        end
        return groups
    end

    function max_bin_in_file(fname)
        if !isfile(fname)
            return 0
        end
        maxbin = 0
        h5open(fname, "r") do f
            for g in collect(keys(f))
                # g is a top-level name
                gpath = "/" * g
                # attempt to iterate children using keys on the group object
                try
                    grp = f[gpath]
                    for n in collect(keys(grp))
                        b = parse_bin(n)
                        if b !== nothing && b > maxbin
                            maxbin = b
                        end
                    end
                catch
                    # ignore non-group entries
                end
            end
        end
        return maxbin
    end

    # Copy attributes from src object to dst object (datasets only)
    function copy_attrs_obj!(srcobj, dstobj)
        for (k, v) in pairs(attrs(srcobj))
            attrs(dstobj)[k] = v
        end
    end

    # read everything from a given group/bin in src and write into dst with optionally renumbered bin index
    function copy_group_bin!(srcfname, dstf, groupname::String, binname::String, new_bin_index::Int)
        src_gpath = "/" * groupname * "/" * binname
        dst_gpath = "/" * groupname
        dst_binname = @sprintf("bin-%d", new_bin_index)
        dst_full = dst_gpath * "/" * dst_binname

        h5open(srcfname, "r") do sf
            if !haskey(sf, src_gpath)
                return
            end

            src_node = sf[src_gpath]

            # try to list children using keys()
            children = try
                collect(keys(src_node))
            catch
                String[]
            end

            if isempty(children)
                data = read(src_node)
                if haskey(dstf, dst_full)
                    delete!(dstf, dst_full)
                end
                # writing to dst_full will create intermediate groups automatically
                dstf[dst_full] = data
                try
                    copy_attrs_obj!(src_node, dstf[dst_full])
                catch
                end
            else
                # src bin is a group with children; copy each child dataset
                for child in children
                    child_src = src_gpath * "/" * child
                    child_dst = dst_full * "/" * child
                    src_child_node = sf[child_src]
                    grand = try
                        collect(keys(src_child_node))
                    catch
                        String[]
                    end

                    if isempty(grand)
                        data = read(src_child_node)
                        if haskey(dstf, child_dst)
                            delete!(dstf, child_dst)
                        end
                        dstf[child_dst] = data
                        try
                            copy_attrs_obj!(src_child_node, dstf[child_dst])
                        catch
                        end
                    else
                        # nested group - copy leaf datasets one level deep
                        for gchild in grand
                            gchild_src = child_src * "/" * gchild
                            gchild_dst = child_dst * "/" * gchild
                            data = read(sf[gchild_src])
                            if haskey(dstf, gchild_dst)
                                delete!(dstf, gchild_dst)
                            end
                            dstf[gchild_dst] = data
                            try
                                copy_attrs_obj!(sf[gchild_src], dstf[gchild_dst])
                            catch
                            end
                        end
                    end
                end
            end
        end
    end

    max_opt = max_bin_in_file(opt_file)
    # @info "max bin index found in opt file: $max_opt"

    out_mode = isfile(out_file) ? "r+" : "w"
    h5open(out_file, out_mode) do outf
        # copy file attrs from opt if present else sim (best-effort)
        if isfile(opt_file)
            h5open(opt_file, "r") do optf
                for (k,v) in pairs(attrs(optf))
                    if !haskey(attrs(outf), k)
                        attrs(outf)[k] = v
                    end
                end
            end
        elseif isfile(sim_file)
            h5open(sim_file, "r") do simf
                for (k,v) in pairs(attrs(simf))
                    if !haskey(attrs(outf), k)
                        attrs(outf)[k] = v
                    end
                end
            end
        end

        # union of top-level groups
        top_groups = union(top_level_groups(opt_file), top_level_groups(sim_file))

        # copy opt file content (original bin numbers)
        if isfile(opt_file)
            h5open(opt_file, "r") do optf
                for g in collect(keys(optf))
                    gpath = "/" * g
                    # iterate children of the group one-level deep
                    grp = try
                        optf[gpath]
                    catch
                        nothing
                    end
                    if grp === nothing
                        continue
                    end
                    for binname in collect(keys(grp))
                        b = parse_bin(binname)
                        if b === nothing
                            # non-bin entry: copy as dataset (overwrite if present)
                            src_path = gpath * "/" * binname
                            dst_path = src_path
                            data = read(optf[src_path])
                            if haskey(outf, dst_path)
                                delete!(outf, dst_path)
                            end
                            outf[dst_path] = data
                            try
                                copy_attrs_obj!(optf[src_path], outf[dst_path])
                            catch
                            end
                        else
                            copy_group_bin!(opt_file, outf, g, binname, b)
                        end
                    end
                end
            end
        end

        # copy sim file content, offsetting bin indices
        if isfile(sim_file)
            h5open(sim_file, "r") do simf
                for g in collect(keys(simf))
                    gpath = "/" * g
                    grp = try
                        simf[gpath]
                    catch
                        nothing
                    end
                    if grp === nothing
                        continue
                    end
                    for binname in collect(keys(grp))
                        b = parse_bin(binname)
                        if b === nothing
                            src_path = gpath * "/" * binname
                            dst_path = src_path
                            data = read(simf[src_path])
                            if haskey(outf, dst_path)
                                delete!(outf, dst_path)
                            end
                            outf[dst_path] = data
                            try
                                copy_attrs_obj!(simf[src_path], outf[dst_path])
                            catch
                            end
                        else
                            new_index = max_opt + b
                            copy_group_bin!(sim_file, outf, g, binname, new_index)
                        end
                    end
                end
            end
        end
    end
    # @info "Merge complete. Output written to $(abspath(out_file))"
        # optionally delete source files
    if clear_src
        if isfile(opt_file)
            rm(opt_file; force=true)
            # @info "Deleted source file: $opt_file"
        end
        if isfile(sim_file)
            rm(sim_file; force=true)
            # @info "Deleted source file: $sim_file"
        end
    end

    return nothing
end
