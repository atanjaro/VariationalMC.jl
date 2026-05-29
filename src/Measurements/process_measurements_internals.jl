# function to process the statistics for a single pID
function _process_measurements(
    folder::String,
    filename::String,
    pIDs::Vector{Int},
    n_bins::Union{Nothing,Int},
    rm_binned_data::Bool,
    process_global_measurements::Bool,
    process_local_measurements::Bool,
    process_all_equal_time_measurements::Bool
)

    # open output HDF5 stats file
    h5_stats_filename = joinpath(folder, filename)
    H5StatsFile = h5open(h5_stats_filename, "w")

    # get directory containing bin folder
    bin_folder = joinpath(folder, "sim_bins")

    # HDF5 bin file names
    h5_bin_filenames = [joinpath(bin_folder, "bins_pID-$(pID).h5") for pID in pIDs]

    # open all the input HDF5 bin files
    H5BinFiles = [h5open(file, "r") for file in h5_bin_filenames]

    # number of HDF5 files containing binned data
    N_pIDs = length(pIDs)

    # get the number of bins per pID
    n_data_bins = read_attribute(H5BinFiles[1], "N_BINS")
    n_bins = isnothing(n_bins) ? n_data_bins : n_bins
    @assert (n_data_bins % n_bins) == 0
    @assert N_pIDs * n_bins > 1 "The total number of data bins is one or smaller, and therefore measurement errors cannot be estimated."

    # calculate total number bins across all pIDs
    N_bins = n_bins * N_pIDs

    # get all equal-time correlations if necessary
    standard_equal_time = nothing
    if process_all_equal_time_measurements
        standard_equal_time = keys(H5BinFiles[1]["CORRELATIONS/EQUAL-TIME"])
    end

    # initialize the HDF5 stats file
    allocate_stats_file!(
        H5StatsFile, H5BinFiles[1],
        process_global_measurements, process_local_measurements,
        process_all_equal_time_measurements, standard_equal_time
    )

    # record pIDs associated with stats HDF5 files
    attributes(H5StatsFile)["PIDS"] = pIDs

    # get system size and numbers of parameters
    N_orbitals = read_attribute(H5BinFiles[1], "N_ORBITALS")
    n_params = read_attribute(H5BinFiles[1], "NUM_PARAMS")
    n_opt_params = read_attribute(H5BinFiles[1], "NUM_OPT_PARAMS")

    # record metadata about stats to computes
    attributes(H5StatsFile)["NUM_PARAMS"] = n_params
    attributes(H5StatsFile)["NUM_OPT_PARAMS"] = n_opt_params
    attributes(H5StatsFile)["N_ORBITALS"] = N_orbitals
    attributes(H5StatsFile)["N_BINS"] = N_bins

    # preallocate arrays for jackknife using first global measurement to determine size/type
    first_global_key = first(keys(H5BinFiles[1]["GLOBAL"]))
    binned_vals_init = vcat((rebin(read(H5BinFiles[n]["GLOBAL"][first_global_key]), n_bins) for n in 1:N_pIDs)...)
    jackknife_sample_means = (similar(binned_vals_init),)
    jackknife_g = similar(binned_vals_init)

    # type of reported mean and standard deviations
    T_mean = eltype(binned_vals_init)
    T_std = real(T_mean)

    # get global measurement stats group
    Global_Out = H5StatsFile["GLOBAL"]

    # iterate over global measurements
    for key in keys(Global_Out)
        # read in global measurement
        binned_vals = vcat((rebin(read(H5BinFiles[n]["GLOBAL"][key]), n_bins) for n in 1:N_pIDs)...)
        binned_vals = rebin(binned_vals, N_bins)
        # calculate mean and error using jackknife
        average, stdev = jackknife(
            identity, binned_vals;
            jackknife_samples = jackknife_sample_means, jackknife_g
        )
        # get global measurement group
        Global_Measurement = Global_Out[key]
        # record final mean and standard deviation
        Global_Measurement["MEAN"] = average
        Global_Measurement["STD"] = stdev
    end

    # iterate over local measurements
    Local_Out = H5StatsFile["LOCAL"]

    for key in keys(Local_Out)

        Measurement_Out = Local_Out[key]
        Measurement_In  = H5BinFiles[1]["LOCAL"][key]

        if Measurement_In isa HDF5.Dataset

            # gather rebinned data from all pIDs
            data = vcat((
                rebin(read(H5BinFiles[n]["LOCAL"][key]), n_bins)
                for n in 1:N_pIDs
            )...)

            # number of IDs/sites
            N_IDs = size(data, 2)

            average = zeros(T_mean, N_IDs)
            stdev   = zeros(T_std, N_IDs)

            for ID in 1:N_IDs

                average[ID], stdev[ID] = jackknife(
                    identity,
                    view(data, :, ID);
                    jackknife_samples = jackknife_sample_means,
                    jackknife_g
                )
            end

            Measurement_Out["MEAN"] = average
            Measurement_Out["STD"]  = stdev

        elseif Measurement_In isa HDF5.Group

            for orbital_key in keys(Measurement_In)

                Orbital_Out = create_group(Measurement_Out, orbital_key)

                # gather rebinned site data
                data = vcat((
                    rebin(
                        read(H5BinFiles[n]["LOCAL"][key][orbital_key]),
                        n_bins
                    )
                    for n in 1:N_pIDs
                )...)

                N_IDs = size(data, 2)

                average = zeros(T_mean, N_IDs)
                stdev   = zeros(T_std, N_IDs)

                for ID in 1:N_IDs

                    average[ID], stdev[ID] = jackknife(
                        identity,
                        view(data, :, ID);
                        jackknife_samples = jackknife_sample_means,
                        jackknife_g
                    )
                end

                Orbital_Out["MEAN"] = average
                Orbital_Out["STD"]  = stdev

                # preserve metadata
                for attr in keys(attributes(Measurement_In[orbital_key]))

                    attributes(Orbital_Out)[attr] =
                        read_attribute(
                            Measurement_In[orbital_key],
                            attr
                        )
                end
            end
        end
    end

    # process standard equal-time correlation measurements
    if process_all_equal_time_measurements
        process_correlations!(
            H5StatsFile, H5BinFiles, "EQUAL-TIME",
            jackknife_sample_means, jackknife_g
        )
    end

    # close all the input HDF5 bin files
    close.(H5BinFiles)

    # close the output HDF5 stats file
    close(H5StatsFile)

    # delete HDF5 files containing binned data
    if rm_binned_data
        rm_bins(folder)
    end

    return h5_stats_filename
end


# process a certain type of correlation
function process_correlations!(
    H5StatsFile::HDF5.File,
    H5BinFiles::Vector{HDF5.File},
    correlation_type::String,
    jackknife_samples::Tuple{Vector{Complex{T}}}, 
    jackknife_g::Vector{Complex{T}}
) where {T<:AbstractFloat}

    # get output correlations
    Correlations_Out = H5StatsFile["CORRELATIONS"][correlation_type]

    # get the number of pIDs
    N_pIDs = length(H5BinFiles)

    # total number of bins
    N_bins = sum(size(H5BinFiles[n]["CORRELATIONS"][correlation_type][first(keys(H5BinFiles[n]["CORRELATIONS"][correlation_type]))]["POSITION"], 1) for n in 1:N_pIDs)

    # number of bins per pID
    n_bins = N_bins ÷ N_pIDs

    for key in keys(Correlations_Out)

        # get the output group for the position space correlation
        Position_Out = Correlations_Out[key]["POSITION"]

        # get the input binned datasets for the position correlations
        Positions_In = tuple((H5BinFiles[n]["CORRELATIONS"][correlation_type][key]["POSITION"] for n in 1:N_pIDs)...)

        # get the dimensionality of the correlation data.
        # first dimension is cut off as that one corresponds to the bins.
        dims = size(first(Positions_In))[2:end]

        # allocate arrays to contain stats
        average = zeros(Complex{T}, dims)
        stdev = zeros(T, dims)

        # iterate over all correlation displacement vectors
        for c in CartesianIndices(dims)
            # concatenate rebinned data for each pID together
            data = vcat((
                # rebin the data associated with single pID
                rebin(
                    # read in the binned data associated with each pID
                    Positions_In[n][:,c.I...],
                    n_bins
                )
                # iterate over pIDs
                for n in 1:N_pIDs
            )...)
            # calculate the stats
            average[c], stdev[c] = jackknife(
            identity, data;
            jackknife_samples, jackknife_g  # ← renamed
        )
        end

        # record the final stats
        Position_Out["MEAN"] = average
        Position_Out["STD"] = stdev

        # get the output group for the momentum space correlation
        Momentum_Out = Correlations_Out[key]["MOMENTUM"]

        # get the input binned datasets for the momentum correlations
        Momentum_In = tuple((H5BinFiles[n]["CORRELATIONS"][correlation_type][key]["MOMENTUM"] for n in 1:N_pIDs)...)

        # iterate over all scattering momentum vectors
        for c in CartesianIndices(dims)
            # concatenate rebinned data for each pID together
            data = vcat((
                # rebin the data associated with single pID
                rebin(
                    # read in the binned data associated with each pID
                    Momentum_In[n][:,c.I...],
                    n_bins
                )
                # iterate over pIDs
                for n in 1:N_pIDs
            )...)
            # calculate the stats
            average[c], stdev[c] = jackknife(
                identity, data;
                jackknife_samples, jackknife_g  
            )
        end

        # record the final stats
        Momentum_Out["MEAN"] = average
        Momentum_Out["STD"] = stdev
    end

    return nothing
end


# allocate HDF5 stats file
function allocate_stats_file!(
    H5StatsFile::HDF5.File,
    H5BinFile::HDF5.File,
    process_global_measurements::Bool,
    process_local_measurements::Bool,
    process_all_equal_time_measurements::Bool,
    standard_equal_time::Union{Vector{String}, Nothing}
)

    # initialize global measurements group
    Global_In = H5BinFile["GLOBAL"]
    Global_Out = create_group(H5StatsFile, "GLOBAL")
    if process_global_measurements
        for key in keys(Global_In)
            create_group(Global_Out, key)
        end
    end

    # initialize local measurements group
    Local_In = H5BinFile["LOCAL"]
    Local_Out = create_group(H5StatsFile, "LOCAL")
    if process_local_measurements
        for key in keys(Local_In)
            Local_Measurement_In = Local_In[key]
            Local_Measurement_Out = create_group(Local_Out, key)
            attributes(Local_Measurement_Out)["ID_TYPE"] = read_attribute(Local_Measurement_In, "ID_TYPE")
        end
    end

    # create groups to contain correlation measurements
    Correlations_In = create_group(H5StatsFile, "CORRELATIONS")
    create_group(Correlations_In, "EQUAL-TIME")

    # initialize standard equal-time correlation measurements
    if process_all_equal_time_measurements
        init_correlation_type_group!(
            H5StatsFile["CORRELATIONS/EQUAL-TIME"],
            H5BinFile["CORRELATIONS/EQUAL-TIME"],
            standard_equal_time
        )
    end

    return nothing
end

# initialize a specific type of correlation measurement
function init_correlation_type_group!(
    CorrelationTypeGroup_Out::HDF5.Group,
    CorrelationTypeGroup_In::HDF5.Group,
    correlations::Vector{String}
)

    for key in correlations
        Correlation_In = CorrelationTypeGroup_In[key]
        Correlation_Out = create_group(CorrelationTypeGroup_Out, key)
        if haskey(attrs(Correlation_In), "ID_TYPE")
            attributes(Correlation_Out)["ID_TYPE"] = read_attribute(Correlation_In, "ID_TYPE")
            attributes(Correlation_Out)["ID_PAIRS"] = read_attribute(Correlation_In, "ID_PAIRS")
        end
        Position = create_group(Correlation_Out, "POSITION")
        Momentum = create_group(Correlation_Out, "MOMENTUM")
        attributes(Position)["DIM_LABELS"] = read_attribute(Correlation_In["POSITION"], "DIM_LABELS")[2:end]
        attributes(Momentum)["DIM_LABELS"] = read_attribute(Correlation_In["MOMENTUM"], "DIM_LABELS")[2:end]
    end

    return nothing
end

