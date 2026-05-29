# function to process the statistics for a single pID
function _process_measurements(
    comm::MPI.Comm,
    folder::String,
    filename::String,
    pIDs::Vector{Int},
    n_bins::Union{Nothing,Int},
    rm_binned_data::Bool,
    process_global_measurements::Bool,
    process_local_measurements::Bool,
    process_all_equal_time_measurements::Bool,
    standard_equal_time::Vector{String}
)
    # make sure number of pIDs matches number of MPI processes
    @assert length(pIDs) == MPI.Comm_size(comm)

    # number of pID processes
    N_pIDs = length(pIDs)

    # get MPI rank
    rank = MPI.Comm_rank(comm)

    # if is the root process
    isroot = iszero(rank)

    # get the process ID
    pID = pIDs[rank+1]

    # open output HDF5 stats file if root process
    h5_stats_filename = joinpath(folder, filename)
    H5StatsFile = isroot ? h5open(h5_stats_filename, "w") : nothing

    # directory containing binned data
    bin_folder = joinpath(folder, "sim_bins")

    # open bin file
    h5_bin_filename = joinpath(bin_folder, "bins_pID-$(pID).h5")
    H5BinFile = h5open(h5_bin_filename, "r")

    # get the number of bins per pID
    n_data_bins = read_attribute(H5BinFile, "N_BINS")
    n_bins = isnothing(n_bins) ? n_data_bins : n_bins
    @assert (n_data_bins % n_bins) == 0
    @assert N_pIDs * n_bins > 1 "The total number of data bins is one or smaller, and therefore measurement errors cannot be estimated."

    # calculate total number bins across all pIDs
    N_bins = n_bins * N_pIDs

    # get all equal-time correlations if necessary
    standard_equal_time = nothing
    if process_all_equal_time_measurements
        standard_equal_time = keys(H5BinFile["CORRELATIONS/EQUAL-TIME"])
    end

    # allocate HDF5 stats file if root process
    if isroot
        allocate_stats_file!(
            H5StatsFile, H5BinFile,
            process_global_measurements, process_local_measurements,
            process_all_equal_time_measurements, standard_equal_time
        )

        # record pIDs associated with stats HDF5 files
        attributes(H5StatsFile)["PIDS"] = pIDs
    end
    MPI.Barrier(comm)

    # get system size and numbers of parameters
    N_orbitals = read_attribute(H5BinFile, "N_ORBITALS")
    n_params = read_attribute(H5BinFile, "NUM_PARAMS")
    n_opt_params = read_attribute(H5BinFile, "NUM_OPT_PARAMS")

    # record metadata about stats
    if isroot
        attributes(H5StatsFile)["NUM_PARAMS"] = n_params
        attributes(H5StatsFile)["NUM_OPT_PARAMS"] = n_opt_params
        attributes(H5StatsFile)["N_ORBITALS"] = N_orbitals
        attributes(H5StatsFile)["N_BINS"] = N_bins
    end

    # if root process, preallocate jackknife arrays using first global measurement
    if isroot
        first_global_key = first(keys(H5BinFile["GLOBAL"]))
        binned_vals_init = rebin(read(H5BinFile["GLOBAL"][first_global_key]), n_bins)
        # note: binned_vals_init here is just one pID's worth; full size is N_bins after gather
        jackknife_sample_means = (zeros(eltype(binned_vals_init), N_bins),)
        jackknife_g = zeros(eltype(binned_vals_init), N_bins)

        # type of reported mean and standard deviations
        T_mean = eltype(binned_vals_init)
        T_std = real(T_mean)

        # get global measurement stats group
        Global_Out = H5StatsFile["GLOBAL"]

    else
        jackknife_sample_means, jackknife_g = nothing, nothing
    end

    # get the input global measurements group
    Global_In = H5BinFile["GLOBAL"]

    # iterate over global measurements
    for key in keys(Global_In)
        # read in the global measurement bins for the current pID and rebin it
        binned_vals = MPI.gather(rebin(read(Global_In[key]), n_bins), comm)
        # if root process
        if isroot
            # concatenate all the data together
            binned_vals = vcat(binned_vals...)
            # calculate mean and error using jackknife
            avg, stdev = jackknife(
                identity, binned_vals;
                jackknife_samples = jackknife_sample_means, jackknife_g
            )
            # record the global measurement stats
            Global_Out[key]["MEAN"] = avg
            Global_Out[key]["STD"] = stdev
        end
    end

    # get the output local measurement group
    Local_Out = isroot ? H5StatsFile["LOCAL"] : nothing

    # get the input local measurements group
    Local_In = H5BinFile["LOCAL"]

    # iterate over local measurements
    for key in keys(Local_In)
        # get the input measurement bins dataset
        Measurement_In = Local_In[key]
        # gather rebinned measurement data
        data = MPI.gather(rebin(read(Measurement_In), n_bins), comm)
        # if root process
        if isroot
            # concatenate all the gathered data across pIDs
            data = vcat(data...)
            # number of IDs associated with local measurement
            N_IDs = size(data, 2)
            # allocate array to contain stats
            average = zeros(T_mean, N_IDs)
            stdev = zeros(T_std, N_IDs)
            # iterate over IDs associated with local measurement
            for ID in 1:N_IDs
                # calculate mean and error for current local measurement ID
                average[ID], stdev[ID] = jackknife(
                    identity, view(data, :, ID);
                    jackknife_samples = jackknife_sample_means, jackknife_g
                )
            end
            # record the local measurement stats
            Measurement_Out = Local_Out[key]
            Measurement_Out["MEAN"] = average
            Measurement_Out["STD"] = stdev
        end
    end

    # process standard equal-time correlations
    if process_all_equal_time_measurements
        process_correlations!(
            comm, H5StatsFile, H5BinFile,
            "EQUAL-TIME",
            standard_equal_time, n_bins,
            jackknife_sample_means, jackknife_g
        )
    end

    # close the HDF5 bin file
    close(H5BinFile)

    # close HDF5 stats file
    if isroot
        close(H5StatsFile)
    end

    # delete HDF5 files containing binned data
    if rm_binned_data
        rm_bins(comm, folder)
    end

    return h5_stats_filename
end


function process_correlations!(
    comm::MPI.Comm,
    H5StatsFile::Union{HDF5.File,Nothing},
    H5BinFile::HDF5.File,
    correlation_group::String,
    correlation_type::String,
    correlations::Vector{String},
    n_bins::Int,
    jackknife_sample_means,
    jackknife_g
)

    # get current MPI rank
    rank = MPI.Comm_rank(comm)
    # get input and output groups containing correlation type
    Correlations_In = H5BinFile["CORRELATIONS"][correlation_group][correlation_type]
    Correlations_Out = isnothing(H5StatsFile) ? nothing : H5StatsFile["CORRELATIONS"][correlation_group][correlation_type]
    # iterate over correlations
    for key in correlations
        # input correlations
        Correlation_In = Correlations_In[key]
        Position_In = Correlation_In["POSITION"]
        Momentum_In = Correlation_In["MOMENTUM"]
        # get dimensions of output correlation group
        dims = size(Position_In)
        # initialize array to contain stats
        T = iszero(rank) ? real(eltype(first(jackknife_sample_means))) : nothing
        average = iszero(rank) ? zeros(Complex{T}, dims[2:end]) : nothing
        stdev = iszero(rank) ? zeros(T, dims[2:end]) : nothing
        # number of spatial dimensions of lattice
        D = length(dims) - 1 - Int("STANDARD" == correlation_group) - Int("TIME-DISPLACED" == correlation_type)
        # get extent of lattice size
        L = dims[2:D+1]
        # get index range associated with lattice size
        Ls = tuple((1:l for l in L)...)
        # iterate over imaginary-time slice and ID pairs when relevant
        for n in CartesianIndices(dims[D+2:end])
            # read in and rebin the position data
            data = rebin(Position_In[:,Ls...,n.I...], n_bins)
            # gather all the data into the root process
            data = MPI.gather(data, comm)
            # if root process
            if iszero(rank)
                # concatenate the bins from each pID together
                data = vcat(data...)
                # iterate over displacement vectors
                for c in CartesianIndices(L)
                    # calculate the stats
                    average[c.I..., n.I...], stdev[c.I..., n.I...] = jackknife(
                        identity, view(data, :, c.I...);
                        jackknife_sample_means, jackknife_g
                    )
                end
            end
        end
        # record the position space correlation stats
        if iszero(rank)
            Correlations_Out[key]["POSITION"]["MEAN"] = average
            Correlations_Out[key]["POSITION"]["STD"] = stdev
        end
        # iterate over imaginary-time slice and ID pairs if relevant
        for n in CartesianIndices(dims[D+2:end])
            # read in and rebin the momentum data
            data = rebin(Momentum_In[:,Ls...,n.I...], n_bins)
            # gather all the data into the root process
            data = MPI.gather(data, comm)
            # if root process
            if iszero(rank)
                # concatenate the bins from each pID together
                data = vcat(data...)
                # iterate over displacement vectors
                for c in CartesianIndices(L)
                    # calculate the stats
                    average[c.I..., n.I...], stdev[c.I..., n.I...] = jackknife(
                        identity, view(data, :, c.I...);
                        jackknife_sample_means, jackknife_g
                    )
                end
            end
        end
        # record the momentum space correlation stats
        if iszero(rank)
            Correlations_Out[key]["MOMENTUM"]["MEAN"] = average
            Correlations_Out[key]["MOMENTUM"]["STD"] = stdev
        end
    end

    return nothing
end