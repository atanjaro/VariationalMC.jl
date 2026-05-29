@doc raw"""

    merge_bins(
        simulation_info::SimulationInfo;
        opt::Bool = false
    )

Merge the separate HDF5 files containing the binned measurements into a single HDF5 file.
This is true even if the HDF5 "files" containing the binned data were [held in memory](https://juliaio.github.io/HDF5.jl/stable/#In-memory-HDF5-files)
during the simulation (`simulation_info.write_bins_concurrent = false`) instead of being actively written to file
during the simulation (`simulation_info.write_bins_concurrent = true`). Set `opt = true` if merging binds from the optimization step.

"""
function merge_bins(
    simulation_info::SimulationInfo;
    opt::Bool = false
)
    (; opt_bin_files, sim_bin_files, write_bins_concurrent, datafolder, pID) = simulation_info
    bin_files = opt ? opt_bin_files : sim_bin_files

    # get the directory containing the HDF5 bin files
    binfolder = opt ? joinpath(datafolder, "opt_bins") : joinpath(datafolder, "sim_bins")

    # construct filename for the merged HDF5 file
    filename = joinpath(binfolder, "bins_pID-$(pID).h5")

    # check HDF5 files still needs to be merged by checking if merged HDF5 already exists
    if !isfile(filename)

        # count the number of bins
        N_bins = length(bin_files)

        # open new HDF5 file to contain all the binned data
        h5open(filename, "w") do fout

            # Record the process ID
            attributes(fout)["PID"] = pID

            # Record the number of bins
            attributes(fout)["N_BINS"] = N_bins

            # open first HDF5 bin file
            if bin_files[1] isa AbstractVector{UInt8}
                fin = h5open(bin_files[1], "r")
            else
                fin = h5open(String(bin_files[1]), "r")
            end

            # initialize/allocate HDF5 file to contain all binned data
            init_hdf5_bins_file(fout, fin, N_bins, opt)

            # copy contents of first bin file over
            copyto_hdf5_bin(fout, fin, 1, opt)

            # close first HDF5 bin file
            close(fin)

            # iterate over remaining bin files
            for bin in 2:N_bins

                # open HDF5 bin file
                if bin_files[bin] isa AbstractVector{UInt8}
                    fin = h5open(bin_files[bin], "r")
                else
                    fin = h5open(String(bin_files[bin]), "r")
                end

                # copy contents of bin file over
                copyto_hdf5_bin(fout, fin, bin, opt)

                # close HDF5 bin file
                close(fin)
            end
        end

        # if HDF5 bin files were written during the simulation
        if write_bins_concurrent

            # get the directory containing the HDF5 bin files for current pID
            pID_binfolder = joinpath(binfolder, "pID-$pID")

            # delete pID bin folder
            rm(pID_binfolder, recursive = true)
        end
    end

    return nothing
end


function init_hdf5_bins_file(
    fout::HDF5.File,
    fin::HDF5.File,
    N_bins::Int,
    opt::Bool
)

    # record some standard meta data
    attributes(fout)["NUM_PARAMS"] = read_attribute(fin, "NUM_PARAMS")
    attributes(fout)["NUM_OPT_PARAMS"] = read_attribute(fin, "NUM_OPT_PARAMS")
    attributes(fout)["N_ORBITALS"] = read_attribute(fin, "N_ORBITALS")

    # initialize global measurements
    Global = create_group(fout, "GLOBAL")
    Global_Bin = fin["GLOBAL"]
    for key in keys(Global_Bin)
        write_dataset(Global, key, zeros(eltype(Global_Bin[key]), N_bins))
        Global_Measurement = Global[key]
        attributes(Global_Measurement)["DIM_LABELS"] = ["BIN"]
    end

    # initialize local measurements
    Local = create_group(fout, "LOCAL")
    Local_Bin = fin["LOCAL"]
    for key in keys(Local_Bin)
        item = Local_Bin[key]
        if item isa HDF5.Dataset
            # simple dataset, handle as before
            write_dataset(Local, key, zeros(eltype(item), (N_bins, size(item)...)))
            Local_Measurement = Local[key]
            attributes(Local_Measurement)["DIM_LABELS"] = ["BIN", "ID"]
            attributes(Local_Measurement)["ID_TYPE"] = LOCAL_MEASUREMENTS[key]
        else
            # parent group indexes orbitals
            Orbital_Group = create_group(Local, key)

            # metadata for orbital indexing
            attributes(Orbital_Group)["DIM_LABELS"] = ["ORBITAL"]
            attributes(Orbital_Group)["ID_TYPE"] = "ORBITAL_ID"

            for orbital_key in keys(item)

                orbital_data = item[orbital_key]

                write_dataset(
                    Orbital_Group,
                    orbital_key,
                    zeros(eltype(orbital_data), (N_bins, size(orbital_data)...))
                )

                Orbital_Measurement = Orbital_Group[orbital_key]

                # metadata for actual site-resolved data
                attributes(Orbital_Measurement)["DIM_LABELS"] = ["BIN", "SITE"]
                attributes(Orbital_Measurement)["ID_TYPE"] = "SITE_ID"
                attributes(Orbital_Measurement)["ORBITAL_ID"] = parse(Int, orbital_key)

            end
        end
    end

    # initialize optimization measurements
    if opt
        Optimization = create_group(fout, "OPTIMIZATION")
        Optimization_Bin = fin["OPTIMIZATION"]
        write_dataset(Optimization, "parameters", zeros(eltype(Optimization_Bin["parameters"]), (N_bins, size(Optimization_Bin["parameters"])...)))
        Optimization_Measurement = Optimization["parameters"]
        attributes(Optimization_Measurement)["DIM_LABELS"] = ["BIN"]
    end

    # initialize group structure of correlations
    Correlations = create_group(fout, "CORRELATIONS")
    EqualTime = create_group(Correlations, "EQUAL-TIME")

    # initialize standard equal-time correlation measurements
    EqualTime_Bin = fin["CORRELATIONS/EQUAL-TIME"]
    for key in keys(EqualTime_Bin)
        Correlation_Bin = EqualTime_Bin[key]
        Correlation = create_group(EqualTime, key)
        attributes(Correlation)["ID_PAIRS"] = read_attribute(Correlation_Bin, "ID_PAIRS")
        attributes(Correlation)["ID_TYPE"] = read_attribute(Correlation_Bin, "ID_TYPE")
        Position_Bin = Correlation_Bin["POSITION"]
        write_dataset(Correlation, "POSITION", zeros(eltype(Position_Bin), (N_bins, size(Position_Bin)...)))
        Position = Correlation["POSITION"]
        Momentum_Bin = Correlation_Bin["MOMENTUM"]
        write_dataset(Correlation, "MOMENTUM", zeros(eltype(Momentum_Bin), (N_bins, size(Momentum_Bin)...)))
        Momentum = Correlation["MOMENTUM"]
        D = ndims(Position_Bin) - 1 # number of spatial dimensions
        attributes(Position)["DIM_LABELS"] = ["BIN", position_column_labels(D)..., "ID_PAIR"]
        attributes(Momentum)["DIM_LABELS"] = ["BIN", momentum_column_labels(D)..., "ID_PAIR"]
    end

    return nothing
end

# copy contents of single bin HDF5 file to HDF5 file containing all binned data
function copyto_hdf5_bin(
    fout::HDF5.File,
    fin::HDF5.File,
    bin::Int,
    opt::Bool
)

    # copy global measurements over
    Global_out = fout["GLOBAL"]
    Global_in = fin["GLOBAL"]

    for key in keys(Global_out)
        Global_out[key][bin] = read(Global_in[key])
    end

    # copy local measurements over
    Local_out = fout["LOCAL"]
    Local_in = fin["LOCAL"]

    for key in keys(Local_out)

        item_out = Local_out[key]
        item_in = Local_in[key]

        # standard local measurement dataset
        if item_in isa HDF5.Dataset

            inds = (bin, ntuple(_ -> Colon(), ndims(item_in))...)
            item_out[inds...] = read(item_in)

        # subgroup containing orbital-resolved datasets
        elseif item_in isa HDF5.Group

            for orbital_key in keys(item_in)

                orbital_out = item_out[orbital_key]
                orbital_in = item_in[orbital_key]

                inds = (bin, ntuple(_ -> Colon(), ndims(orbital_in))...)
                orbital_out[inds...] = read(orbital_in)

            end
        end
    end

    # copy optimization measurements over
    if opt
        Optimization_out = fout["OPTIMIZATION"]
        Optimization_in = fin["OPTIMIZATION"]

        inds = (bin, ntuple(_ -> Colon(), ndims(Optimization_in["parameters"]))...)
        Optimization_out["parameters"][inds...] = read(Optimization_in["parameters"])
    end

    # copy equal-time correlation measurements
    copy_correlation_bins(
        fout["CORRELATIONS/EQUAL-TIME"],
        fin["CORRELATIONS/EQUAL-TIME"],
        bin
    )

    return nothing
end


# copy binned correlation table over
function copy_correlation_bins(
    Correlations_out::HDF5.Group,
    Correlations_in::HDF5.Group,
    bin::Int
)

    # iterate over correlation measurements
    for key in keys(Correlations_out)

        # copy the specified bin over
        Cout = Correlations_out[key]
        Cin = Correlations_in[key]
        Position_out = Cout["POSITION"]
        Position_in = Cin["POSITION"]
        Momentum_out = Cout["MOMENTUM"]
        Momentum_in = Cin["MOMENTUM"]
        slice = tuple((1:d for d in size(Position_in))...)
        Position_out[bin, slice...] = read(Position_in)
        Momentum_out[bin, slice...] = read(Momentum_in)
    end

    return nothing
end

# axis/column labels for displacement vectors in position and momentum space
position_column_labels(D::Int) = ((@sprintf("R_%d", d) for d in 1:D)...,)
momentum_column_labels(D::Int) = ((@sprintf("K_%d", d) for d in 1:D)...,)


@doc raw"""
    rm_bins(
        comm::MPI.Comm,
        datafolder::String
    )

    rm_bins(
        datafolder::String
    )

Delete the binned data stored in the `datafolder` directory.
This function essentially deletes the directory `datafolder/bins` and its contents.
"""
function rm_bins(
    comm::MPI.Comm,
    datafolder::String
)
    
    pID = MPI.Comm_rank(comm)
    MPI.Barrier(comm)
    rm(joinpath(datafolder, "bins", "pID-$(pID)"), recursive = true, force = true)
    rm(joinpath(datafolder, "bins", "bins_pID-$(pID).h5"), recursive = true, force = true)
    MPI.Barrier(comm)
    if iszero(pID)
        rm(joinpath(datafolder, "bins"), recursive = true, force = true)
    end
    MPI.Barrier(comm)

    return nothing
end

function rm_bins(
    datafolder::String
)
    rm(joinpath(datafolder,"bins"), recursive = true, force = true)
    return nothing
end


@doc raw"""

    merge_opt_sim_bins(
        # ARGUMENTS
        simulation_info::SimulationInfo
    )

Merge the two separately merged HDF5 bin files — one from the optimization step (`opt_bins/bins_pID-\$(pID).h5`) and one from the simulation step
(`sim_bins/bins_pID-\$(pID).h5`) — into a single file named `vmc_bins_pID-\$(pID).h5` located in `datafolder`.

"""
function merge_opt_sim_bins(
    # ARGUMENTS
    simulation_info::SimulationInfo
)
    (; datafolder, pID) = simulation_info

    # paths to the two already-merged bin files
    opt_filename = joinpath(datafolder, "opt_bins", "bins_pID-$(pID).h5")
    sim_filename = joinpath(datafolder, "sim_bins", "bins_pID-$(pID).h5")

    # path for the combined output file
    out_filename = joinpath(datafolder, "vmc_bins_pID-$(pID).h5")

    # skip if already merged
    isfile(out_filename) && return nothing

    h5open(opt_filename, "r") do f_opt
        h5open(sim_filename, "r") do f_sim
            h5open(out_filename, "w") do fout

                # --- top-level attributes ---
                N_opt_bins = read_attribute(f_opt, "N_BINS")
                N_sim_bins = read_attribute(f_sim, "N_BINS")
                attributes(fout)["PID"]          = pID
                attributes(fout)["N_OPT_BINS"]   = N_opt_bins
                attributes(fout)["N_SIM_BINS"]   = N_sim_bins

                # --- create top-level groups ---
                OPT = create_group(fout, "OPT")
                SIM = create_group(fout, "SIM")

                # copy opt metadata attributes into OPT group
                for attr in ("NUM_PARAMS", "NUM_OPT_PARAMS", "N_ORBITALS")
                    attributes(OPT)[attr] = read_attribute(f_opt, attr)
                end

                # copy sim metadata attributes into SIM group
                for attr in ("NUM_PARAMS", "NUM_OPT_PARAMS", "N_ORBITALS")
                    attributes(SIM)[attr] = read_attribute(f_sim, attr)
                end

                # recursively copy opt and sim bin file contents into their groups
                copy_hdf5_group(f_opt, OPT)
                copy_hdf5_group(f_sim, SIM)
            end
        end
    end

    return nothing
end

# copies all datasets and groups from `src` into `dst`, preserving attributes and dimension labels.
function copy_hdf5_group(src::Union{HDF5.File, HDF5.Group}, dst::HDF5.Group)
    for key in keys(src)
        obj = src[key]
        if obj isa HDF5.Group
            grp = create_group(dst, key)
            # copy group-level attributes
            for attr in keys(attributes(obj))
                attributes(grp)[attr] = read_attribute(obj, attr)
            end
            copy_hdf5_group(obj, grp)
        elseif obj isa HDF5.Dataset
            data = read(obj)
            write_dataset(dst, key, data)
            ds = dst[key]
            # copy dataset-level attributes (e.g. DIM_LABELS, ID_TYPE, ID_PAIRS)
            for attr in keys(attributes(obj))
                attributes(ds)[attr] = read_attribute(obj, attr)
            end
        end
    end
    return nothing
end