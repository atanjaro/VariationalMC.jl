@doc raw"""
    SimulationInfo

Contains identification information about simulation, including the location data is written to,
the simulation ID, and MPI process ID, and whether this simulation started a new simulation or resumed
a previous simulation.

# Fields

- `filepath::String`: File path to where data folder lives.
- `datafolder_prefix`: Prefix for the data folder name.
- `datafolder_name::String`: The data folder name, given by `$(datafolder_prefix)_$(sID)`.
- `datafolder::String`: The data folder, including filepath, given by `joinpath(filepath, datafolder_name)`.
- `pID::Int`: MPI process ID, defaults to 0 if MPI not being used.
- `sID::Int`: Simulation ID.
- `bin_files::Vector{Vector{UInt8}}`: Represents the HDF5 files containing the binned data.
- `resuming::Bool`: Whether current simulation is resuming a previous simulation (`true`) or starting a new one (`false`).
- `variationalmc_version::VersionNumber`: Version of [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl) used in simulation.

"""
mutable struct SimulationInfo

    # filepath to where data directory will live
    filepath::String

    # prefix of data directory name
    datafolder_prefix::String

    # data directory name
    datafolder_name::String

    # data directory including filepath
    datafolder::String

    # process ID number (for MPI)
    pID::Int

    # simulation ID number
    sID::Int

    # if binned data will be held in memory or written to file during the simulation
    write_bins_concurrent::Bool

    # HDF5 optimization bin files
    opt_bin_files::Vector{Vector{UInt8}}

    # HDF5 simulation bin files
    sim_bin_files::Vector{Vector{UInt8}}

    # whether previous simulation is being resumed or a new one is begininning
    resuming::Bool

    # version of package (= VARIATIONALMC_VERSION)
    varmc_version::VersionNumber
end


@doc raw"""

    SimulationInfo(;
        # KEYWORD ARGUMENTS
        datafolder_prefix::String,
        filepath::String = ".",
        sID::Int=0,
        pID::Int=0
    )

Initialize and return in instance of the type [`SimulationInfo`](@ref).
"""
function SimulationInfo(;
    # KEYWORD ARGUMENTS
    datafolder_prefix::String,
    filepath::String = ".",
    write_bins_concurrent::Bool = true,
    sID::Int=0,
    pID::Int=0
)
    # initialize data folder names
    datafolder_name = @sprintf "%s-%d" datafolder_prefix sID
    datafolder = joinpath(filepath, datafolder_name)
    complete_datafolder_name = @sprintf "complete_%s-%d" datafolder_prefix sID
    complete_datafolder = joinpath(filepath, complete_datafolder_name)

    # if null data folder id given, determine data name and id
    if sID==0
        while isdir(datafolder) || isdir(complete_datafolder) || sID==0
            sID += 1
            datafolder_name = @sprintf "%s-%d" datafolder_prefix sID
            datafolder = joinpath(filepath, datafolder_name)
            complete_datafolder_name = @sprintf "complete_%s-%d" datafolder_prefix sID
            complete_datafolder = joinpath(filepath, complete_datafolder_name)
        end
    end

    # if directory already exists then must be resuming simulation
    resuming = isdir(datafolder)

    # initialize bin files
    opt_bin_files = Vector{UInt8}[]
    sim_bin_files = Vector{UInt8}[]

    return SimulationInfo(filepath, datafolder_prefix, datafolder_name, datafolder, pID, sID, write_bins_concurrent, opt_bin_files, sim_bin_files, resuming, VARIATIONALMC_VERSION)
end

# print struct info as TOML format
function Base.show(io::IO, ::MIME"text/plain", sim_info::SimulationInfo)

    @printf io "[simulation_info]\n\n"
    @printf io "name = \"%s\"\n" sim_info.datafolder_prefix
    @printf io "sID = %d\n" sim_info.sID
    @printf io "pID = %d\n" sim_info.pID
    @printf io "variationalmc_version = \"%s\"\n" sim_info.varmc_version
    @printf io "julia_version = \"%s\"" VERSION

    return nothing
end


@doc raw"""

    save_simulation_info(
        # ARGUMENTS
        sim_info::SimulationInfo,
        metadata = nothing;
        # KEYWORD ARGUMENTS
        filename = @sprintf "simulation_info_sID-%d_pID-%d.toml" sim_info.sID sim_info.pID
    )

Save the contents `sim_info` to a TOML file, and add an optional additional table to the
output file based on the contents of a dictionary `metadata`.
"""
function save_simulation_info(
    # ARGUMENTS
    sim_info::SimulationInfo,
    metadata = nothing;
    # KEYWORD ARGUMENTS
    filename = @sprintf "simulation_info_sID-%d_pID-%d.toml" sim_info.sID sim_info.pID
)

    (; datafolder ) = sim_info
    open(joinpath(datafolder, filename), "w") do fout
        show(fout, "text/plain", sim_info)
        if !isnothing(metadata)
            @printf fout "\n\n"
            TOML.print(fout, Dict("metadata" => metadata), sorted = true)
        end
    end

    return nothing
end


@doc raw"""
    initialize_datafolder(comm::MPI.Comm, sim_info::SimulationInfo)

    initialize_datafolder(sim_info::SimulationInfo)

Initialize `sim_info.datafolder` directory if it does not already exist.
If `comm::MPI.Comm` is passed as the first argument, this this function will synchronize
all the MPI processes, ensuring that none proceed beyond this function call until
the data folder that results will be written to is successfully initialized.
"""
function initialize_datafolder(comm::MPI.Comm, sim_info::SimulationInfo)

    MPI.Barrier(comm)
    initialize_datafolder(sim_info)

    return nothing
end

function initialize_datafolder(sim_info::SimulationInfo)

    (; pID, datafolder, write_bins_concurrent) = sim_info

    # make subdirectories for binned data to be written to
    opt_bin_dir = mkpath(joinpath(datafolder, "opt_bins"))
    sim_bin_dir = mkpath(joinpath(datafolder, "sim_bins"))

    if write_bins_concurrent
        mkpath(joinpath(opt_bin_dir, "pID-$(pID)"))
        mkpath(joinpath(sim_bin_dir, "pID-$(pID)"))
    end

    return nothing
end
















































# @doc raw"""

#     SimulationInfo{S<:AbstractString, I<:Int}

# A type containing information about the simulation, including where data is written,
# the simulation ID, and MPI process ID, and whether this simulation started a new 
# simulation or resumed a previous simulation.

# - `filepath::S`: path to data directory.
# - `datafolder_prefix::S`: prefix of the data directory name.
# - `datafolder_name::S`: data directory name.
# - `datafolder::S`: data directory including filepath.
# - `pID::I`: processor ID number/MPI rank.
# - `sID::I`: simulation ID number.
# - `resuming::Bool`: whether a previous simulation is being resumed (`true`) or starting a new one (`false`).
# - `variationalmc_version::VersionNumber`: Version of [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl) used in simulation.
#     # version of package (= VARIATIONALMC_VERSION)
#     variationalmc_version::VersionNumber
    
# """
# mutable struct SimulationInfo{S<:AbstractString, I<:Int}

#     # filepath to where data directory will live
#     filepath::S

#     # prefix of data directory name
#     datafolder_prefix::S

#     # data directory name
#     datafolder_name::S

#     # data directory including filepath
#     datafolder::S

#     # process ID number/MPI rank
#     pID::I

#     # simulation ID number
#     sID::I

#     # whether previous simulation is being resumed or a new one is begininning
#     resuming::Bool
# end


# @doc raw"""

#     SimulationInfo( datafolder_prefix::String 
#                     filepath::S = ".", 
#                     sID::I = 0, 
#                     pID::I = 0 ) where {S<:AbstractString, I<:Integer}

# Initialize and return an instance of the SimulationInfo type.

# """
# function SimulationInfo(; 
#     datafolder_prefix::S,
#     filepath::S = ".",
#     sID::I = 0,
#     pID::I = 0
# ) where {S<:AbstractString, I<:Integer}
#     # same body as before
#     datafolder_name = @sprintf("%s-%d", datafolder_prefix, sID)
#     datafolder = joinpath(filepath, datafolder_name)

#     if sID == 0
#         while isdir(datafolder) || sID == 0
#             sID += 1
#             datafolder_name = @sprintf("%s-%d", datafolder_prefix, sID)
#             datafolder = joinpath(filepath, datafolder_name)
#         end
#     end

#     resuming = isdir(datafolder)

#     return SimulationInfo(filepath, datafolder_prefix, datafolder_name, datafolder, pID, sID, resuming) #VARIATIONALMC_VERSION
# end


# # print struct info as TOML format
# function Base.show(io::IO, ::MIME"text/plain", sim_info::SimulationInfo)

#     @printf io "[simulation_info]\n\n"
#     @printf io "name  = \"%s\"\n" sim_info.datafolder_prefix
#     @printf io "sID   = %d\n" sim_info.sID
#     @printf io "pID   = %d\n" sim_info.pID
#     # @printf io "variationalmc_version = \"%s\"\n" sim_info.variationalmc_version
#     @printf io "julia_version     = \"%s\"" VERSION

#     return nothing
# end


# @doc raw"""

#     save_simulation_info( sim_info::SimulationInfo, 
#                           metadata::Dict{Any, Any} = nothing )

# Save the contents `sim_info` to a TOML file, and add an optional additional table to the
# output file based on the contents of a dictionary `metadata`.

# """
# function save_simulation_info(
#     sim_info::SimulationInfo, 
#     metadata::Dict{Any, Any} = nothing
# )
#     (; datafolder, pID, sID) = sim_info
#     fn = @sprintf "simulation_info_sID-%d_rank-%d.toml" sID pID
#     open(joinpath(datafolder, fn), "w") do fout
#         show(fout, "text/plain", sim_info)
#         if !isnothing(metadata)
#             @printf fout "\n\n"
#             TOML.print(fout, Dict("metadata" => metadata), sorted = true)
#         end
#     end

#     return nothing
# end


# @doc raw"""

#     initialize_datafolder( comm::MPI.Comm, 
#                            sim_info::SimulationInfo )

#     initialize_datafolder( sim_info::SimulationInfo )

# Initalize `sim_info.datafolder` directory if it does not already exist.
# If `comm::MPI.Comm` is passed as the first argument, this this function will synchronize
# all the MPI processes, ensuring that none proceed beyond this function call until
# the data folder that results will be written to is successfully initialized.
# """
# function initialize_datafolder(comm::MPI.Comm, sim_info::SimulationInfo)

#     (; pID, datafolder, resuming) = sim_info

#     # synchrnize MPI walkers
#     MPI.Barrier(comm)

#     # if main process and starting new simulation (not resuming an existing simulation)
#     if iszero(pID) && !resuming

#         # make data folder directory
#         mkdir(datafolder)
#     end

#     # synchrnize MPI walkers
#     MPI.Barrier(comm)

#     return nothing
# end

# function initialize_datafolder(
#     sim_info::SimulationInfo
# )
#     (; pID, datafolder, resuming) = sim_info

#     # if main process and starting new simulation (not resuming an existing simulation)
#     if iszero(pID) && !resuming

#         # make data folder diretory
#         mkdir(datafolder)
#     end

#     return nothing
# end


# @doc raw"""

#     create_datafolder_prefix( optimize::NamedTuple, 
#                               df_prefix::S ) where {S<:AbstractString}

# Check the optimization fields and appends parameter names to the end of the foldername.
# Returns the datafolder prefix.

# - `optimize::NamedTuple,`: field of optimization flags.
# - `df_prefix::S`: datafolder prefix.

# """
# function create_datafolder_prefix(
#     optimize::NamedTuple, 
#     df_prefix::S
# ) where {S<:AbstractString}
#     opt_keys = keys(optimize)

#     enabled_opts = [begin
#         k_str = string(k)
#         k_str_clean = replace(k_str, "_" => "")  # Remove underscores

#         if k_str_clean == "Δ0"
#             "swave"
#         elseif k_str_clean == "Δd"
#             "dwave"
#         elseif k_str_clean == "μ"
#             "mu"
#         elseif startswith(k_str_clean, "Δ")
#             k_str_clean[3:end]  # drop "Δ" = 2-byte Unicode, start from index 3
#         else
#             k_str_clean
#         end
#     end for k in opt_keys if optimize[k]]

#     opt_suffix = isempty(enabled_opts) ? "_none" : "_" * join(enabled_opts, "_")
#     return df_prefix * opt_suffix
# end