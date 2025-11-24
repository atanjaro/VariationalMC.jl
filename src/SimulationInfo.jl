@doc raw"""

    SimulationInfo{S<:AbstractString, I<:Int}

A type containing information about the simulation, including where data is written,
the simulation ID, and MPI process ID, and whether this simulation started a new 
simulation or resumed a previous simulation.

- `filepath::S`: path to data directory.
- `datafolder_prefix::S`: prefix of the data directory name.
- `datafolder_name::S`: data directory name.
- `datafolder::S`: data directory including filepath.
- `pID::I`: processor ID number.
- `sID::I`: simulation ID number.
- `resuming::Bool`: whether a previous simulation is being resumed.

"""
mutable struct SimulationInfo{S<:AbstractString, I<:Int}

    # filepath to where data directory will live
    filepath::S

    # prefix of data directory name
    datafolder_prefix::S

    # data directory name
    datafolder_name::S

    # data directory including filepath
    datafolder::S

    # process ID number (for MPI)
    pID::I

    # simulation ID number
    sID::I

    # whether previous simulation is being resumed or a new one is begininning
    resuming::Bool
end


@doc raw"""

    SimulationInfo( datafolder_prefix::String 
                    filepath::S = ".", 
                    sID::I = 0, 
                    pID::I = 0 ) where {S<:AbstractString, I<:Integer}

Creates an instance of the SimulationInfo type.

"""
function SimulationInfo(; 
    datafolder_prefix::S,
    filepath::S = ".",
    sID::I = 0,
    pID::I = 0
) where {S<:AbstractString, I<:Integer}
    # same body as before
    datafolder_name = @sprintf("%s-%d", datafolder_prefix, sID)
    datafolder = joinpath(filepath, datafolder_name)

    if sID == 0
        while isdir(datafolder) || sID == 0
            sID += 1
            datafolder_name = @sprintf("%s-%d", datafolder_prefix, sID)
            datafolder = joinpath(filepath, datafolder_name)
        end
    end

    resuming = isdir(datafolder)

    return SimulationInfo(filepath, datafolder_prefix, datafolder_name, datafolder, pID, sID, resuming)
end


# print struct info as TOML format
function Base.show(io::IO, ::MIME"text/plain", sim_info::SimulationInfo)

    @printf io "[simulation_info]\n\n"
    @printf io "name  = \"%s\"\n" sim_info.datafolder_prefix
    @printf io "sID   = %d\n" sim_info.sID
    @printf io "pID   = %d\n" sim_info.pID
    @printf io "julia_version     = \"%s\"" VERSION

    return nothing
end


@doc raw"""

    save_simulation_info( sim_info::SimulationInfo, 
                          additional_info = nothing )

Save the contents `sim_info` to a TOML file, and add an optional additional table to the
output file based on the contents of a dictionary `additional_info`.

"""
function save_simulation_info(
    sim_info::SimulationInfo, 
    additional_info = nothing
)
    (; datafolder, pID, sID) = sim_info
    fn = @sprintf "simulation_info_pID%d_sID%d.toml" pID sID
    open(joinpath(datafolder, fn), "w") do fout
        show(fout, "text/plain", sim_info)
        if !isnothing(additional_info)
            @printf fout "\n\n"
            TOML.print(fout, Dict("additional_info" => additional_info), sorted = true)
        end
    end

    return nothing
end


@doc raw"""

    initialize_datafolder( sim_info::SimulationInfo )

Initalize `sim_info.datafolder` directory if it does not already exist.

"""
function initialize_datafolder(
    sim_info::SimulationInfo
)
    (; pID, datafolder, resuming) = sim_info

    # if main process and starting new simulation (not resuming an existing simulation)
    if iszero(pID) && !resuming

        # make data folder diretory
        mkdir(datafolder)
    end

    return nothing
end


@doc raw"""

    create_datafolder_prefix( optimize::NamedTuple, 
                              df_prefix::S ) where {S<:AbstractString}

Check the optimization fields and appends parameter names to the end of the foldername.
Returns the datafolder prefix.

- `optimize::NamedTuple,`: field of optimization flags.
- `df_prefix::S`: datafolder prefix.

"""
function create_datafolder_prefix(
    optimize::NamedTuple, 
    df_prefix::S
) where {S<:AbstractString}
    opt_keys = keys(optimize)

    enabled_opts = [begin
        k_str = string(k)
        k_str_clean = replace(k_str, "_" => "")  # Remove underscores

        if k_str_clean == "Δ0"
            "swave"
        elseif k_str_clean == "Δd"
            "dwave"
        elseif k_str_clean == "μ"
            "mu"
        elseif startswith(k_str_clean, "Δ")
            k_str_clean[3:end]  # drop "Δ" = 2-byte Unicode, start from index 3
        else
            k_str_clean
        end
    end for k in opt_keys if optimize[k]]

    opt_suffix = isempty(enabled_opts) ? "_none" : "_" * join(enabled_opts, "_")
    return df_prefix * opt_suffix
end