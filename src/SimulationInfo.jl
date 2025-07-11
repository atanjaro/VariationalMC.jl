@doc raw"""

    SimulationInfo( filepath::String, 
                    datafolder_prefix::String,
                    datafolder_name::String,
                    datafolder::String,
                    pID::Int,
                    sID::Int,
                    resuming::Bool )

A type containing information about the simulation, including the location data is written to,
the simulation ID, and MPI process ID, and whether this simulation started a new simulation or resumed
a previous simulation.

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

    # whether previous simulation is being resumed or a new one is begininning
    resuming::Bool
end


@doc raw"""

    SimulationInfo( datafolder_prefix::String, 
                    filepath::String = ".", 
                    sID::Int=0, 
                    pID::Int=0 )::SimulationInfo

Creates an instance of the SimulationInfo type.

"""
function SimulationInfo( 
    datafolder_prefix::String, 
    filepath::String = ".", 
    sID::Int=0, 
    pID::Int=0
)::SimulationInfo
    # initialize data folder names
    datafolder_name = @sprintf "%s-%d" datafolder_prefix sID
    datafolder = joinpath(filepath, datafolder_name)

    # if null data folder id given, determine data name and id
    if sID==0
        while isdir(datafolder) || sID==0
            sID += 1
            datafolder_name = @sprintf "%s-%d" datafolder_prefix sID
            datafolder = joinpath(filepath, datafolder_name)
        end
    end

    # if directory already exists then must be resuming simulation
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
                          additional_info = nothing )::Nothing

Save the contents `sim_info` to a TOML file, and add an optional additional table to the
output file based on the contents of a dictionary `additional_info`.

"""
function save_simulation_info(
    sim_info::SimulationInfo, 
    additional_info = nothing
)::Nothing
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

    initialize_datafolder( sim_info::SimulationInfo )::Nothing

Initalize `sim_info.datafolder` directory if it does not already exist.

"""
function initialize_datafolder(sim_info::SimulationInfo)::Nothing

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
                              df_prefix::String )::String

Check the optimization fields and appends parameter names to the end of the foldername.
Returns the datafolder prefix.

"""
function create_datafolder_prefix(
    optimize::NamedTuple, 
    df_prefix::String
)::String
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