@doc raw"""

    write_measurements!( bin::Int, 
                         step::Int, 
                         measurement_container::NamedTuple, 
                         simulation_info::SimulationInfo )::Nothing

Writes optimization and simulation measurements in the current bin to file. 
Files created are in JLD2 format. 

- `bin::Int`: current bin.
- `step::Int`: current step within the bin.
- `measurement_container::NamedTuple`: container where measurements are stored.
- `simulation_info`: contains datafolder names.

""" 
# TODO: this function will be modifed to write to HDF5 files instead of JLD2.
function write_measurements!(
    bin::Int, 
    step::Int, 
    measurement_container::NamedTuple, 
    simulation_info::SimulationInfo
)::Nothing
    (; datafolder, pID) = simulation_info
    (; optimization_measurements, simulation_measurements) = measurement_container

    fn = @sprintf "bin-%d_pID-%d.jld2" bin pID  

    file_path_dblocc  = joinpath(datafolder, "simulation", "double_occ", fn)
    file_path_den     = joinpath(datafolder, "simulation", "density", fn)
    file_path_configs = joinpath(datafolder, "simulation", "configurations", fn)
    file_path_energy  = joinpath(datafolder, "simulation", "energy", fn)

    # Define step key string
    step_key = "step_$step"

    # Extract latest measurements
    dblocc_measurement  = simulation_measurements["double_occ"][2]
    den_measurement     = simulation_measurements["density"][2]
    config_measurement  = simulation_measurements["pconfig"][2]
    energy_measurement  = simulation_measurements["energy"][2]

    # Use jldopen to save with dynamic keys
    jldopen(file_path_dblocc, "a") do file
        file[step_key] = dblocc_measurement
    end

    jldopen(file_path_den, "a") do file
        file[step_key] = den_measurement
    end

    jldopen(file_path_configs, "a") do file
        file[step_key] = config_measurement
    end

    jldopen(file_path_energy, "a") do file
        file[step_key] = energy_measurement
    end

    # Reset for next measurement
    reset_measurements!(simulation_measurements)
    reset_measurements!(optimization_measurements)

    return nothing
end


@doc raw"""

    write_measurements!( measurement_container::NamedTuple, 
                         energy_bin::Vector{Any}, 
                         dblocc_bin::Vector{Any}, 
                         param_bin::Vector{Any} )::Nothing

DEBUG version of the write_measurements!() method. Will write binned energies, double occupancy, and parameters
to specified vectors.  

- `measurement_container::NamedTuple`: container where measurements are stored.
- `energy_bin::Vector{Float64}`: externally specified vector for storing energy measurements.
- `dblocc_bin::Vector{Float64}`: externally specified vector for storing double occupancy measurements.
- `param_bin::Vector{Any}`: externally specified vector for storing parameter measurements.

""" 
function write_measurements!(
    measurement_container::NamedTuple, 
    energy_bin::Vector{Float64}, 
    dblocc_bin::Vector{Float64}, 
    param_bin::Vector{Any}
)::Nothing
    # extract container info
    simulation_measurements = measurement_container.simulation_measurements
    optimization_measurements = measurement_container.optimization_measurements

    # append accumulated values to the storage vectors
    push!(energy_bin, simulation_measurements["energy"][1])
    push!(dblocc_bin, simulation_measurements["double_occ"][1])
    push!(param_bin, optimization_measurements["parameters"][1])
    
    # reset all measurements
    reset_measurements!(simulation_measurements)
    reset_measurements!(optimization_measurements)

    return nothing
end