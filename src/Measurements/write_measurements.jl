@doc raw"""

    write_measurements!( bin::Int, 
                         step::Int, 
                         measurement_container::NamedTuple, 
                         simulation_info::SimulationInfo )::Nothing

Writes optimization and simulation measurements in the current bin to file. 
Files created are in HDF5 format. 

- `bin::Int`: current bin.
- `step::Int`: current step within the bin.
- `measurement_container::NamedTuple`: container where measurements are stored.
- `simulation_info`: contains datafolder names.

""" 
function write_measurements!(
    bin::Int, 
    step::Int, 
    measurement_container::NamedTuple, 
    simulation_info::SimulationInfo
)::Nothing
    (; datafolder, pID) = simulation_info
    (; simulation_measurements, optimization_measurements, correlation_measurements) = measurement_container

    # build path to simulation HDF5 bin file
    sim_file_path = joinpath(datafolder, "simulation", @sprintf("bin-%d.h5", bin))
    opt_file_path = joinpath(datafolder, "optimization", @sprintf("bin-%d.h5", bin))
    corr_file_path = joinpath(datafolder, "correlation", @sprintf("bin-%d.h5", bin))

    # prepare step key
    step_key = "step_$step"

    # extract simulation measurements
    density_measurement = simulation_measurements["density"][2]
    double_occ_measurement = simulation_measurements["double_occ"][2]
    energy_measurement = simulation_measurements["energy"][2]
    pconfig_measurement = simulation_measurements["pconfig"][2]
    # TODO: check for site-dependent density and spin measurements

    # extract correlation measurements, if measuring
    if !isempty(correlation_measurements)
        if !isempty(correlation_measurements["density"])
            dcorr_measurement = correlation_measurements["density"][2]
        elseif !isempty(correlation_measurements["spin"])
            scorr_measurement = correlation_measurements["spin"][2]
        end
    end

    # extract optimization measurements, if storing
    # Δk_measurement = optimization_measurements["Δk"][2]
    # ΔkΔkp_measurement = optimization_measurements["ΔkΔkp"][2]
    # ΔkE_measurement = optimization_measurements["ΔkE"][2]

    # write to HDF5
    h5open(sim_file_path, "a") do file
        file["density/$step_key"] = density_measurement
        file["double_occ/$step_key"] = double_occ_measurement
        file["energy/$step_key"] = energy_measurement
        file["pconfig/$step_key"] = pconfig_measurement
        # file["density-density"] = dcorr_measurement
        # file["spin-spin"] = scorr_measurement
        # file["sd-density"] = sdden_measurement
        # file["sd-spin"] = sdspin_measurement
    end

    # reset all measurements
    reset_measurements!(simulation_measurements)
    reset_measurements!(optimization_measurements)
    # reset_measurements!(correlation_measurements)

    return nothing
end


@doc raw"""

    write_measurements!( measurement_container::NamedTuple, 
                         energy_bin::Vector{Any}, 
                         dblocc_bin::Vector{Any}, 
                         param_bin::Vector{Any} )::Nothing

DEBUG version of the write_measurements!() method. Will write individuallly binned energies, double occupancy, and parameters
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