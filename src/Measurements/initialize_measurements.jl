@doc raw"""

    initialize_measurement_container( N_opts::Int, 
                                      opt_bin_size::Int, 
                                      N_bins::Int, 
                                      bin_size::Int,
                                      determinantal_parameters::DeterminantalParameters, 
                                      model_geometry::ModelGeometry )::NamedTuple

Initializes a set of dictionaries containing generic arrays for storing measurements. Each dictionary in the
container has [keys => values] of the form: [observable name => (local value(s), current bin value(s))] 
    
- `N_opts::Int`: number of optimization updates.
- `opt_bin_size::Int`: length of an optimization bin.
- `N_bins::Int`: number of simulation bins.
- `bin_size::Int`: length of a simulation bin. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_measurement_container(
    N_opts::Int, 
    opt_bin_size::Int, 
    N_bins::Int, 
    bin_size::Int,
    determinantal_parameters::DeterminantalParameters, 
    model_geometry::ModelGeometry
)::NamedTuple
    # total number of lattice sites
    N = model_geometry.lattice.N
    
    # one side of the lattice
    L = model_geometry.lattice.L

    # number of determinantal parameters
    num_det_pars = determinantal_parameters.num_det_pars

    # number of variational parameters to be optimized
    num_vpars = determinantal_parameters.num_det_opts

    # initial parameters
    init_vpars = collect(values(determinantal_parameters.det_pars))

    # container to store optimization measurements
    optimization_measurements = Dict{String, Any}([
        ("parameters", (init_vpars, init_vpars)),                    
        ("Δk", (zeros(num_det_pars), zeros(num_det_pars),[])),                         
        ("ΔkΔkp", (zeros(num_det_pars, num_det_pars), zeros(num_det_pars, num_det_pars),[])),         
        ("ΔkE", (zeros(num_det_pars), zeros(num_det_pars),[])),                        
    ])      

    # dictionary to store simulation measurements
    simulation_measurements = Dict{String, Any}([
        ("density", (0.0,  0.0)),         
        ("double_occ", (0.0,  0.0)),       
        ("energy", (0.0,  0.0)),           
        ("pconfig", zeros(N))              
    ])                     

    # TODO: dictionary to store correlation measurements
    correlation_measurements = Dict{String, Any}()

    # create container
    measurement_container = (
        simulation_measurements   = simulation_measurements,
        optimization_measurements = optimization_measurements,         
        correlation_measurements  = correlation_measurements,                       
        L                         = L,
        N                         = N,
        N_opts                    = N_opts,
        opt_bin_size              = opt_bin_size,
        N_updates                 = N_updates,
        N_bins                    = N_bins,
        bin_size                  = bin_size,
        num_vpars                 = num_vpars,
        num_detpars               = num_det_pars                 
    )

    return measurement_container
end


@doc raw"""

    initialize_measurement_container( N_opts::Int, 
                                      opt_bin_size::Int, 
                                      N_bins::Int, 
                                      bin_size::Int,
                                      determinantal_parameters::DeterminantalParameters, 
                                      jastrow_parameters::JastrowParameters,
                                      model_geometry::ModelGeometry )::NamedTuple

Initializes a set of dictionaries containing generic arrays for storing measurements. Each dictionary in the
container has [keys => values] of the form: [observable name => (local value(s), current bin value(s))] 
    
- `N_opts::Int`: number of optimization updates.
- `opt_bin_size::Int`: length of an optimization bin.
- `N_bins::Int`: number of simulation bins.
- `bin_size::Int`: length of a simulation bin. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_measurement_container(
    N_opts::Int, 
    opt_bin_size::Int, 
    N_bins::Int, 
    bin_size::Int,
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters, 
    model_geometry::ModelGeometry
)::NamedTuple
    # total number of lattice sites
    N = model_geometry.lattice.N
    
    # one side of the lattice
    L = model_geometry.lattice.L

    # number of determinantal_parameters
    num_detpars = determinantal_parameters.num_det_pars

    # number of Jastrow parameters
    num_jpars = jastrow_parameters.num_jpars

    # number of variational parameters to be optimized
    num_vpars = determinantal_parameters.num_det_opts + jastrow_parameters.num_jpar_opts

    # container to store optimization measurements
    optimization_measurements = Dict{String, Any}([
        ("parameters", (zeros(num_vpars), zeros(num_vpars))),                     
        ("Δk", (zeros(num_vpars), zeros(num_vpars),[])),                         
        ("ΔkΔkp", (zeros(num_vpars,num_vpars), zeros(num_vpars,num_vpars),[])),       
        ("ΔkE", (zeros(num_vpars), zeros(num_vpars),[])),                         
    ])      

    # dictionary to store simulation measurements
    simulation_measurements = Dict{String, Any}([
        ("density", (0.0,  0.0)),          
        ("double_occ", (0.0,  0.0)),      
        ("energy", (0.0,  0.0)),           
        ("pconfig", zeros(N))              
    ])                     

    # TODO: dictionary to store correlation measurements
    correlation_measurements = Dict{String, Any}()

    # create container
    measurement_container = (
        simulation_measurements   = simulation_measurements,
        optimization_measurements = optimization_measurements,         
        correlation_measurements  = correlation_measurements,                       
        L                         = L,
        N                         = N,
        N_opts                    = N_opts,
        opt_bin_size              = opt_bin_size,
        N_updates                 = N_updates,
        N_bins                    = N_bins,
        bin_size                  = bin_size,
        num_vpars                 = num_vpars,
        num_detpars               = num_detpars,
        num_jpars                 = num_jpars                      
    )

    return measurement_container
end


"""

    initialize_measurement_directories( simulation_info::SimulationInfo, 
                                        measurement_container::NamedTuple )::Nothing

Creates file directories and subdirectories for storing measurements. 

- `simulation_info::SimulationInfo`: contains datafolder information.
- `measurement_container::NamedTuple`: container where measurements are contained.

"""
function initialize_measurement_directories(
    simulation_info::SimulationInfo, 
    measurement_container::NamedTuple
)::Nothing
    (; datafolder, resuming, pID) = simulation_info
    (; optimization_measurements, simulation_measurements, correlation_measurements) = measurement_container

    # only initialize folders if pID = 0
    if iszero(pID) && !resuming

        # make optimization measurements directory
        optimization_directory = joinpath(datafolder, "optimization")
        mkdir(optimization_directory)

        # make simulation measurements directory
        simulation_directory = joinpath(datafolder, "simulation")
        mkdir(simulation_directory)

        # # make global measurements directory
        # global_directory = joinpath(datafolder, "global")
        # mkdir(global_directory)

        # make directories for each parameter measurement
        detpars_directory = joinpath(optimization_directory, "determinantal")
        mkdir(detpars_directory)
        jpars_directory = joinpath(optimization_directory, "Jastrow")
        mkdir(jpars_directory)

        # make energy measurement directory
        energy_directory = joinpath(simulation_directory, "energy")
        mkdir(energy_directory)

        # make configuration measurement directory
        config_directory = joinpath(simulation_directory, "configurations")
        mkdir(config_directory)

        # make double occupancy measurement directory
        dblocc_directory = joinpath(simulation_directory, "double_occ")
        mkdir(dblocc_directory)

        # make density measurement directory
        density_directory = joinpath(simulation_directory, "density")
        mkdir(density_directory)

        # TODO: add correlation measurement directory initialization
    end

    return nothing
end