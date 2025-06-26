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
    
    # dimensions of the lattice
    L = model_geometry.lattice.L

    # number of determinantal parameters
    num_det_pars = determinantal_parameters.num_det_pars

    # number of variational parameters to be optimized
    num_vpars = determinantal_parameters.num_det_opts

    # initial parameters
    init_vpars = collect(values(determinantal_parameters.det_pars))

    # container to store optimization measurements
    optimization_measurements = Dict{String, Any}([
        ("parameters", init_vpars),                    
        ("Δk", zeros(num_det_pars)),                         
        ("ΔkΔkp", zeros(num_det_pars)),         
        ("ΔkE", zeros(num_det_pars)),                        
    ])      

    # dictionary to store simulation measurements
    simulation_measurements = Dict{String, Any}([
        ("local_density", 0.0),         
        ("double_occ", 0.0),       
        ("local_energy", ComplexF64(0.0)),           
        ("configuration", zeros(Int, N))              
    ])                     

    # dictionary to store correlation measurements
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

    # initial parameters
    init_vpars = collect(values(determinantal_parameters.det_pars))

    # container to store optimization measurements
    optimization_measurements = Dict{String, Any}([
        ("parameters", init_vpars),                     
        ("Δk", zeros(num_vpars)),                         
        ("ΔkΔkp", zeros(num_vpars,num_vpars)),       
        ("ΔkE", zeros(num_vpars)),                         
    ])      

    # dictionary to store simulation measurements
    simulation_measurements = Dict{String, Any}([
        ("global_density", 0.0),          
        ("double_occ", 0.0),      
        ("local_energy", ComplexF64(0.0)),           
        ("configuration", zeros(Int, N))              
    ])                     

    # dictionary to store correlation measurements
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


@doc raw"""

    initialize_correlation_measurement!( correlation_type::String,
                                         measurement_container::NamedTuple,
                                         model_geometry::ModelGeometry)::Nothing

Initializes a specified equal-time correlation measurement. 

- `correlation_type::String`: "density" or "spin" correlations.
- `measurement_container::NamedTuple`: container where mesurements are stored.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_correlation_measurement!(
    correlation_type::String,
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry
)::Nothing
    @assert correlation_type == "density" || correlation_type == "spin"

    # number of lattice sites
    N = model_geometry.lattice.n

    # add to measurement container
    measurement_container.correlation_measurements[correlation_type] = zeros(N, N)

    return nothing
end


@doc raw"""

    initialize_simulation_measurement!( type::String,
                                        observable::String,
                                        measurement_container::NamedTuple,
                                        model_geometry::ModelGeometry )::Nothing

Initializes a specified simulation observable measurement.

- `type::String`: "local" or "site-dependent"
- `observable::String`: "density" or "spin".
- `measurement_container::NamedTuple`: container where mesurements are stored.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_simulation_measurement!(
    type::String,
    observable::String,
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry
)::Nothing
    @assert type == "local" || type == "site-dependent"
    @assert observable == "density" || observable == "spin"

    # number of lattice sites
    N = model_geometry.lattice.N

    if type == "local"
        measurement_container.simulation_measurements[type * "_" * observable] = 0.0
    elseif type == "site-dependent"
        measurement_container.simulation_measurements[type * "_" * observable] = zeros(N)
    end

    return nothing
end


"""

    initialize_measurement_directories( simulation_info::SimulationInfo, 
                                        measurement_container::NamedTuple )::Nothing

Creates file directories and for storing measurements. 

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

        # make correlation measurements directory
        if !isempty(correlation_measurements)
            correlation_directory = joinpath(datafolder, "correlation")
            mkdir(correlation_directory)
        end
    end

    return nothing
end