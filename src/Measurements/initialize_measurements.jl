@doc raw"""

    initialize_measurement_container( N_opt::I, 
                                      opt_bin_size::I, 
                                      N_sim::I, 
                                      sim_bin_size::I,
                                      determinantal_parameters::DeterminantalParameters{I}, 
                                      model_geometry::ModelGeometry ) where {I<:Integer}

Initializes a set of dictionaries containing generic arrays for storing measurements.
    
- `N_opt::I`: number of optimization updates.
- `opt_bin_size::I`: length of an optimization bin.
- `N_sim::I`: number of simulation bins.
- `sim_bin_size::I`: length of a simulation bin. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_measurement_container(
    N_opt::I, 
    opt_bin_size::I, 
    N_sim::I, 
    sim_bin_size::I,
    determinantal_parameters::DeterminantalParameters{I}, 
    model_geometry::ModelGeometry
) where {I<:Integer}
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
        ("ΔkΔkp", zeros(num_det_pars, num_det_pars)),         
        ("ΔkE", zeros(num_det_pars)),                        
    ])      

    # dictionary to store simulation measurements
    simulation_measurements = Dict{String, Any}([
        ("global_density", 0.0),         
        ("double_occ", 0.0),       
        ("local_energy", ComplexF64(0.0)),           
        ("pconfig", zeros(Int, 2*N))              
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
        N_opt                     = N_opt,
        opt_bin_size              = opt_bin_size,
        N_sim                     = N_sim,
        sim_bin_size              = sim_bin_size,
        num_vpars                 = num_vpars,
        num_detpars               = num_det_pars                 
    )

    return measurement_container
end


@doc raw"""

    initialize_measurement_container( N_opt::I, 
                                      opt_bin_size::I, 
                                      N_sim::I, 
                                      sim_bin_size::I,
                                      determinantal_parameters::DeterminantalParameters{I}, 
                                      jastrow_parameters::JastrowParameters{S, K, V, I},
                                      model_geometry::ModelGeometry ) where {I<:Integer, S<:AbstractString, K, V}

Initializes a set of dictionaries containing generic arrays for storing measurements.
    
- `N_opt::I`: number of optimization updates.
- `opt_bin_size::I`: length of an optimization bin.
- `N_sim::I`: number of simulation bins.
- `sim_bin_size::I`: length of a simulation bin. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_measurement_container(
    N_opt::I, 
    opt_bin_size::I, 
    N_sim::I, 
    sim_bin_size::I,
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I}, 
    model_geometry::ModelGeometry
) where {I<:Integer, S<:AbstractString, K, V}
    # total number of lattice sites
    N = model_geometry.lattice.N
    
    # one side of the lattice
    L = model_geometry.lattice.L

    # number of determinantal_parameters
    num_detpars = determinantal_parameters.num_det_pars

    # number of Jastrow parameters
    num_jpars = jastrow_parameters.num_jpars

    # total number of variational parameters 
    num_vpars = num_detpars + num_jpars

    # initial parameters
    init_vpars = collect_parameters(determinantal_parameters, jastrow_parameters)

    # container to store optimization measurements
    optimization_measurements = Dict{String, Any}([
        ("parameters", init_vpars),                     
        ("Δk", zeros(num_vpars)),                         
        ("ΔkΔkp", zeros(num_vpars, num_vpars)),       
        ("ΔkE", zeros(num_vpars)),                         
    ])      

    # dictionary to store simulation measurements
    simulation_measurements = Dict{String, Any}([
        ("global_density", 0.0),          
        ("double_occ", 0.0),      
        ("local_energy", ComplexF64(0.0)),           
        ("pconfig", zeros(Int, 2*N))              
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
        N_opt                     = N_opt,
        opt_bin_size              = opt_bin_size,
        N_sim                     = N_sim,
        sim_bin_size              = sim_bin_size,
        num_vpars                 = num_vpars,
        num_detpars               = num_detpars,
        num_jpars                 = num_jpars                      
    )

    return measurement_container
end


@doc raw"""

    initialize_measurement_container( N_opt::I, 
                                      opt_bin_size::I, 
                                      N_sim::I, 
                                      sim_bin_size::I,
                                      determinantal_parameters::DeterminantalParameters{I}, 
                                      jastrow_parameters_1::JastrowParameters{S, K, V, I},
                                      jastrow_parameters_2::JastrowParameters{S, K, V, I},
                                      model_geometry::ModelGeometry ) where {I<:Integer, S<:AbstractString, K, V}

Initializes a set of dictionaries containing generic arrays for storing measurements.
    
- `N_opt::I`: number of optimization updates.
- `opt_bin_size::I`: length of an optimization bin.
- `N_sim::I`: number of simulation bins.
- `sim_bin_size::I`: length of a simulation bin. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_measurement_container(
    N_opt::I, 
    opt_bin_size::I, 
    N_sim::I, 
    sim_bin_size::I,
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I}, 
    jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
    model_geometry::ModelGeometry
) where {I<:Integer, S<:AbstractString, K, V}
    # total number of lattice sites
    N = model_geometry.lattice.N
    
    # one side of the lattice
    L = model_geometry.lattice.L

    # number of determinantal_parameters
    num_detpars = determinantal_parameters.num_det_pars

    # number of Jastrow parameters
    num_jpars = jastrow_parameters_1.num_jpars + jastrow_parameters_2.num_jpars 

    # number of variational parameters to be optimized
    num_vpars = num_detpars + num_jpars 

    # initial parameters
    init_vpars = collect_parameters(determinantal_parameters, jastrow_parameters_1, jastrow_parameters_2)

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
        ("pconfig", zeros(Int, N))              
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
        N_opt                     = N_opt,
        opt_bin_size              = opt_bin_size,
        N_sim                     = N_sim,
        sim_bin_size              = sim_bin_size,
        num_vpars                 = num_vpars,
        num_detpars               = num_detpars,
        num_jpars                 = num_jpars                      
    )

    return measurement_container
end


@doc raw"""

    initialize_correlation_measurement!( correlation_type::S,
                                         measurement_container::NamedTuple,
                                         model_geometry::ModelGeometry ) where {S<:AbstractString}

Initializes a specified equal-time correlation measurement. 

- `correlation_type::S`: "density" or "spin".
- `measurement_container::NamedTuple`: container where mesurements are stored.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_correlation_measurement!(
    correlation_type::S,
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry
) where {S<:AbstractString}
    @assert correlation_type == "density" || correlation_type == "spin"

    # number of lattice sites
    N = model_geometry.lattice.N

    # add to measurement container
    measurement_container.correlation_measurements[correlation_type] = zeros(N, N)

    return nothing
end


@doc raw"""

    initialize_simulation_measurement!( type::S,
                                        observable::S,
                                        measurement_container::NamedTuple,
                                        model_geometry::ModelGeometry ) where {S<:AbstractString}

Initializes a specified simulation observable measurement.

- `type::S`: "local" or "site-dependent".
- `observable::S`: "spin-z", "double_occ", "density" or "spin".
- `measurement_container::NamedTuple`: container where mesurements are stored.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function initialize_simulation_measurement!(
    type::S,
    observable::S,
    measurement_container::NamedTuple,
    model_geometry::ModelGeometry
) where {S<:AbstractString}
    @assert type == "local" || type == "site-dependent"

    # number of lattice sites
    N = model_geometry.lattice.N

    if type == "local"
        @assert observable == "spin-z" || "density" || observable == "spin"
        measurement_container.simulation_measurements[type * "_" * observable] = 0.0
    elseif type == "site-dependent"
        @assert observable == "density" || observable == "spin"
        measurement_container.simulation_measurements[type * "_" * observable] = zeros(N)
    end

    return nothing
end


@doc raw"""

    initialize_measurement_directories(;
        simulation_info::SimulationInfo,
        measurement_container::NamedTuple
    )

    initialize_measurement_directories(
            comm::MPI.Comm;
            simulation_info::SimulationInfo,
            measurement_container::NamedTuple
    )

    initialize_measurement_directories(
        simulation_info::SimulationInfo,
        measurement_container::NamedTuple
    )

    initialize_measurement_directories(
            comm::MPI.Comm,
            simulation_info::SimulationInfo,
            measurement_container::NamedTuple
    )

Initialize the measurement directories for simulation. If using MPI and a `comm::MPI.Comm` object is passed
as the first argument, then none of the MPI processes will proceed beyond this function call until the measurement
directories have been initialized.

"""
function initialize_measurement_directories(
        comm::MPI.Comm;
        simulation_info::SimulationInfo,
        measurement_container::NamedTuple
)

    initialize_measurement_directories(simulation_info, measurement_container)
    MPI.Barrier(comm)
    return nothing
end

function initialize_measurement_directories(
        comm::MPI.Comm,
        simulation_info::SimulationInfo,
        measurement_container::NamedTuple
)

    initialize_measurement_directories(simulation_info, measurement_container)
    MPI.Barrier(comm)
    return nothing
end

function initialize_measurement_directories(;
    simulation_info::SimulationInfo,
    measurement_container::NamedTuple
)
    initialize_measurement_directories(simulation_info, measurement_container)
    return nothing
end

function initialize_measurement_directories(
    simulation_info::SimulationInfo, 
    measurement_container::NamedTuple
)
    (; datafolder, resuming, pID) = simulation_info
    (; correlation_measurements) = measurement_container

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
