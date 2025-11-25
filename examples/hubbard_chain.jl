# Import the required packages
using LinearAlgebra
using Random
using Printf

using LatticeUtilities
using VariationalMC 

# # Uncomment this to trigger output of debug statements
# global_logger(ConsoleLogger(stderr, Logging.Debug))

# We define a top-level function for running the VMC simulation.
# Note that the arguments of this function correspond to the command line
# arguments used to run this script.
function run_hubbard_chain_simulation(
    sID, 
    L,
    U,
    nup,
    ndn,
    N_equil,
    N_opt,
    N_opt_bins,
    N_sim,
    N_sim_bins; 
    filepath="."
)
    # Select which parameters in the variational wavefunction will be optimized.
    optimize = (
        # spin-x (in-plane magnetization)
        Δ_sx = false,
        # spin-z (out-of-plane magnetization)
        Δ_sz = true,
        # (BCS) chemical potential
        μ = false,
        # uniform charge density
        Δ_cdw = false,
        # density-density Jastrow 
        density_J = false
    )

    # Specify whether the model will be particle-hole transformed.
    # Note that this is required if adding pair fields to the wavefunction.
    pht = false

    # Construct the foldername the data will be written.
    df_prefix = @sprintf("hubbard_chain_U%.2f_nup%.2f_ndn%.2f_L%d_opt", U, nup, ndn, L)

    # Append optimized parameter names to the foldername.
    datafolder_prefix = create_datafolder_prefix(optimize, df_prefix)

    # Initialize an instance of the SimulationInfo type.
    # This type tracks of where the data is written, as well as 
    # which version of Julia and VariationalMC are used in the script. 
    simulation_info = SimulationInfo(
        filepath = filepath, 
        datafolder_prefix = datafolder_prefix,
        sID = sID
    )

    # Initialize the directory the data will be written.
    initialize_datafolder(simulation_info)

    # Initialize a random number generator that will be used throughout the simulation.
    # We seed this function with a randomly sampled number for the
    # global random number generator.
    seed = abs(rand(Int))
    rng = Xoshiro(seed)

    # Set the optimization rate for the VMC simulation.
    dt = 0.1

    # Set the stabilization factor used in parameter optimization. 
    η = 1e-4    

    # Set the frequency in Monte Carlo steps with which the equal-time Green's function
    # matrix will be recomputed using the numerical stabilization procedure.
    n_stab_W = 50

    # Specify the maximum allowed error in the equal-time Green's function matrix that
    # is corrected by numerical stabilization.
    δW = 1e-3

    # Calculate optimization bins size.
    # The bin size is the number of measurements that are averaged over each time data is written
    # to file during optimization.
    opt_bin_size = div(N_opt, N_opt_bins) 

    # Calculate simulation bins size.
    # The bin size is the number of measurements that are averaged over each time data is written
    # to file during the simulation.
    sim_bin_size = div(N_sim, N_sim_bins)

    # Initialize a dictionary to store additional information about the simulation.
    # This is a sort of "notebook" for tracking extraneous parameters during the VMC simulation.
    metadata = Dict()
    metadata["seed"] = seed
    metadata["pht"] = pht
    metadata["N_equil"] = N_equil
    metadata["N_opt"] = N_opt
    metadata["N_sim"] = N_sim
    metadata["N_opt_bins"] = N_opt_bins
    metadata["N_sim_bins"] = N_sim_bins
    metadata["δW"] = δW
    metadata["n_stab_W"] = n_stab_W
    metadata["dt"] = dt 
    metadata["opt_flags"] = optimize 
    metadata["acceptance_rate"] = 0.0
    metadata["opt_time"] = 0.0
    metadata["sim_time"] = 0.0
    metadata["vmc_time"] = 0.0

    #######################
    ### DEFINE THE MODEL ##
    #######################

    # Initialize an instance of the UnitCell type.
    # This struct defines the unit cell.
    unit_cell = UnitCell(
        lattice_vecs = [[1.0]],
        basis_vecs   = [[0.0]]
    )

    # Initialize an instance of the Lattice type.
    # The struct describes the size of the finite periodic lattice to be simulated.
    # Note that the current version of LatticeUtilities requires 
    # periodic boundary conditions be used.
    lattice = Lattice(
        [L], 
        [true]
    )

    # Define the nearest neighbor bonds.
    bond_x = Bond(
        orbitals = (1,1), 
        displacement = [1]
    )

    # Define next-nearest neighbor bonds.
    bond_xp = Bond(
        orbitals = (1,1), 
        displacement = [2]
    )

    # Collect all bond definitions into a single vector.
    # Note that this has the structure [[nearest],[next-nearest]].
    bonds = [[bond_x], [bond_xp]]

    # Initialize an instance of the ModelGeometry type.
    # This type helps keep track of all the relevant features of the lattice
    # geometry being simulated, including the defintion of the unit cell,
    # the size of the finite periodic lattice, and all the relevant
    # bond defintions that may arise in the model.
    model_geometry = ModelGeometry(
        unit_cell, 
        lattice, 
        bonds
    )

    ############################
    ### SET-UP VMC SIMULATION ##
    ############################

    # Determine the total particle density in the canonical ensemble. 
    (density, Np, Ne, nup, ndn) = get_particle_density(nup, ndn, model_geometry, pht) 

    # Define the nearest neighbor hopping amplitude, setting the energy scale of the system. 
    t = 1.0

    # Define the next-nearest neighbor hopping amplitude.
    tp = 0.0

    # Define the third-nearest neighbor hopping amplitude.
    tpd = 0.0

    # Define the non-interacting tight binding model.
    tight_binding_model = TightBindingModel(t, tp, tpd)

    # Initialize determinantal variational parameters.
    determinantal_parameters = DeterminantalParameters(
        optimize, 
        tight_binding_model, 
        model_geometry, 
        Ne, 
        pht
    )

    # Initialize the (fermionic) particle configuration.
    # Will start with a random initial configuration unless provided a starting configuration.
    pconfig = Int[]

    ##############################
    ### INITIALIZE MEASUREMENTS ##
    ##############################

    # Initialize the container that measurements will be accumulated into.
    measurement_container = initialize_measurement_container(
        N_opt, 
        opt_bin_size, 
        N_sim, 
        sim_bin_size,
        determinantal_parameters,
        model_geometry
    )

    # Initialize the sub-directories to which the various measurements will be written.
    initialize_measurement_directories(
        simulation_info, 
        measurement_container
    )

    ###########################
    ### OPTIMIZATION UPDATES ##
    ###########################

    # Start time for optimization.
    opt_start = time()

    # Iterate over optimization bins.
    for bin in 1:N_opt_bins

        # Initialize the determinantal wavefunction.
        detwf = get_determinantal_wavefunction(
            tight_binding_model, 
            determinantal_parameters, 
            optimize, 
            Np, 
            nup, 
            ndn, 
            model_geometry, 
            rng,
            pht,
            pconfig
        )  

        # Iterate over optimization bin length
        for n in 1:opt_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf) = local_fermion_update!(
                    detwf, 
                    Ne, 
                    model_geometry, 
                    pht,
                    n_stab_W,
                    δW, 
                    rng
                )

                # Record acceptance rate.
                metadata["acceptance_rate"] += acceptance_rate
            end

            # Make measurements, with results being recorded in the measurement container.
            make_measurements!(
                measurement_container, 
                detwf, 
                tight_binding_model, 
                determinantal_parameters, 
                optimize,
                model_geometry, 
                U,
                Np, 
                pht
            )
        end

        # Record the last particle configuration used for the start of the next bin.
        pconfig = detwf.pconfig

        # Attempt to update the variational parameters using the Stochastic Reconfiguration procedure. 
        optimize_parameters!( 
            measurement_container,  
            determinantal_parameters, 
            η, 
            dt, 
            opt_bin_size
        )  

        # Write measurement for the current bin to file.
        write_measurements!(
            "opt",
            bin, 
            opt_bin_size,
            measurement_container, 
            simulation_info,
            write_parameters=true
        )
    end

    # End time for optimization.
    opt_end = time()

    # Record the total time for optimization.
    metadata["opt_time"] += opt_end - opt_start

    #########################
    ### SIMULATION UPDATES ##
    #########################

    # Start time for simulation.
    sim_start = time()

    # Iterate over simulation bins.
    for bin in 1:N_sim_bins

        # Initialize the determinantal wavefunction.
        detwf = get_determinantal_wavefunction(
            tight_binding_model, 
            determinantal_parameters, 
            optimize, 
            Np, 
            nup, 
            ndn, 
            model_geometry, 
            rng,
            pht,
            pconfig
        )  

        # Iterate over optimization bin length
        for n in 1:sim_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf) = local_fermion_update!(
                    detwf, 
                    Ne, 
                    model_geometry, 
                    pht,
                    n_stab_W,
                    δW, 
                    rng
                )

                # Record acceptance rate.
                metadata["acceptance_rate"] += acceptance_rate
            end

            # Make measurements, with results being recorded in the measurement container.
            make_measurements!(
                measurement_container, 
                detwf, 
                tight_binding_model, 
                determinantal_parameters, 
                model_geometry, 
                U,
                Np, 
                pht
            )
        end

        # Record the last particle configuration used for the start of the next bin.
        pconfig = detwf.pconfig

        # Write measurement for the current bin to file.
        write_measurements!(
            bin, 
            opt_bin_size,
            measurement_container, 
            simulation_info
        )
    end

    # End time for simulation.
    sim_end = time()

    # Record the total simulation time.
    metadata["sim_time"] += sim_end - sim_start

    # Record the total VMC time.
    metadata["vmc_time"] += metadata["opt_time"] + metadata["sim_time"]

    # Process all optimization and simulation measurements.
    # Each observable will be written to CSV files for later processing.
    process_measurements(
        measurement_container, 
        simulation_info, 
        determinantal_parameters
    )

    # Write model summary to file.
    model_summary(simulation_info, metadata)

    return nothing
end

# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Read in the command line arguments.
    sID = parse(Int, ARGS[1])           # simulation ID
    L = parse(Int, ARGS[2])
    U = parse(Float64, ARGS[3])
    nup = parse(Int, ARGS[4])
    ndn = parse(Int, ARGS[5])
    N_equil = parse(Int, ARGS[6])
    N_opt = parse(Int, ARGS[7])
    N_opt_bins = parse(Int, ARGS[8])
    N_sim = parse(Int, ARGS[9])
    N_sim_bins = parse(Int, ARGS[10])

    # Run the simulation.
    run_hubbard_chain_simulation(sID, L, U, nup, ndn, N_equil, N_opt, N_opt_bins, N_sim, N_sim_bins)
end