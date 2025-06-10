# Import the required packages
using LinearAlgebra
using Random
using Printf

using LatticeUtilities
using VariationalMC

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
    N_opts,
    N_updates,
    N_bins; 
    filepath="."
)
    # DEBUG flag which writes information to terminal during the simulation.
    # For efficeincy, this should always be turned off. 
    debug = false

    # Select which parameters in the variational wavefunction will be optimized.
    optimize = (
        # (BCS) chemical potential
        μ = false,
        # onsite (s-wave) pairing
        Δ_0 = false,
        # d-wave pairing
        Δ_d = false,
        # spin-z (Néel antiferromagnet)
        Δ_afm = true,
        # charge density 
        Δ_cdw = false,
        # site-dependent charge (stripe)
        Δsdc = false,
        # site-dependent spin (stripe)
        Δ_sds = false,
        # density-density Jastrow 
        djastrow = false,
        # spin-spin Jastrow
        sjastrow = false
    )

    # Specify whether the model will be particle-hole transformed.
    # Note that this is required if adding pair fields to the wavefunction.
    pht = false

    # Construct the foldername the data will be written.
    df_prefix = @sprintf "hubbard_chain_U%.2f_nup%.2f_ndn%.2f_L%d_opt" U, nup, ndn, L

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
    dt = 0.03     

    # Set the stabilization factor used in parameter optimization. 
    η = 1e-4    

    # Set the frequency in Monte Carlo steps with which the equal-time Green's function
    # matrix will be recomputed using the numerical stabilization procedure.
    n_stab_W = 50

    # Specify the maximum allowed error in the equal-time Green's function matrix that
    # is corrected by numerical stabilization.
    δW = 1e-3

    # Specify the optimization bin size.
    # The optimization bin size if the number of measurments that are averaged over each time data
    # is written to file during optimization. 
    opt_bin_size = 3000

    # Calculate the simulation bin size.
    # The bin size is the number of measurements that are averaged over each time data is written
    # to file during the simulation.
    bin_size = div(N_updates, N_bins)

    # Initialize a dictionary to store additional information about the simulation.
    # This is a sort of "notebook" for tracking extraneous parameters during the VMC simulation.
    metadata = Dict()
    metadata["N_equil"] = N_equil
    metadata["N_opts"] = N_opts
    metadata["N_updates"] = N_updates
    metadata["N_bins"] = N_bins
    metadata["seed"] = seed
    metadata["pht"] = true
    metadata["δW"] = δW
    metadata["dt"] = dt 
    metadata["opt_flags"] = optimize 
    metadata["avg_acceptance_rate"] = 0.0
    metadata["opt_time"] = 0.0
    metadata["sim_time"] = 0.0
    metadata["total_time"] = 0.0

    #######################
    ### DEFINE THE MODEL ##
    #######################

    # Initialize an instance of the UnitCell type.
    # This struct defines the UnitCell.
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

    # Define the nearest neighbor bond for a 1D chain.
    bond_x = Bond(
        orbitals = (1,1), 
        displacement = [1]
    )

    # Define the next-nearest neighbor bonds for a 1D chain.
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

    ##############################
    ### INITIALIZE MEASUREMENTS ##
    ##############################

    # Initialize the container that measurements will be accumulated into.
    measurement_container = initialize_measurement_container(
        N_opts, 
        opt_bin_size, 
        N_bins, 
        bin_size,
        determinantal_parameters,
        model_geometry
    )

    # Initialize the sub-directories to which the various measurements will be written.
    initialize_measurement_directories(
        simulation_info, 
        measurement_container
    )

    ############################
    ### SET-UP VMC SIMULATION ##
    ############################

    # Determine the total particle density in the canonical ensemble. 
    (density, Np, Ne) = get_particle_density(nup, ndn)

    # Define the nearest neighbor hopping amplitude, setting the energy scale of the system. 
    t = 1.0;

    # Define the next-nearest neighbor hopping amplitude.
    tp = 0.0;

    # Define the non-interacting tight binding model.
    tight_binding_model = TightBindingModel(t, tp)

    # Specify the minimum value of each variational parameter.
    # Note that this is done to avoid open shell issues.
    minabs_vpar = 1e-4;

    # Initialize determinantal variational parameters.
    determinantal_parameters = DeterminantalParameters(
        optimize, 
        tight_binding_model, 
        model_geometry, 
        minabs_vpar, 
        Ne, 
        pht
    )

    # Initialize the (fermionic) particle configuration cache.
    pconfig_cache = nothing

    ###########################
    ### OPTIMIZATION UPDATES ##
    ###########################

    # Record start time for optimization. 
    opt_start_time = time()

    # Iterate over optimization updates.
    for bin in 1:N_opts

        # Initialize the determinantal wavefunction.
        detwf = get_determinantal_wavefunction(
            tight_binding_model, 
            determinantal_parameters, 
            optimize, 
            Ne, 
            nup, 
            ndn, 
            model_geometry, 
            rng,
            pconfig_cache
        )   

        # Iterate over equilibration/thermalization updates.
        for step in 1:N_equil 

            # Attempt to update the fermionic particle configuration.
            (acceptance_rate, detwf) = local_fermion_update!(
                detwf, 
                Ne, 
                model_geometry, 
                n_stab_W,
                δW, 
                rng
            )

            # Record the acceptance rate for attempted local updates of the particle configuration.                                                         
            metadata["avg_acceptance_rate"] += acceptance_rate;
        end

        # Iterate over the number of optimization bins.
        for n in 1:opt_bin_size

            # Attempt to update the fermionic particle configuration.
            (acceptance_rate, detwf) = local_fermion_update!(
                detwf, 
                Ne, 
                model_geometry, 
                n_stab_W,
                δW, 
                rng
            ) 
                                                                                                        
            # Record the acceptance rate for attempted local updates of the particle configuration                                                       
            metadata["avg_acceptance_rate"] += acceptance_rate;
            
            # Make measurements, with results being recorded in the measurement container. 
            make_measurements!(
                measurement_container, 
                detwf, 
                tight_binding_model, 
                determinantal_parameters, 
                model_geometry, 
                Ne, 
                pht
            )

            # Write measurement for the current bin to file.
            write_measurements!(
                bin, 
                n, 
                measurement_container, 
                simulation_info
            )
        end

        # Record the last particle configuration used in the current bin. 
        pconfig_cache = detwf.pconfig

        # Process the measurement results, calculating error bars for all measurements. 
        # process_measurements(simulation_info, opt_bin_size)

        # Attempt an update to the variational parameters using the Stochastic Reconfiguration procedure. 
        stochastic_reconfiguration!( 
            measurement_container,  
            determinantal_parameters, 
            η, 
            dt, 
            bin, 
            opt_bin_size
        )  
    end     

    # Record end time for optimization.
    opt_end_time = time()

    # Calculate the total time for optimization. 
    metadata["opt_time"] += opt_end_time - opt_start_time

    #########################
    ### SIMULATION UPDATES ##
    #########################

    # Record start time for the simulation.
    sim_start_time = time()

    # Iterate over simulation updates.
    for bin in 1:N_updates

        # Initialize the determinantal wavefunction.
        detwf = get_determinantal_wavefunction(
            tight_binding_model, 
            determinantal_parameters, 
            optimize, 
            Ne, 
            nup, 
            ndn, 
            model_geometry, 
            rng,
            pconfig_cache
        )   

        # Iterate over equilibration/thermalization updates.
        for step in 1:N_equil 

            # Attempt to update the fermionic particle configuration.
            (acceptance_rate, detwf) = local_fermion_update!(
                detwf, 
                Ne, 
                model_geometry, 
                n_stab_W,
                δW, 
                rng
            )

            # Record the acceptance rate for attempted local updates of the particle configuration.                                                         
            metadata["avg_acceptance_rate"] += acceptance_rate;
        end

        # Iterate over the number of simulation bins.
        for n in 1:bin_size

            # Attempt to update the fermionic particle configuration.
            (acceptance_rate, detwf) = local_fermion_update!(
                detwf, 
                Ne, 
                model_geometry, 
                n_stab_W,
                δW, 
                rng
            ) 
                                                                                                        
            # Record the acceptance rate for attempted local updates of the particle configuration                                                       
            metadata["avg_acceptance_rate"] += acceptance_rate;
            
            # Make measurements, with results being recorded in the measurement container. 
            make_measurements!(
                measurement_container, 
                detwf, 
                tight_binding_model, 
                determinantal_parameters, 
                model_geometry, 
                Ne, 
                pht
            )

            # Write measurement for the current bin to file.
            write_measurements!(
                bin, 
                n, 
                measurement_container, 
                simulation_info
            )
        end

        # Record the last particle configuration used in the current bin. 
        pconfig_cache = detwf.pconfig

        # Process the measurement results, calculating error bars for all measurements. 
        # process_measurements(simulation_info, bin_size)
    end     

    # Record end time for the simulation. 
    sim_end_time = time()

    # Calculate total time for the simulation. 
    metadata["sim_time"] += sim_end_time - sim_start_time

    # Record the total runtime.
    metadata["total_time"] += metadata["opt_time"] + metadata["sim_time"]

    # # Write simulation summary TOML file.
    # save_simulation_info(simulation_info, metadata)

    return nothing
end 

# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Read in the command line arguments.
    sID = parse(Int, ARGS[1]) # simulation ID
    L = parse(Int, ARGS[2])
    U = parse(Float64, ARGS[3])
    nup = parse(Int, ARGS[4])
    ndn = parse(Int, ARGS[5])
    N_equil = parse(Int, ARGS[6])
    N_opts = parse(Int, ARGS[7])
    N_updates = parse(Int, ARGS[8])
    N_bins = parse(Int, ARGS[9])

    # Run the simulation.
    run_hubbard_chain_simulation(sID, L, U, nup, ndn, N_equil, N_opts, N_updates, N_bins)
end

