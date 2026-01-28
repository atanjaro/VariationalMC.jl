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
function run_hubbard_square_simulation(;
    sID,                    # Simulation ID.
    Lx,                     # System size x.
    Ly,                     # System size y.
    U,                      # Hubbard interaction.
    density,                # Electron density.
    pht,                    # Whether model is particle-hole transformed. 
    N_equil,                # Number of equilibration/thermalization updates.
    N_opt,                  # Number of optimization steps.
    N_opt_bins,             # Number of times bin-averaged measurements are written to file during optimization step.
    N_sim,                  # Number of simulation steps.
    N_sim_bins,             # Number of times bin-averaged measurements are written to file during simulation step.
    dt = 0.1,              # Optimization rate.
    dt_J = 1.0,             # Optional boost in the Jastrow optimization rate.
    η = 1e-4,               # Optimization stablity factor.
    n_stab_W = 50,          # Green's function stabilization frequency.
    δW = 1e-3,              # Maximum allowed error in the Green's function. 
    n_stab_T = 50,          # Jastrow factor stabilization frequency.
    δT = 1e-3,              # Maximum allowed error in the Jastrow factor.           
    seed = abs(rand(Int)),  # Seed for random number generator.
    filepath="."            # Filepath to where data folder will be created.
)
    # Select which parameters in the variational wavefunction will be optimized.
    optimize = (
        # local s-wave pairing
        Δ_0 = false,
        # site-dependent s-wave pairing (Larkin-Ovchinnikov-type)
        Δ_slo = false,
        # site-dependent s-wave pairing (Fulde-Ferrell-type)
        Δ_sff = false,
        # local d-wave pairing
        Δ_d = false,
        # site-dependent d-wave pairing (Larkin-Ovchinnikov-type)
        Δ_dlo = false,
        # site-dependent d-wave pairing (Fulde-Ferrell-type)
        Δ_dff = false,     
        # spin-x (in-plane magnetization)
        Δ_sx = false,
        # spin-z (out-of-plane magnetization)
        Δ_sz = true,
        # site-dependent spin density
        Δ_ssd = false,
        # (BCS) chemical potential
        μ = false,
        # uniform charge density 
        Δ_cdw = false,
        # site-dependent charge density
        Δ_csd = false,
        # density-density Jastrow 
        density_J = true,
        # spin-spin Jastrow 
        spin_J = true
    )

    # Construct the foldername the data will be written.
    df_prefix = @sprintf("hubbard_rectangle_U%.2f_density%.3f_Lx%d_Ly%d_opt", U, density, Lx, Ly)

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
    rng = Xoshiro(seed)

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
    metadata["N_equil"] = N_equil
    metadata["N_opt"] = N_opt
    metadata["N_sim"] = N_sim
    metadata["N_opt_bins"] = N_opt_bins
    metadata["N_sim_bins"] = N_sim_bins
    metadata["δW"] = δW
    metadata["δT"] = δT
    metadata["n_stab_W"] = n_stab_W
    metadata["n_stab_T"] = n_stab_T
    metadata["dt"] = dt 
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
        lattice_vecs = [[1.0, 0.0], [0.0, 1.0]],
        basis_vecs   = [[0.0, 0.0]]
    )

    # Initialize an instance of the Lattice type.
    # The struct describes the size of the finite periodic lattice to be simulated.
    # Note that the current version of LatticeUtilities requires 
    # periodic boundary conditions be used.
    lattice = Lattice(
        [Lx, Ly], 
        [true, true]
    )

    # Define the nearest neighbor x-bond for a square lattice.
    bond_x = Bond(
        orbitals = (1,1), 
        displacement = [1,0]
    )

    # Define the nearest neighbor y-bond for a square lattice.
    bond_y = Bond(
        orbitals = (1,1), 
        displacement = [0,1]
    )

    # Define the next-nearest neighbor bonds for a square lattice.
    bond_xy = Bond(
        orbitals = (1,1), 
        displacement = [1,1]
    )

    # Define the next-nearest neighbor bonds for a square lattice.
    bond_yx = Bond(
        orbitals = (1,1), 
        displacement = [1,-1]
    )

    # Collect all bond definitions into a single vector.
    # Note that this has the structure [[nearest],[next-nearest]].
    bonds = [[bond_x, bond_y], [bond_xy, bond_yx]]

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
    (density, Np, Ne, nup, ndn) = get_particle_density(density, model_geometry, pht) 

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

    # Initialize density-density Jastrow variational parameters.
    density_J_parameters = JastrowParameters(
        "e-den-den",
        optimize, 
        model_geometry,
        rng
    )

    # Initialize spin-spin Jastrow variational parameters.
    spin_J_parameters = JastrowParameters(
        "e-spn-spn",
        optimize, 
        model_geometry,
        rng
    )

    # Write model summary TOML file specifying the Hamiltonian that will be simulated.
    model_summary(
        simulation_info, 
        determinantal_parameters, 
        density_J_parameters, 
        spin_J_parameters,
        pht, 
        model_geometry, 
        tight_binding_model, 
        U
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
        density_J_parameters,
        spin_J_parameters,
        model_geometry
    )

    # add measurement of site-dependent Sz
    initialize_simulation_measurement!(
        "site-dependent",
        "spin-z",
        measurement_container,
        model_geometry
    )

    # add measurement of site-dependent n
    initialize_simulation_measurement!(
        "site-dependent",
        "density",
        measurement_container,
        model_geometry
    )

    # Add density-density correlation measurements.
    initialize_correlation_measurement!(
        "density", 
        measurement_container, 
        model_geometry
    )

    # Add spin-spin correlation measurements.
    initialize_correlation_measurement!(
        "spin", 
        measurement_container, 
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

    # Record start time for optimization.
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

        # Initialize density-density Jastrow factor.
        density_J_factor = get_jastrow_factor(
            density_J_parameters,
            detwf,
            model_geometry,
            pht
        )

        # Initialize spin-spin Jastrow factor.
        spin_J_factor = get_jastrow_factor(
            spin_J_parameters,
            detwf,
            model_geometry,
            pht
        )

        # Iterate over optimization bin length
        for n in 1:opt_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf, density_J_factor, spin_J_factor) = local_fermion_update!(
                    detwf, 
                    density_J_factor,
                    spin_J_factor,
                    density_J_parameters,
                    spin_J_parameters,
                    Np, 
                    model_geometry, 
                    pht,
                    n_stab_W,
                    n_stab_T,
                    δW, 
                    δT,
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
                density_J_parameters,
                spin_J_parameters,
                density_J_factor,
                spin_J_factor,
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
            density_J_parameters,
            spin_J_parameters,
            η, 
            dt, 
            dt_J,
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

    # Record end time for optimization.
    opt_end = time()

    # Record the total time for optimization.
    metadata["opt_time"] += opt_end - opt_start

    #########################
    ### SIMULATION UPDATES ##
    #########################

    # Record start time for simulation.
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

        # Initialize density-density Jastrow factor.
        density_J_factor = get_jastrow_factor(
            density_J_parameters,
            detwf,
            model_geometry,
            pht
        )

        # Initialize spin-spin Jastrow factor.
        spin_J_factor = get_jastrow_factor(
            spin_J_parameters,
            detwf,
            model_geometry,
            pht
        )

        # Iterate over optimization bin length
        for n in 1:sim_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf, density_J_factor, spin_J_factor) = local_fermion_update!(
                    detwf, 
                    density_J_factor,
                    spin_J_factor,
                    density_J_parameters,
                    spin_J_parameters,
                    Np, 
                    model_geometry, 
                    pht,
                    n_stab_W,
                    n_stab_T,
                    δW, 
                    δT,
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
                density_J_parameters,
                spin_J_parameters,
                density_J_factor,
                spin_J_factor,
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
            "sim",
            bin, 
            sim_bin_size,
            measurement_container, 
            simulation_info
        )
    end

    # Record end time for simulation.
    sim_end = time()

    # Record the total simulation time.
    metadata["sim_time"] += sim_end - sim_start

    # Record the total VMC time.
    metadata["vmc_time"] += metadata["opt_time"] + metadata["sim_time"]

    # Calculate the average local acceptance rate.
    metadata["acceptance_rate"] /= ((N_opt + N_sim) * N_equil)

    # Write simulation summary TOML file.
    save_simulation_info(simulation_info, metadata)

    # Process all optimization and simulation measurements.
    # Each observable will be written to CSV files for later processing.
    process_measurements(
        measurement_container, 
        simulation_info, 
        determinantal_parameters, 
        density_J_parameters,
        spin_J_parameters,
        model_geometry
    )
    
    return nothing
end

# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Run the simulation.
    run_hubbard_square_simulation(;
        sID         = parse(Int,     ARGS[1]), 
        Lx          = parse(Int,     ARGS[2]), 
        Ly          = parse(Int,     ARGS[3]),       
        U           = parse(Float64, ARGS[4]), 
        density     = parse(Float64, ARGS[5]), 
        pht         = parse(Bool,    ARGS[6]),
        N_equil     = parse(Int,     ARGS[7]), 
        N_opt       = parse(Int,     ARGS[8]), 
        N_opt_bins  = parse(Int,     ARGS[9]), 
        N_sim       = parse(Int,     ARGS[10]), 
        N_sim_bins  = parse(Int,     ARGS[11])
    )
end
