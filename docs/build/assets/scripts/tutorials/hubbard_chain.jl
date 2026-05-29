using VariationalMC
import VariationalMC.LatticeUtilities as lu

using Random
using Printf

# Top-level function to run the simulation
function run_simulation(;
    # KEYWORD ARGUMENTS
    sID,                    # Simulation ID.
    L,                      # System size.
    U,                      # Hubbard interaction.
    nup,                    # Number of spin-up Fermions.
    ndn,                    # Number of spin-down Fermions.
    ph_transform,           # Whether the model is particle-hole transformed.
    N_equil,                # Number of equilibration/thermalization updates.
    N_opt,                  # Number of optimization steps.
    N_sim,                  # Number of simulation steps.
    N_opt_bins,             # Number of optimization bins.
    N_sim_bins,             # Number of simulation bins.
    dt = 0.1,               # Optimization rate.
    η = 1e-4,               # Optimization stablity factor.
    n_stab_W = 50,          # Green's function stabilization frequency.
    δW = 1e-3,              # Maximum allowed error in the Green's function.
    n_stab_T = 50,          # Jastrow factor stabilization frequency.
    δT = 1e-3,              # Maximum allowed error in the Jastrow factor.
    seed = abs(rand(Int)),  # Seed for random number generator.
    filepath="."            # Filepath to where data folder will be created.
)

    # Construct the foldername the data will be written to.
    datafolder_prefix = @sprintf "hubbard_chain_nup%d_ndn%d_U%.2f_L%d" nup ndn U L

    # Initialize simulation info.
    simulation_info = SimulationInfo(
        filepath = filepath,
        datafolder_prefix = datafolder_prefix,
        sID = sID
    )

    # Initialize the directory the data will be written to.
    initialize_datafolder(simulation_info)

    # Initialize a random number generator.
    rng = Xoshiro(seed)

    # Initialize the metadata dictionary.
    metadata = Dict()

    # Record simulation parameters.
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
    metadata["seed"] = seed

    # Define the unit cell.
    unit_cell = lu.UnitCell(
        lattice_vecs = [[1.0]],
        basis_vecs   = [[0.0]]
    )

    # Define a finite lattice with periodic boundary conditions.
    lattice = lu.Lattice(
        [L],
        [true]
    )

    # Initialize model geometry.
    model_geometry = ModelGeometry(
        unit_cell, lattice
    )

    # Define the nearest-neighbor bond in the x-direction.
    bond_x = lu.Bond(
        orbitals = (1,1),
        displacement = [1]
    )

    # Add this bond definition to the model by adding it to the `model_geometry`.
    add_bond!(model_geometry, bond_x)

    # Determine the density and initialize the configuration of Fermions.
    particle_configuration = ParticleConfiguration(
        nup = nup, ndn = ndn,
        model_geometry = model_geometry,
        ph_transform = ph_transform
    )

    # Set nearest-neighbor hopping amplitude to unity,
    # setting the energy scale in the model.
    t = 1.0

    # Define the non-interacting tight binding model.
    tight_binding_model = TightBindingModel(
        model_geometry = model_geometry,
        μ       =  0.,           # set an initial estimate for the chemical potential.
        t_bonds = [bond_x],      # defines hopping.
        t_mean  = [t]            # defines corresponding mean hopping amplitude.
    )

    # Define the Hubbard interaction in the model.
    hubbard_model = HubbardModel(
        U_orbital   = [1],
        U_mean      = [U]
    )

    # Initialize tight-binding parameters.
    tight_binding_parameters = TightBindingParameters(
        tight_binding_model     = tight_binding_model,
        particle_configuration  = particle_configuration,
        model_geometry          = model_geometry
    )

    # Initialize Hubbard interaction parameters.
    hubbard_parameters = HubbardParameters(
        hubbard_model   = hubbard_model,
        model_geometry  = model_geometry
    )

    # Initialize variational parameters in the determinantal-part of the wavefunction.
    determinantal_parameters = DeterminantalParameters(
        tight_binding_parameters    = tight_binding_parameters,
        model_geometry              = model_geometry,
        rng                         = rng
    )

    # Add the spin-z variational parameter for optimization.
    add_parameter!(
        determinantal_parameters,
        param_name      = "spin-z",
        order_type      = "spin",
        symmetry        = "uniform"
    )

    # Initialize variational parameters in the form of Jastrow pseudopotentials.
    jastrow_parameters = JastrowParameters(
        particle_pair   = "electron-electron",
        order_pair      = "density-density",
        orbitals        = [1],
        optimize        = true,
        model_geometry  = model_geometry,
        rng             = rng
    )

    # Write model summary TOML file specifying Hamiltonian that will be simulated.
    model_summary(
        simulation_info         = simulation_info,
        model_geometry          = model_geometry,
        particle_configuration  = particle_configuration,
        tight_binding_model     = tight_binding_model,
        interactions            = (hubbard_model,),
        parameters              = (determinantal_parameters, jastrow_parameters,)
    )

    # Initialize the container where measurements will be accumulated.
    measurement_container = initialize_measurement_container(
        determinantal_parameters,
        model_geometry,
        (jastrow_parameters,)
    )

    # Initialize the tight-binding model related measurements, namely the hopping energy.
    initialize_measurements!(measurement_container, tight_binding_model)

    # Initialize the Hubbard interaction related measurements.
    initialize_measurements!(measurement_container, hubbard_model)

    # Initialize the optimzer that will be used for Stochastic Reconfiguration.
    optimizer = Optimizer(
        determinantal_parameters = determinantal_parameters,
        η = η, dt = dt,
        jas_parameters = (jastrow_parameters,)
    )

    # Calculate the size of each optimization bin.
    opt_bin_size = div(N_opt, N_opt_bins)

    # Iterate over the optimization bins.
    for bin in 1:N_opt_bins

        # Construct the auxiliary/mean-field Hamiltonian.
        hamiltonian = Hamiltonian(
            tight_binding_parameters    = tight_binding_parameters,
            determinantal_parameters    = determinantal_parameters,
            particle_configuration      = particle_configuration,
            model_geometry              = model_geometry
        )

        # Construct the determinantal/Slater wavefunction.
        detwf = DeterminantalWavefunction(
            determinantal_parameters    = determinantal_parameters,
            hamiltonian                 = hamiltonian,
            particle_configuration      = particle_configuration,
            model_geometry              = model_geometry,
            rng                         = rng
        )

        # Construct the Jastrow factor.
        jfac = JastrowFactor(
            jastrow_parameters      = jastrow_parameters,
            particle_configuration  = particle_configuration,
            model_geometry          = model_geometry
        )

        # Iterate over the length of the bin.
        for n in 1:opt_bin_size
            # Equilibrate/thermalize the system
            for equil in 1:N_equil
                # Attempt to update the fermion configuration.
                accepted = local_fermion_update!(
                    detwf = detwf, jas_factors = (jfac,),
                    tight_binding_parameters = tight_binding_parameters, jas_parameters = (jastrow_parameters,),
                    particle_configuration = particle_configuration,
                    model_geometry = model_geometry,
                    δW = δW, δT = δT,
                    n_stab_W = n_stab_W, n_stab_T = n_stab_T, rng = rng
                )

                # Record the acceptance rate.
                metadata["acceptance_rate"] += accepted
            end

            # Make measurements
            make_measurements!(
                measurement_container,
                detwf = detwf, jas_factors = (jfac,),
                tight_binding_parameters = tight_binding_parameters,
                determinantal_parameters = determinantal_parameters,
                jas_parameters = (jastrow_parameters,),
                coupling_parameters = (hubbard_parameters,),
                particle_configuration = particle_configuration,
                model_geometry = model_geometry, opt_step = true
            )
        end

        # Attempt to update the variational parameters.
        update_optimizer!(
            optimizer = optimizer,
            measurement_container = measurement_container,
            bin_size = opt_bin_size,
            determinantal_parameters = determinantal_parameters,
            jas_parameters = (jastrow_parameters,)
        )

        # Write the measurements.
        write_measurements!(
            measurement_container = measurement_container,
            simulation_info = simulation_info,
            model_geometry = model_geometry,
            bin = bin,
            bin_size = opt_bin_size,
            opt_step = true
        )
    end

    # Calculate the size of each simulation bin.
    sim_bin_size = div(N_sim, N_sim_bins)

    # Iterate over the simulation bins.
    for bin in 1:N_sim_bins

        # Construct the auxiliary/mean-field Hamiltonian.
        hamiltonian = Hamiltonian(
            tight_binding_parameters    = tight_binding_parameters,
            determinantal_parameters    = determinantal_parameters,
            particle_configuration      = particle_configuration,
            model_geometry              = model_geometry
        )

        # Construct the determinantal/Slater wavefunction.
        detwf = DeterminantalWavefunction(
            determinantal_parameters    = determinantal_parameters,
            hamiltonian                 = hamiltonian,
            particle_configuration      = particle_configuration,
            model_geometry              = model_geometry,
            rng                         = rng
        )

        # Construct the Jastrow factor.
        jfac = JastrowFactor(
            jastrow_parameters      = jastrow_parameters,
            particle_configuration  = particle_configuration,
            model_geometry          = model_geometry
        )

        # Iterate over the length of the bin.
        for n in 1:opt_bin_size
            # Equilibrate/thermalize the system
            for equil in 1:N_equil
                # Attempt to update the fermion configuration.
                accepted = local_fermion_update!(
                    detwf = detwf, jas_factors = (jfac,),
                    tight_binding_parameters = tight_binding_parameters, jas_parameters = (jastrow_parameters,),
                    particle_configuration = particle_configuration,
                    model_geometry = model_geometry,
                    δW = δW, δT = δT,
                    n_stab_W = n_stab_W, n_stab_T = n_stab_T, rng = rng
                )

                # Record the acceptance rate.
                metadata["acceptance_rate"] += accepted
            end

            # Make measurements
            make_measurements!(
                measurement_container,
                detwf = detwf, jas_factors = (jfac,),
                tight_binding_parameters = tight_binding_parameters,
                determinantal_parameters = determinantal_parameters,
                jas_parameters = (jastrow_parameters,),
                coupling_parameters = (hubbard_parameters,),
                particle_configuration = particle_configuration,
                model_geometry = model_geometry
            )
        end

        # Write the measurements.
        write_measurements!(
            measurement_container = measurement_container,
            simulation_info = simulation_info,
            model_geometry = model_geometry,
            bin = bin,
            bin_size = sim_bin_size
        )
    end

    # Merge binned optimization data into a single HDF5 file.
    merge_bins(simulation_info, opt = true)

    # Merge binned simulation data into a single HDF5 file.
    merge_bins(simulation_info)

    # Normalize the acceptance rate.
    metadata["acceptance_rate"] /= N_equil * (N_opt + N_sim)

    # Write the simulation summary TOML file.
    save_simulation_info(simulation_info, metadata)

    # Process the simulation results, calculating final error bars for all measurements.
    # writing final statistics to CSV files.
    process_measurements(
        datafolder = simulation_info.datafolder,
        n_bins = N_sim_bins,
        export_to_csv = true,
        scientific_notation = false,
        decimals = 7,
        delimiter = " "
    )

    return nothing
end # end of run_simulation function

# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Run the simulation.
    run_simulation(;
        sID             = parse(Int,     ARGS[1]),  # Simulation ID.
        L               = parse(Int,     ARGS[2]),  # System size.
        U               = parse(Float64, ARGS[3]),  # Hubbard interaction.
        nup             = parse(Int,     ARGS[4]),  # Number of spin-up Fermions.
        ndn             = parse(Int,     ARGS[5]),  # Number of spin-down Fermions.
        ph_transform    = parse(Bool,    ARGS[6]),  # Whether the model is particle-hole transformed.
        N_equil         = parse(Int,     ARGS[7]),  # Number of equilibration/thermalization updates.
        N_opt           = parse(Int,     ARGS[8]),  # Number of optimization steps.
        N_sim           = parse(Int,     ARGS[9]),  # Number of simulation steps.
        N_opt_bins      = parse(Int,     ARGS[10]), # Number of optimization bins.
        N_sim_bins      = parse(Int,     ARGS[11])  # Number of simulation bins.
    )
end
