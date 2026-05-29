```@meta
EditURL = "../../../examples/hubbard_ladders_mpi.jl"
```

# Interacting Hubbard Ladders with MPI Parallelization
Download this example as a [Julia script](../assets/scripts/examples/hubbard_ladders_mpi.jl).

This is a ready-to-run script for simulating a two-dimensional (2D) Hubbard model
on a rectangular lattice goemetry (ladders) with a stripe wavefunction using MPI parallelization.

````julia
# Import the required packages
using VariationalMC
import VariationalMC.LatticeUtilities as lu

using Random
using Printf
using MPI

# Top-level function to run the simulation
function run_simulation(
    comm::MPI.Comm;         # MPI communicator.
    # KEYWORD ARGUMENTS
    sID,                    # Simulation ID.
    Lx,                     # System size in the x-direction.
    Ly,                     # System size in the y-direction.
    t₁′,                    # Next-nearest neighbor hopping amplitude in ladder 1.
    t₂′,                    # Next-nearest neighbor hopping amplitude in ladder 2.
    U,                      # Hubbard interaction.
    V,                      # Extended Hubbard interaction.
    density,                # Density of Fermions.
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
    datafolder_prefix = @sprintf "hubbard_ladder_n%.3f_U%.2f_Lx%d_Ly%d" density U Lx Ly

    # Get MPI process ID.
    pID = MPI.Comm_rank(comm)

    # Initialize simulation info.
    simulation_info = SimulationInfo(
        filepath = filepath,
        datafolder_prefix = datafolder_prefix,
        sID = sID,
        pID = pID
    )

    # Initialize the directory the data will be written to.
    initialize_datafolder(simulation_info)

    # Initialize a random number generator.
    rng = Xoshiro(seed)
````

#Initialize the metadata dictionary.

````julia
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
        lattice_vecs = [[1.0, 0.0], [0.0, 1.0]],
        basis_vecs   = [[0.0, 0.0], [1.0, 0.0]]
    )

    # Define a finite lattice with periodic boundary conditions in the x-direction.
    lattice = lu.Lattice(
        [Lx, Ly],
        [true, false]       # TODO: check that doing this is OK
    )

    # Initialize model geometry.
    model_geometry = ModelGeometry(
        unit_cell, lattice
    )

    # Define the nearest-neighbor x-bond in ladder 1.
    bond_x1 = lu.Bond(
        orbitals = (1,1),
        displacement = [1,0]
    )
    bond_x1_id = add_bond!(model_geometry, bond_x1)

    # Define the nearest-neighbor y-bond in ladder 1.
    bond_y1 = lu.Bond(
        orbitals = (1,1),
        displacement = [0,1]
    )
    bond_y1_id = add_bond!(model_geometry, bond_y1)

    # Define the postive next-nearest neighbor in ladder 1.
    bond_xy1 = lu.Bond(
        orbitals = (1,1),
        displacement = [1,1]
    )
    bond_xy1_id = add_bond!(model_geometry, bond_xy1)

    # Define the negative next-nearest neighbor in ladder 1.
    bond_yx1 = lu.Bond(
        orbitals = (1,1),
        displacement = [1,-1]
    )
    bond_yx1_id = add_bond!(model_geometry, bond_yx1)

    # Define the nearest-neighbor x-bond in ladder 2.
    bond_x2 = lu.Bond(
        orbitals = (2,2),
        displacement = [1,0]
    )
    bond_x2_id = add_bond!(model_geometry, bond_x2)

    # Define the nearest-neighbor y-bond in ladder 2.
    bond_y2 = lu.Bond(
        orbitals = (2,2),
        displacement = [0,1]
    )
    bond_y2_id = add_bond!(model_geometry, bond_y2)

    # Define the postive next-nearest neighbor in ladder 2.
    bond_xy2 = lu.Bond(
        orbitals = (2,2),
        displacement = [1,1]
    )
    bond_xy2_id = add_bond!(model_geometry, bond_xy2)

    # Define the negative next-nearest neighbor in ladder 2.
    bond_yx2 = lu.Bond(
        orbitals = (2,2),
        displacement = [1,-1]
    )
    bond_yx2_id = add_bond!(model_geometry, bond_yx2)

    # Define a nearest-neighbor interaction bond between the ladders 1 and 2.
    bond_int1 = lu.Bond(
        orbitals = (1,2),
        displacement = [-1,-1]
    )
    bond_int1_id = add_bond!(model_geometry, bond_int1)

    # Define a nearest-neighbor interaction bond between the ladders 1 and 2.
    bond_int2 = lu.Bond(
        orbitals = (1,2),
        displacement = [0,-1]
    )
    bond_int2_id = add_bond!(model_geometry, bond_int2)

    # Determine the density and initialize the configuration of Fermions.
    particle_configuration = ParticleConfiguration(
        density,
        model_geometry = model_geometry,
        ph_transform = ph_transform
    )

    # Nearest-neighbor hopping amplitude for the first ladder.
    t₁ = 1.0

    # Nearest-neighbor hopping ampplitude for the second ladder.
    t₂ = 1.0

    # Define the non-interacting tight binding model.
    tight_binding_model = TightBindingModel(
        model_geometry = model_geometry,
        μ       =  0., # set initial estimate for the chemical potential
        t_bonds = [bond_x1, bond_y1, bond_x2, bond_y2, bond_xy1, bond_yx1, bond_xy2, bond_yx2], # defines hopping
        t_mean  = [t₁, t₁, t₂, t₂, t₁′, t₁′, t₂′, t₂′] # defines corresponding mean hopping amplitude
    )

    # Define the Hubbard interaction in the model.
    hubbard_model = HubbardModel(
        U_orbital   = [1, 2],
        U_mean      = [U, U]
    )

    # Define the Hubbard interaction between the ladders.
    extended_hubbard_model = ExtendedHubbardModel(
        model_geometry = model_geometry,
        V_bond = [bond_int1, bond_int2],
        V_mean = [V, V]
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

    # Initialize extended Hubbard interaction parameters.
    extended_hubbard_parameters = ExtendedHubbardParameters(
        extended_hubbard_model = extended_hubbard_model,
        model_geometry = model_geometry
    )

    # Initialize variational parameters in the determinantal-part of the wavefunction.
    determinantal_parameters = DeterminantalParameters(
        tight_binding_parameters    = tight_binding_parameters,
        model_geometry              = model_geometry,
        rng                         = rng
    )

    # Add charge stripe parameters for optimization.
    add_parameter!(
        determinantal_parameters,
        param_name  = "density",
        order_type  = "charge",
        symmetry    = "site-dependent",
        optimize    = true,
        orbital     = [1, 2],
        model_geometry = model_geometry
    )

    # Add spin stripe parameters for optimization.
    add_parameter!(
        determinantal_parameters,
        param_name  = "spin-z",
        order_type  = "spin",
        symmetry    = "site-dependent",
        optimize    = true,
        orbital     = [1, 2],
        model_geometry = model_geometry
    )

    # Initialize variational parameters in the form of Jastrow pseudopotentials.
    jastrow_parameters = JastrowParameters(
        particle_pair   = "electron-electron",
        order_pair      = "spin-spin",
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
        interactions            = (hubbard_model, extended_hubbard_model,),
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

    # Initialize the extended Hubbard interaction related measurements.
    initialize_measurements!(measurement_container, extended_hubbard_model)

    # Initialize equal-time charge correlation measurements.
    initialize_correlation_measurements!(
        measurement_container = measurement_container,
        model_geometry = model_geometry,
        correlation = "density",
        pairs = [(1,2)]
    )

    # Initialize equal-time spin correlation measurements.
    initialize_correlation_measurements!(
        measurement_container = measurement_container,
        model_geometry = model_geometry,
        correlation = "spin-z",
        pairs = [(1,2)]
    )

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

            # Make measurements.
            make_measurements!(
                measurement_container,
                detwf = detwf, jas_factors = (jfac,),
                tight_binding_parameters = tight_binding_parameters,
                determinantal_parameters = determinantal_parameters,
                jas_parameters = (jastrow_parameters,),
                coupling_parameters = (hubbard_parameters, extended_hubbard_parameters,),
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
            # Equilibrate/thermalize the system.
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
                coupling_parameters = (hubbard_parameters, extended_hubbard_parameters,),
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
end # end of `run_simulation` function


# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Initialize MPI
    MPI.Init()

    # Initialize the MPI communicator.
    comm = MPI.COMM_WORLD

    # Run the simulation, reading in command line arguments.
    run_simulation(
        comm;
        sID             = parse(Int,     ARGS[1]),
        Lx              = parse(Int,     ARGS[2]),
        Ly              = parse(Int,     ARGS[3]),
        t₁′             = parse(Float64, ARGS[4]),
        t₂′             = parse(Float64, ARGS[5]),
        U               = parse(Float64, ARGS[6]),
        V               = parse(Float64, ARGS[7]),
        density         = parse(Float64, ARGS[8]),
        ph_transform    = parse(Bool,    ARGS[9]),
        N_equil         = parse(Int,     ARGS[10]),
        N_opt           = parse(Int,     ARGS[11]),
        N_sim           = parse(Int,     ARGS[12]),
        N_opt_bins      = parse(Int,     ARGS[13]),
        N_sim_bins      = parse(Int,     ARGS[14])
    )

    # Finalize MPI.
    MPI.Finalize()
end
````

