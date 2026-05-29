```@meta
EditURL = "../../../tutorials/hubbard_square_mpi.jl"
```

# 2b) Square Hubbard Model with MPI Parallelization
Download this example as a [Julia script](../assets/scripts/tutorials/hubbard_square_mpi.jl).

This tutorial will build on the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial,
demonstrating how to add parallelization with MPI using the [MPI.jl](https://github.com/JuliaParallel/MPI.jl.git)
package. By this we mean that each MPI process will act as independent walker, running it's own independent VMC simulation,
with the final reported estimates for measured quantities being the average across all walkers.

The exposition in this tutorial will focus on the changes that need to be made to the [2a) Square Hubbard Model with a Pairing Wavefunction](@ref)
tutorial to introduce MPI parallelization, omitting a more comprehensive discussion of other parts of the code that
were included in the previous tutorial.

## Import Packages
We now need to import the [MPI.jl](https://github.com/JuliaParallel/MPI.jl.git) package as well.

````julia
using VariationalMC
import VariationalMC.LatticeUtilities as lu

using Random
using Printf
using MPI
````

## Specify simulation parameters
Here we have introduced the `comm` argument to the `run_simulation` function, which is a type exported by the
[MPI.jl](https://github.com/JuliaParallel/MPI.jl.git) package to facilitate communication and synchronization
between the different MPI processes.

````julia
# Top-level function to run simulation.
function run_simulation(
    comm::MPI.Comm;         # MPI communicator.
    # KEYWORD ARGUMENTS
    sID,                    # Simulation ID.
    L,                      # System size.
    U,                      # Hubbard interaction.
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
````

## Initialize simulation
Now when initializing the [`SimulationInfo`](@ref) type, we also need to include the
MPI process ID `pID`, which can be retrieved using the
[`MPI.Comm_rank`](https://juliaparallel.org/MPI.jl/stable/reference/comm/#MPI.Comm_rank)
function.

We also the [`initialize_datafolder`](@ref) function such that it takes the `comm` as the
first argument. This ensures that all the MPI processes remained synchronized, and none
try proceeding beyond this point until the data folder has been initialized.

````julia
    # Construct the foldername the data will be written to.
    datafolder_prefix = @sprintf "hubbard_square_n%.3f_U%.2f_L%d" density U L

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
````

## Initialize simulation metadata
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
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
````

## Initialize model
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
    # Define the unit cell.
    unit_cell = lu.UnitCell(
        lattice_vecs = [[1.0, 0.0], [0.0, 1.0]],
        basis_vecs   = [[0.0, 0.0]]
    )

    # Define a finite lattice with periodic boundary conditions.
    lattice = lu.Lattice(
        [L, L],
        [true, true]
    )

    # Initialize model geometry.
    model_geometry = ModelGeometry(
        unit_cell, lattice
    )

    # Define the nearest-neighbor bond in the x-direction.
    bond_x = lu.Bond(
        orbitals = (1,1),
        displacement = [1,0]
    )

    # Add this bond definition to the model by adding it to the `model_geometry`.
    bond_x_id = add_bond!(model_geometry, bond_x)

    # Define the nearest-neighbor bond in the y-direction.
    bond_y = lu.Bond(
        orbitals = (1,1),
        displacement = [0,1]
    )

    # Add this bond definition to the model by adding it to the `model_geometry`.
    bond_y_id = add_bond!(model_geometry, bond_y)

    # Determine the density and initialize the configuration of Fermions.
    particle_configuration = ParticleConfiguration(
        density,
        model_geometry = model_geometry,
        ph_transform = ph_transform
    )

    # Set nearest-neighbor hopping amplitude to unity,
    # setting the energy scale in the model.
    t = 1.0

    # Define the non-interacting tight binding model.
    tight_binding_model = TightBindingModel(
        model_geometry = model_geometry,
        μ       =  0.,                    # set an initial estimate for the chemical potential.
        t_bonds = [bond_x, bond_y],       # defines hopping.
        t_mean  = [t, t]                  # defines corresponding mean hopping amplitude.
    )

    # Define a Hubbard model
    hubbard_model = HubbardModel(
        U_orbital   = [1],
        U_mean      = [U]
    )
````

## Initialize model parameters
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
    # Initialize tight binding parameters.
    tight_binding_parameters = TightBindingParameters(
        tight_binding_model     = tight_binding_model,
        particle_configuration  = particle_configuration,
        model_geometry          = model_geometry
    )

    # Initialize Hubbard parameters.
    hubbard_parameters = HubbardParameters(
        hubbard_model   = hubbard_model,
        model_geometry  = model_geometry
    )

    # Initialize determinantal parameters.
    determinantal_parameters = DeterminantalParameters(
        tight_binding_parameters    = tight_binding_parameters,
        model_geometry              = model_geometry,
        rng                         = rng
    )

    # Add the onsite s-wave variational parameter for optimization.
    add_parameter!(
        determinantal_parameters,
        param_name      = "s-wave",
        order_type      = "pair",
        symmetry        = "uniform",
        optimize        = true,
        orbital         = [1]
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
````

## Initialize measurements
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
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
````

## Set up VMC simulation
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
    # Initialize the optimzer that will be used for Stochastic Reconfiguration.
    optimizer = Optimizer(
        determinantal_parameters = determinantal_parameters,
        η = η, dt = dt,
        jas_parameters = (jastrow_parameters,)
    )
````

## Optimize the parameters
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
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
````

## Simulate the system
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
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

            # Make measurements.
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
````

## Merge binned data
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
    # Merge binned optimization data into a single HDF5 file.
    merge_bins(simulation_info, opt = true)

    # Merge binned simulation data into a single HDF5 file.
    merge_bins(simulation_info)
````

## Record simulation metadata
No changes need to made to this section of the code from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial.

````julia
    # Normalize the acceptance rate.
    metadata["acceptance_rate"] /= N_equil * (N_opt + N_sim)

    # Write the simulation summary TOML file.
    save_simulation_info(simulation_info, metadata)
````

## Post-process results
The main change we need to make from the previous [2a) Square Hubbard Model with a Pairing Wavefunction](@ref) tutorial is to call
the [`process_measurements`](@ref) and functions
such that the first argument is the `comm` object, thereby ensuring a parallelized version of each method is called.

Process the simulation results, calculating final error bars for all measurements.
writing final statistics to CSV files.

````julia
    process_measurements(
        comm;
        datafolder = simulation_info.datafolder,
        n_bins = N_sim_bins,
        export_to_csv = true,
        scientific_notation = false,
        decimals = 7,
        delimiter = " "
    )

    return nothing
end # end of run_simulation function
````

## Execute script
Here we first need to initialize MPI using the
[`MPI.Init`](https://juliaparallel.org/MPI.jl/stable/reference/environment/#MPI.Init) command.
Then, we need to make sure to pass the `comm = MPI.COMM_WORLD` to the `run_simulation` function.
At the very end of simulation it is good practice to run the `MPI.Finalize()` function even though
it is typically not strictly required.

````julia
# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Initialize MPI
    MPI.Init()

    # Initialize the MPI communicator.
    comm = MPI.COMM_WORLD

    # Run the simulation, reading in command line arguments.
    run_simulation(
        comm;
        sID             = parse(Int,     ARGS[1]),  # Simulation ID.
        L               = parse(Int,     ARGS[2]),  # System size.
        U               = parse(Float64, ARGS[3]),  # Hubbard interaction.
        density         = parse(Float64, ARGS[4]),  # Density of Fermions.
        ph_transform    = parse(Bool,    ARGS[5]),  # Whether the model is particle-hole transformed.
        N_equil         = parse(Int,     ARGS[6]),  # Number of equilibration/thermalization updates.
        N_opt           = parse(Int,     ARGS[7]),  # Number of optimization steps.
        N_sim           = parse(Int,     ARGS[8]),  # Number of simulation steps.
        N_opt_bins      = parse(Int,     ARGS[9]),  # Number of optimization bins.
        N_sim_bins      = parse(Int,     ARGS[10])  # Number of simulation bins.
    )

    # Finalize MPI.
    MPI.Finalize()
end
````

Here is an example of what the command to run this script might look like:
```bash
mpiexecjl -n 8 julia hubbard_square_mpi.jl 1 4 2.0 0.875 true 200 3000 6000 300 600
```
This will initialize 8 MPI processes, each running and independent simulation using a different random seed
the final results arrived at by averaging over all 16 walkers.
Here `mpiexecjl` is the MPI executable that can be easily install using the directions
found [here](https://juliaparallel.org/MPI.jl/stable/usage/#Julia-wrapper-for-mpiexec) in the
[MPI.jl](https://github.com/JuliaParallel/MPI.jl) documentation. However, you can substitute a
different MPI executable here if one is already configured on your system.

Also, when submitting jobs via [SLURM](https://slurm.schedmd.com/documentation.html)
on a High-Performance Computing (HPC) cluster, if a default MPI executable
is already configured on the system, as is frequently the case, then the script can likely be run inside the
`*.sh` job file using the [`srun`](https://slurm.schedmd.com/srun.html) command:
```bash
srun julia hubbard_square_mpi.jl 1 4 2.0 0.875 true 200 3000 6000 300 600
```
The `srun` command should automatically detect the number of available cores requested by the job and run
the script using the MPI executable with the appropriate number of processes.

