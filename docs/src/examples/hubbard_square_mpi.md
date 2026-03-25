# 3) Hubbard model on a square lattice with MPI Parallelization

The script that follows this example can be found in the `/example/` directory in the repo. 

This Example will use the same base script as the previous Hubbard model on a square lattice (Example 2), demonstrating how to add parallelization with MPI using the [MPI.jl](https://github.com/JuliaParallel/MPI.jl) package. By this we mean that each MPI process will act as independent walker, running it's own independent VMC simulation.

## Import packages
We now need to import the [MPI.jl](https://github.com/JuliaParallel/MPI.jl) package as well. 

````julia
using LinearAlgebra
using Random
using Printf

# Import MPI
using MPI

using LatticeUtilities
using VariationalMC 
````

## Specify simulation parameters

Here we have introduced the comm argument to the run_simulation function, which is a type exported by the [MPI.jl](https://github.com/JuliaParallel/MPI.jl) package to facilitate communication and synchronization between the different MPI processes.

````julia
# We define a top-level function for running the VMC simulation.
function run_hubbard_square_simulation(
    comm::MPI.Comm;                 # MPI communicator.
    # KEYWORD ARGUMENTS
    sID,                            # Simulation ID.
    L,                              # System size.
    U,                              # Hubbard interaction.
    density,                        # Electron density.
    pht,                            # Whether model is particle-hole transformed. 
    N_equil,                        # Number of equilibration/thermalization updates.
    N_opt,                          # Number of optimization steps.
    N_opt_bins,                     # Number of times bin-averaged measurements are written to file during optimization step.
    N_sim,                          # Number of simulation steps.
    N_sim_bins,                     # Number of times bin-averaged measurements are written to file during simulation step.
    dt          = 0.1,              # Optimization rate.
    dt_J        = 1.0,              # Optional boost in the Jastrow optimization rate.
    η           = 1e-4,             # Optimization stablity factor.
    n_stab_W    = 50,               # Green's function stabilization frequency.
    δW          = 1e-3,             # Maximum allowed error in the Green's function. 
    n_stab_T    = 50,               # Jastrow factor stabilization frequency.
    δT          = 1e-3,             # Maximum allowed error in the Jastrow factor.           
    seed        = abs(rand(Int)),   # Seed for random number generator.
    filepath    ="."                # Filepath to where data folder will be created.
)
````

## Initialize simulation

While choosing which variational parameters to optimize is the same as before, Now when initializing the `SimulationInfo` type, we also need to include the MPI process ID `pID`, which can be retrieved using the `MPI.Comm_rank` function.

We also the initialize_datafolder function such that it takes the comm as the first argument. This ensures that all the MPI processes remained synchronized, and none try proceeding beyond this point until the data folder has been initialized.

````julia
    # Construct the foldername the data will be written.
    df_prefix = @sprintf("hubbard_square_U%.2f_density%.2f_Lx%d_Ly%d_opt", U, density, L, L)

    # Append optimized parameter names to the foldername.
    datafolder_prefix = create_datafolder_prefix(optimize, df_prefix)

    # Get the MPI comm rank, which fixes the processor ID (pID).
    pID = MPI.Comm_rank(comm)

    # Initialize an instance of the SimulationInfo type.
    # This type tracks of where the data is written, as well as 
    # which version of Julia and VariationalMC are used in the script. 
    simulation_info = SimulationInfo(
        filepath = filepath, 
        datafolder_prefix = datafolder_prefix,
        sID = sID,
        pID = pID
    )

    # Initialize the directory the data will be written.
    initialize_datafolder(comm, simulation_info)
````

There are no changes to the simulation script in terms of initializing the simulation metadata, model, model parameters, and measurements. There are also no changes in performing the optimization and simulation steps. Please refer to Example 2 for details. 

## Execute script
Here we first need to initialize MPI using the `MPI.Init` command. Then, we need to make sure to pass the `comm = MPI.COMM_WORLD` to the run_simulation function. At the very end of simulation it is good practice to run the `MPI.Finalize()` function even though it is typically not strictly required.

Only excute if the script is run directly from the command line.

````julia
# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__
    # Initialize MPI.
    MPI.Init()

    # Initialize the MPI communicator.
    comm = MPI.COMM_WORLD

    # Run the simulation.
    run_hubbard_square_simulation(
        comm; 
        sID        = parse(Int,     ARGS[1]), 
        L          = parse(Int,     ARGS[2]), 
        U          = parse(Float64, ARGS[3]),
        density    = parse(Float64, ARGS[4]),
        pht        = parse(Bool,    ARGS[5]),
        N_equil    = parse(Int,     ARGS[6]), 
        N_opt      = parse(Int,     ARGS[7]), 
        N_opt_bins = parse(Int,     ARGS[8]), 
        N_sim      = parse(Int,     ARGS[9]),
        N_sim_bins = parse(Int,     ARGS[10])
    )

    # Finalize MPI.
    MPI.Finalize()
end
````