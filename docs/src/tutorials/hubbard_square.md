```@meta
EditURL = "../../../tutorials/hubbard_square.jl"
```

# 2a) Square Hubbard Model with a Pairing Wavefunction
Download this example as a [Julia script](../assets/scripts/tutorials/hubbard_square.jl).

In this tutorial, we will work through simulating a repulsive Hubbard model on a square lattice
using a variational wavefunction with $s$-wave pairing.
The Hubbard Hamiltonian for a square lattice given by
```math
\begin{align}
\hat{H} = & -t \sum_{\langle i, j \rangle, \sigma} (\hat{c}^{\dagger}_{i, \sigma}, \hat{c}^{\phantom \dagger}_{j, \sigma} + {\rm H.c.})\\
& + U \sum_i \hat{n}_{i,\uparrow}\hat{n}_{i,\downarrow}\\
& - \mu \sum_{i,\sigma} \hat{n}_{i,\sigma},
\end{align}
```
where $\hat{c}^\dagger_{i,\sigma} \ (\hat{c}^{\phantom \dagger}_{i,\sigma})$ creates (annihilates) a spin $\sigma$
electron on site $i$ in the lattice, and $\hat{n}_{i,\sigma} = \hat{c}^\dagger_{i,\sigma} \hat{c}^{\phantom \dagger}_{i,\sigma}$
is the spin-$\sigma$ electron number operator for site $i$. In the above Hamiltonian $U > 0$ controls the strength of the on-site Hubbard repulsion.

## [Import packages](@id hubbard_square_import_packages)
Let us begin by importing [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl), and its relevant submodules.

````julia
using VariationalMC
import VariationalMC.LatticeUtilities as lu
````

The [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl) package re-exports the [LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl.git)
from [SmoQyDQMC](https://github.com/SmoQySuite/SmoQyDQMC.jl.git), which we will use to define the lattice geometry for our model.

We will also  use the Standard Library packages [Random](https://docs.julialang.org/en/v1/stdlib/Random/)
and [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/) for random number generation and C-style string
formatting, respectively.

````julia
using Random
using Printf
````

## Specify simulation parameters

The entire main body of the simulation we will wrapped in a top-level function named `run_simulation`
that will take as keyword arguments various model and simulation parameters that we may want to change.
The function arguments with default values are ones that are typically left unchanged between simulations.
The specific meaning of each argument will be discussed in more detail in later sections of the tutorial.

````julia
# Top-level function to run simulation.
function run_simulation(;
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

## [Initialize simulation](@id hubbard_square_initialize_simulation)
In this first part of the script we name and initialize our simulation, creating the data folder our simulation results will be written to.
This is done by initializing an instances of the [`SimulationInfo`](@ref) type, and then calling the [`initialize_datafolder`](@ref) function.

Note that the `write_bins_concurrent` keyword arguments controls whether or not binned simulation measurement data
is written to HDF5 file during the simulation, or held in memory and only written to file once the simulation is complete.
Here we decide how to set `write_bins_concurrent` based on the system size being simulated.
This is because when performing simulations of small systems that do not take very long, writing data to file too frequently can
sometimes cause network latency problems on clusters and HPC systems. However, for larger systems that take longer to simulate,
you are not limited by file IO frequency but rather by available memory, so writing data to file more frequently is preferred
in these cases.

````julia
    # Construct the foldername the data will be written to.
    datafolder_prefix = @sprintf "hubbard_square_n%.3f_U%.2f_L%d" density U L
````

Initialize simulation info.

````julia
    simulation_info = SimulationInfo(
        filepath = filepath,
        datafolder_prefix = datafolder_prefix,
        write_bins_concurrent = (L > 10),
        sID = sID
    )

    # Initialize the directory the data will be written to.
    initialize_datafolder(simulation_info)
````

## Initialize simulation metadata
In this section of the code we record important metadata about the simulation, including initializing the random number
generator that will be used throughout the simulation.
The important metadata within the simulation will be recorded in the `metadata` dictionary.

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

In the above, `sID` stands for simulation ID, which is used to distinguish simulations that would otherwise be identical i.e. to
distinguish simulations that use the same parameters and are only different in the random seed used to initialize the simulation.
A valid `sID` is any positive integer greater than zero, and is used when naming the data folder the simulation results will be written to.
Specifically, the actual data folder created above will be `"$(filepath)/$(datafolder_prefix)-$(sID)"`.
Note that if you set `sID = 0`, then it will instead be assigned smallest previously unused integer value. For instance, suppose the directory
`"$(filepath)/$(datafolder_prefix)-1"` already exits. Then if you pass `sID = 0` to [`SimulationInfo`](@ref), then the simulation ID
`sID = 2` will be used instead, and a directory `"$(filepath)/$(datafolder_prefix)-2"` will be created.

Another useful resource in the documentation is the [Simulation Output Overview](@ref) page, which describes the output written to the
data folder generated during a [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl) simulation.

## Initialize model

The next step is define the model we wish to simulate.
In this example the relevant model parameters are the Hubbard interaction strength ``U`` (`U`) and lattice size ``L`` (`L`).

First we define the lattice geometry for our model, relying on the
[LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl.git) package to do so.
We define a the unit cell and size of our finite lattice using the [`UnitCell`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.UnitCell)
and [`Lattice`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.Lattice) types, respectively.
Lastly, we define various instances of the [`Bond`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.Bond) type to represent the
the nearest-neighbor bonds.
All of this information regarding the lattice geometry is then stored in an instance of the [`ModelGeometry`](@ref) type.
Further documentation, with usage examples, for [LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl.git) package
can be found [here](https://smoqysuite.github.io/LatticeUtilities.jl/stable/).

We also initialize the Canonical Ensemble particle configuration that lives on the lattice. In this case, given the total Fermion density that we want to live on the
lattice, [`ParticleConfiguration`](@ref) will determine the correct total number of particles and the correct number of spin-up and spin-down particles.

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
````

Next we specify the non-interacting tight-binding term in our Hamiltonian with the [`TightBindingModel`](@ref) type.

````julia
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
````

Then, we define the Hubbard interaction with the [`HubbardModel`](@ref) type.

````julia
    # Define the Hubbard interaction in the model.
    hubbard_model = HubbardModel(
        U_orbital   = [1],
        U_mean      = [U]
    )
````

## Initialize model parameters
The next step is to initialize our model parameters given the size of our finite lattice.
To clarify, both the [`TightBindingModel`](@ref) and [`HubbardModel`](@ref) types are agnostic to the size of the lattice being simulated,
defining the model in a translationally invariant way.
To do so we need to initialize an instance of the [`TightBindingParameters`](@ref) and [`HubbardParameters`](@ref) types.

Note that the initial guess for the chemical potential given in [`TightBindingModel`](@ref) will be overriden here with the exact one.

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
````

Next, we add variational parameters to our model. This first set represents parameters in the determinantal or Slater part of the wavefunction.
By initializing the [`DeterminantalParameters`](@ref) type, a base number of parameters are automatically added to the model, which represent
the minimal possible orders. These include charge order in the form of `"mu"` and `"density"`, and spin order in the form of `"spin-x"` and
`"spin-z"`. To flag parameters for optimization, they will need to be added using the [`add_parameter!`](@ref) method. In this case, we will
add an additional parameter in the form of $s$-wave pairing.

We also add variational parameters in the form of Jastrow pseudopotentials, in this case of the `"density-density` type.

````julia
    # Initialize variational parameters in the determinantal-part of the wavefunction.
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
````

Lastly, the [`model_summary`](@ref) function is used to write a `model_summary.toml` file,
specifying the parameters which will be used in the Hamiltonian that will be simulated.

````julia
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

## [Initialize measurements](@id hubbard_square_initialize_measurements)
Having initialized both our model and the corresponding model parameters,
the next step is to initialize the various measurements we want to make during our VMC simulation.

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

## [Set up VMC simulation](@id hubbard_square_setup_vmc)
This section of the code sets up the VMC simulation by initializing the optimizer that will be used during Stochastic Reconfiguration.

````julia
    # Initialize the optimzer that will be used for Stochastic Reconfiguration.
    optimizer = Optimizer(
        determinantal_parameters = determinantal_parameters,
        η = η, dt = dt,
        jas_parameters = (jastrow_parameters,)
    )
````

## [Optimize the parameters](@id hubbard_square_optimize)
This next section of the code performs updates to optimize the variational parameters.
The structure of this function should be fairly intuitive, mainly consisting of a loop,
inside which we perform updates to the particle configuration.

Before the start of the each update loop, we initialize the auxiliary Hamiltonian which is then used to
construct the variational wavefunction using the update parametrs from the previous loop.

In the block of code below, `N_equil` is the number of local updates performed to equilibrate/thermalize the system.
Unlike a DQMC simulation, where thermalization updates can be done before any measurements are made, VMC requires the system be
thermalized/equilibrated before each measurement. This is because the VMC optimizer is more prone to getting stuck in a
local minimum and by equilibrating before each measurement, this reduces this risk.

The total number of optimization iterations `N_opt` and the total number of bins `N_opt_bins` is automatically
translated into the correct bin lengths. The bin-averaged measurements are written to file once `opt_bin_size` measurements are
accumulated using the [`write_measurements!`](@ref) function.

Note that we set `opt_step = true` in the [`make_measurements!`](@ref) and [`write_measurements!`](@ref) functions to ensure that the
logarithmic derivatives are being measured that to let the simulation know that the optimizer will handle normalization during this stage
of the VMC simulation.

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

            # Make measurements.
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

## [Simulate the system](@id hubbard_square_simulate_system)
Now that the variational parameters have been optimzied, we simulate the system using these parameters to obtain
clean measurements of any observables. The structure of this section is similar to the optimization block, sans
the actual optimizer.

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
At this point the simulation is essentially complete, with all updates and measurements having been performed.
However, the binned measurement data resides in many separate HDF5 files currently.
Here we will merge these separate HDF5 files into a single file containing all the binned data
using the [`merge_bins`](@ref) function. This must be done for both the optimization and simulation bins.

````julia
    # Merge binned optimization data into a single HDF5 file.
    merge_bins(simulation_info, opt = true)

    # Merge binned simulation data into a single HDF5 file.
    merge_bins(simulation_info)
````

## Record simulation metadata
Next, we want to calculate the final acceptance rate for the Monte Carlo updates,
as well as write the simulation metadata to file, including the contents of the `metadata` dictionary.
This is done using the [`save_simulation_info`](@ref) function.

````julia
    # Normalize the acceptance rate.
    metadata["acceptance_rate"] /= N_equil * (N_opt + N_sim)

    # Write the simulation summary TOML file.
    save_simulation_info(simulation_info, metadata)
````

## [Post-process results](@id hubbard_square_process_results)
In this final section of code we post-process the binned data.
This includes calculating the final estimates for the mean and error of all measured observables,
which will be written to an HDF5 file using the [`process_measurements`](@ref) function.
Inside this function the binned data gets further re-binned into `n_bins`,
where `n_bins` is any positive integer satisfying the constraints `(N_bins ≥ n_bin)` and `(N_bins % n_bins == 0)`.
Note that the [`process_measurements`](@ref) function has many additional keyword arguments that can be used to control the output.
For instance, in this example in addition to writing the statistics to an HDF5 file, we also export the statistics to CSV files
by setting `export_to_csv = true`, with additional keyword arguments controlling the formatting of the CSV files.
Again, for more information on how to interpret the output refer the [Simulation Output Overview](@ref) page.

````julia
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
````

## Execute script

VMC simulations are typically run from the command line as jobs on a computing cluster.
With this in mind, the following block of code only executes if the Julia script is run from the command line,
also reading in additional command line arguments.

````julia
# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Run the simulation, reading in command line arguments.
    run_simulation(;
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
end
````

For instance, the command
```
> julia hubbard_square.jl 1 4 2.0 0.875 true 200 3000 6000 300 600
```
runs a VMC simulation of a $N = 4 \times 4$ square Hubbard model at $1/8$-doping
with interaction strength $U = 2.0$. Becuase we added a pairing wavefunction, we ensure that `ph_transform = true`.
In the VMC simulation, ``200`` sweeps through the lattice are be performed to equilibrate the system,
while sweeping through ``3,000`` optimiztion iterations and ``6,000`` simulation iterations.
Measurements are then written to file ``300`` times during the optimization and ``600`` times
during simulation, with each bin containing the average of ``3,000/300 = 10`` (``6,000/600 = 10``)
sequential measurements.

