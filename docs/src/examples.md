# 1) Hubbard model on a 1D chain (without Jastrow factor)

Download this example as a [Julia script](https://github.com/atanjaro/VariationalMC.jl/tree/e07eb401cbe7123b875da159f0372c9001c6cf41/examples).

In this example, we will work through simulating a repulsive Hubbard model on a one dimensional (1D) chain.
The Hubbard hamiltonian in 1D is given by
```math
\hat{H} = -t \sum_{i, \sigma} (\hat{c}^{\dagger}_{i, \sigma}, \hat{c}^{\phantom \dagger}_{i+1, \sigma} + {\rm H.c.})
-t^{\prime} \sum_{i, \sigma} (\hat{c}^{\dagger}_{i, \sigma}, \hat{c}^{\phantom \dagger}_{i+2, \sigma} + {\rm h.c.})
+ U \sum_i \hat{n}_{i, \uparrow}\hat{n}_{i, \downarrow} 
- \mu \sum_{i, \sigma} \hat{n}_{i, \sigma},
```
where ``\hat{c}^\dagger_{i, \sigma} \ (\hat{c}^{\phantom \dagger}_{i, \sigma})`` creates (annihilates) a spin ``\sigma``
electron on site ``i`` in the lattice, and ``\hat{n}_{i, \sigma} = \hat{c}^\dagger_{i, \sigma} \hat{c}^{\phantom \dagger}_{i, \sigma}``
is the spin-``\sigma`` electron number operator for site ``i``. In the above Hamiltonian ``(t^{\prime}) \ t`` is the (next-) nearest-neighbor hopping amplitude and ``U > 0`` controls the strength of the on-site Hubbard repulsion.
Lastly, we note the system is half-filled and particle-hole symmetric when the next-nearest-neighbor hopping amplitude
and the chemical potential is zero ``(t^{\prime} = \mu = 0.0)``.

## Import packages
We begin by importing the necessary packages, including the Standard Library packages [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/), [Random](https://docs.julialang.org/en/v1/stdlib/Random/), and [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/) for perfroming linear algebra operations, random number generation, and C-style string formatting, respectively.

````julia
using Linear Algebra
using Random
using Printf
````

Next, we use [LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl) which exports a suite of types and methods useful for defining arbitrary lattice geometries, and the construction of neighbor tables. Finally, we import the [VariationalMC](https://github.com/atanjaro/VariationalMC.jl) package.

````julia
using LatticeUtilities
using VariationalMC
````

## Specify simulation parameters

The main body of the simulation is wrapped in a top-level function named `run_hubbard_chain_simulation`that will take as keyword arguments various model and simulation parameters that we may want to change. The specific meaning of each argument will be discussed in more detail in later sections of the tutorial.

````julia
# We define a top-level function for running the VMC simulation.
function run_hubbard_chain_simulation(;
    # KEYWORD ARGUMENTS
    sID,                    # Simulation ID.
    L,                      # System size.
    U,                      # Hubbard interaction.
    nup,                    # number of spin-up electrons.
    ndn,                    # number of spin-down electrons.
    pht,                    # Whether model is particle-hole transformed. 
    N_equil,                # Number of equilibration/thermalization updates.
    N_opt,                  # Number of optimization steps.
    N_opt_bins,             # Number of times bin-averaged measurements are written to file during optimization step.
    N_sim,                  # Number of simulation steps.
    N_sim_bins,             # Number of times bin-averaged measurements are written to file during simulation step.
    dt = 0.03,              # Optimization rate.
    η = 1e-4,               # Optimization stablity factor.
    n_stab_W = 50,          # Green's function stabilization frequency.
    δW = 1e-3,              # Maximum allowed error in the Green's function.          
    seed = abs(rand(Int)),  # Seed for random number generator.
    filepath="."            # Filepath to where data folder will be created.
)
````

## Initialize simulation
In this initial part of the script, we will specify what parameters in the trial wavefunction will be optimized and name our simulation. The `NamedTuple` called `optimized` will contain all possible parameters that are automatically initialized in the wavefunction. Setting a flag to `true` will cause that parameter value (or values) to change during the optimization. We also specify whether the model is particle-hole transformed and is required to be `true` if pairing symmetery is being added to the wavefunction. 

````julia
    # Select which parameters in the variational wavefunction will be optimized.
    optimize = (
        # Spin-x or in-plane magnetization
        Δ_sx = false,
        # Spin-z or out-of-plane magnetization
        Δ_sz = true,
        # (BCS) Chemical potential
        μ = false,
        # Uniform charge density or Charge Density Wave
        Δ_cdw = false,
        # Density-density Jastrow factor
        density_J = false
    )
````

The datafolder is created by initializing an instances of the `SimulationInfo` type, and then calling the `initialize_datafolder` function. Iin addition, we initialize the random number generator that will be used throughout the simulation. 


````julia
    # Construct the foldername the data will be written.
    df_prefix = @sprintf("hubbard_chain_U%.2f_nup%.2f_ndn%.2f_L%d_opt", U, nup, ndn, L)

    # Append optimized parameter names to the foldername.
    datafolder_prefix = create_datafolder_prefix(optimize, df_prefix)

    # Initialize an instance of the SimulationInfo type.
    simulation_info = SimulationInfo(
        filepath = filepath, 
        datafolder_prefix = datafolder_prefix,
        sID = sID
    )

    # Initialize the directory the data will be written.
    initialize_datafolder(simulation_info)

    # Initialize a random number generator that will be used throughout the simulation.
    rng = Xoshiro(seed)
````

In addition, we can calculate the length of each bin by dividing the number of iterations/steps by the number of bins. The bin size is the number of measurements that are averaged over each time data is written during either the optimization or simulation steps.

````julia
    # Calculate optimization bins size.
    opt_bin_size = div(N_opt, N_opt_bins) 

    # Calculate simulation bins size.
    sim_bin_size = div(N_sim, N_sim_bins)
````

## Initialize simulation metadata
In this section of the code we record important metadata about the simulation, including initializing the random number
generator that will be used throughout the simulation.
The important metadata within the simulation will be recorded in the `metadata` dictionary.

````julia
    # Record simulation metadata.
    metadata = Dict()
    metadata["N_equil"] = N_equil
    metadata["N_opt"] = N_opt
    metadata["N_sim"] = N_sim
    metadata["N_opt_bins"] = N_opt_bins
    metadata["N_sim_bins"] = N_sim_bins
    metadata["δW"] = δW
    metadata["n_stab_W"] = n_stab_W
    metadata["dt"] = dt 
    metadata["acceptance_rate"] = 0.0
    metadata["opt_time"] = 0.0
    metadata["sim_time"] = 0.0
    metadata["vmc_time"] = 0.0
````

In the above, `sID` stands for simulation ID, which is used to distinguish simulations that would otherwise be identical i.e. to
distinguish simulations that use the same parameters and are only different in the random seed used to initialize the simulation.
A valid `sID` is any positive integer greater than zero, and is used when naming the data folder the simulation results will be written to.
Specifically, the actual data folder created above will be `"$(filepath)/$(datafolder_prefix)-$(sID)"`.
Note that if you set `sID = 0`, then it will instead be assigned smallest previously unused integer value. For instance, suppose the directory
`"$(filepath)/$(datafolder_prefix)-1"` already exits. Then if you pass `sID = 0` to `SimulationInfo`, then the simulation ID
`sID = 2` will be used instead, and a directory `"$(filepath)/$(datafolder_prefix)-2"` will be created.

## Initialize model

The next part of the script defines the model that we will simulate.
First we define the lattice geometry for our model, relying on the
[LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl.git) package to do so.
We define a the unit cell and size of our finite lattice using the [`UnitCell`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.UnitCell)
and [`Lattice`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.Lattice) types, respectively.
Lastly, we define various instances of the [`Bond`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.Bond) type to represent the
the nearest-neighbor and next-nearest-neighbor bonds.
All of this information regarding the lattice geometry is then stored in an instance of the `ModelGeometry` type.

````julia
    # Initialize an instance of the UnitCell type.
    unit_cell = UnitCell(
        lattice_vecs = [[1.0]],
        basis_vecs   = [[0.0]]
    )

    # Initialize an instance of the Lattice type.
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
    model_geometry = ModelGeometry(
        unit_cell, 
        lattice, 
        bonds
    )
````

## Initialize model parameters

The next step is to initialize our model parameters, which includes calculating the particle density in the canonical ensemble. The `get_particle_density` method includes mutliple ways of specifying the particle information. Here, we have chosen to pass the number of spin-up and spin-down particles. The function then computes the corresponding total density, particle number, and electron number.

````julia
    # Determine the total particle density in the canonical ensemble. 
    (density, Np, Ne, nup, ndn) = get_particle_density(nup, ndn, model_geometry, pht) 
````

We then specify parameters for our tight binding model and initializes all of the parameters that are in the determinantal part of the trial wavefunction. The `model_summary` function is used to write a `model_summary.toml` file, completely specifying the Hamiltonian that will be simulated. Lastly, we can also set an initial particle configuration, if we have one. If an empty array is provided instead, a random configuration will be given at the start of the simulation. 

````julia
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

    # Write model summary TOML file specifying the Hamiltonian that will be simulated.
    model_summary(
        simulation_info, 
        determinantal_parameters, 
        pht, 
        model_geometry, 
        tight_binding_model, 
        U
    )

    # Initialize the (fermionic) particle configuration.
    pconfig = Int[]
````

## Initialize measurements
Finally, we initialize the mesaurement container, which accumulates the sums of measurements made in a bin. The standard measurements are of the local energy `local_energy`, double occupancy `double_occ`, particle confiugurations `pconfig`, and variational parameters `parameters`. There is also the option to add additional observables to measure through the `initialize_simulation_measurement!`and `initialize_correlation_measurement!` methods. 

````julia
    # Initialize the container that measurements will be accumulated into.
    measurement_container = initialize_measurement_container(
        N_opt, 
        opt_bin_size, 
        N_sim, 
        sim_bin_size,
        determinantal_parameters,
        model_geometry
    )
````

Here, we have added measurements of the z-component of the spin ``S_z``. After this, the `initialize_measurement_directories` can now be used used to initialize the various subdirectories in the data folder that the measurements will be written to.
For more information, please refer to the Simulation Output Overview page. 

````julia
    # Add local Sz measurements.
    initialize_simulation_measurement!(
        "local",
        "spin-z",
        measurement_container,
        model_geometry
    )

    # Initialize the sub-directories to which the various measurements will be written.
    initialize_measurement_directories(
        simulation_info, 
        measurement_container
    )
````

## Optimize the variational parameters
Now that we have set-up the VMC simulation, we can begin optimizing the variational parameters. At the start of each bin, we initialize the variational wavefunction. On the first iteration, if no initial configuration is specified, a random one will be generated. In this example, we are neglecting the inclusion of the Jastrow correlation factor here for demonstration purposes but it can be included (see Example 2).

````julia
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
````

This next section of the code equilibrates/thermalizes the system before making measurements. Within the equilibration loop, the structure is fairly simple: the `local_fermion_update!` function is called to sweep over the particle configuration, attempting to update particle positions. Here, the number of updates performed before measurement is `N_equil` and `opt_bin_size` refers to the length of each bin from `N_opt_bins`.

The quantities `n_stab_W` and `δW` are passed to the `local_fermion_update!` function and controls the stability of the equal-time Green's function matrix `W`. After `n_stab_W` fermionic updates, there is a deviation check performed and if it exceeds the threshold `δW`, then the Green's function is recomputed from scratch.

Finally, the number of measurements that are averaged over per bin is given by `opt_bin_size = N_opt ÷ N_opt_bins`.
The bin-averaged measurements are written to file once `bin_size` measurements are accumulated using the `write_measurements!` function.


````julia
        # Iterate over optimization bin length
        for n in 1:opt_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf) = local_fermion_update!(
                    detwf, 
                    Np, 
                    model_geometry, 
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
````

After the last measurement is made, the last particle configuration is recorded globally for use in the initialization of the determinantal wavefunction for the next bin. 


````julia
        # Record the last particle configuration used for the start of the next bin.
        pconfig = detwf.pconfig
````

The primary part of this step of the simulation is the optimization of the variational parameters. Since all measurements for this bin have already been accumulated in the `measurement_container`, the `optimize_parameters!` function need only apply the Stochastic Reconfiguration (SR) procedure. Once an SR iteration is complete, all measurements are written to file. The `write_measurements!` method has an optional argument to write the parameter values to file, which we want for the optimization step.


````julia
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
````

## Simulate the system with optimized parameters
In this next section, we continue to sample particle configurations using `local_fermion_update!` function but without the SR optimization. This is mainly done to ensure proper statistics in calculating observables like the local energy. 


````julia
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
                    Np, 
                    model_geometry, 
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
````

## Record simulation metadata
Now that the optimization and simulation of the system are complete, we calculate the total time of the VMC simulation and the average final acceptance rate. Such information is saved to file using the `save_simulation_info` function. 

````julia
    # Record the total VMC time.
    metadata["vmc_time"] += metadata["opt_time"] + metadata["sim_time"]

    # Normalize acceptance rate.
    metadata["acceptance_rate"] /=  (N_opt + N_sim)

    # Write simulation summary TOML file.
    save_simulation_info(simulation_info, metadata)
````

## Post-processing
During the simulation, all measurments are written to file in HDF5 format for speed and portability; however, for analyzing data, having CSV files are mush more convenient. Calling the `process_measurements` function accomoplishes this. It will then up to the user to determine final processing and statistics. It should be noted that the next version of the code will have a convergence detection module and plotting function. 

````julia
    # Process all optimization and simulation measurements.
    process_measurements(
        measurement_container, 
        simulation_info, 
        determinantal_parameters,
        model_geometry
    )
````


## Execute script
VMC simulations are typically run from the command line as jobs on a computing cluster.
With this in mind, the following block of code only executes if the Julia script is run from the command line,
also reading in additional command line arguments.


````julia
# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Run the simulation.
    run_hubbard_chain_simulation(;
        sID         = parse(Int,     ARGS[1]), 
        L           = parse(Int,     ARGS[2]), 
        U           = parse(Float64, ARGS[3]), 
        nup         = parse(Float64, ARGS[4]), 
        ndn         = parse(Float64, ARGS[5]),
        pht         = parse(Bool,    ARGS[6]),
        N_equil     = parse(Int,     ARGS[7]), 
        N_opt       = parse(Int,     ARGS[8]), 
        N_opt_bins  = parse(Int,     ARGS[9]), 
        N_sim       = parse(Int,     ARGS[10]), 
        N_sim_bins  = parse(Int,     ARGS[11])
    )
end
````

For instance, the command
```
> julia hubbard_chain.jl 1 4 2.0 2 2 false 200 3000 100 6000 100
```
runs a VMC simulation of a ``N = 4`` half-filled 1D Hubbard model with interaction strength ``U = 2.0``.
In the VMC simulation, ``200`` sweeps through the lattice are be performed to thermalize the system, with ``3000`` optimization steps and ``6000`` simulation steps. During the simulation, bin-averaged measurements are written to file ``100`` times, with each bin of data containing the average of ``3,000/100 = 30`` sequential optimization measurements and ``6,000/100 = 60`` simulation measurements.


# 2) Hubbard model on a square lattice (with Jastrow factor)

Download this example as a [Julia script](https://github.com/atanjaro/VariationalMC.jl/blob/14b1c221029e963f09efcc8fe7b6be34dce7e126/examples/hubbard_square.jl).

In this example, we will work through simulating a repulsive Hubbard model on a two dimensional (2D) square lattice.
The Hubbard hamiltonian in 2D is given by
```math
\hat{H} = & -t \sum_{\sigma,\langle i, j \rangle} (\hat{c}^{\dagger}_{i, \sigma}, \hat{c}^{\phantom \dagger}_{j, \sigma} + {\rm H.c.})
+ \sum_i \hat{n}_{i, \uparrow}\hat{n}_{i, \downarrow}
- \mu \sum_{i, \sigma}\hat{n}_{i,\sigma}
```
where ``\hat{c}^\dagger_{i, \sigma} \ (\hat{c}^{\phantom \dagger}_{i, \sigma})`` creates (annihilates) a spin ``\sigma``
electron on site ``i`` in the lattice, and ``\hat{n}_{i, \sigma} = \hat{c}^\dagger_{i, \sigma} \hat{c}^{\phantom \dagger}_{i, \sigma}``
is the spin-``\sigma`` electron number operator for site ``i``. The nearest-neighbor hopping amplitude is ``t`` and ``\mu`` is the
chemical potential. The strength of the repulsive Hubbard interaction is controlled by ``U>0``. 

## Import packages
We begin by importing the necessary packages, including the Standard Library packages [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/), [Random](https://docs.julialang.org/en/v1/stdlib/Random/), and [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/) for perfroming linear algebra operations, random number generation, and C-style string formatting, respectively.

````julia
using Linear Algebra
using Random
using Printf
````

Next, we use [LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl) which exports a suite of types and methods useful for defining arbitrary lattice geometries, and the construction of neighbor tables. Finally, we import the [VariationalMC](https://github.com/atanjaro/VariationalMC.jl) package.

````julia
using LatticeUtilities
using VariationalMC
````

## Specify simulation parameters

The main body of the simulation is wrapped in a top-level function named `run_hubbard_square_simulation`that will take as keyword arguments various model and simulation parameters that we may want to change. The specific meaning of each argument will be discussed in more detail in later sections of the tutorial.

````julia
# We define a top-level function for running the VMC simulation.
function run_hubbard_square_simulation(;
    # KEYWORD ARGUMENTS
    sID,                    # Simulation ID.
    L,                      # System size.
    U,                      # Hubbard interaction.
    density,                # Electron density.
    pht,                    # Whether model is particle-hole transformed. 
    N_equil,                # Number of equilibration/thermalization updates.
    N_opt,                  # Number of optimization steps.
    N_opt_bins,             # Number of times bin-averaged measurements are written to file during optimization step.
    N_sim,                  # Number of simulation steps.
    N_sim_bins,             # Number of times bin-averaged measurements are written to file during simulation step.
    dt = 0.03,              # Optimization rate.
    dt_J = 1.0,             # Optional boost in the Jastrow optimization rate.
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
In this initial part of the script, we will specify what parameters in the trial wavefunction will be optimized and name our simulation. The `NamedTuple` called `optimized` will contain all possible parameters that are automatically initialized in the wavefunction. Setting a flag to `true` will cause that parameter value (or values) to change during the optimization. We also specify whether the model is particle-hole transformed and is required to be `true` if pairing symmetery is being added to the wavefunction. 

````julia
    # Select which parameters in the variational wavefunction will be optimized.
    optimize = (
        # local s-wave pairing
        Δ_0 = false,
        # site-dependent s-wave pairing  
        Δ_spd = false,
        # local d-wave pairing
        Δ_d = false,
        # site-dependent d-wave pairing 
        Δ_dpd = false,          
        # pairing momentum
        q_p = false,
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
    )
````

The datafolder is created by initializing an instances of the `SimulationInfo` type, and then calling the `initialize_datafolder` function. We also initialize the random number generator that will be used throughout the simulation. 


````julia
    # Construct the foldername the data will be written.
    df_prefix = @sprintf("hubbard_square_U%.2f_density%.2f_Lx%d_Ly%d_opt", U, density, L, L)

    # Append optimized parameter names to the foldername.
    datafolder_prefix = create_datafolder_prefix(optimize, df_prefix)

    # Initialize an instance of the SimulationInfo type.
    simulation_info = SimulationInfo(
        filepath = filepath, 
        datafolder_prefix = datafolder_prefix,
        sID = sID
    )

    # Initialize the directory the data will be written.
    initialize_datafolder(simulation_info)

    # Initialize a random number generator that will be used throughout the simulation.
    rng = Xoshiro(seed)
````

In addition, we can calculate the length of each bin by dividing the number of iterations/steps by the number of bins. The bin size is the number of measurements that are averaged over each time data is written during either the optimization or simulation steps.

````julia
    # Calculate optimization bins size.
    opt_bin_size = div(N_opt, N_opt_bins) 

    # Calculate simulation bins size.
    sim_bin_size = div(N_sim, N_sim_bins)
````

## Initialize simulation metadata
In this section of the code we record important metadata about the simulation, including initializing the random number
generator that will be used throughout the simulation.
The important metadata within the simulation will be recorded in the `metadata` dictionary.

````julia
    # Record simulation metadata.
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
````

In the above, `sID` stands for simulation ID, which is used to distinguish simulations that would otherwise be identical i.e. to
distinguish simulations that use the same parameters and are only different in the random seed used to initialize the simulation.
A valid `sID` is any positive integer greater than zero, and is used when naming the data folder the simulation results will be written to.
Specifically, the actual data folder created above will be `"$(filepath)/$(datafolder_prefix)-$(sID)"`.
Note that if you set `sID = 0`, then it will instead be assigned smallest previously unused integer value. For instance, suppose the directory
`"$(filepath)/$(datafolder_prefix)-1"` already exits. Then if you pass `sID = 0` to `SimulationInfo`, then the simulation ID
`sID = 2` will be used instead, and a directory `"$(filepath)/$(datafolder_prefix)-2"` will be created.

## Initialize model

The next part of the script defines the model that we will simulate.
First we define the lattice geometry for our model, relying on the
[LatticeUtilities](https://github.com/SmoQySuite/LatticeUtilities.jl.git) package to do so.
We define a the unit cell and size of our finite lattice using the [`UnitCell`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.UnitCell)
and [`Lattice`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.Lattice) types, respectively.
Lastly, we define various instances of the [`Bond`](https://smoqysuite.github.io/LatticeUtilities.jl/stable/api/#LatticeUtilities.Bond) type to represent the
the nearest-neighbor and next-nearest-neighbor bonds.
All of this information regarding the lattice geometry is then stored in an instance of the `ModelGeometry` type.

````julia
    # Initialize an instance of the UnitCell type.
    unit_cell = UnitCell(
        lattice_vecs = [[1.0, 0.0], [0.0, 1.0]],
        basis_vecs   = [[0.0, 0.0]]
    )

    # Initialize an instance of the Lattice type.
    lattice = Lattice(
        [L, L], 
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
    model_geometry = ModelGeometry(
        unit_cell, 
        lattice, 
        bonds
    )
````

## Initialize model parameters

The next step is to initialize our model parameters, which includes calculating the particle density in the canonical ensemble. The `get_particle_density` method includes mutliple ways of specifying the particle information. Here, we have chosen to pass the total particle density. The function then computes the corresponding total density, particle number, and electron number.

````julia
    # Determine the total particle density in the canonical ensemble. 
    (density, Np, Ne, nup, ndn) = get_particle_density(density, model_geometry, pht) 
````

We then specify parameters for our tight binding model and initializes all of the parameters that are in the determinantal part of the trial wavefunction. Parameters for the density-desnity Jastrow factor are initialized seperately. The `model_summary` function is used to write a `model_summary.toml` file, completely specifying the Hamiltonian that will be simulated. Lastly, we can also set an initial particle configuration, if we have one. If an empty array is provided instead, a random configuration will be given at the start of the simulation. 

````julia
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

    # Write model summary TOML file specifying the Hamiltonian that will be simulated.
    model_summary(
        simulation_info, 
        determinantal_parameters, 
        density_J_parameters, 
        pht, 
        model_geometry, 
        tight_binding_model, 
        U
    )

    # Initialize the (fermionic) particle configuration.
    pconfig = Int[]
````

## Initialize measurements
Finally, we initialize the mesaurement container, which accumulates the sums of measurements made in a bin. The standard measurements are of the local energy `local_energy`, double occupancy `double_occ`, particle confiugurations `pconfig`, and variational parameters `parameters`. There is also the option to add additional observables to measure through the `initialize_simulation_measurement!`and `initialize_correlation_measurement!` methods. Here, we have added density-density correlation measurements.

````julia
    # Initialize the container that measurements will be accumulated into.
    measurement_container = initialize_measurement_container(
        N_opt, 
        opt_bin_size, 
        N_sim, 
        sim_bin_size,
        determinantal_parameters,
        model_geometry
    )

    # Add density-density correlation measurements.
    initialize_correlation_measurement!(
        "density", 
        measurement_container, 
        model_geometry
    )
````

The `initialize_measurement_directories` can now be used used to initialize the various subdirectories in the data folder that the measurements will be written to.
For more information, please refer to the Simulation Output Overview page.

````julia
    # Initialize the sub-directories to which the various measurements will be written.
    initialize_measurement_directories(
        simulation_info, 
        measurement_container
    )
````

## Optimize the variational parameters
Now that we have set-up the VMC simulation, we can begin optimizing the variational parameters. At the start of each bin, we initialize the variational wavefunction which includes the determinantal part and the density-density Jastrow factor. On the first iteration, if no initial configuration is specified, a random one will be generated.

````julia
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
````

This next section of the code equilibrates/thermalizes the system before making measurements. Within the equilibration loop, the structure is fairly simple: the `local_fermion_update!` function is called to sweep over the particle configuration, attempting to update particle positions. Here, the number of updates performed before measurement is `N_equil` and `opt_bin_size` refers to the length of each bin from `N_opt_bins`.

The quantities `n_stab_W` and `δW` are passed to the `local_fermion_update!` function and controls the stability of the equal-time Green's function matrix `W`. After `n_stab_W` fermionic updates, there is a deviation check performed and if it exceeds the threshold `δW`, then the Green's function is recomputed from scratch.

Finally, the number of measurements that are averaged over per bin is given by `opt_bin_size = N_opt ÷ N_opt_bins`.
The bin-averaged measurements are written to file once `bin_size` measurements are accumulated using the `write_measurements!` function.


````julia
        # Iterate over optimization bin length
        for n in 1:opt_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf, density_J_factor) = local_fermion_update!(
                    detwf, 
                    density_J_factor,
                    density_J_parameters,
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
                density_J_factor,
                optimize,
                model_geometry, 
                U,
                Np, 
                pht
            )
        end
````

After the last measurement is made, the last particle configuration is recorded globally for use in the initialization of the determinantal wavefunction for the next bin. 


````julia
        # Record the last particle configuration used for the start of the next bin.
        pconfig = detwf.pconfig
````

The primary part of this step of the simulation is the optimization of the variational parameters. Since all measurements for this bin have already been accumulated in the `measurement_container`, the `optimize_parameters!` function need only apply the Stochastic Reconfiguration (SR) procedure. Once an SR iteration is complete, all measurements are written to file. The `write_measurements!` method has an optional argument to write the parameter values to file, which we want for the optimization step.


````julia
        # Attempt to update the variational parameters using the Stochastic Reconfiguration procedure. 
        optimize_parameters!( 
            measurement_container,  
            determinantal_parameters, 
            density_J_parameters,
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
````

## Simulate the system with optimized parameters
In this next section, we continue to sample particle configurations using `local_fermion_update!` function but without the SR optimization. This is mainly done to ensure proper statistics in calculating observables like the local energy. 


````julia
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

        # Iterate over optimization bin length
        for n in 1:sim_bin_size

            # Iterate over equilibration/thermalization updates
            for equil in 1:N_equil
                (acceptance_rate, detwf, density_J_factor) = local_fermion_update!(
                    detwf, 
                    density_J_factor,
                    density_J_parameters,
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
                density_J_factor,
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
````

## Record simulation metadata
Now that the optimization and simulation of the system are complete, we calculate the total time of the VMC simulation and the average final acceptance rate. Such information is saved to file using the `msave_simulation_info` function. 

````julia
    # Record the total VMC time.
    metadata["vmc_time"] += metadata["opt_time"] + metadata["sim_time"]

    # Normalize acceptance rate.
    metadata["acceptance_rate"] /=  (N_opt + N_sim)

    # Write simulation summary TOML file.
    save_simulation_info(simulation_info, metadata)
````

## Post-processing
During the simulation, all measurments are written to file in HDF5 format for speed and portability; however, for analyzing data, having CSV files are mush more convenient. Calling the `process_measurements` function accomoplishes this. It will then up to the user to determine final processing and statistics. It should be noted that the next version of the code will have a convergence detection module and plotting function. 

````julia
    # Process all optimization and simulation measurements.
    process_measurements(
        measurement_container, 
        simulation_info, 
        determinantal_parameters, 
        density_J_parameters,
        model_geometry
    )
````


## Execute script
VMC simulations are typically run from the command line as jobs on a computing cluster.
With this in mind, the following block of code only executes if the Julia script is run from the command line,
also reading in additional command line arguments.


````julia
# Only execute if the script is run directly from the command line.
if abspath(PROGRAM_FILE) == @__FILE__

    # Run the simulation.
    run_hubbard_square_simulation(;
        sID         = parse(Int,     ARGS[1]), 
        L           = parse(Int,     ARGS[2]), 
        U           = parse(Float64, ARGS[3]), 
        density     = parse(Float64, ARGS[4]), 
        pht         = parse(Bool,    ARGS[5]),
        N_equil     = parse(Int,     ARGS[6]), 
        N_opt       = parse(Int,     ARGS[7]), 
        N_opt_bins  = parse(Int,     ARGS[8]), 
        N_sim       = parse(Int,     ARGS[9]), 
        N_sim_bins  = parse(Int,     ARGS[10])
    )
end

````

For instance, the command
```
> julia hubbard_square.jl 1 4 2.0 0.5 1000 3000 100 6000 100
```
runs a VMC simulation of a ``N = 4 \times 4`` quarter-filled 2D Hubbard model with interaction strength ``U = 2.0``.
In the VMC simulation, ``1000`` sweeps through the lattice are be performed to thermalize the system, with ``3000`` optimization steps and ``6000`` simulation steps. During the simulation, bin-averaged measurements are written to file ``100`` times, with each bin of data containing the average of ``3,000/100 = 30`` sequential optimization measurements and ``6,000/100 = 60`` simulation measurements.


# 3) Hubbard model on a square lattice with MPI Parallelization

Download this example as a [Julia script](https://github.com/atanjaro/VariationalMC.jl/blob/cc49763acfbd71c9dae5340603e4dc5f96d76e5e/examples/hubbard_square_mpi.jl).

This example will build on the previous 2) Hubbard model on a square lattice (withh Jastrow factor) example, demonstrating how to add parallelization with MPI using the [MPI.jl](https://github.com/JuliaParallel/MPI.jl) package. By this we mean that each MPI process will act as independent walker, running it's own independent VMC simulation.

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
    comm::MPI.Comm;         # MPI communicator.
    # KEYWORD ARGUMENTS
    sID,                    # Simulation ID.
    L,                      # System size.
    U,                      # Hubbard interaction.
    density,                # Electron density.
    pht,                    # Whether model is particle-hole transformed. 
    N_equil,                # Number of equilibration/thermalization updates.
    N_opt,                  # Number of optimization steps.
    N_opt_bins,             # Number of times bin-averaged measurements are written to file during optimization step.
    N_sim,                  # Number of simulation steps.
    N_sim_bins,             # Number of times bin-averaged measurements are written to file during simulation step.
    dt = 0.03,              # Optimization rate.
    dt_J = 1.0,             # Optional boost in the Jastrow optimization rate.
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