# 4) Hubbard model on a rectangular lattice with initial parameters

The script that follows this example can be found in the `/example/` directory in the repo. 

This Example uses the same base script as the previous Hubbard model on a square lattice (Example 2) with the geometry altered to be rectuangular and our trial wavefunction will have a set of initial values instead of being randomly initialized. 

## Specify simulation parameters

The top-level function named `run_hubbard_rectangle_simulation` is similar to the previous `run_hubbard_square_simulation` that we defined in Example 2 but we have additional arguments corresponding to the system size in the ``x`` and ``y`` directions. 

````julia
# We define a top-level function for running the VMC simulation.
function run_hubbard_rectangle_simulation(;
    sID,                                # Simulation ID.
    Lx,                                 # System size x. 
    Ly,                                 # System size y.
    U,                                  # Hubbard interaction.
    density,                            # Electron density.
    pht,                                # Whether model is particle-hole transformed. 
    N_equil,                            # Number of equilibration/thermalization updates.
    N_opt,                              # Number of optimization steps.
    N_opt_bins,                         # Number of times bin-averaged measurements are written to file during optimization step.
    N_sim,                              # Number of simulation steps.
    N_sim_bins,                         # Number of times bin-averaged measurements are written to file during simulation step.
    dt              = 0.1,              # Optimization rate.
    dt_J            = 1.0,              # Optional boost in the Jastrow optimization rate.
    η               = 1e-4,             # Optimization stablity factor.
    n_stab_W        = 50,               # Green's function stabilization frequency.
    δW              = 1e-3,             # Maximum allowed error in the Green's function. 
    n_stab_T        = 50,               # Jastrow factor stabilization frequency.
    δT              = 1e-3,             # Maximum allowed error in the Jastrow factor.           
    seed            = abs(rand(Int)),   # Seed for random number generator.
    filepath        ="."                # Filepath to where data folder will be created.
)
````

## Initialize simulation
For a recatngle, we will define a charge and spin stripe wavefunction in which the site-dependent spin order parameter ``\Delta_{ssd}`` and the site-dependent charge order parameter ``\Delta_{csd}`` are optimized as well as the pseudopotentials in both the density-density and spin-spin Jastrow factors, ``\mathcal{J}_d`` and ``\mathcal{J}_s``, respectively.

````julia
    # Select which parameters in the variational wavefunction will be optimized.
    optimize = (
        # Uniform s-wave pairing
        Δ_0         = false,
        # Site-dependent s-wave pairing (Larkin-Ovchinnikov)
        Δ_slo       = false,
        # Site-dependent s-wave pairing (Fulde-Ferrell)
        Δ_sff       = false,
        # Uniform d-wave pairing
        Δ_d         = false,
        # Site-dependent d-wave pairing (Larkin-Ovchinnikov)
        Δ_dlo       = false,
        # Site-dependent d-wave pairing (Fulde-Ferrell)
        Δ_dff       = false,     
        # In-plane magnetization
        Δ_sx        = false,
        # Out-of-plane magnetization
        Δ_sz        = false,
        # Site-dependent magnetization
        Δ_ssd       = true,
        # Chemical potential
        μ           = false,
        # Charge density wave
        Δ_cdw       = false,
        # Site-dependent charge density
        Δ_csd       = false,
        # Density-density Jastrow pseudopotentials
        density_J   = true,
        # Spin-spin Jastrow pseudopotentials
        spin_J      = true
    )
````

To setup the VMC simulation, we define all of the parameters that will be used to build our variational wavefunction. In this case, we introduce a finite value of the next nearest neighbor hopping amplitide ``t^\prime``. We create an initialize an instance of determinantal parameters before overriding them with a set of initial values using the `manual_input_parameters!` method. Since we want to use both types of Jastrow factor in our wavefunction, we initialize the density-density and spin-spin parameters seperately. Note that we also add both types of Jastrow parameters to the `model_summary`. The rest of the simulation set-up remains unchanged from Example 2. 

```julia
    # Define the nearest neighbor hopping amplitude, setting the energy scale of the system. 
    t = 1.0

    # Define the next-nearest neighbor hopping amplitude.
    tp = -0.25

    # Define the third-nearest neighbor hopping amplitude.
    tpd = 0.0

    # Define the non-interacting tight binding model.
    tight_binding_model = TightBindingModel(t, tp, tpd)

    # Define a Hubbard model.
    hubbard_model = HubbardModel(U, 0.0)

    # Initialize determinantal variational parameters.
    determinantal_parameters = DeterminantalParameters(
        tight_binding_model,
        model_geometry, 
        optimize,
        Ne, 
        pht
    )

    # Override random determinantal parameters with some initial values.
    init_det_pars = [0.001,                                                                                     # In-plane magnetization
                     0.0001,                                                                                    # Out-of-plane magnetization
                     0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, 0.0, 0.0,      # Site-dependent magnetization
                     -0.15224093497742808,                                                                      # Exact chemical potential
                     0.001,                                                                                     # Charge density wave
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3]            # Site-dependent charge density
    manual_input_parameters!(init_det_pars, determinantal_parameters)

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
        tight_binding_model,
        hubbard_model,
        determinantal_parameters, 
        density_J_parameters, 
        spin_J_parameters, 
        model_geometry, 
        pht
    )

    # Initialize the (fermionic) particle configuration.
    # Will start with a random initial configuration unless provided a starting configuration.
    pconfig = Int[]
```

Both types of Jastrow parameters need to be added to the measurement container. The rest of the measurement set-up remains the same as Example 2. 

```julia
    # Initialize the container that measurements will be accumulated into.
    measurement_container = initialize_measurement_container(
        determinantal_parameters,
        density_J_parameters,
        spin_J_parameters,
        model_geometry,
        N_opt, 
        opt_bin_size, 
        N_sim, 
        sim_bin_size
    )
```

In addition to the standard measurements, we add site-dependent spin and density measurements.

```julia
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
```

The only modifications to the optimization steps are the seperate initialization of each Jastrow factor with the wavefunction at the start of every bin as well as the addition of both types of Jastrow factor in the arguments of the updating functions. 

```julia
    # Iterate over optimization bins.
    for bin in 1:N_opt_bins

        # Initialize the determinantal wavefunction.
        detwf = get_determinantal_wavefunction(
            tight_binding_model, 
            determinantal_parameters, 
            model_geometry,
            optimize, 
            pconfig,
            Np, 
            nup, 
            ndn, 
            rng,
            pht
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
                    model_geometry,
                    density_J_parameters,
                    spin_J_parameters,
                    Np, 
                    δW, 
                    δT,
                    n_stab_W,
                    n_stab_T,
                    rng,
                    pht
                )

                # Record acceptance rate.
                metadata["acceptance_rate"] += acceptance_rate
            end

            # Make measurements, with results being recorded in the measurement container.
            make_measurements!(
                measurement_container, 
                detwf, 
                density_J_factor,
                spin_J_factor,
                tight_binding_model, 
                hubbard_model,
                determinantal_parameters, 
                density_J_parameters,
                spin_J_parameters,
                model_geometry,
                optimize,
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
```

This is done similarly for the simulation steps.

```julia
    # Iterate over simulation bins.
    for bin in 1:N_sim_bins

        # Initialize the determinantal wavefunction.
        detwf = get_determinantal_wavefunction(
            tight_binding_model, 
            determinantal_parameters, 
            model_geometry,
            optimize, 
            pconfig,
            Np, 
            nup, 
            ndn, 
            rng,
            pht
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
                    model_geometry,
                    density_J_parameters,
                    spin_J_parameters,
                    Np, 
                    δW, 
                    δT,
                    n_stab_W,
                    n_stab_T,
                    rng,
                    pht
                )

                # Record acceptance rate.
                metadata["acceptance_rate"] += acceptance_rate
            end

            # Make measurements, with results being recorded in the measurement container.
            make_measurements!(
                measurement_container, 
                detwf, 
                density_J_factor,
                spin_J_factor,
                tight_binding_model, 
                hubbard_model,
                density_J_parameters,
                spin_J_parameters,
                model_geometry,
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
```

Finally, ensure that the measurement processing also accounts for both type of Jastrow factors.

```julia
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
```

All other aspects of the main simulation function are carried over from Example 2. 


## Execute script
This part of the script is the same as Example 2, save for the arguments `Lx` and `Ly` for a rectangular lattice geometry.

```julia
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
```