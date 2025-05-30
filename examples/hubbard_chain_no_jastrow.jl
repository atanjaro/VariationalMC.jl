using LatticeUtilities
using Random
using LinearAlgebra
using DelimitedFiles
using Profile
using OrderedCollections
using CSV
using DataFrames
using DataStructures
using Printf
using JLD2
using Revise

# files to include
include("ModelGeometry.jl")
include("Parameters.jl")
include("Hamiltonian.jl");
include("ParticleConfiguration.jl");
include("Wavefunction.jl");
include("Greens.jl")
include("Jastrow.jl");
include("Markov.jl");
include("Optimizer.jl");
include("SimulationInfo.jl");
include("Measurements/initialize_measurements.jl")
include("Measurements/make_measurements.jl")
include("Measurements/optimization_measurements.jl")
include("Measurements/scalar_measurements.jl")
include("Measurements/simulation_measurements.jl")
include("Measurements/write_measurements.jl")


###########################################
##          LATTICE PARAMETERS           ##
###########################################
# define the size of the lattice
Lx = 4;
Ly = 1;

# chain unit cell
unit_cell = UnitCell(
    lattice_vecs = [[1.0]], # lattice vectors
    basis_vecs   = [[0.0]]  # basis vectors
);

# build the lattice
lattice = Lattice(
    [Lx], 
    [true]
);

# define nearest neighbor bonds
bond_x = Bond(
    orbitals = (1,1), 
    displacement = [1]
);

# define next nearest neighbor bonds
bond_xp = Bond(
    orbitals = (1,1), 
    displacement = [2]
);

# collect all bond definitions: [[nearest],[next-nearest]]
bonds = [[bond_x], [bond_xp]];      

# define model geometry
model_geometry = ModelGeometry(
    unit_cell, 
    lattice, 
    bonds
);


#########################################
##          MODEL PARAMETERS           ##
#########################################
# nearest neighbor hopping amplitude
t = 1.0;

# next nearest neighbor hopping amplitude
tp = 0.0;

# onsite Hubbard repulsion
U = 1.0;

# # OPTIONAL: use to read-in initial parameters from file
# path_to_parameter_file = "/Users/xzt/Documents/VariationalMC/src/vpar_out.dat"

# minimum value of each variational parameter (this is to avoid open shell issues)
minabs_vpar = 1e-4;

# select which parameters will be optimized
optimize = (
    # (BCS) chemical potential
    μ = false,
    # onsite (s-wave)
    Δ_0 = false,
    # d-wave
    Δ_d = false,
    # spin-z (AFM)
    Δ_afm = true,
    # charge density 
    Δ_cdw = false,
    # site-dependent charge (stripe)
    Δsdc = false,
    # site-dependent spin (stripe)
    Δ_sds = false,
    # density Jastrow 
    djastrow = false,
    # spin Jastrow
    sjastrow = false
)

# whether model is particle-hole transformed
pht = false;

# define electron density
n̄ = 1.0;

# # define electron numbers     
# nup = 2
# ndn = 2

# Get particle numbers 
(Np, Ne, nup, ndn) = get_particle_numbers(n̄);       # Use this if an initial density is specified

# # Get particle density 
# (density, Np, Ne) = get_particle_density(nup, ndn);   # Use this if initial particle numbers are specified


######################################
##          FILE HANDLING           ##
######################################
# specify filepath
filepath = ".";

# simulation ID
sID = 1;

# construct the foldername the data will be written 
df_prefix = @sprintf "hubbard_chain_U%.2f_n%.2f_Lx%d_Ly%d_opt" U n̄ Lx Ly;

# append parameters to the foldername
datafolder_prefix = create_datafolder_prefix(optimize, df_prefix)

# initialize an instance of the SimulationInfo type
simulation_info = SimulationInfo(
    filepath = filepath, 
    datafolder_prefix = datafolder_prefix,
    sID = sID
);

# initialize the directory the data will be written to
initialize_datafolder(simulation_info);


##############################################
##          SIMULATION PARAMETERS           ##
##############################################
# random seed
seed = abs(rand(Int)) 
println("seed = ", seed)

# initialize random number generator
rng = Xoshiro(seed);

# number of equilibration/thermalization steps 
N_equil = 1000;

# number of minimization/optimization updates
N_opts = 1000;

# optimization bin size
opt_bin_size = 100;

# number of simulation updates 
N_updates = 1000;

# number of simulation bins
N_bins = 100;

# simulation bin size
bin_size = div(N_updates, N_bins);

# number of steps until numerical stabilization is performed 
n_stab_W = 50;
n_stab_T = 50;

# maximum allowed error in the equal-time Green's function
δW = 1e-3;

# maximum allowed error in the Jastrow T vector
δT = 1e-3;

# stabilization factor for Stochastic Reconfiguration
η = 1e-4;    

# optimization rate for Stochastic Reconfiguration
dt = 0.03;       

# whether debug statements are printed (for speed, this should be turned off)
debug = false;

# initialize metadata dictionary
metadata = Dict()

# record simulation parameters
metadata["N_equil"] = N_equil
metadata["N_opts"] = N_opts
metadata["N_updates"] = N_updates
metadata["N_bins"] = N_bins
metadata["seed"] = seed
metadata["pht"] = true
metadata["δW"] = δW
metadata["δT"] = δT
metadata["dt"] = dt 
metadata["opt_flags"] = optimize 
metadata["avg_acceptance_rate"] = 0.0


##############################################
##          SET-UP VMC SIMULATION           ##
##############################################

# define non-interacting tight binding model
tight_binding_model = TightBindingModel(t, tp);

# initialize determinantal parameters
determinantal_parameters = DeterminantalParameters(
    optimize, 
    tight_binding_model, 
    model_geometry, 
    minabs_vpar, 
    Ne, 
    pht
);

# # OPTIONAL: initialize determinantal parameters from file
# determinantal_parameters =  DeterminantalParameters(
#     optimize, 
#     model_geometry, pht, 
#     path_to_parameter_file
# );

# initialize particle configuration cache
pconfig_cache = nothing

# initialize measurement container for VMC measurements
measurement_container = initialize_measurement_container(
    N_opts, 
    opt_bin_size, 
    N_bins, 
    bin_size,
    determinantal_parameters,
    model_geometry
);

# initialize the sub-directories to which the various measurements will be written
initialize_measurement_directories(
    simulation_info, 
    measurement_container
);


#############################################
##          OPTIMIZATION UPDATES           ##
#############################################
opt_start_time = time();
for bin in 1:N_opts
    # initialize determinantal wavefunction
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
    );   

    # equilibrate/thermalize the system   
    for step in 1:N_equil 
        # perform local update to electronic degrees of freedom
        (acceptance_rate, detwf) = local_fermion_update!(
            detwf, 
            Ne, 
            model_geometry, 
            n_stab_W,
            δW, 
            rng
        );

        # record acceptance rate                                                        
        metadata["avg_acceptance_rate"] += acceptance_rate;
    end

    # perform measurements for optimization
    for n in 1:opt_bin_size
        println("Starting Monte Carlo step = ", n)

        # perform local update to electronic degrees of freedom
        (acceptance_rate, detwf) = local_fermion_update!(
            detwf, 
            Ne, 
            model_geometry, 
            n_stab_W,
            δW, 
            rng); 
                                                                                                    
        # record acceptance rate                                                        
        metadata["avg_acceptance_rate"] += acceptance_rate;
        
        # make basic measurements
        make_measurements!(
            measurement_container, 
            detwf, 
            tight_binding_model, 
            determinantal_parameters, 
            model_geometry, 
            Ne, 
            pht
        );

        # write measurements (to file)
        write_measurements!(
            bin, 
            n, 
            measurement_container, 
            simulation_info
        );
    end

    # save particle configuration from the current bin
    pconfig_cache = detwf.pconfig

    # process measurements

    # perform update to variational parameters
    stochastic_reconfiguration!( 
        measurement_container,  
        determinantal_parameters, 
        η, 
        dt, 
        bin, 
        opt_bin_size
    );  
end     
opt_end_time = time();

# time for optimization 
opt_time = opt_end_time - opt_start_time;


###########################################
##          SIMULATION UPDATES           ##
###########################################
sim_start_time = time();
for bin in 1:N_opts
    # initialize determinantal wavefunction
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
    );   

    # equilibrate/thermalize the system   
    for step in 1:N_equil 
        # perform local update to electronic degrees of freedom
        (acceptance_rate, detwf) = local_fermion_update!(
            detwf, 
            Ne, 
            model_geometry, 
            n_stab_W,
            δW, 
            rng
        );

        # record acceptance rate                                                        
        metadata["avg_acceptance_rate"] += acceptance_rate;
    end

    # perform measurements for optimization
    for n in 1:opt_bin_size
        println("Starting Monte Carlo step = ", n)

        # perform local update to electronic degrees of freedom
        (acceptance_rate, detwf) = local_fermion_update!(
            detwf, 
            Ne, 
            model_geometry, 
            n_stab_W,
            δW, 
            rng); 
                                                                                                    
        # record acceptance rate                                                        
        metadata["avg_acceptance_rate"] += acceptance_rate;
        
        # make basic measurements
        make_measurements!(
            measurement_container, 
            detwf, 
            tight_binding_model, 
            determinantal_parameters, 
            model_geometry, 
            Ne, 
            pht
        );

        # write measurements (to file)
        write_measurements!(
            bin, 
            n, 
            measurement_container, 
            simulation_info
        );
    end

    # save particle configuration from the current bin
    pconfig_cache = detwf.pconfig

     # process measurements
     
end     
sim_end_time = time();