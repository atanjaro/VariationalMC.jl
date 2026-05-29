module VariationalMC

# import external dependencies
using LinearAlgebra
using FFTW
using Random
using Statistics
using Printf
using JLD2   
using Reexport
using PkgVersion
using TOML
using MPI
using Format
using HDF5
using OffsetArrays

# import `SmoQy` packages
using LatticeUtilities

# import methods from to overload
import LinearAlgebra: mul!, lmul!, rmul!

# re-export external package/modules as sub-modules
@reexport import LatticeUtilities

# get and set package version number as global constant
const VARIATIONALMC_VERSION = PkgVersion.@Version 0

##################################
## GENERAL SIMULATION FRAMEWORK ##
##################################

# Contains functions that are of general use.
include("utilities.jl")

# Define SimulationInfo struct for tracking things like simulation ID, process ID,
# and data folder location.
include("SimulationInfo.jl")
export SimulationInfo, save_simulation_info, initialize_datafolder

# Defines all aspects of the geometry of the model, including the `UnitCell`,
# `Lattice`, and a list of `Bond` definitions as defined in the `LatticeUtilities` package.
include("ModelGeometry.jl")
export ModelGeometry, add_bond!, get_bond_id

# Defines all aspects of initializing and tracking of particles in the Canonical Ensemble.
include("ParticleConfiguration.jl")
export ParticleConfiguration

# Define TightBindingModel.
include("Models/TightBindingModel.jl")
export TightBindingModel

# Define TightBindingParameters
include("../src/Parameters/TightBindingParameters.jl")
export TightBindingParameters

###################
## HUBBARD MODEL ##
###################

# Define HubbardModel and ExtendedHubbardModel.
include("Models/HubbardModel.jl")
export HubbardModel, ExtendedHubbardModel

# Define HubbardParameters and ExtendedHubbardParameters.
include("../src/Parameters/HubbardParameters.jl")
export HubbardParameters, ExtendedHubbardParameters

##############################
## VARIATIONAL WAVEFUNCTION ##
##############################

# Define DeterminantalParameters.
include("Parameters/DeterminantalParameters.jl")
export DeterminantalParameters, add_parameter!

# Define JastrowParameters.
include("Parameters/JastrowParameters.jl")
export JastrowParameters

# Defines dictionaries as global variables that contain the names of all
# allowed parameters in the variational wavefunction.
include("Parameters/parameter_names.jl")
export DETERMINANTAL_PARAMETERS
export JASTROW_PARAMETERS

# Defines all aspects of constructing an auxiliary Hamiltonian.
include("Hamiltonian.jl")
export Hamiltonian

# Define DeterminantalWavefunction.
include("DeterminantalWavefunction.jl")
export DeterminantalWavefunction

# Contains all methods for initializing and updating the equal-time Green's function.
include("greens_function.jl")

# Define JastrowFactor.
include("JastrowFactor.jl")
export AbstractJastrowFactor, FermionJastrowFactor, JastrowFactor

##############################
## MARKOV CHAIN MONTE CARLO ##
##############################

# Defines all methods used to perform Monte Carlo updates.
include("../src/MarkovMove.jl")
export local_fermion_update!

# Define Optimizer.
include("../src/Optimizer.jl")
export Optimizer, update_optimizer!

##########################################
## MEASUREMENTS, DATA ANALYSIS & OUTPUT ##
##########################################

# Method to write model summary/definition to TOML file.
include("../src/model_summary.jl")
export model_summary

# Tight-binding specific measurements.
include("Measurements/tight_binding_measurements.jl")
export measure_hopping_energy

# Hubbard specific measurements.
include("Measurements/hubbard_model_measurements.jl")
export measure_hubbard_energy, measure_ext_hubbard_energy

# Various local and global observable measurements
include("Measurements/scalar_measurements.jl")
export measure_double_occ
export measure_n
export measure_spinz
export measure_variational_energy
export measure_coupling_energy

# Opimization specific measurements.
include("Measurements/optimization_measurements.jl")
export measure_Dk

# Defines dictionaries as global variables that contain the names of all
# local measurements and correlation measurements that can be made, and the
# type ID type they are reported in terms of.
include("Measurements/measurement_names.jl")
export GLOBAL_MEASUREMENTS
export LOCAL_MEASUREMENTS
export CORRELATION_FUNCTIONS

# Define CorrelationContainer struct to store correlation measurements in.
include("Measurements/CorrelationContainer.jl")

# Define methods for performing fast Fourier transforms.
include("../src/Measurements/fourier_transform.jl")

# Initialize measurement container
include("Measurements/initialize_measurements.jl")
export initialize_measurement_container
export initialize_measurements!
export initialize_site_dependent_measurements!
export initialize_correlation_measurements!
export initialize_correlation_measurement!

# Make measurements
include("Measurements/make_measurements.jl")
export make_measurements!

# Write measurements to file.
include("Measurements/write_measurements.jl")
export write_measurements!

# Implements function to merge HDF5 bin files for a given pID (process ID)
# into a single HDF5 file. Also includes functions to delete all binned data.
include("Measurements/merge_bins.jl")
export merge_bins, rm_bins

# Implements utility function for converting numbers to string.
include("Measurements/num_to_string_formatter.jl")

# Function for exporting global measurement stats to csv file.
include("Measurements/export_global_stats_to_csv.jl")
export export_global_stats_to_csv

# Function for exporting global measurement stats to csv file.
include("Measurements/export_local_stats_to_csv.jl")
export export_local_stats_to_csv

# Function for exporting correlation measurement stats to csv file.
include("Measurements/export_correlation_stats_to_csv.jl")
export export_correlation_stats_to_csv

# Function for performing Jackknife resampling of measurement data.
include("Measurements/jackknife.jl")

# Internal functions for processing the binned data to calculate final
# statistics using a single process.
include("Measurements/process_measurements_internals.jl")

# Internal functions for processing the binned data to calculate final
# statistics using MPI parallelization to accelerate the computation.
include("Measurements/process_measurements_internals_mpi.jl")

# public api functions for processing measurements
include("Measurements/process_measurements.jl")
export process_measurements


############################
## PACKAGE INITIALIZATION ##
############################

# set number of threads for BLAS to 1.
# we assume for now the default OpenBLAS that ships with Julia is used.
# this behavior will in general need to be changed if we want to use other BLAS/LAPACK libraries.
function __init__()

    BLAS.set_num_threads(1)
    FFTW.set_num_threads(1)
    return nothing
end

end # VARIATIONALMC.JL
