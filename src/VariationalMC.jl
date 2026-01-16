module VariationalMC

using LatticeUtilities
using Random
using LinearAlgebra
using SparseArrays
using OrderedCollections
using Printf
using CSV
using HDF5
using TOML
using JLD2
using Profile
using Revise
using DelimitedFiles
using DataFrames
using DataStructures
using MPI

# # get and set package version number as global constant
# const VARIATIONALMC_VERSION = PkgVersion.@Version 0

include("ModelGeometry.jl")
export ModelGeometry

include("Parameters.jl")
export TightBindingModel, SpinModel, DeterminantalParameters, JastrowParameters

include("Hamiltonian.jl")

include("ParticleConfiguration.jl")
export get_particle_numbers, get_particle_density

include("DeterminantalWavefunction.jl")
export DeterminantalWavefunction, DeterminantalWavefunctionTABC, get_determinantal_wavefunction

include("Greens.jl")

include("Jastrow.jl")
export JastrowFactor, get_jastrow_factor

include("Markov.jl")
export local_fermion_update!

include("Optimizer.jl")
export optimize_parameters!

include("SimulationInfo.jl")
export SimulationInfo, save_simulation_info, initialize_datafolder, create_datafolder_prefix 

include("model_summary.jl")
export model_summary

include("Measurements/initialize_measurements.jl")
export initialize_measurement_container, initialize_measurement_directories, initialize_simulation_measurement!, initialize_correlation_measurement!

include("Measurements/make_measurements.jl")
export make_measurements!

include("Measurements/optimization_measurements.jl")

include("Measurements/scalar_measurements.jl")

include("Measurements/simulation_measurements.jl")

include("Measurements/correlation_measurements.jl")

include("Measurements/write_measurements.jl")
export write_measurements!

include("Measurements/process_measurements.jl")
export process_measurements

############################
## PACKAGE INITIALIZATION ##
############################

# set number of threads for BLAS to 1.
# we assume for now the default OpenBLAS that ships with Julia is used.
# this behavior will in general need to be changed if we want to use other BLAS/LAPACK libraries.
# function __init__()

#     BLAS.set_num_threads(1)
#     return nothing
# end

end
