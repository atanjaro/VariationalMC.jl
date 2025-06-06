# API

<!-- ## SimulationInfo Type and Methods

- [`SimulationInfo`](@ref)
- [`SimulationInfo(;)`](@ref)
- [`save_simulation_info`](@ref)
- [`initialize_datafolder`](@ref)
- [`create_datafolder_prefix`](@ref)
- [`model_summary`](@ref)

```@docs
SimulationInfo
SimulationInfo(;)
save_simulation_info
initialize_datafolder
create_datafolder_prefix
model_summary
``` -->

## ModelGeometry Type and Methods

- [`ModelGeometry`](@ref)
- [`reduce_index_2d`](@ref)
- [`reduce_index_1d`](@ref)
- [`max_dist`](@ref)
- [`x`](@ref)
- [`y`](@ref)

```@docs
ModelGeometry
reduce_index_2d
reduce_index_1d
max_dist
x
y
d
```

## Parameters Types and Methods

- [`TightBindingModel`](@ref)
- [`SpinModel`](@ref)
- [`DeterminantalParameters`](@ref)
- [`DeterminantalParameters(;)`](@ref)
- [`JastrowParameters`](@ref)
- [`JastrowParameters(;)`](@ref)
- [`collect_parameters`](@ref)
- [`update_parameters!`](@ref)
- [`readin_parameters`](@ref)

```@docs
TightBindingModel
SpinModel
DeterminantalParameters
DeterminantalParameters(;)
JastrowParameters
JastrowParameters(;)
collect_parameters
update_parameters!
readin_parameters
```

## Hamiltonian Methods

- [`build_auxiliary_hamiltonian`](@ref)
- [`build_tight_binding_hamiltonian`](@ref)
- [`build_variational_hamiltonian`](@ref)
- [`diagonalize`](@ref)
- [`is_openshell`](@ref)
- [`get_variational_matrices`](@ref)
- [`get_tb_chem_pot`](@ref)

```@docs
build_auxiliary_hamiltonian
build_tight_binding_hamiltonian
build_variational_hamiltonian
diagonalize
is_openshell
get_variational_matrices
get_tb_chem_pot
```

### Adding different symmetries to the Hamiltonian

- [`add_pairing_symmetry!`](@ref)
- [`add_spin_order!`](@ref)
- [`add_charge_order!`](@ref)
- [`add_chemical_potential!`](@ref)

```@docs
add_pairing_symmetry!
add_spin_order!
add_charge_order!
add_chemical_potential!
```

## DeterminantalWavefunction Type and Methods

- [`DeterminantalWavefunction`](@ref)
- [`get_determinantal_wavefunction`](@ref)

```@docs
DeterminantalWavefunction
get_determinantal_wavefunction
```

## Jastrow Types and Methods

- [`JastrowFactor`](@ref)
- [`get_jastrow_factor`](@ref)
- [`get_fermionic_Tvec`](@ref)
- [`update_fermionic_Tvec!`](@ref)
- [`get_fermionic_jastrow_ratio`](@ref)
- [`map_jastrow_parameters`](@ref)
- [`check_deviation`](@ref)

```@docs
JastrowFactor
get_jastrow_factor
get_fermionic_Tvec
update_fermionic_Tvec!
get_fermionic_jastrow_ratio
map_jastrow_parameters
check_deviation
```

## Markov Methods

- [`local_fermion_update!`](@ref)
- [`metropolis_step`](@ref)


```@docs
local_fermion_update!
metropolis_step
```

## ParticleConfiguration Types and Methods

- [`MarkovMove`](@ref)
- [`propose_random_move`](@ref)
- [`hop!`](@ref)
- [`exchange!`](@ref)
- [`generate_initial_fermion_configuration`](@ref)
- [`get_onsite_fermion_occupation`](@ref)
- [`get_particle_numbers`](@ref)
- [`get_particle_density`](@ref)
- [`get_spindex_type`](@ref)
- [`get_index_from_spindex`](@ref)
- [`get_spindices_from_index`](@ref)
- [`get_linked_spindex`](@ref)

```@docs
MarkovMove
propose_random_move
hop!
exchange!
generate_initial_fermion_configuration
get_onsite_fermion_occupation
get_particle_numbers
get_particle_density
get_spindex_type
get_index_from_spindex
get_spindices_from_index
get_linked_spindex
```

## Greens Methods

- [`initialize_equal_time_greens`](@ref)
- [`update_equal_time_greens!`](@ref)
- [`rank1_update!`](@ref)
- [`recalculate_equal_time_greens`](@ref)
- [`check_deviation`](@ref)

```@docs
initialize_equal_time_greens
update_equal_time_greens!
rank1_update!
recalculate_equal_time_greens
check_deviation
```

## Optimizer Methods

- [`stochastic_reconfiguration!`](@ref)
- [`get_Δk`](@ref)
- [`get_covariance_matrix`](@ref)
- [`get_force_vector`](@ref)

```@docs
stochastic_reconfiguration!
get_Δk
get_covariance_matrix
get_force_vector
```

<!-- # Measurement Methods

### Intitialize Measurements

- [`initialize_measurement_container`](@ref)
- [`initialize_measurement_directories`](@ref)


```@docs
initialize_measurement_container
initialize_measurement_directories
```

### Scalar Measurements

- [`get_local_energy`](@ref)
- [`get_local_kinetic_energy`](@ref)
- [`get_local_hubbard_energy`](@ref)
- [`get_double_occ`](@ref)
- [`get_n`](@ref)

```@docs
get_local_energy
get_local_kinetic_energy
get_local_hubbard_energy
get_double_occ
get_n
```

### Optimization Measurements

- [`measure_parameters!`](@ref)
- [`measure_Δk!`](@ref)
- [`measure_ΔkΔkp!`](@ref)
- [`measure_ΔkE!`](@ref)

```@docs
measure_parameters!
measure_Δk!
measure_ΔkΔkp!
measure_ΔkE!
```

### Simulation Measurements

- [`measure_local_energy!`](@ref)
- [`measure_double_occ!`](@ref)
- [`measure_n!`](@ref)

```@docs
measure_local_energy!
measure_double_occ!
measure_n!
```

### Make Measurements

- [`make_measurements!`](@ref)
- [`reset_measurements!`](@ref)

```@docs
make_measurements!
reset_measurements!
```

### Write Measurements

- [`write_measurements!`](@ref)

```@docs
write_measurements!
```

### Process Measurements

<!-- - [``](@ref)

```@docs

``` --> -->






