# API

<!-- ## SimulationInfo Type and Methods

- [`SimulationInfo`](@ref)
- [`save_simulation_info`](@ref)
- [`initialize_datafolder`](@ref)
- [`create_datafolder_prefix`](@ref)
- [`model_summary`](@ref)

```@docs
SimulationInfo
VariationalMC.save_simulation_info
initialize_datafolder
create_datafolder_prefix
VariationalMC.model_summary
``` -->

## ModelGeometry Type and Methods

- [`ModelGeometry`](@ref)


```@docs
ModelGeometry
```

## Parameters Types and Methods

- [`TightBindingModel`](@ref)
- [`SpinModel`](@ref)
- [`DeterminantalParameters`](@ref)
- [`JastrowParameters`](@ref)


```@docs
TightBindingModel
SpinModel
DeterminantalParameters
JastrowParameters
```

<!-- ## Hamiltonian Methods

- [`build_auxiliary_hamiltonian`](@ref)
- [`build_tight_binding_hamiltonian`](@ref)
- [`build_variational_hamiltonian`](@ref)
- [`diagonalize`](@ref)
- [`is_openshell`](@ref)
- [`get_variational_matrices`](@ref)
- [`get_tb_chem_pot`](@ref)

```@docs
VariationalMC.build_auxiliary_hamiltonian
VariationalMC.build_tight_binding_hamiltonian
VariationalMC.build_variational_hamiltonian
VariationalMC.diagonalize
VariationalMC.is_openshell
VariationalMC.get_variational_matrices
VariationalMC.get_tb_chem_pot
``` -->

<!-- ### Adding different symmetries to the Hamiltonian

- [`add_pairing_symmetry!`](@ref)
- [`add_spin_order!`](@ref)
- [`add_charge_order!`](@ref)
- [`add_chemical_potential!`](@ref)

```@docs
VariationalMC.add_pairing_symmetry!
VariationalMC.add_spin_order!
VariationalMC.add_charge_order!
VariationalMC.add_chemical_potential!
``` -->

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
VariationalMC.get_fermionic_Tvec
VariationalMC.update_fermionic_Tvec!
VariationalMC.get_fermionic_jastrow_ratio
VariationalMC.map_jastrow_parameters
VariationalMC.check_deviation
```

## Markov Methods

- [`local_fermion_update!`](@ref)

```@docs
local_fermion_update!
```

## ParticleConfiguration Types and Methods

- [`get_particle_numbers`](@ref)
- [`get_particle_density`](@ref)

```@docs
get_particle_numbers
get_particle_density
```

<!-- ## Greens Methods

- [`initialize_equal_time_greens`](@ref)
- [`update_equal_time_greens!`](@ref)
- [`rank1_update!`](@ref)
- [`recalculate_equal_time_greens`](@ref)
- [`check_deviation`](@ref)

```@docs
VariationalMC.initialize_equal_time_greens
VariationalMC.update_equal_time_greens!
VariationalMC.rank1_update!
VariationalMC.recalculate_equal_time_greens
VariationalMC.check_deviation
``` -->

## Optimizer Methods

- [`stochastic_reconfiguration!`](@ref)

```@docs
stochastic_reconfiguration!
```

## Measurement Methods

### Intitialize Measurements

- [`initialize_measurement_container`](@ref)
- [`initialize_measurement_directories`](@ref)


```@docs
initialize_measurement_container
initialize_measurement_directories
```

<!-- ### Scalar Measurements

- [`get_local_energy`](@ref)
- [`get_local_kinetic_energy`](@ref)
- [`get_local_hubbard_energy`](@ref)
- [`get_double_occ`](@ref)
- [`get_n`](@ref)

```@docs
VariationalMC.get_local_energy
VariationalMC.get_local_kinetic_energy
VariationalMC.get_local_hubbard_energy
VariationalMC.get_double_occ
VariationalMC.get_n
``` -->

<!-- ### Optimization Measurements

- [`measure_parameters!`](@ref)
- [`measure_Δk!`](@ref)
- [`measure_ΔkΔkp!`](@ref)
- [`measure_ΔkE!`](@ref)

```@docs
VariationalMC.measure_parameters!
VariationalMC.measure_Δk!
VariationalMC.measure_ΔkΔkp!
VariationalMC.measure_ΔkE!
``` -->

<!-- ### Simulation Measurements

- [`measure_local_energy!`](@ref)
- [`measure_double_occ!`](@ref)
- [`measure_n!`](@ref)

```@docs
VariationalMC.measure_local_energy!
VariationalMC.measure_double_occ!
VariationalMC.measure_n!
``` -->

### Make Measurements

- [`make_measurements!`](@ref)

```@docs
make_measurements!
```

### Write Measurements

- [`write_measurements!`](@ref)

```@docs
write_measurements!
```

### Process Measurements








