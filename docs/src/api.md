# API

## ModelGeometry Type and Methods

- [`ModelGeometry`](@ref)

```@docs
ModelGeometry
```

### Internal Methods

- [`VariationalMC.x`](@ref)
- [`VariationalMC.y`](@ref)
- [`VariationalMC.d`](@ref)
- [`VariationalMC.reduce_index_2d`](@ref)
- [`VariationalMC.reduce_index_1d`](@ref)
- [`VariationalMC.max_dist`](@ref)

```@docs
VariationalMC.x
VariationalMC.y
VariationalMC.d
VariationalMC.reduce_index_2d
VariationalMC.reduce_index_1d
VariationalMC.max_dist
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

### Internal Methods

- [`VariationalMC.collect_parameters`](@ref)
- [`VariationalMC.update_parameters!`](@ref)
- [`VariationalMC.readin_parameters`](@ref)

```@docs
VariationalMC.collect_parameters
VariationalMC.update_parameters!
VariationalMC.readin_parameters
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

```@docs
JastrowFactor
get_jastrow_factor
```

### Internal Methods

- [`VariationalMC.get_fermionic_Tvec`](@ref)
- [`VariationalMC.update_fermionic_Tvec!`](@ref)
- [`VariationalMC.get_fermionic_jastrow_ratio`](@ref)
- [`VariationalMC.map_jastrow_parameters`](@ref)

```@docs
VariationalMC.get_fermionic_Tvec
VariationalMC.update_fermionic_Tvec!
VariationalMC.get_fermionic_jastrow_ratio
VariationalMC.map_jastrow_parameters
```

## Markov Methods

- [`local_fermion_update!`](@ref)

```@docs
local_fermion_update!
```

### Internal Methods

- [`VariationalMC.metropolis_step`](@ref)

```@docs
VariationalMC.metropolis_step
```

## ParticleConfiguration Types and Methods

- [`get_particle_numbers`](@ref)
- [`get_particle_density`](@ref)

```@docs
get_particle_numbers
get_particle_density
```

### Internal Types and Methods

- [`VariationalMC.MarkovMove`](@ref)
- [`VariationalMC.propose_random_move`](@ref)
- [`VariationalMC.hop!`](@ref)
- [`VariationalMC.exchange!`](@ref)
- [`VariationalMC.generate_initial_fermion_configuration`](@ref)
- [`VariationalMC.get_onsite_fermion_occupation`](@ref)
- [`VariationalMC.get_spindex_type`](@ref)
- [`VariationalMC.get_index_from_spindex`](@ref)
- [`VariationalMC.get_spindices_from_index`](@ref)
- [`VariationalMC.get_linked_spindex`](@ref)

```@docs
VariationalMC.MarkovMove
VariationalMC.propose_random_move
VariationalMC.hop!
VariationalMC.exchange!
VariationalMC.generate_initial_fermion_configuration
VariationalMC.get_onsite_fermion_occupation
VariationalMC.get_spindex_type
VariationalMC.get_index_from_spindex
VariationalMC.get_spindices_from_index
VariationalMC.get_linked_spindex
```

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








