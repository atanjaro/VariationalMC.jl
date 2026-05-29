```@meta
CollapsedDocStrings = true
```

# API

## Simulation Information Type and Methods

```@docs
SimulationInfo
SimulationInfo(;)
initialize_datafolder
model_summary
```

## Model Geometry Type and Methods

```@docs
ModelGeometry
ModelGeometry(::UnitCell{D}, ::Lattice{D}) where {D}
add_bond!
get_bond_id
```

## Particle Configuration Type and Methods

```@docs
ParticleConfiguration
ParticleConfiguration(;)
```

## Tight-Binding Model and Parameters

```@docs
TightBindingModel
TightBindingModel(;)
TightBindingParameters
TightBindingParameters(;)
measure_hopping_energy
```

## Hubbard Model and Parameters

```@docs
HubbardModel
HubbardModel(;)
HubbardParameters
HubbardParameters(;)
ExtendedHubbardModel
ExtendedHubbardModel(;)
ExtendedHubbardParameters
ExtendedHubbardParameters(;)
measure_hubbard_energy
measure_ext_hubbard_energy
```

## Variational Parameters

- [Variational Parameter Names](@ref)
- [Determinantal Parameters](@ref)
- [Jastrow Parameters](@ref)

### Variational Parameter Names

```@docs
DETERMINANTAL_PARAMETERS
JASTROW_PARAMETERS
```

### Determinantal Parameters
```@docs
DeterminantalParameters
DeterminantalParameters(;)
add_parameter!
```

### Jastrow Parameters
```@docs
JastrowParameters
JastrowParameters(;)
```

## Monte Carlo

- [Variational Wavefunction](@ref)
- [Metropolis Updater](@ref)
- [Stochastic Reconfiguration](@ref)

### Variational Wavefunction
```@docs
DeterminantalWavefunction
DeterminantalWavefunction(;)
AbstractJastrowFactor
FermionJastrowFactor
JastrowFactor(;)
```

### Metropolis Updater
```@docs
local_fermion_update!
```

### Stochastic Reconfiguration 
```@docs
Optimizer
update_optimizer!
```

## Measurement Methods

- [Measurement Names](@ref)
- [Initialize Measurements](@ref)
- [Make Measurements](@ref)
- [Write Measurements](@ref)
- [Process Measurements](@ref)
- [Export Measurements](@ref)

### Measurement Names

```@docs
GLOBAL_MEASUREMENTS
LOCAL_MEASUREMENTS
CORRELATION_FUNCTIONS
```

### Initialize Measurements
```@docs
initialize_measurement_container
initialize_measurements!
initialize_correlation_measurements!
```

### Make Measurements
```@docs
make_measurements!
```

### Write Measurements
```@docs
write_measurements!
merge_bins
rm_bins
```

### Process Measurements
```@docs
save_simulation_info
process_measurements
```

### Export Measurements
```@docs
export_global_stats_to_csv
export_local_stats_to_csv
export_correlation_stats_to_csv
```