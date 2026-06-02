![VariationalMC logo](images/variationalmc-logo-dark.png)

Documentation for `VariationalMC.jl`. This package implements the Variational Monte Carlo (VMC) method for the Hubbard model and models of electron-phonon interaction (coming in version 1.1.0). 

**This code is currently the experimental stage of development. Use with caution.**

## Funding 

The development of this code was supported by the U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists, Office of Science Graduate Student Research (SCGSR) program. The SCGSR program is administered by the Oak Ridge Institute for Science and Education (ORISE) for the DOE. ORISE is managed by ORAU under contract number DE-SC0014664.

## Installation

**Julia version > 1.12.0 is required**

To install [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl),
simply open the Julia REPL and run the commands
```
julia> ]
pkg> add https://github.com/atanjaro/VariationalMC.jl
```
or equivalently via `Pkg` do
```
julia> using Pkg; Pkg.add(url="https://github.com/atanjaro/VariationalMC.jl")
```

## Documentation

- [`DEV`](https://atanjaro.github.io/VariationalMC.jl/dev/): Documentation associated with most recent commit to the main branch.

## Notable Package Dependencies

This section reviews some notable package dependencies.

### Re-exported Packages

[VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl) re-exports certain packages from [SmoQyDQMC.jl](https://github.com/SmoQySuite/SmoQyDQMC.jl) using the [Reexport.jl](https://github.com/simonster/Reexport.jl.git) package in order to simplify the installation process.

- [LatticeUtilties.jl](https://github.com/SmoQySuite/LatticeUtilities.jl.git): Used to represent arbitrary lattice geometries.

### External Dependencies

- [HDF5.jl](https://github.com/JuliaIO/HDF5.jl): Used for handling binary files.
- [JLD2.jl](https://github.com/JuliaIO/JLD2.jl.git): Used to write data to binary files in an HDF5 compatible format.

## Contact Us

For questions and comments regarding this package, please email Andy Tanjaroon Ly at [atanjaro@vols.utk.edu](mailto:atanjaro@vols.utk.edu) or [atanjaly82@gmail.com](mailto:atanjaly82@gmail.com).


