# VariationalMC.jl

Documentation for `VariationalMC.jl`. This package impelments the variational Monte Carlo (VMC) method for Hubbard and electron-phonon interactions (coming soon).

This code is currently the experimental stage of development. Use with caution. 

## Funding 

The development of this code was supported by the U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists, Office of Science Graduate Student Research (SCGSR) program. The SCGSR program is administered by the Oak Ridge Institute for Science and Education (ORISE) for the DOE. ORISE is managed by ORAU under contract number DE-SC0014664.

## Installation

To install the [VariationalMC.jl](https://github.com/atanjaro/VariationalMC.jl),
simply open the Julia REPL and run the commands
```
julia> ]
pkg> add VariationalMC
```
or equivalently via `Pkg` do
```
julia> using Pkg; Pkg.add("VariationalMC")
```

## Notable External Package Dependencies

This section reviews some notable package dependencies.

- [LatticeUtilties.jl](https://github.com/SmoQySuite/LatticeUtilities.jl.git): Package that is used to represent arbitrary lattice geometries.
- [OrderedCollections.jl](https://github.com/JuliaCollections/OrderedCollections.jl): Package that implements associative containers that preserve the order of insertion.
- [JLD2.jl](https://github.com/JuliaIO/JLD2.jl.git): Package used to write data to binary files in an HDF5 compatible format. 
- [CSV.jl](https://github.com/JuliaData/CSV.jl): Package used for handling delimited text data.

## Contact

For questions and comments regarding this package, please email Andy Tanjaroon Ly at [atanjaro@vols.utk.edu](mailto:atanjaro@vols.utk.edu).

