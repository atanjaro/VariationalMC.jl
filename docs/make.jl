using Pkg
Pkg.activate("..")  # Activate the package environment
Pkg.instantiate()   # Ensure all dependencies are installed

using Documenter
using VariationalMC

makedocs(
    sitename = "VariationalMC.jl",
    modules = [VariationalMC],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
)
