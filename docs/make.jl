using Pkg
Pkg.activate(@__DIR__)  # Activate the package environment
Pkg.instantiate()       # Ensure all dependencies are installed
push!(LOAD_PATH, "../src")

using Documenter
using VariationalMC

makedocs(
    sitename = "VariationalMC.jl",
    modules = [VariationalMC],
    format = Documenter.HTML(),
    checkdocs = :none,    
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Hubbard Chain" => "examples/hubbard_chain.md",
            "Hubbard Square" => "examples/hubbard_square.md",
            "Hubbard Square with MPI" => "examples/hubbard_square_mpi.md",
            "Hubbard Rectangle" => "examples/hubbard_rectangle.md"
    ],
        "Simulation Output Overview" => "simulation_output.md",
        "API" => "api.md"
    ],
)

deploydocs(
    repo = "github.com/atanjaro/VariationalMC.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
