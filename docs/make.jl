using Pkg
Pkg.activate(@__DIR__)  # Activate the package environment
Pkg.instantiate()       # Ensure all dependencies are installed

using Documenter
using VariationalMC

makedocs(
    sitename = "VariationalMC.jl",
    modules = [VariationalMC],
    format = Documenter.HTML(),
    checkdocs = :none,    
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
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
