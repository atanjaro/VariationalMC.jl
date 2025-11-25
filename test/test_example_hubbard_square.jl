# integration test that runs a VMC simulation of a Hubbard model on a square lattice at half-filling
@testitem "Hubbard Square Example" begin
    include("../examples/hubbard_square.jl")
    sID = abs(rand(Int))
    L = 4
    U = 4.0
    density = 1.0
    N_equil = 2
    N_opt = 2
    N_opt_bins= 2
    N_sim = 2
    N_sim_bins = 2
    @test isnothing(run_hubbard_square_simulation(sID, L, U, density, N_equil, N_opt, N_opt_bins, N_sim, N_sim_bins, filepath = tempdir()))
end
