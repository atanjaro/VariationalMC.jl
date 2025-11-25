# integration test that runs a VMC simulation of a Hubbard model on a chain at half-filling
@testitem "Hubbard Chain Example" begin
    include("../examples/hubbard_chain.jl")
    sID = abs(rand(Int))
    L = 4
    U = 2.0
    nup = 2
    ndn = 2
    N_equil = 2
    N_opt = 2
    N_opt_bins = 2
    N_sim = 2
    N_sim_bins = 2
    @test isnothing(run_hubbard_chain_simulation(sID, L, U, nup, ndn, N_equil, N_opt, N_opt_bins, N_sim, N_sim_bins, filepath = tempdir()))
end