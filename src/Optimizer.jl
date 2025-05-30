@doc raw"""

    stochastic_reconfiguration!( measurement_container::NamedTuple, 
                                 determinantal_parameters::DeterminantalParameters, 
                                 η::Float64, 
                                 dt::Float64, 
                                 bin::Int,
                                 bin_size::Int64 )::Nothing

Updates variational parameters through stochastic optimization.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `η::Float64`: optimization stabilization factor. 
- `dt::Float64`: optimization rate. 
- `bin::Int`: current bin number.
- `bin_size::Int64`: length of the current bin. 

"""
function stochastic_reconfiguration!(
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters, 
    η::Float64, 
    dt::Float64, 
    bin::Int, 
    opt_bin_size::Int64
)::Nothing
    debug && println("Optimizer::stochastic_reconfiguration!() : ")
    debug && println("Start of optimization")

    # get S matrix
    S = get_covariance_matrix(
        measurement_container, 
        opt_bin_size
    ) 

    # get f vector
    f = get_force_vector(
        measurement_container, 
        opt_bin_size
    )

    # solve for variation in the parameters
    δvpars = (S + η * I(size(S,1))) \ f  

    # new varitaional parameters
    new_vpars = reduce(vcat, (isa(x, AbstractVector) ? x : [x] for x in values(determinantal_parameters.det_pars))) 
    new_vpars .+= dt .* δvpars                                     

    debug && println("Optimizer::stochastic_reconfiguration!() : ")
    debug && println("Parameters have been updated")
    debug && println("new parameters = ", new_vpars)

    # push back determinantal parameters
    update_parameters!(
        new_vpars, 
        determinantal_parameters
    )

    # extract container info
    optimization_measurements = measurement_container.optimization_measurements

    (; datafolder, pID) = simulation_info

    # extract other container components
    (; optimization_measurements) = measurement_container

    # construct filenames
    fn = @sprintf "bin-%d_pID-%d.jld2" bin pID  

    # path to parameter directory
    file_path_det_parameters = joinpath(datafolder, "optimization", "determinantal", fn)

    # append parameter measurements to file
    det_parameter_measurements = optimization_measurements["parameters"][2][1:end-num_jpars] 
    JLD2.@save file_path_det_parameters det_parameter_measurements append=true

    debug && println("Optimizer::stochastic_reconfiguration!() : ")
    debug && println("End of optimization")

    return nothing
end


@doc raw"""

    stochastic_reconfiguration!( measurement_container::NamedTuple, 
                                 determinantal_parameters::DeterminantalParameters, 
                                 jastrow_parameters::JastrowParameters, 
                                 η::Float64, 
                                 dt::Float64, 
                                 bin::Int, 
                                 bin_size::Int64 )::Nothing

Updates variational parameters through stochastic optimization.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: set of Jastrow variational parameters. 
- `η::Float64`: optimization stabilization factor. 
- `dt::Float64`: optimization rate. 
- `bin::Int`: current bin number.
- `bin_size::Int64`: length of the current bin. 

"""
function stochastic_reconfiguration!( 
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters, 
    η::Float64, 
    dt::Float64, 
    bin::Int, 
    opt_bin_size::Int64
)::Nothing
    debug && println("Optimizer::stochastic_reconfiguration!() : ")
    debug && println("Start of optimization")

    # get S matrix
    S = get_covariance_matrix(
        measurement_container, 
        opt_bin_size
    ) 

    # get f vector
    f = get_force_vector(
        measurement_container, 
        opt_bin_size
    )

    # solve for variation in the parameters
    δvpars = (S + η * I(size(S,1))) \ f  

    # new varitaional parameters
    vpars = collect_parameters(determinantal_parameters, jastrow_parameters)
    vpars += dt * δvpars

    debug && println("Optimizer::stochastic_reconfiguration!() : ")
    debug && println("Parameters have been updated")
    debug && println("new parameters = ", vpars)

    # push back determinantal_parameters
    update_parameters!(
        new_vpars, 
        determinantal_parameters,
        jastrow_parameters
    )

    # extract container info
    optimization_measurements = measurement_container.optimization_measurements

    (; datafolder, pID) = simulation_info

    # extract other container components
    (; optimization_measurements) = measurement_container

    # construct filenames
    fn = @sprintf "bin-%d_pID-%d.jld2" bin pID  

    # path to parameter directories
    file_path_det_parameters = joinpath(datafolder, "optimization", "determinantal", fn)
    file_path_jas_parameters = joinpath(datafolder, "optimization", "Jastrow", fn)

    # seperate determinantal and Jastrow parameters  
    det_parameter_measurements = optimization_measurements["parameters"][2][1:end-num_jpars] 
    jas_parameter_measurements = optimization_measurements["parameters"][2][end-num_jpars+1:end]

    # append parameter measurements to file
    JLD2.@save file_path_det_parameters det_parameter_measurements append=true
    JLD2.@save file_path_jas_parameters jas_parameter_measurements append=true

    debug && println("Optimizer::stochastic_reconfiguration!() : ")
    debug && println("End of optimization")

    return nothing
end


@doc raw"""

    get_Δk( optimize::NamedTuple, 
            determinantal_parameters::DeterminantalParameters, 
            detwf::DeterminantalWavefunction, 
            model_geometry::ModelGeometry, 
            Ne::Int )::Vector{Float64}

Calculates the local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect to the kth variational parameter αₖ,
in the determinantal part of the wavefunction. Returns a vector of derivatives.

- `optimize::NamedTuple`: field of optimization flags.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `detwf::DeterminantalWavefunction`: current variational wavefunction. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int`: total number of electrons. 

"""
function get_Δk(
    optimize::NamedTuple, 
    determinantal_parameters::DeterminantalParameters, 
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry, 
    Ne::Int
)::Vector{Float64}
    # number of lattice sites
    N = model_geometry.unit_cell.n * model_geometry.lattice.N

    # current determinantal_parameters
    det_pars = determinantal_parameters.det_pars

    # total number of parameters in the model
    num_det_pars = determinantal_parameters.num_det_pars

    # # number of parameters being optimized
    # num_det_opts = determinantal_parameters.num_det_opts

    # resultant derivatives
    result = zeros(Float64, num_det_pars)

    # perform derivatives
    @assert any(values(optimize))
    G = zeros(Complex, 2*N, 2*N)
    for β in 1:Ne
        k = findfirst(x -> x == β, detwf.pconfig)

        G[k, :] .= detwf.W[:, β]
    end

    opt_idx = 1
    for (i, pname) in enumerate(keys(det_pars))
        if getfield(optimize, pname)
            result[i] = sum(detwf.A[opt_idx] .* G)
            opt_idx += 1
        end
    end

    return result
end


"""

    get_Δk( optimize::NamedTuple, 
            determinantal_parameters::DeterminantalParameters, 
            jastrow_parameters::JastrowParameters,
            detwf::DeterminantalWavefunction, 
            model_geometry::ModelGeometry, 
            Ne::Int,
            pht::Bool )::Vector{Float64}

Calculates the local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect to the kth variational parameter αₖ,
in the determinantal part of the wavefunction. Returns a vector of derivatives.

- `optimize::NamedTuple`: field of optimization flags.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: current set of Jastrow variational parameters. 
- `detwf::DeterminantalWavefunction`: current variational wavefunction. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Ne::Int`: total number of electrons. 
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_Δk(
    optimize::NamedTuple, 
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters,
    detwf::DeterminantalWavefunction, 
    model_geometry::ModelGeometry, 
    Ne::Int,
    pht::Bool
)::Vector{Float64}
    # number of lattice sites
    N = model_geometry.unit_cell.n * model_geometry.lattice.N

    # current determinantal_parameters
    det_pars = determinantal_parameters.det_pars

    # total number of parameters in the model
    num_det_pars = determinantal_parameters.num_det_pars

    # # number of determinantal parameters being optimized
    # num_det_opts = determinantal_parameters.num_det_opts

    # type of Jastrow parameters
    jastrow_type = jastrow_parameters.jastrow_type

    # map of current Jastrow parameters
    jpar_map = jastrow_parameters.jpar_map

    # # total number of Jastrow parameters
    # num_jpars = jastrow_parameters.num_jpars

    # number of Jastrow parametes being optimized
    num_jpars = jastrow_parameters.num_jpar_opts

    # resultant derivatives
    result = zeros(Float64, num_det_pars)

    # perform determinantal derivatives
    @assert any(values(optimize)) 
    G = zeros(Complex, 2*N, 2*N)
    for β in 1:Ne
        k = findfirst(x -> x == β, detwf.pconfig)

        G[k, :] .= detwf.W[:, β]
    end

    opt_idx = 1
    for (i, pname) in enumerate(keys(det_pars))
        if getfield(optimize, pname)
            result[i] = sum(detwf.A[opt_idx] .* G)
            opt_idx += 1
        end
    end

    # perform Jastrow derivatives
    @assert optimize.djastrow == true || optimize.sjastrow == true

    irr_indices = collect(keys(jpar_map))
    if jastrow_type == "e-den-den"
        for num in 1:num_jpars
            indices, _ = jpar_map[irr_indices[num]]

            for idx in indices
                i = idx[1]
                j = idx[2]

                # double counting correction
                dbc_corr = (j==i) ? 0.5 : 1.0

                # get occupations
                nup_i = get_onsite_fermion_occupation(i+1, detwf.pconfig)[1]
                ndn_i = get_onsite_fermion_occupation(i+1, detwf.pconfig)[2]
                nup_j = get_onsite_fermion_occupation(j+1, detwf.pconfig)[1]
                ndn_j = get_onsite_fermion_occupation(j+1, detwf.pconfig)[2]
                if pht
                    result[num_det_pars + num] += -dbc_corr * (nup_i - ndn_i) * (nup_j - ndn_j)
                else
                    result[num_det_pars + num] += -dbc_corr * (nup_i + ndn_i) * (nup_j + ndn_j)
                end
            end
        end
    elseif jastrow_type == "e-spn-spn"
        for num in 1:num_jpars
            indices, _ = jpar_map[irr_indices[num]]

            for idx in indices
                i = idx[1]
                j = idx[2]

                # double counting correction
                dbc_corr = (j==i) ? 0.5 : 1.0

                # get occupations
                nup_i = get_onsite_fermion_occupation(i+1, detwf.pconfig)[1]
                ndn_i = get_onsite_fermion_occupation(i+1, detwf.pconfig)[2]
                nup_j = get_onsite_fermion_occupation(j+1, detwf.pconfig)[1]
                ndn_j = get_onsite_fermion_occupation(j+1, detwf.pconfig)[2]
                if pht
                    result[num_det_pars + num] += -0.5 * dbc_corr * (nup_i - ndn_i) * (nup_j - ndn_j)
                else
                    result[num_det_pars + num] += -0.5 * dbc_corr * (nup_i + ndn_i) * (nup_j + ndn_j)
                end
            end
        end
    end

    return result
end


@doc raw"""

    get_covariance_matrix( measurement_container::NamedTuple, 
                           opt_bin_size::Int64 )

Calculates the covariance matrix S, for Stochastic Reconfiguration, with elements
S_kk' = <Δ_kΔk'> - <Δ_k><Δ_k'>.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `opt_bin_size::Int64`: length of the current bin.

"""
function get_covariance_matrix(
    measurement_container::NamedTuple, 
    opt_bin_size::Int64
)
    # measure local parameters derivatives ⟨Δₖ⟩ for the current bin
    Δk = measurement_container.optimization_measurements["Δk"][1]/opt_bin_size
    
    # measure the product of local derivatives ⟨ΔₖΔₖ'⟩ for the current bin
    ΔkΔkp = measurement_container.optimization_measurements["ΔkΔkp"][1]/opt_bin_size
    
    # calculate the product of local derivatives ⟨Δₖ⟩⟨Δₖ'⟩
    ΔkΔk = Δk * Δk'  

    # get S
    S = ΔkΔkp - ΔkΔk
    
    return S
end


@doc raw"""

    get_force_vector( measurement_container::NamedTuple, 
                      opt_bin_size::Int64 )

Generates the force vector f, for Stochastic Reconfiguration with elements 
f_k = <Δ_k><H> - <Δ_kH>.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `opt_bin_size::Int64`: length of the current bin.

"""
function get_force_vector(
    measurement_container::NamedTuple, 
    opt_bin_size::Int64
)
    # initialize force vector
    f = Float64[]

    # measure local parameters derivatives ⟨Δₖ⟩ for the current bin
    Δk = measurement_container.optimization_measurements["Δk"][1]/opt_bin_size

    # measure local energy E = ⟨H⟩ for the current bin
    E = measurement_container.simulation_measurements["energy"][1]/opt_bin_size

    # measure product of local derivatives with energy ⟨ΔkE⟩ for the current bin
    ΔkE = measurement_container.optimization_measurements["ΔkE"][1]/opt_bin_size

    # calculate product of local derivative with the local energy ⟨Δk⟩⟨H⟩
    ΔktE = Δk * E
    
    for (i, j) in zip(ΔktE, ΔkE)
        fk = i - j
        push!(f, fk)
    end

    return f 
end


