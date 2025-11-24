@doc raw"""

    optimize_parameters!( measurement_container::NamedTuple, 
                          determinantal_parameters::DeterminantalParameters{I}, 
                          η::E, 
                          dt::E, 
                          bin_size::I ) where {I<:Integer, E<:AbstractFloat}

Optimizes and updates variational parameters using the Stochastic Reconfiguration method, in
which the variation in parameter ``\alpha_k`` can be found by solving ``\delta\alpha_k = (S + \eta\mathbb{I})f^{-1}``,
where ``S`` is the covariance matrix, ``f`` is the force vector, and ``\eta`` is a stabilization parameter.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `η::E`: optimization stabilization factor. 
- `dt::E`: optimization rate. 
- `bin_size::I`: length of the bins. 

"""
function optimize_parameters!(
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters{I}, 
    η::E, 
    dt::E, 
    opt_bin_size::I
) where {I<:Integer, E<:AbstractFloat}
    @debug """
    Optimizer::optimize_parameters!() :
    Start of optimization!
    """

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
    δvpars = (S + η * LinearAlgebra.I(size(S,1))) \ f  

    # new variational parameters
    new_vpars = reduce(vcat, [x isa AbstractVector ? x : [x] for x in measurement_container.optimization_measurements["parameters"]])
    new_vpars += dt .* δvpars               
    
    @debug """
    Optimizer::optimize_parameters!() :
    Parameters have been updated =>
    new parameters = $(new_vpars)
    """

    # push back determinantal parameters
    update_parameters!(
        new_vpars, 
        determinantal_parameters
    )

    @debug """
    Optimizer::optimize_parameters!() :
    End of optimization!
    """

    return nothing
end


@doc raw"""

    optimize_parameters!( measurement_container::NamedTuple, 
                          determinantal_parameters::DeterminantalParameters{I}, 
                          jastrow_parameters::JastrowParameters{T, K, V, I}, 
                          η::E, 
                          dt::E, 
                          dt_J::E,
                          bin_size::I ) where {I<:Integer, T<:AbstractString, K, V, E<:AbstractFloat}

Optimizes and updates variational parameters using the Stochastic Reconfiguration method, in
which the variation in parameter ``\alpha_k`` can be found by solving ``\delta\alpha_k = (S + \eta\mathbb{I})f^{-1}``,
where ``S`` is the covariance matrix, ``f`` is the force vector, and ``\eta`` is a stabilization parameter.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{T, K, V, I}`: set of Jastrow variational parameters. 
- `η::E`: optimization stabilization factor. 
- `dt::E`: optimization rate. 
- `dt_J::E`: boost in the optimization rate of the Jastrow parameters.
- `bin_size::I`: length of the bins. 

"""
function optimize_parameters!( 
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{T, K, V, I}, 
    η::E, 
    dt::E, 
    dt_J::E,
    opt_bin_size::I
) where {I<:Integer, T<:AbstractString, K, V, E<:AbstractFloat}
    @debug """
    Optimizer::optimize_parameters!() :
    Start of optimization!
    """

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
    δvpars = (S + η * LinearAlgebra.I(size(S,1))) \ f  

    # new variational parameters
    new_vpars = reduce(vcat, [x isa AbstractVector ? x : [x] for x in measurement_container.optimization_measurements["parameters"]])

    # apply Jastrow convergence boost
    new_vpars[determinantal_parameters.num_det_pars+1:end] *= dt_J

    # update parameters
    new_vpars += dt * δvpars

    @debug """
    Optimizer::optimize_parameters!() :
    Parameters have been updated =>
    new parameters = $(new_vpars)
    """

    # push back determinantal_parameters
    update_parameters!(
        new_vpars, 
        determinantal_parameters,
        jastrow_parameters
    )

    @debug """
    Optimizer::optimize_parameters!() :
    End of optimization!
    """

    return nothing
end


@doc raw"""

    optimize_parameters!( measurement_container::NamedTuple, 
                          determinantal_parameters::DeterminantalParameters{I}, 
                          jastrow_parameters_1::JastrowParameters{T, K, V, I},
                          jastrow_parameters_2::JastrowParameters{T, K, V, I}, 
                          η::E, 
                          dt::E, 
                          bin_size::I ) where {I<:Integer, T<:AbstractString, K, V, E<:AbstractFloat}

Optimizes and updates variational parameters using the Stochastic Reconfiguration method, in
which the variation in parameter ``\alpha_k`` can be found by solving ``\delta\alpha_k = (S + \eta\mathbb{I})f^{-1}``,
where ``S`` is the covariance matrix, ``f`` is the force vector, and ``\eta`` is a stabilization parameter.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `jastrow_parameters_1::JastrowParameters{T, K, V, I}`: set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{T, K, V, I}`: set of Jastrow variational parameters. 
- `η::E`: optimization stabilization factor. 
- `dt::E`: optimization rate. 
- `bin_size::I`: length of the bins. 

"""
function optimize_parameters!( 
    measurement_container::NamedTuple, 
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{T, K, V, I},
    jastrow_parameters_2::JastrowParameters{T, K, V, I}, 
    η::E, 
    dt::E, 
    opt_bin_size::I
) where {I<:Integer, T<:AbstractString, K, V, E<:AbstractFloat}
    @debug """
    Optimizer::optimize_parameters!() :
    Start of optimization!
    """

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
    δvpars = (S + η * LinearAlgebra.I(size(S,1))) \ f  

    # new varitaional parameters
    new_vpars = reduce(vcat, [x isa AbstractVector ? x : [x] for x in measurement_container.optimization_measurements["parameters"]])
    new_vpars += dt * δvpars 

    @debug """
    Optimizer::optimize_parameters!() :
    Parameters have been updated =>
    new parameters = $(new_vpars)
    """

    # push back determinantal_parameters
    update_parameters!(
        new_vpars, 
        determinantal_parameters,
        jastrow_parameters_1,
        jastrow_parameters_2
    )

    @debug """
    Optimizer::optimize_parameters!() :
    End of optimization!
    """

    return nothing
end


@doc raw"""

    get_Δk( optimize::NamedTuple, 
            determinantal_parameters::DeterminantalParameters{I}, 
            detwf::DeterminantalWavefunction{T, V, E, I}, 
            model_geometry::ModelGeometry, 
            Np::I ) where {I<:Integer, T<:Number, V, E<:AbstractFloat}

Calculates the local logarithmic derivative ``\Delta_k = \frac{\partial\ln\Psi_{T}}{\partial\alpha_k}``, 
for each determinantal variational parameter ``\alpha_k`` and returns a vector of derivatives.

- `optimize::NamedTuple`: field of optimization flags.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `detwf::DeterminantalWavefunction{T, V, E, I}`: current variational wavefunction. 
- `model_geometry::ModelGeometry: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.

"""
function get_Δk(
    optimize::NamedTuple, 
    determinantal_parameters::DeterminantalParameters{I}, 
    detwf::DeterminantalWavefunction{T, V, E, I}, 
    model_geometry::ModelGeometry, 
    Np::I
) where {I<:Integer, T<:Number, V, E<:AbstractFloat}
    # number of lattice sites
    N = model_geometry.unit_cell.n * model_geometry.lattice.N

    # current determinantal_parameters
    det_pars = determinantal_parameters.det_pars

    # number of determinantal parameters 
    num_det_pars = determinantal_parameters.num_det_pars

    # resultant derivatives
    result = zeros(Float64, num_det_pars)

    @assert any(values(optimize))
    G = zeros(Complex, 2*N, 2*N)
    for β in 1:Np
        k = findfirst(x -> x == β, detwf.pconfig)

        G[k, :] .= detwf.W[:, β]
    end

    # perform derivatives    
    opt_idx = 1
    res_idx = 1
    for (pname, pval) in pairs(det_pars)
        if getfield(optimize, pname)
            if isa(pval, AbstractVector)
                for _ in eachindex(pval)
                    result[res_idx] = sum(detwf.A[opt_idx] .* G)
                    opt_idx += 1
                    res_idx += 1
                end
            else
                result[res_idx] = sum(detwf.A[opt_idx] .* G)
                opt_idx += 1
                res_idx += 1
            end
        else
            # skip over these but still count them
            res_idx += isa(pval, AbstractVector) ? length(pval) : 1
        end
    end
    
    return result
end


@doc raw"""

    get_Δk( optimize::NamedTuple, 
            jastrow_parameters::JastrowParameters{T, K, V, I},
            detwf::DeterminantalWavefunction{Q, F, E, I}, 
            pht::Bool ) where {T<:AbstractString, K, V, I<:Integer, Q<:Number, F, E<:AbstractFloat}

Calculates the local logarithmic derivative ``\Delta_k = \frac{\partial\ln\Psi_{T}}{\partial\alpha_k}``, 
for each determinantal variational parameter ``\alpha_k`` and returns a vector of derivatives.

- `optimize::NamedTuple`: field of optimization flags.
- `jastrow_parameters::JastrowParameters{T, K, V, I}`: current set of Jastrow variational parameters. 
- `detwf::DeterminantalWavefunction{Q, F, E, I}`: current variational wavefunction. 
- `pht::Bool`: whether model is particle-hole transformed.

"""
function get_Δk(
    optimize::NamedTuple, 
    jastrow_parameters::JastrowParameters{T, K, V, I},
    detwf::DeterminantalWavefunction{Q, F, E, I}, 
    pht::Bool
) where {T<:AbstractString, K, V, I<:Integer, Q<:Number, F, E<:AbstractFloat}
    # type of Jastrow parameters
    jastrow_type = jastrow_parameters.jastrow_type

    # map of current Jastrow parameters
    jpar_map = jastrow_parameters.jpar_map

    # number of Jastrow parameters
    num_jpars = jastrow_parameters.num_jpars

    # resultant derivatives
    result = zeros(num_jpars)

    jpar = 0
    for (irr_index, _) in jpar_map
        jpar += 1
        if jpar < num_jpars
            for (i,j) in jpar_map[irr_index][1]
                # current site occupations
                n_i_up, n_i_dn = get_onsite_fermion_occupation(i+1, detwf.pconfig)[1:2]
                n_j_up, n_j_dn = get_onsite_fermion_occupation(j+1, detwf.pconfig)[1:2]

                if jastrow_type == "e-den-den"
                    # calculate derivative
                    if pht
                        result[jpar] += 0.5 * (n_i_up - n_i_dn) * (n_j_up - n_j_dn)
                    else
                        result[jpar] += 0.5 * (n_i_up + n_i_dn) * (n_j_up + n_j_dn)
                    end
                elseif jastrow_type == "e-spn-spn"
                    # calculate derivative
                    if pht
                        result[jpar] += 0.25 * (n_i_up - n_i_dn) * (n_j_up - n_j_dn)
                    else
                        result[jpar] += 0.25 * (n_i_up + n_i_dn) * (n_j_up + n_j_dn)
                    end
                end
            end
        end
    end

    return result
end


@doc raw"""

    get_covariance_matrix( measurement_container::NamedTuple, 
                           opt_bin_size::I ) where {I<:Integer}

Calculates the covariance matrix ``S`` with elements ``S_{kk}' = \langle \Delta_{k}\Delta_{k^\prime}\rangle - \langle \Delta_{k} \rangle\langle \Delta_{k^\prime} \rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `opt_bin_size::I`: length of the current bin.

"""
function get_covariance_matrix(
    measurement_container::NamedTuple, 
    opt_bin_size::I
) where {I<:Integer}
    # measure local parameters derivatives ⟨Δₖ⟩ for the current bin
    Δk = measurement_container.optimization_measurements["Δk"]/opt_bin_size
    
    # measure the product of local derivatives ⟨ΔₖΔₖ'⟩ for the current bin
    ΔkΔkp = measurement_container.optimization_measurements["ΔkΔkp"]/opt_bin_size
    
    # calculate the product of local derivatives ⟨Δₖ⟩⟨Δₖ'⟩
    ΔkΔk = Δk * Δk'  

    # get S
    S = ΔkΔkp - ΔkΔk
    
    return S
end


@doc raw"""

    get_force_vector( measurement_container::NamedTuple, 
                      opt_bin_size::I ) where {I<:Integer}

Constructs force vector ``f`` with elements ``f_k = \langle \Delta_{k} \rangle\langle H\rangle - \langle \Delta_{k}H\rangle``.

- `measurement_container::NamedTuple`: container where measurements are stored. 
- `opt_bin_size::I`: length of the current bin.

"""
function get_force_vector(
    measurement_container::NamedTuple, 
    opt_bin_size::I
) where {I<:Integer}
    # initialize force vector
    f = Float64[]

    # measure local parameters derivatives ⟨Δₖ⟩ for the current bin
    Δk = measurement_container.optimization_measurements["Δk"]/opt_bin_size

    # measure local energy E = ⟨H⟩ for the current bin
    E = measurement_container.simulation_measurements["local_energy"]/opt_bin_size 

    # measure product of local derivatives with energy ⟨ΔkE⟩ for the current bin
    ΔkE = measurement_container.optimization_measurements["ΔkE"]/opt_bin_size

    # calculate product of local derivative with the local energy ⟨Δk⟩⟨H⟩
    ΔktE = Δk * E
    
    for (i, j) in zip(ΔktE, ΔkE)
        fk = i - j
        push!(f, fk)
    end

    return f 
end

