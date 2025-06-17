@doc raw"""

    DeterminantalParameters( pars::Vector{AbstractString}, 
                             vals::Vector{AbstractFloat}, 
                             num_detpars::Int )

A type defining a set of variational parameters for the determinantal wavefunction.

"""
mutable struct DeterminantalParameters
    # determinantal parameters and their values
    det_pars::NamedTuple

    # total number of determinantal parameters
    num_det_pars::Int

    # total number of determinantal parameters being optimized
    num_det_opts::Int
end


@doc raw"""

    TightBindingModel( t₀::Float64,
                       t₁::Float64 )

A type defining a non-interacting tight binding model.

"""
struct TightBindingModel    
    # nearest neighbor hopping amplitude
    t₀::Float64

    # next nearest neighbor hopping amplitude
    t₁::Float64
end


@doc raw"""

    SpinModel( J₁::Float64,
               J₂::Float64,
               J₃::Float64 )

A type defining a spin model.

"""
struct SpinModel
    # nearest neighbor spin exchange coupling 
    J₁::Float64

    # next nearest neighbor spin exchange coupling
    J₂::Float64

    # next next nearest neighbor spin exchange coupling
    J₃::Float64
end


@doc raw"""

    DeterminantalParameters( optimize::NamedTuple, 
                             tight_binding_model::TightBindingModel, 
                             model_geometry::ModelGeometry, 
                             minabs_vpar::Float64, 
                             Ne::Int, 
                             pht::Bool )::DeterminantalParameters

Given an intial set of parameters and optimization flags, generates a set of variational parameters.

- `optimize::NamedTuple`: field of optimization flags.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice qunatities.
- `minabs_vpar::Float64`: minimum value of initialized variational parameters.
- `Ne::Int`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function DeterminantalParameters(
    optimize::NamedTuple, 
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry, 
    minabs_vpar::Float64, 
    Ne::Int, 
    pht::Bool
)::DeterminantalParameters
    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # x-dimension
    Lx = model_geometry.lattice.L[1]

    # number of lattice sites
    N = model_geometry.lattice.N

    if dims > 1
        if pht
            det_pars = (
                Δ_0 = minabs_vpar,
                Δ_spd = fill(0.0, N),
                Δ_d = 0.0,
                Δ_dpd = fill(0.0, N),
                q_p = fill(0.0, dims),
                Δ_sx = 0.0,
                Δ_sz = minabs_vpar,
                Δ_ssd = fill(0.0, Lx),
                μ = 3.0,
                Δ_cdw = 0.0,
                Δ_csd = fill(0.0, Lx)
            )
        else
            det_pars = (
                Δ_sx = 0.0,
                Δ_sz = minabs_vpar,
                Δ_ssd = fill(0.0, Lx),
                μ = get_tb_chem_pot(Ne, tight_binding_model, model_geometry),
                Δ_cdw = 0.0,
                Δ_csd = fill(0.0, Lx),
            )
        end
    else
        if pht
            det_pars = (
                Δ_0 = minabs_vpar,
                Δ_sx = 0.0,
                Δ_sz = minabs_vpar,
                μ = 0.0,
                Δ_cdw = 0.0,
            )
        else
            det_pars = (
                Δ_sx = 0.0,
                Δ_sz = minabs_vpar,
                μ = get_tb_chem_pot(Ne, tight_binding_model, model_geometry),
                Δ_cdw = 0.0,
            )
        end
    end

    # determine total number of determinantal parameters being added to the model
    num_det_pars = sum(x -> isa(x, AbstractArray) ? length(x) : 1, values(det_pars))

    # determine the number of determinantal parameters being optimized
    opt_keys = intersect(keys(optimize), keys(det_pars))
    num_det_opts = sum(opt_keys) do key
        opt = getfield(optimize, key)
        val = getfield(det_pars, key)
        opt ? (isa(val, AbstractArray) ? length(val) : 1) : 0
    end


    debug && println("Hamiltonian::DeterminantalParameters() : ")
    debug && println("Number of determinantal parameters = $num_det_pars")
    debug && println("Number of determinantal parameters to be optimized = $num_det_opts")

    return DeterminantalParameters(det_pars, num_det_pars, num_det_opts)
end


@doc raw"""

    DeterminantalParameters( optimize::NamedTuple, 
                             model_geometry::ModelGeometry, 
                             pht::Bool, 
                             path_to_parameter_file::String )::DeterminantalParameters

Given an intial set of parameters and set of optimization flags, generates a set of variational parameters
from an intial parameter file.

- `optimize::NamedTuple`: field of optimization flags
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.
- `path_to_parameter_file::String`: filepath to initial parameter file.


"""
function DeterminantalParameters(
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    pht::Bool, 
    path_to_parameter_file::String
)::DeterminantalParameters
    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # get parameters from file
    vpar_dict = readin_parameters(path_to_parameter_file)

    if dims > 1
        if pht
            det_pars = (
                Δ_0 = vpar_dict[:pairing0],
                Δ_spd = vpar_dict[:spd],
                Δ_d = vpar_dict[:pairingd], 
                Δ_dpd = vpar_dict[:dpd],
                q_p = vpar_dict[:q_p],
                Δ_sx = vpar_dict[:sx],
                Δ_sz = vpar_dict[:sz],
                Δ_ssd = vpar_dict[:ssd],
                μ = vpar_dict[:chemical_potential],
                Δ_cdw = vpar_dict[:cdw],
                Δ_csd = vpar_dict[:csd]
            )
        else
            det_pars = (
                Δ_sx = vpar_dict[:sx],
                Δ_sz = vpar_dict[:sz],
                Δ_ssd = vpar_dict[:ssd],
                μ = vpar_dict[:chemical_potential],
                Δ_cdw = vpar_dict[:cdw],
                Δ_csd = vpar_dict[:csd]
            )
        end
    else
        if pht
            det_pars = (
                Δ_0 = vpar_dict[:pairing0],
                Δ_sx = vpar_dict[:sx],
                Δ_sz = vpar_dict[:sz],
                μ = vpar_dict[:chemical_potential],
                Δ_cdw = vpar_dict[:cdw],
            )
        else
            det_pars = (
                Δ_sx = vpar_dict[:sx],
                Δ_sz = vpar_dict[:sz],
                μ = vpar_dict[:chemical_potential],
                Δ_cdw = vpar_dict[:cdw],
            )
        end
    end

    # determine total number of determinantal parameters being added to the model
    num_det_pars = sum(x -> isa(x, AbstractArray) ? length(x) : 1, values(det_pars))

    # determine the number of determinantal parameters being optimized
    opt_keys = intersect(keys(optimize), keys(det_pars))
    num_det_opts = sum(opt_keys) do key
        opt = getfield(optimize, key)
        val = getfield(det_pars, key)
        opt ? (isa(val, AbstractArray) ? length(val) : 1) : 0
    end

    debug && println("Hamiltonian::DeterminantalParameters() : ")
    debug && println("Number of determinantal parameters = $num_det_pars")
    debug && println("Number of determinantal parameters to be optimized = $num_det_opts")

    return DeterminantalParameters(det_pars, num_det_pars, num_det_opts)
end


@doc raw"""

    JastrowParameters( jastrow_type::String,
                       jpar_map::OrderedDict{Any, Any},
                       num_jpars::Int,
                       num_jpar_opts::Int )

A type defining quantities related to Jastrow variational parameters.

"""
mutable struct JastrowParameters
    # type of Jastrow parameter
    jastrow_type::String

    # map of Jastrow parameters
    jpar_map::OrderedDict{Any, Any}

    # total number of Jastrow parameters
    num_jpars::Int

    # number of Jastrow parameters to be optimized
    num_jpar_opts::Int
end


@doc raw"""

    JastrowParameters( jastrow_type::String, 
                       optimize::NamedTuple,
                       model_geometry::ModelGeometry,
                       rng::Xoshiro )::JastrowParameters

Given a type of Jastrow factor and set of optimization flags, generates a random initial set of 
Jastrow parameters.

- `jastrow_type::String`: type of Jastrow factor: "e-den-den", "e-spn-spn". TBA: "eph-den-den", "ph-den-den"
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `rng::Xoshiro`: random number generator.

"""
function JastrowParameters(
    jastrow_type::String, 
    optimize::NamedTuple,
    model_geometry::ModelGeometry,
    rng::Xoshiro
)::JastrowParameters
    # create map of Jastrow parameters
    jpar_map = map_jastrow_parameters(
        model_geometry, 
        rng
    )

    # get total number of Jastrow parameters
    num_jpars = length(jpar_map)

    if optimize.djastrow && jastrow_type == "e-den-den"
        num_jpar_opts = num_jpars - 1
    elseif optimize.sjastrow && jastrow_type == "e-spn-spn"
        num_jpar_opts = num_jpars - 1
    # elseif optimize.ephjastrow && jastrow_type == "eph-den-den"
    # elseif optimize.phjastrow && jastrow_type == "ph-den-den"
    else
        num_jpar_opts = 0
    end

    debug && println("Jastrow::JastrowParameters() : type: ", jastrow_type)
    debug && println("Number of Jastrow parameters = ", num_jpars)
    debug && println("Number of Jastrow parameters to be optimized = ", num_jpar_opts)

    return JastrowParameters(jastrow_type, jpar_map, num_jpars, num_jpar_opts)
end


@doc raw"""

    JastrowParameters( jastrow_type::String, 
                       optimize::NamedTuple,
                       model_geometry::ModelGeometry,
                       path_to_parameter_file::String )::JastrowParameters

Given a type of Jastrow factor and set of optimization flags, generates an initial set of 
Jastrow parameters from a specfied file. 

- `jastrow_type::String`: type of Jastrow factor: "e-den-den", "e-spn-spn". TBA: "eph-den-den", "ph-den-den"
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `path_to_parameter_file::String`: filepath to initial parameters.

"""
function JastrowParameters(
    jastrow_type::String, 
    optimize::NamedTuple,
    model_geometry::ModelGeometry,
    path_to_parameter_file::String
)::JastrowParameters
    # get parameters from file
    vpar_dict = readin_parameters(path_to_parameter_file)

    if jastrow_type == "e-den-den"
        init_jpars = vpar_dict[:density_jastrow]
    elseif jastrow_type == "e-spn-spn"
        init_jpars = vpar_dict[:spin_jastrow]
    end

    # create map of Jastrow parameters 
    jpar_map = map_jastrow_parameters(
        model_geometry, 
        init_jpars
    ) 

    # get total number of Jastrow parameters
    num_jpars = length(jpar_map)

    if optimize.djastrow && jastrow_type == "e-den-den"
        num_jpar_opts = num_jpars - 1
    elseif optimize.sjastrow && jastrow_type == "e-spn-spn"
        num_jpar_opts = num_jpars - 1
    # elseif optimize.ephjastrow && jastrow_type == "eph-den-den"
    # elseif optimize.phjastrow && jastrow_type == "ph-den-den"
    else
        num_jpar_opts = 0
    end

    debug && println("Jastrow::JastrowParameters() : type: ", jastrow_type)
    debug && println("Number of Jastrow parameters = ", num_jpars)
    debug && println("Number of Jastrow parameters to be optimized = ", num_jpar_opts)

    return JastrowParameters(jastrow_type, jpar_map, num_jpars, num_jpar_opts)
end


@doc raw"""

    collect_parameters( determinantal_parameters::DeterminantalParameters, 
                        jastrow_parameters::JastrowParameters )::Vector{AbstractFloat}

Concatenates all values of determinantal and Jastrow parameters into a single vector.

"""
function collect_parameters(
    determinantal_parameters::DeterminantalParameters, 
    jastrow_parameters::JastrowParameters
)::Vector{AbstractFloat}
    # determinantal parameters
    det_pars = collect(values(determinantal_parameters.det_pars))

    # Jastrow parameters
    keys_sorted = sort(collect(keys(jastrow_parameters.jpar_map)))
    jpars = [jastrow_parameters.jpar_map[k][2] for k in keys_sorted[1:end-1]]

    return vcat(det_pars, jpars)
end


@doc raw"""

    update_parameters!( new_vpars::AbstractVector, 
                        determinantal_parameters::DeterminantalParameters )::Nothing

Updates variational (determinantal) parameters after Stochastic Reconfiguration.

"""
function update_parameters!(
    new_vpars::AbstractVector, 
    determinantal_parameters::DeterminantalParameters
)::Nothing
    # extract current parameters
    current_pars = determinantal_parameters.det_pars
    @assert length(new_vpars) != length(current_pars)

    # preserve parameter names from current_pars
    param_names = keys(current_pars)

    # build updated NamedTuple with same keys and new values
    new_det_pars = NamedTuple{Tuple(param_names)}(Tuple(new_vpars))

    # update the struct
    determinantal_parameters.det_pars = new_det_pars

    return nothing
end


@doc raw"""

    update_parameters!( new_vpars::AbstractVector, 
                        determinantal_parameters::DeterminantalParameters, 
                        jastrow_parameters::JastrowParameters )::Nothing

Updates variational (determinantal) parameters after Stochastic Reconfiguration.

"""
function update_parameters!(
    new_vpars::AbstractVector, 
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters::JastrowParameters
)::Nothing
    # extract current parameters from their containers
    current_det_pars = determinantal_parameters.det_pars
    current_jpar_map = jastrow_parameters.jpar_map

    # seperate the new determinantal and Jastrow parameters
    num_jpars = jastrow_parameters.num_jpar_opts

    new_det_pars = new_vpars[1:end-num_jpars]
    new_jpars = new_vpars[end-num_jpars+1:end]

    @assert length(new_vpars) != length(current_det_pars) + num_jpars
    
    # update determinantal parameters
    param_names = keys(current_det_pars)
    new_det_pars = NamedTuple{Tuple(param_names)}(Tuple(new_vpars))
    determinantal_parameters.det_pars = new_det_pars

    # update Jastrow parameters
    irr_indices = collect(keys(current_jpar_map))
    for i in 1:num_jpars
        indices, _ = current_jpar_map[irr_indices[i]]
    
        current_jpar_map[irr_indices[i]] = (indices, new_jpars[i])
    end

    return nothing
end


@doc raw"""

    readin_parameters( filename::String )

Parses TOML file containing initial variational parameters. 

- `filename::String`: name of parameter summary file in TOML format.

"""
function readin_parameters(
    filename::String
)
    toml_data = TOML.parsefile(filename)

    det_dict = get(toml_data, "DeterminantalParameters", Dict())
    jastrow_dict = get(toml_data, "JastrowParameters", Dict())
    det_pars_raw = get(det_dict, "det_pars", Dict())
    jpar_map = get(jastrow_dict, "jpar_map", Dict())

    # Mapping from TOML keys to internal symbols
    keymap = Dict(
        "Δ_sx" => :sx,
        "Δ_sz" => :sz,
        "Δ_0" => :pairing0,
        "Δ_d" => :pairingd,
        "Δ_spd" => :spd,
        "Δ_dpd" => :dpd,
        "Δ_ssd" => :ssd,
        "Δ_csd" => :csd,
        "μ" => :chemical_potential,
        "Δ_cdw" => :cdw,
        "q_p" => :q_p
    )

    # Construct `vpar_dict` from whatever exists in the TOML file
    vpar_dict = Dict{Symbol, Any}()

    for (toml_key, sym_key) in keymap
        if haskey(det_pars_raw, toml_key)
            val = det_pars_raw[toml_key]
            if startswith(toml_key, "Δ_") && toml_key == "Δ_0"
                # Store pairing as first element in vector
                vpar_dict[:pairing] = get(vpar_dict, :pairing, [])[1:0]
                push!(vpar_dict[:pairing], val)
            elseif toml_key == "Δ_d"
                vpar_dict[:pairing] = get(vpar_dict, :pairing, [])[1:1]
                push!(vpar_dict[:pairing], val)
            else
                vpar_dict[sym_key] = val
            end
        end
    end

    # # Also include jastrow if present
    # if !isempty(jpar_map)
    #     vpar_dict[:density_jastrow] = get(jpar_map, "(1, 2)", [])
    #     vpar_dict[:spin_jastrow] = get(jpar_map, "(3, 4)", [])
    # end

    return vpar_dict
end





# """

#     apply_tabc!( H::Matrix{ComplexF64}, Lx::Int, Ly::Int, θ::Tuple{Float64, Float64} )

# Given a Hamiltonian matrix H, phase θ, and lattice dimensions Lx and Ly, applies twisted boundary conditions. 

# """
# function apply_tbc!(H::Matrix{ComplexF64}, θ::Tuple{Float64, Float64}, Lx::Int, Ly::Int)
#     θx, θy = θ
#     N = Lx * Ly

#     # Apply the twist in the x-direction for spin-up and spin-down sectors
#     for y in 1:Lx
#         idx1 = Lx * (y - 1) + 1  # First element of row y
#         idx2 = Lx * y            # Last element of row y

#         # Spin-up sector
#         H[idx1, idx2] *= cis(θx) #exp(1im * θx)
#         H[idx2, idx1] *= cis(-θx) #exp(-1im * θx)

#         # Spin-down sector
#         H[idx1 + N, idx2 + N] *= cis(-θx) #exp(1im * -θx)
#         H[idx2 + N, idx1 + N] *= cis(θx) #exp(-1im * -θx)
#     end

#     # Apply the twist in the y-direction for spin-up and spin-down sectors
#     for x in 1:Ly
#         idx1 = x                  # Top element of column x
#         idx2 = Ly * (Ly - 1) + x  # Bottom element of column x

#         # Spin-up sector
#         H[idx1, idx2] *= cis(θy) #exp(1im * θy)
#         H[idx2, idx1] *= cis(-θy) #exp(-1im * θy)

#         # Spin-down sector
#         H[idx1 + N, idx2 + N] *= cis(-θy) #exp(1im * -θy)
#         H[idx2 + N, idx1 + N] *= cis(θy) #exp(-1im * -θy)
#     end
# end


# # Apply twist phases to the Hamiltonian
# θ = (rand(rng) * 2 * π, rand(rng) * 2 * π)  
# apply_tbc!(H, θ, Lx, Ly)