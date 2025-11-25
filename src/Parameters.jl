@doc raw"""

    DeterminantalParameters{I<:Integer} 

A type defining a set of variational parameters obtained from the auxiliary Hamiltonian, 
otherwise known as the "determinantal" part of the variational wavefunction.

- `det_pars::NamedTuple`: contains parameter names and their values.
- `num_det_pars::I`: total number of determinantal parameters.
- `num_det_opts::I`: number of determinantal parameters being optimized.

"""
mutable struct DeterminantalParameters{I<:Integer}
    # determinantal parameters and their values
    det_pars::NamedTuple

    # total number of determinantal parameters
    num_det_pars::I

    # total number of determinantal parameters being optimized
    num_det_opts::I
end


@doc raw"""

    TightBindingModel{E<:AbstractFloat}

A type defining a non-interacting tight binding model with nearest, next 
nearest, and third nearest hopping parameters.

- `t₀::E`: zeroth hopping parameter (nearest neighbor)
- `t₁::E`: first hopping parameter (next nearest neighbor)
- `t₂::E`: second hopping parameter (third nearest neighbor)

"""
struct TightBindingModel{E<:AbstractFloat}    
    # nearest neighbor hopping amplitude
    t₀::E

    # next nearest neighbor hopping amplitude
    t₁::E

    # third nearest neighbor hopping amplitude
    t₂::E
end


@doc raw"""

    SpinModel{E<:AbstractFloat}

A type defining a non-interacting spin model with nearest, next nearest, and 
third nearest neighbor exchange coupling parameters.

"""
struct SpinModel{E<:AbstractFloat}
    # nearest neighbor spin exchange coupling 
    J₁::E

    # next nearest neighbor spin exchange coupling
    J₂::E

    # third nearest neighbor spin exchange coupling
    J₃::E
end


@doc raw"""

    DeterminantalParameters( optimize::NamedTuple, 
                             tight_binding_model::TightBindingModel{E}, 
                             model_geometry::ModelGeometry, 
                             Ne::I, 
                             pht::Bool;
                             vpar_overrides::NT = NamedTuple() ) where {E<:AbstractFloat, I<:Integer}

Given an intial set of tight-binding and determinantal parameters, and optimization flags, 
generates a set of determinantal variational parameters.

- `optimize::NamedTuple`: field of optimization flags.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice qunatities.
- `Ne::I`: total number of electrons.
- `pht::Bool`: whether model is particle-hole transformed. 
- `vpar_overrides::NamedTuple = NamedTuple()`: user-defined overrides for initial variational parameters values.

"""
function DeterminantalParameters(
    optimize::NamedTuple, 
    tight_binding_model::TightBindingModel{E}, 
    model_geometry::ModelGeometry, 
    Ne::I, 
    pht::Bool;
    vpar_overrides::NamedTuple= NamedTuple()
) where {E<:AbstractFloat, I<:Integer}
    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # x-dimension
    Lx = model_geometry.lattice.L[1]

    # number of lattice sites
    N = model_geometry.lattice.N

    if dims > 1
        if pht
            det_pars = (
                Δ_0 = 0.01,
                Δ_spd = fill(0.01, N),
                Δ_d = 0.01,
                Δ_dpd = fill(0.01, N),
                q_p = fill(0.0, dims),
                Δ_sx = 0.01,
                Δ_sz = 0.01,
                Δ_ssd = fill(0.01, Lx),
                μ = get_tb_chem_pot(Ne, tight_binding_model, model_geometry),
                Δ_cdw = 0.01,
                Δ_csd = fill(0.01, Lx)
            )
        else
            det_pars = (
                Δ_sx = 0.01,
                Δ_sz = 0.01,
                Δ_ssd = fill(0.01, Lx),
                μ = get_tb_chem_pot(Ne, tight_binding_model, model_geometry),
                Δ_cdw = 0.01,
                Δ_csd = fill(0.01, Lx),
            )
        end
    else
        if pht 
            det_pars = (
                Δ_0 = 0.01, 
                Δ_sx = 0.01,
                Δ_sz = 0.01,
                μ = get_tb_chem_pot(Ne, tight_binding_model, model_geometry),
                Δ_cdw = 0.01,
            )
        else
            det_pars = (
                Δ_sx = 0.01,
                Δ_sz = 0.01,
                μ = get_tb_chem_pot(Ne, tight_binding_model, model_geometry),
                Δ_cdw = 0.01,
            )
        end
    end

    # override starting values (if provdided)
    det_pars = merge(det_pars, vpar_overrides)

    # determine total number of determinantal parameters being added to the model
    num_det_pars = sum(x -> isa(x, AbstractArray) ? length(x) : 1, values(det_pars))

    # determine the number of determinantal parameters being optimized
    opt_keys = intersect(keys(optimize), keys(det_pars))
    num_det_opts = sum(opt_keys) do key
        opt = getfield(optimize, key)
        val = getfield(det_pars, key)
        opt ? (isa(val, AbstractArray) ? length(val) : 1) : 0
    end

    @debug """
    Parameters::DeterminantalParameters() :
    Number of determinantal parameters = $(num_det_pars)
    Number of determinantal parameters to be optimized = $(num_det_opts)
    """

    return DeterminantalParameters(det_pars, num_det_pars, num_det_opts)
end


@doc raw"""

    DeterminantalParameters( optimize::NamedTuple, 
                             model_geometry::ModelGeometry, 
                             pht::Bool, 
                             path_to_parameter_file::S ) where {S<:AbstractString}

Given an intial set of determinantal parameters from file, and optimization flags,
generates a set of determinantal variational parameters

- `optimize::NamedTuple`: field of optimization flags
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.
- `path_to_parameter_file::S`: filepath to initial parameter file.

"""
function DeterminantalParameters(
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    pht::Bool, 
    path_to_parameter_file::S
) where {S<:AbstractString}
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

    @debug """
    Parameters::DeterminantalParameters() :
    Number of determinantal parameters = $(num_det_pars)
    Number of determinantal parameters to be optimized = $(num_det_opts)
    """

    return DeterminantalParameters(det_pars, num_det_pars, num_det_opts)
end


@doc raw"""

    JastrowParameters{S<:AbstractString, K, V, I<:Integer}

A type defining quantities related to Jastrow variational parameters.

- `jastrow_type::S`: type of Jastrow parameters: "e-den-den", "e-spn-spn"
- `jpar_map::OrderedDict{K, V}`: dictionary of irreducible indices to their correct index pairs and Jastrow parameter.
- `num_jpars::I`: total number of Jastrow parameters.
- `num_jpar_opts::I`: total number of Jastrow parameters being optimized.

"""
mutable struct JastrowParameters{S<:AbstractString, K, V, I<:Integer}
    # type of Jastrow parameter
    jastrow_type::S

    # map of Jastrow parameters
    jpar_map::OrderedDict{K, V}

    # total number of Jastrow parameters
    num_jpars::I

    # number of Jastrow parameters to be optimized
    num_jpar_opts::I
end


@doc raw"""

    JastrowParameters( jastrow_type::S, 
                       optimize::NamedTuple,
                       model_geometry::ModelGeometry,
                       rng::AbstractRNG ) where {S<:AbstractString}

Given a type of Jastrow factor and optimization flags, generates a random initial set of 
Jastrow variational parameters.

- `jastrow_type::S`: type of Jastrow factor: "e-den-den", "e-spn-spn". TBA: "eph-den-den", "ph-den-den"
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `rng:AbstractRNG`: random number generator.

"""
function JastrowParameters(
    jastrow_type::S, 
    optimize::NamedTuple,
    model_geometry::ModelGeometry,
    rng::AbstractRNG
) where {S<:AbstractString}
    # create map of Jastrow parameters
    jpar_map = map_jastrow_parameters(
        model_geometry, 
        rng
    )

    # get total number of Jastrow parameters
    num_jpars = length(jpar_map)

    if optimize.density_J && jastrow_type == "e-den-den"
        num_jpar_opts = num_jpars - 1
    elseif optimize.spin_J && jastrow_type == "e-spn-spn"
        num_jpar_opts = num_jpars - 1
    else
        num_jpar_opts = 0
    end

    @debug """
    Parameters::JastrowParameters() : 
    Type: $(jastrow_type)
    Total number of Jastrow parameters = $(num_jpars)
    Number of Jastrow parameters to be optimized = $(num_jpar_opts)
    """

    return JastrowParameters(jastrow_type, jpar_map, num_jpars, num_jpar_opts)
end


@doc raw"""

    JastrowParameters( jastrow_type::S, 
                       optimize::NamedTuple,
                       model_geometry::ModelGeometry,
                       path_to_parameter_file::S ) where {S<:AbstractString}

Given a type of Jastrow factor, an initial set of Jastrow parameters from file, and optimization 
flags, generates an initial set of Jastrow variational.

- `jastrow_type::S`: type of Jastrow factor: "e-den-den", "e-spn-spn".
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `path_to_parameter_file::S`: filepath to initial parameters.

"""
function JastrowParameters(
    jastrow_type::S, 
    optimize::NamedTuple,
    model_geometry::ModelGeometry,
    path_to_parameter_file::S
) where {S<:AbstractString}
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

    if optimize.density_J && jastrow_type == "e-den-den"
        num_jpar_opts = num_jpars - 1
    elseif optimize.spin_J && jastrow_type == "e-spn-spn"
        num_jpar_opts = num_jpars - 1
    else
        num_jpar_opts = 0
    end

    @debug """
    Parameters::JastrowParameters() : 
    Type: $(jastrow_type)
    Total number of Jastrow parameters = $(num_jpars)
    Number of Jastrow parameters to be optimized = $(num_jpar_opts)
    """

    return JastrowParameters(jastrow_type, jpar_map, num_jpars, num_jpar_opts)
end


@doc raw"""

    collect_parameters( determinantal_parameters::DeterminantalParameters{I}, 
                        jastrow_parameters::JastrowParameters{S, K, V, I} ) where {S<:AbstractString, K, V, I<:Integer}

Concatenates all values of determinantal and Jastrow parameters into a single vector of
variational parameters.

- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow parameters.

"""
function collect_parameters(
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters::JastrowParameters{S, K, V, I}
) where {S<:AbstractString, K, V, I<:Integer}
    # determinantal parameters
    det_pars = collect(values(determinantal_parameters.det_pars))
    det_pars = reduce(vcat, (isa(x, AbstractVector) ? x : [x] for x in values(det_pars)))

    # Jastrow parameters
    keys_sorted = sort(collect(keys(jastrow_parameters.jpar_map)))
    jpars = [jastrow_parameters.jpar_map[k][2] for k in keys_sorted]

    return vcat(det_pars, jpars)
end


@doc raw"""

    collect_parameters( determinantal_parameters::DeterminantalParameters{I}, 
                        jastrow_parameters_1::JastrowParameters{S, K, V, I},
                        jastrow_parameters_2::JastrowParameters{S, K, V, I} ) where {S<:AbstractString, K, V, I<:Integer}

Concatenates all values of determinantal and Jastrow parameters into a single vector of 
variational parameters.

- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal parameters.
- `jastrow_parameters_1::JastrowParameters{S, I, V}`: first set of Jastrow parameters.
- `jastrow_parameters_2::JastrowParameters{S, I, V}`: second set of Jastrow parameters.

"""
function collect_parameters(
    determinantal_parameters::DeterminantalParameters{I}, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I}
) where {S<:AbstractString, K, V, I<:Integer}
    # determinantal parameters
    det_pars = collect(values(determinantal_parameters.det_pars))

    # Jastrow parameters
    keys_sorted_1 = sort(collect(keys(jastrow_parameters_1.jpar_map)))
    jpars_1 = [jastrow_parameters_1.jpar_map[k][2] for k in keys_sorted_1[1:end-1]]

    keys_sorted_2 = sort(collect(keys(jastrow_parameters_2.jpar_map)))
    jpars_2 = [jastrow_parameters_2.jpar_map[k][2] for k in keys_sorted_2[1:end-1]]

    return vcat(det_pars, jpars_1, jpars_2)
end



@doc raw"""

    update_parameters!( measurement_container::NamedTuple,
                        new_vpars::AbstractVector, 
                        determinantal_parameters::DeterminantalParameters{I} ) where {I<:Integer}

Updates variational parameters after Stochastic Reconfiguration.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `new_vpars::AbstractVector`: updated variational parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.

"""
function update_parameters!(
    measurement_container::NamedTuple,
    new_vpars::AbstractVector,
    determinantal_parameters::DeterminantalParameters{I}
) where {I<:Integer}
    # extract current parameters (a NamedTuple)
    current_det_pars = determinantal_parameters.det_pars

    # tuple of parameter names (Symbols) in the same order as current_pars
    param_names = Tuple(keys(current_det_pars))

    # Prepare to collect reconstructed parameter *values*
    reconstructed_vals = Vector{Any}(undef, length(param_names))
    pos = 1

    for (i, pname) in enumerate(param_names)
        old_val = current_det_pars[pname]

        if old_val isa AbstractVector
            n = length(old_val)
            # bounds-check
            if pos + n - 1 > length(new_vpars)
                throw(ArgumentError("new_vpars too short for parameter $(pname): need $n elements starting at $pos"))
            end
            slice = new_vpars[pos:pos + n - 1]

            # create a vector of the same shape/type as the old one when reasonable.
            # If old_val is a plain Vector, collect(slice) is fine.
            # For other custom vector types this may still produce a plain Vector; adapt if needed.
            reconstructed_vals[i] = collect(slice)

            pos += n
        else
            # scalar parameter
            if pos > length(new_vpars)
                throw(ArgumentError("new_vpars too short for parameter $(pname) (expecting a scalar at position $pos)"))
            end
            reconstructed_vals[i] = new_vpars[pos]
            pos += 1
        end
    end

    # ensure we've consumed exactly all entries of new_vpars
    if pos - 1 != length(new_vpars)
        throw(ArgumentError("Length mismatch: consumed $(pos - 1) elements but new_vpars has length $(length(new_vpars))."))
    end

    # build NamedTuple with the same field names in the same order
    recon_det_pars = NamedTuple{param_names}(Tuple(reconstructed_vals))

    # update determinantal_parameters and measurement container
    determinantal_parameters.det_pars = recon_det_pars
    measurement_container.optimization_measurements["parameters"] = new_vpars

    return nothing
end



@doc raw"""

    update_parameters!( measurement_container::NamedTuple,
                        new_vpars::AbstractVector, 
                        determinantal_parameters::DeterminantalParameters{I}, 
                        jastrow_parameters::JastrowParameters{S, K, V, I} ) where {S<:AbstractString, K, V, I<:Integer}

Updates variational parameters after Stochastic Reconfiguration.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `new_vpars::AbstractVector`: updated variational parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: set of Jastrow variational parameters.

"""
function update_parameters!(
    measurement_container::NamedTuple,
    new_vpars::AbstractVector, 
    determinantal_parameters::DeterminantalParameters{I},
    jastrow_parameters::JastrowParameters{S, K, V, I}
) where {S<:AbstractString, K, V, I<:Integer}
    # extract current parameters from their containers
    current_det_pars = determinantal_parameters.det_pars
    current_jpar_map = jastrow_parameters.jpar_map

    # seperate the new determinantal and Jastrow parameters
    num_jpars = jastrow_parameters.num_jpars

    new_det_pars = new_vpars[1:end-num_jpars]
    new_jpars = new_vpars[end-num_jpars+1:end]
    
    # update determinantal parameters
    # tuple of parameter names (Symbols) in the same order as current_pars
    param_names = Tuple(keys(current_det_pars))

    # Prepare to collect reconstructed parameter *values*
    reconstructed_vals = Vector{Any}(undef, length(param_names))
    pos = 1

    for (i, pname) in enumerate(param_names)
        old_val = current_det_pars[pname]

        if old_val isa AbstractVector
            n = length(old_val)
            # bounds-check
            if pos + n - 1 > length(new_det_pars)
                throw(ArgumentError("new_vpars too short for parameter $(pname): need $n elements starting at $pos"))
            end
            slice = new_det_pars[pos:pos + n - 1]

            # create a vector of the same shape/type as the old one when reasonable.
            # If old_val is a plain Vector, collect(slice) is fine.
            # For other custom vector types this may still produce a plain Vector; adapt if needed.
            reconstructed_vals[i] = collect(slice)

            pos += n
        else
            # scalar parameter
            if pos > length(new_det_pars)
                throw(ArgumentError("new_det_pars too short for parameter $(pname) (expecting a scalar at position $pos)"))
            end
            reconstructed_vals[i] = new_det_pars[pos]
            pos += 1
        end
    end

    # ensure we've consumed exactly all entries of new_vpars
    if pos - 1 != length(new_det_pars)
        throw(ArgumentError("Length mismatch: consumed $(pos - 1) elements but new_det_pars has length $(length(new_det_pars))."))
    end

    # build NamedTuple with the same field names in the same order
    recon_det_pars = NamedTuple{param_names}(Tuple(reconstructed_vals))

    # update determinantal_parameters
    determinantal_parameters.det_pars = recon_det_pars

    # update Jastrow parameters
    irr_indices = collect(keys(current_jpar_map))
    for i in 1:num_jpars
        indices, _ = current_jpar_map[irr_indices[i]]
    
        current_jpar_map[irr_indices[i]] = (indices, new_jpars[i])
    end

    # update the measurement container
    measurement_container.optimization_measurements["parameters"] = new_vpars

    return nothing
end


@doc raw"""

    update_parameters!( measurement_container::NamedTuple,
                        new_vpars::AbstractVector, 
                        determinantal_parameters::DeterminantalParameters{I}, 
                        jastrow_parameters_1::JastrowParameters{S, K, V, I},
                        jastrow_parameters_2::JastrowParameters{S, K, V, I} ) where {S<:AbstractString, K, V, I<:Integer}

Updates variational parameters after Stochastic Reconfiguration.

- `measurement_container::NamedTuple`: container where measurements are stored.
- `new_vpars::AbstractVector`: updated variational parameters.
- `determinantal_parameters::DeterminantalParameters{I}`: current set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.

"""
# TODO: need to finish Jastrow updating method
function update_parameters!(
    measurement_container::NamedTuple,
    new_vpars::AbstractVector, 
    determinantal_parameters::DeterminantalParameters{I},
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I}
) where {S<:AbstractString, K, V, I<:Integer}
    # extract current parameters from their containers
    current_det_pars = determinantal_parameters.det_pars
    current_jpar_map_1 = jastrow_parameters_1.jpar_map
    current_jpar_map_2 = jastrow_parameters_2.jpar_map

    # seperate the new determinantal and Jastrow parameters
    num_jpars_1 = jastrow_parameters_1.num_jpar_opts
    num_jpars_2 = jastrow_parameters_2.num_jpar_opts
    num_jpars = num_jpars_1 + num_jpars_2

    new_det_pars = new_vpars[1:end-num_jpars]
    new_jpars = new_vpars[end-num_jpars+1:end]

    @assert length(new_vpars) != length(current_det_pars) + num_jpars
    
    # update determinantal parameters
    param_names = keys(current_det_pars)
    new_det_pars = NamedTuple{Tuple(param_names)}(Tuple(new_vpars))
    determinantal_parameters.det_pars = new_det_pars

    # update Jastrow parameters
    irr_indices_1 = collect(keys(current_jpar_map_1))
    for i in 1:num_jpars_1
        indices, _ = current_jpar_map_1[irr_indices[i]]
    
        current_jpar_map_1[irr_indices[i]] = (indices, new_jpars[i])
    end

    irr_indices_2 = collect(keys(current_jpar_map_2))
    for i in 1:num_jpars_2
        indices, _ = current_jpar_map_2[irr_indices[i]]
    
        current_jpar_map_2[irr_indices[i]] = (indices, new_jpars[i])
    end

    # update the measurement container
    measurement_container.optimization_measurements["parameters"] = new_vpars

    return nothing
end


@doc raw"""

    readin_parameters( filename::S ) where {S<:AbstractString}

Parses TOML file containing initial variational parameters. 

- `filename::S`: name of parameter summary file in TOML format.

"""
function readin_parameters(
    filename::S
) where {S<:AbstractString}
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
