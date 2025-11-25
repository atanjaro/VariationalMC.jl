@doc raw"""

    process_measurements( measurement_container::NamedTuple,
                          simulation_info::SimulationInfo,
                          determinantal_parameters::DeterminantalParameters )

Processes all simulation and optimization measurements by organinzing them into CSV files.

- `measurement_container::NamedTuple`: contains measurement quantities.
- `simulation_info::SimulationInfo`: contains all simulation info.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.

"""
function process_measurements(
    measurement_container, 
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters
)
    (; datafolder, pID) = simulation_info
    (; N_opt, N_sim) = measurement_container

    # merge `opt` and `sim` bins
    merge_bin_measurements!(datafolder * "/simulation/")

    # process all scalar measurements
    process_scalar_measurements(
        datafolder, 
        "local_energy"
    )
    process_scalar_measurements(
        datafolder, 
        "double_occ"
    )
    process_scalar_measurements(
        datafolder, 
        "global_density"
    )
    process_scalar_measurements(
        datafolder, 
        "pconfig"
    )

   # process optimization measurements
    process_optimization_measurements(
        datafolder, 
        determinantal_parameters
    )

   return nothing
end


@doc raw"""

    process_measurements( measurement_container::NamedTuple,
                          simulation_info::SimulationInfo,
                          determinantal_parameters::DeterminantalParameters,
                          jastrow_parameters::JastrowParameters )

Processes all simulation and optimization measurements by organinzing them into CSV files.

- `measurement_container::NamedTuple`: contains measurement quantities.
- `simulation_info::SimulationInfo`: contains all simulation info.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: set of Jastrow parameters.

"""
function process_measurements(
    measurement_container, 
    simulation_info::SimulationInfo,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters::JastrowParameters
)
    (; datafolder, pID) = simulation_info
    (; N_opt, N_sim) = measurement_container

    # merge `opt` and `sim` bins
    merge_bin_measurements!(datafolder * "/simulation/")

    # process all scalar measurements
    process_scalar_measurements(
        datafolder, 
        "local_energy"
    )
    process_scalar_measurements(
        datafolder, 
        "double_occ"
    )
    process_scalar_measurements(
        datafolder, 
        "global_density"
    )
    process_scalar_measurements(
        datafolder, 
        "pconfig"
    )

   # process optimization measurements
    process_optimization_measurements(
        datafolder, 
        determinantal_parameters,
        jastrow_parameters
    )

   return nothing
end


@doc raw"""

    process_scalar_measurements( datafolder::T,
                                 measurement::T;
                                 N_bins::Union{I, Nothing}=nothing) where {T<:AbstractString, I<:Integer}

Write binned simulation measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `measurement::T`: `local_energy`, `double_occ`, `global_density`, or `pconfig`.
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.

"""
function process_scalar_measurements(
    datafolder::T,
    measurement::T,
    N_bins::Union{I, Nothing}=nothing
) where {T<:AbstractString, I<:Integer}
    sim_file = joinpath(datafolder, "simulation" ,"bin_measurements.h5")
    @assert isfile(sim_file) "HDF5 file not found: $sim_file"
    # sim_file = joinpath(datafolder, "simulation", "simulation_measurements.h5")
    
    allowed = Set(["local_energy", "double_occ", "global_density", "pconfig"])
    results = Vector{Any}()

    h5open(sim_file, "r") do f
        # helper to form the dataset path
        path_for_bin(bin_num) = "/$measurement/bin-$bin_num"

        if N_bins === nothing
            # probe bins until a bin is missing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    # Post-process results into a DataFrame.
    bins = collect(1:length(results))

    if measurement == "local_energy"
        # Separate real and imaginary components
        mean_r = [real(x) for x in results]
        mean_i = [imag(x) for x in results]
        df = DataFrame(BIN = bins, MEAN_R = mean_r, MEAN_I = mean_i)
    else
        mean_r = [x for x in results]
        df = DataFrame(BIN = bins, MEAN_R = mean_r)
    end

    output_csv = joinpath(datafolder, "simulation", measurement*"_stats.csv")

    CSV.write(output_csv, df)

    return nothing
end


@doc raw"""

    process_optimization_measurements( datafolder::T,
                                       determinantal_parameters::DeterminantalParameters;
                                       N_bins::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

Writes binned optimization measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.

"""
function process_optimization_measurements(
    datafolder::T,
    determinantal_parameters::DeterminantalParameters;
    N_bins::Union{I, Nothing}=nothing
    
) where {T<:AbstractString, I<:Integer}
    opt_file = joinpath(datafolder, "optimization", "opt_bin_measurements.h5")
    @assert isfile(opt_file) "HDF5 file not found: $opt_file"

    # helper function to convert filenames to ASCII
    convert_name(name::AbstractString) = replace(name, "Δ" => "delta", "μ" => "mu", " " => "_")

    # extract parameter information 
    current_det_pars = determinantal_parameters.det_pars
    param_names = Tuple(keys(current_det_pars))  

    # helper function to reconstruct a parameter NamedTuple
    function reconstruct_from_vpars(vpars::AbstractVector)
        reconstructed = Dict{Symbol, Any}()
        pos = 1
        for pname in param_names
            old_val = current_det_pars[pname]
            if old_val isa AbstractVector
                n = length(old_val)
                if pos + n - 1 > length(vpars)
                    throw(ArgumentError("vpars too short for parameter $(pname): need $n elements starting at $pos"))
                end
                slice = vpars[pos:pos + n - 1]
                reconstructed[pname] = collect(slice) 
                pos += n
            else
                if pos > length(vpars)
                    throw(ArgumentError("vpars too short for scalar parameter $(pname) at position $pos"))
                end
                reconstructed[pname] = vpars[pos]
                pos += 1
            end
        end
        if pos - 1 != length(vpars)
            throw(ArgumentError("Length mismatch: consumed $(pos - 1) elements but vpars has length $(length(vpars))."))
        end
        return reconstructed
    end

    results = Vector{Any}()

    h5open(opt_file, "r") do f
        path_for_bin(bin_num) = "/parameters/bin-$bin_num"

        if N_bins === nothing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    bins = collect(1:length(results))

    # reconstruct parameter dictionary
    per_bin_params = Vector{Union{Dict{Symbol,Any}, Missing}}(undef, length(results))
    for (i, entry) in enumerate(results)
        if entry === missing
            per_bin_params[i] = missing
        else
            vpars = collect(entry)
            per_bin_params[i] = reconstruct_from_vpars(vpars)
        end
    end

    outdir = joinpath(datafolder, "optimization")
    isdir(outdir) || mkpath(outdir)

    for pname in param_names
        converted = convert_name(String(pname))
        outcsv = joinpath(outdir, "parameter_$(converted)_stats.csv")

        template_val = current_det_pars[pname]
        if template_val isa AbstractVector
            # for vector parameters
            n = length(template_val)
            el_is_complex = Base.eltype(template_val) <: Complex || any(x->x isa Complex, template_val)

            df_cols = Dict{Symbol, Any}()
            df_cols[:BIN] = bins
            for j in 1:n
                df_cols[Symbol("MEAN_R_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
            end
            if el_is_complex
                for j in 1:n
                    df_cols[Symbol("MEAN_I_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
                end
            end

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    for j in 1:n
                        df_cols[Symbol("MEAN_R_$(j)")][i] = missing
                    end
                    if el_is_complex
                        for j in 1:n
                            df_cols[Symbol("MEAN_I_$(j)")][i] = missing
                        end
                    end
                else
                    val = entry[pname]
                    if !(val isa AbstractVector)
                        throw(ArgumentError("Expected vector for parameter $(pname) but got scalar in bin $i"))
                    end
                    for j in 1:n
                        comp = val[j]
                        if comp isa Complex
                            df_cols[Symbol("MEAN_R_$(j)")][i] = real(comp)
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = imag(comp)
                            end
                        else
                            df_cols[Symbol("MEAN_R_$(j)")][i] = comp
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = 0.0
                            end
                        end
                    end
                end
            end

            df = DataFrame(df_cols)
            CSV.write(outcsv, df)

        else
            # for scalar parameters
            el_is_complex = (template_val isa Complex) ||
                            any(b -> (b !== missing && (b[pname] isa Complex)), per_bin_params)

            MEAN_R = Vector{Union{Float64, Missing}}(undef, length(bins))
            MEAN_I = el_is_complex ? Vector{Union{Float64, Missing}}(undef, length(bins)) : nothing

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    MEAN_R[i] = missing
                    if el_is_complex
                        MEAN_I[i] = missing
                    end
                else
                    val = entry[pname]
                    if val isa Complex
                        MEAN_R[i] = real(val)
                        if el_is_complex
                            MEAN_I[i] = imag(val)
                        end
                    else
                        MEAN_R[i] = val
                        if el_is_complex
                            MEAN_I[i] = 0.0
                        end
                    end
                end
            end

            if el_is_complex
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R, MEAN_I = MEAN_I)
            else
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R)
            end
            CSV.write(outcsv, df)
        end
    end

    return nothing
end


@doc raw"""

    process_optimization_measurements( datafolder::T,
                                       determinantal_parameters::DeterminantalParameters,
                                       jastrow_parameters::JastrowParameters;
                                       N_bins::Union{I, Nothing}=nothing ) where {T<:AbstractString, I<:Integer}

Writes binned optimization measurements to CSV.

- `datafolder::T`: path to folder where simulation files are written.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `jastrow_parameters::JastrowParameters`: set of Jastrow parameters.
- `N_bins::Union{I, Nothing}=nothing`: (optional) total number of bins.

"""
function process_optimization_measurements(
    datafolder::T,
    determinantal_parameters::DeterminantalParameters,
    jastrow_parameters::JastrowParameters;
    N_bins::Union{I, Nothing}=nothing
) where {T<:AbstractString, I<:Integer}
    opt_file = joinpath(datafolder, "optimization", "opt_bin_measurements.h5")
    @assert isfile(opt_file) "HDF5 file not found: $opt_file"

    # helper function to convert filenames to ASCII
    convert_name(name::AbstractString) = replace(name, "Δ" => "delta", "μ" => "mu", " " => "_")

    # extract parameter information 
    current_det_pars = determinantal_parameters.det_pars
    param_names = Tuple(keys(current_det_pars))

    # helper function to reconstruct a parameter NamedTuple and Jastrow parameters
    function reconstruct_from_vpars(vpars::AbstractVector)
        reconstructed = Dict{Symbol, Any}()
        pos = 1
        for pname in param_names
            old_val = current_det_pars[pname]
            if old_val isa AbstractVector
                n = length(old_val)
                if pos + n - 1 > length(vpars)
                    throw(ArgumentError("vpars too short for parameter $(pname): need $n elements starting at $pos"))
                end
                slice = vpars[pos:pos + n - 1]
                reconstructed[pname] = collect(slice)
                pos += n
            else
                if pos > length(vpars)
                    throw(ArgumentError("vpars too short for scalar parameter $(pname) at position $pos"))
                end
                reconstructed[pname] = vpars[pos]
                pos += 1
            end
        end
        # If there are extra elements after the determinantal parameters, treat them as jastrow tail
        if pos <= length(vpars)
            jtail = collect(vpars[pos:end])
        else
            jtail = missing
        end
        return reconstructed, jtail, (pos - 1)
    end

    results = Vector{Any}()

    h5open(opt_file, "r") do f
        path_for_bin(bin_num) = "/parameters/bin-$bin_num"

        if N_bins === nothing
            bin_num = 1
            while true
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                    bin_num += 1
                else
                    break
                end
            end
        else
            for bin_num in 1:N_bins
                path = path_for_bin(bin_num)
                if haskey(f, path)
                    push!(results, read(f[path]))
                else
                    push!(results, missing)
                end
            end
        end
    end

    bins = collect(1:length(results))

    # reconstruct parameter dictionary and collect Jastrow parameters
    per_bin_params = Vector{Union{Dict{Symbol,Any}, Missing}}(undef, length(results))
    per_bin_jpars  = Vector{Union{Vector{Float64}, Vector{Any}, Missing}}(undef, length(results))
    consumed_lengths = Vector{Int}(undef, length(results))

    for (i, entry) in enumerate(results)
        if entry === missing
            per_bin_params[i] = missing
            per_bin_jpars[i]  = missing
            consumed_lengths[i] = 0
        else
            vpars = collect(entry)
            det_dict, jtail, consumed = reconstruct_from_vpars(vpars)
            per_bin_params[i] = det_dict
            per_bin_jpars[i] = jtail
            consumed_lengths[i] = consumed
        end
    end

    # get number of Jastrow parameters
    n_jpars = nothing
    try
        if hasproperty(jastrow_parameters, :num_jpar_opts)
            n_jpars = Int(getproperty(jastrow_parameters, :num_jpar_opts))
        end
    catch
        n_jpars = nothing
    end

    if n_jpars === nothing
        try
            if hasproperty(jastrow_parameters, :jpars)
                n_jpars = length(getproperty(jastrow_parameters, :jpars))
            elseif hasproperty(jastrow_parameters, :vpars)
                n_jpars = length(getproperty(jastrow_parameters, :vpars))
            elseif hasproperty(jastrow_parameters, :parameters)
                n_jpars = length(getproperty(jastrow_parameters, :parameters))
            end
        catch
            n_jpars = nothing
        end
    end

    if n_jpars === nothing
        for jtail in per_bin_jpars
            if jtail !== missing
                n_jpars = length(jtail)
                break
            end
        end
    end
    n_jpars = n_jpars === nothing ? 0 : n_jpars

    outdir = joinpath(datafolder, "optimization")
    isdir(outdir) || mkpath(outdir)

    # write determinantal parameters
    for pname in param_names
        converted = convert_name(String(pname))
        outcsv = joinpath(outdir, "parameter_$(converted)_stats.csv")

        template_val = current_det_pars[pname]
        if template_val isa AbstractVector
            # for vector parameters
            n = length(template_val)
            el_is_complex = Base.eltype(template_val) <: Complex || any(x->x isa Complex, template_val)

            df_cols = Dict{Symbol, Any}()
            df_cols[:BIN] = bins
            for j in 1:n
                df_cols[Symbol("MEAN_R_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
            end
            if el_is_complex
                for j in 1:n
                    df_cols[Symbol("MEAN_I_$(j)")] = Vector{Union{Float64, Missing}}(undef, length(bins))
                end
            end

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    for j in 1:n
                        df_cols[Symbol("MEAN_R_$(j)")][i] = missing
                    end
                    if el_is_complex
                        for j in 1:n
                            df_cols[Symbol("MEAN_I_$(j)")][i] = missing
                        end
                    end
                else
                    val = entry[pname]
                    if !(val isa AbstractVector)
                        throw(ArgumentError("Expected vector for parameter $(pname) but got scalar in bin $i"))
                    end
                    for j in 1:n
                        comp = val[j]
                        if comp isa Complex
                            df_cols[Symbol("MEAN_R_$(j)")][i] = real(comp)
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = imag(comp)
                            end
                        else
                            df_cols[Symbol("MEAN_R_$(j)")][i] = comp
                            if el_is_complex
                                df_cols[Symbol("MEAN_I_$(j)")][i] = 0.0
                            end
                        end
                    end
                end
            end

            df = DataFrame(df_cols)
            CSV.write(outcsv, df)

        else
            # for scalar parameters
            el_is_complex = (template_val isa Complex) ||
                            any(b -> (b !== missing && (b[pname] isa Complex)), per_bin_params)

            MEAN_R = Vector{Union{Float64, Missing}}(undef, length(bins))
            MEAN_I = el_is_complex ? Vector{Union{Float64, Missing}}(undef, length(bins)) : nothing

            for (i, entry) in enumerate(per_bin_params)
                if entry === missing
                    MEAN_R[i] = missing
                    if el_is_complex
                        MEAN_I[i] = missing
                    end
                else
                    val = entry[pname]
                    if val isa Complex
                        MEAN_R[i] = real(val)
                        if el_is_complex
                            MEAN_I[i] = imag(val)
                        end
                    else
                        MEAN_R[i] = val
                        if el_is_complex
                            MEAN_I[i] = 0.0
                        end
                    end
                end
            end

            if el_is_complex
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R, MEAN_I = MEAN_I)
            else
                df = DataFrame(BIN = bins, MEAN_R = MEAN_R)
            end
            CSV.write(outcsv, df)
        end
    end

    # write Jastrow parameters
    jtype = try
        string(getproperty(jastrow_parameters, :jastrow_type))
    catch
        "jastrow"
    end
    jfilename = joinpath(outdir, "$(jtype)_jastrow_parameters.csv")

    if n_jpars > 0
        jcolnames = [Symbol("MEAN_V_$(k)") for k in 1:n_jpars]
        jcols = Dict{Symbol, Any}()
        jcols[:BIN] = bins
        for sym in jcolnames
            jcols[sym] = Vector{Union{Float64, Missing}}(undef, length(bins))
        end

        for i in 1:length(bins)
            jtail = per_bin_jpars[i]
            if jtail === missing
                for sym in jcolnames
                    jcols[sym][i] = missing
                end
            else
                for k in 1:n_jpars
                    if k <= length(jtail)
                        val = jtail[k]
                        if val isa Complex
                            jcols[jcolnames[k]][i] = real(val)
                        else
                            jcols[jcolnames[k]][i] = Float64(val)
                        end
                    else
                        jcols[jcolnames[k]][i] = missing
                    end
                end
            end
        end

        jdf = DataFrame(jcols)
        CSV.write(jfilename, jdf)
    else
        CSV.write(jfilename, DataFrame(BIN = bins))
    end

    return nothing
end




@doc raw"""

    process_correlation_measurements()

For either density-density or spin-spin correlation data, calculates the static structure factor.

"""
function process_correlation_measurements(correlation_type::String, model_geometry::ModelGeometry)
    # calculate all momentum points that will used in the FT
    q_points = calc_k_points(model_geometry.unit_cell, model_geometry.lattice)

    if correlation_type == "density"
        # read in the site-dependent density data 

    elseif correlation_type == "spin"
        # read in the site-dependent spin data

    end

    # calculate_structure_factor()

    # write structure factor info to file
end


@doc raw"""

    calculate_structure_factor()

Calculates the static structure factor from either density-density or spin-spin correlations.

"""
function calculate_structure_factor(correlation_type::String, model_geometry::ModelGeometry)
    # datafolder PATH
    # correlation/density OR correlation/spin-z

    # check that there are indeed the required correlation measurements

    # collect all q-points for the lattice
    uc = model_geometry.unit_cell
    lat = model_geometry.lattice
    q_points = calc_k_points(uc, lat)

    # read-in the correlation measurements

    # store them in a matrix

    # open csv file where final results will be written

    # calculate onsite part of the sum
    N = lat.N
    onsite = diag(nn_corr)/N       # OR onsite = diag(ss_corr)/lat.N for spin-spin

    for q in q_points
        sum = 0.0
        for l in 1:N
            for k in (l+1):N
                # if include_r(l, k)      # this is a consideration for if we have a bilayer model
                # end
                r_diff = r(1,k) .- r(1,l)
                sum += cos(dot(q, r_diff)) * nn_corr[l, k]
            end
        end
        Nq = 2.0 / N * sum + onsite     # or Sq
    end
end


"""

    merge_bin_measurements!( datafolder::AbstractString;
                             opt_name::AbstractString="opt_bin_measurements.h5",
                             sim_name::AbstractString="sim_bin_measurements.h5",
                             out_name::AbstractString="bin_measurements.h5",
                             clear_src::Bool=false )

Merge the HDF5 files `opt_name` and `sim_name` (located inside `datafolder`)
into a single HDF5 file named `out_name` (also inside `datafolder`).

Bins from the sim file are offset by the maximum bin index found in the opt file,
so the result has contiguous, non-overlapping `bin-%d` subgroup names.

- `datafolder::AbstractString`: path to data files.
- `opt_name::AbstractString="opt_bin_measurements.h5"`: name of opt_bin file.
- `sim_name::AbstractString="sim_bin_measurements.h5"`: name of sim_bin file.
- `out_name::AbstractString="bin_measurements.h5"`: output file name.
- `clear_src::Bool=false`: optionallyy clear source bin files after merger.

"""
function merge_bin_measurements!(
    datafolder::AbstractString;
    opt_name::AbstractString="opt_bin_measurements.h5",
    sim_name::AbstractString="sim_bin_measurements.h5",
    out_name::AbstractString="bin_measurements.h5",
    clear_src::Bool=false
)
    opt_file = joinpath(datafolder, opt_name)
    sim_file = joinpath(datafolder, sim_name)
    out_file = joinpath(datafolder, out_name)

    bin_regex = r"^bin-(\d+)$"
    parse_bin(s::AbstractString) = begin
        m = match(bin_regex, s)
        m === nothing ? nothing : parse(Int, m.captures[1])
    end

    function top_level_groups(fname)
        if !isfile(fname)
            @warn "File not found: $fname"
            return String[]
        end
        groups = String[]
        h5open(fname, "r") do f
            append!(groups, collect(keys(f)))   # use keys() instead of names()
        end
        return groups
    end

    function max_bin_in_file(fname)
        if !isfile(fname)
            return 0
        end
        maxbin = 0
        h5open(fname, "r") do f
            for g in collect(keys(f))
                # g is a top-level name
                gpath = "/" * g
                # attempt to iterate children using keys on the group object
                try
                    grp = f[gpath]
                    for n in collect(keys(grp))
                        b = parse_bin(n)
                        if b !== nothing && b > maxbin
                            maxbin = b
                        end
                    end
                catch
                    # ignore non-group entries
                end
            end
        end
        return maxbin
    end

    # Copy attributes from src object to dst object (datasets only)
    function copy_attrs_obj!(srcobj, dstobj)
        for (k, v) in pairs(attrs(srcobj))
            attrs(dstobj)[k] = v
        end
    end

    # read everything from a given group/bin in src and write into dst with optionally renumbered bin index
    function copy_group_bin!(srcfname, dstf, groupname::String, binname::String, new_bin_index::Int)
        src_gpath = "/" * groupname * "/" * binname
        dst_gpath = "/" * groupname
        dst_binname = @sprintf("bin-%d", new_bin_index)
        dst_full = dst_gpath * "/" * dst_binname

        h5open(srcfname, "r") do sf
            if !haskey(sf, src_gpath)
                return
            end

            src_node = sf[src_gpath]

            # try to list children using keys()
            children = try
                collect(keys(src_node))
            catch
                String[]
            end

            if isempty(children)
                data = read(src_node)
                if haskey(dstf, dst_full)
                    delete!(dstf, dst_full)
                end
                # writing to dst_full will create intermediate groups automatically
                dstf[dst_full] = data
                try
                    copy_attrs_obj!(src_node, dstf[dst_full])
                catch
                end
            else
                # src bin is a group with children; copy each child dataset
                for child in children
                    child_src = src_gpath * "/" * child
                    child_dst = dst_full * "/" * child
                    src_child_node = sf[child_src]
                    grand = try
                        collect(keys(src_child_node))
                    catch
                        String[]
                    end

                    if isempty(grand)
                        data = read(src_child_node)
                        if haskey(dstf, child_dst)
                            delete!(dstf, child_dst)
                        end
                        dstf[child_dst] = data
                        try
                            copy_attrs_obj!(src_child_node, dstf[child_dst])
                        catch
                        end
                    else
                        # nested group - copy leaf datasets one level deep
                        for gchild in grand
                            gchild_src = child_src * "/" * gchild
                            gchild_dst = child_dst * "/" * gchild
                            data = read(sf[gchild_src])
                            if haskey(dstf, gchild_dst)
                                delete!(dstf, gchild_dst)
                            end
                            dstf[gchild_dst] = data
                            try
                                copy_attrs_obj!(sf[gchild_src], dstf[gchild_dst])
                            catch
                            end
                        end
                    end
                end
            end
        end
    end

    max_opt = max_bin_in_file(opt_file)
    # @info "max bin index found in opt file: $max_opt"

    out_mode = isfile(out_file) ? "r+" : "w"
    h5open(out_file, out_mode) do outf
        # copy file attrs from opt if present else sim (best-effort)
        if isfile(opt_file)
            h5open(opt_file, "r") do optf
                for (k,v) in pairs(attrs(optf))
                    if !haskey(attrs(outf), k)
                        attrs(outf)[k] = v
                    end
                end
            end
        elseif isfile(sim_file)
            h5open(sim_file, "r") do simf
                for (k,v) in pairs(attrs(simf))
                    if !haskey(attrs(outf), k)
                        attrs(outf)[k] = v
                    end
                end
            end
        end

        # union of top-level groups
        top_groups = union(top_level_groups(opt_file), top_level_groups(sim_file))

        # copy opt file content (original bin numbers)
        if isfile(opt_file)
            h5open(opt_file, "r") do optf
                for g in collect(keys(optf))
                    gpath = "/" * g
                    # iterate children of the group one-level deep
                    grp = try
                        optf[gpath]
                    catch
                        nothing
                    end
                    if grp === nothing
                        continue
                    end
                    for binname in collect(keys(grp))
                        b = parse_bin(binname)
                        if b === nothing
                            # non-bin entry: copy as dataset (overwrite if present)
                            src_path = gpath * "/" * binname
                            dst_path = src_path
                            data = read(optf[src_path])
                            if haskey(outf, dst_path)
                                delete!(outf, dst_path)
                            end
                            outf[dst_path] = data
                            try
                                copy_attrs_obj!(optf[src_path], outf[dst_path])
                            catch
                            end
                        else
                            copy_group_bin!(opt_file, outf, g, binname, b)
                        end
                    end
                end
            end
        end

        # copy sim file content, offsetting bin indices
        if isfile(sim_file)
            h5open(sim_file, "r") do simf
                for g in collect(keys(simf))
                    gpath = "/" * g
                    grp = try
                        simf[gpath]
                    catch
                        nothing
                    end
                    if grp === nothing
                        continue
                    end
                    for binname in collect(keys(grp))
                        b = parse_bin(binname)
                        if b === nothing
                            src_path = gpath * "/" * binname
                            dst_path = src_path
                            data = read(simf[src_path])
                            if haskey(outf, dst_path)
                                delete!(outf, dst_path)
                            end
                            outf[dst_path] = data
                            try
                                copy_attrs_obj!(simf[src_path], outf[dst_path])
                            catch
                            end
                        else
                            new_index = max_opt + b
                            copy_group_bin!(sim_file, outf, g, binname, new_index)
                        end
                    end
                end
            end
        end
    end
    # @info "Merge complete. Output written to $(abspath(out_file))"
        # optionally delete source files
    if clear_src
        if isfile(opt_file)
            rm(opt_file; force=true)
            # @info "Deleted source file: $opt_file"
        end
        if isfile(sim_file)
            rm(sim_file; force=true)
            # @info "Deleted source file: $sim_file"
        end
    end

    return nothing
end
