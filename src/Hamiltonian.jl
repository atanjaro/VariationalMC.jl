@doc raw"""

    build_auxiliary_hamiltonian( tight_binding_model::TightBindingModel{E}, 
                                 determinantal_parameters::DeterminantalParameters{I}, 
                                 optimize::NamedTuple, 
                                 model_geometry::ModelGeometry, 
                                 pht::Bool ) where {E<:AbstractFloat, I<:Integer}

Constructs an auxiliary Hamiltonian matrix ``H_{\mathrm{aux}}`` by combining the non-interacting
hopping matrix ``H_t`` with matrices of variational terms ``H_{\mathrm{var}}``.

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hol transformed.

"""
function build_auxiliary_hamiltonian(
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    pht::Bool
) where {E<:AbstractFloat, I<:Integer}
    # hopping matrix
    H_tb = build_tight_binding_hamiltonian(
        tight_binding_model, 
        model_geometry, 
        pht
    )

    # variational matrices and operators
    H_var, V = build_variational_hamiltonian(
        determinantal_parameters, 
        model_geometry,
        optimize, 
        pht
    )

    return H_tb + H_var, V
end


@doc raw"""

    build_auxiliary_hamiltonian( tight_binding_model::TightBindingModel{E}, 
                                 determinantal_parameters::DeterminantalParameters{I}, 
                                 optimize::NamedTuple, 
                                 model_geometry::ModelGeometry, 
                                 twist_angles::AbstractRange{E},
                                 pht::Bool ) where {E<:AbstractFloat, I<:Integer}

Constructs an auxiliary Hamiltonian matrix ``H_{\mathrm{aux}}^{\theta}`` for ``N_\theta`` twist angles
by combining the non-interacting hopping matrix ``H_t`` with matrices of variational terms ``H_{\mathrm{var}}``.

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `twist_angles::AbstractRange{E}`: set of twist angles.
- `pht::Bool`: whether model is particle-hol transformed.

"""
function build_auxiliary_hamiltonian(
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    twist_angles::AbstractRange{E},
    pht::Bool
) where {E<:AbstractFloat, I<:Integer}
    H_ts = []

    for θ_twist in twist_angles
        # untwisted hopping matrix
        H_tb = build_tight_binding_hamiltonian(
            tight_binding_model, 
            model_geometry, 
            pht
        )

        # apply twist phase
        apply_twist_angle!(H_tb, θ_twist, model_geometry)
        push!(H_ts, H_tb)
    end
    
    # variational matrices and operators
    H_var, V = build_variational_hamiltonian(
        determinantal_parameters, 
        model_geometry,
        optimize, 
        pht
    )

    for H_t in H_ts
        H_t .+= H_var
    end 

    return H_ts, V
end


@doc raw"""

    build_tight_binding_hamiltonian( tight_binding_model::TightBindingModel{E},
                                     model_geometry::ModelGeometry,
                                     pht::Bool ) where {E<:AbstractFloat}

Constructs a ``2N`` by ``2N`` tight-binding hopping matrix where ``N`` is the number of lattice sites, 
given hopping parameters ``t``, ``t^{\prime}``, and ``t^{\prime\prime}``.

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function build_tight_binding_hamiltonian(
    tight_binding_model::TightBindingModel{E}, 
    model_geometry::ModelGeometry, 
    pht::Bool
) where {E<:AbstractFloat}
    dims = length(model_geometry.lattice.L)
    N    = model_geometry.unit_cell.n * model_geometry.lattice.N 
    Lx   = model_geometry.lattice.L[1]
    if dims != 1
        Ly = model_geometry.lattice.L[2]
    else
        Ly = 1
    end

    # neighbor table
    bonds = model_geometry.bond
    nbr0 = build_neighbor_table(
        bonds[1], 
        model_geometry.unit_cell, 
        model_geometry.lattice
    )

    # remove double counting if any lattice dimension equals 2
    if Lx == 2 || Ly == 2
        keep = trues(size(nbr0, 2))
        seen = Set{Tuple{Int,Int}}()

        for (j, col) in enumerate(eachcol(nbr0))
            key = Tuple(sort(col))
            if key in seen
                keep[j] = false
            else
                push!(seen, key)
            end
        end

        nbr0 = nbr0[:, keep]
    end

    # allocate Hamiltonian once
    H_t = Matrix{Complex}(undef, 2N, 2N)
    fill!(H_t, 0)

    # hopping parameters
    t₀, t₁, t₂ = tight_binding_model.t₀, tight_binding_model.t₁, tight_binding_model.t₂

    if pht
        # nearest neighbors
        for col in axes(nbr0, 2)
            i, j = nbr0[1, col], nbr0[2, col]
            H_t[i, j] += -t₀
            H_t[j, i] += -t₀

            ip, jp = i + N, j + N
            H_t[ip, jp] += t₀
            H_t[jp, ip] += t₀
        end
        # next-nearest neighbors
        if t₁ != 0
            nbr1 = build_neighbor_table(
                bonds[2], 
                model_geometry.unit_cell, 
                model_geometry.lattice
            )
            for col in axes(nbr1, 2)
                i, j = nbr1[1, col], nbr1[2, col]
                H_t[i, j]   += t₁
                H_t[j, i]   += t₁
                ip, jp = i + N, j + N
                H_t[ip, jp] += -t₁
                H_t[jp, ip] += -t₁
            end
        end
        # third nearest neighbors
        if t₂ != 0
            nbr2 = build_neighbor_table(
                bonds[3], 
                model_geometry.unit_cell, 
                model_geometry.lattice
            )
            for col in axes(nbr2, 2)
                i, j = nbr2[1, col], nbr2[2, col]
                H_t[i, j]   += t₂
                H_t[j, i]   += t₂
                ip, jp = i + N, j + N
                H_t[ip, jp] += -t₂
                H_t[jp, ip] += -t₂
            end
        end
    else
        # nearest neighbors
        for col in axes(nbr0, 2)
            i, j = nbr0[1, col], nbr0[2, col]
            H_t[i, j] += -t₀
            H_t[j, i] += -t₀

            ip, jp = i + N, j + N
            H_t[ip, jp] += -t₀
            H_t[jp, ip] += -t₀
        end
        # next-nearest neighbors
        if t₁ != 0
            nbr1 = build_neighbor_table(
                bonds[2], 
                model_geometry.unit_cell, 
                model_geometry.lattice
            )
            for col in axes(nbr1, 2)
                i, j = nbr1[1, col], nbr1[2, col]
                H_t[i, j]   += t₁
                H_t[j, i]   += t₁
                ip, jp = i + N, j + N
                H_t[ip, jp] += t₁
                H_t[jp, ip] += t₁
            end
        end
        # third nearest neighbors
        if t₂ != 0
            nbr2 = build_neighbor_table(
                bonds[3], 
                model_geometry.unit_cell, 
                model_geometry.lattice
            )
            for col in axes(nbr2, 2)
                i, j = nbr2[1, col], nbr2[2, col]
                H_t[i, j]   += t₂
                H_t[j, i]   += t₂
                ip, jp = i + N, j + N
                H_t[ip, jp] += t₂
                H_t[jp, ip] += t₂
            end
        end
    end

    return H_t
end


@doc raw"""

    build_variational_hamiltonian( determinantal_parameters::DeterminantalParameters{I}, 
                                   model_geometry::ModelGeometry,
                                   optimize::NamedTuple, 
                                   pht::Bool ) where {I<:Integer}

Constructs a set of ``2N`` by ``2N`` matrices for each variational parameter, where `N` is the number of 
lattice sites. Returns a total variational Hamiltonian matrix ``H_{\mathrm{var}}`` as well has a vector 
of operators ``V``.

- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`: field of optimization flags.
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function build_variational_hamiltonian(
    determinantal_parameters::DeterminantalParameters{I}, 
    model_geometry::ModelGeometry,
    optimize::NamedTuple, 
    pht::Bool
) where {I<:Integer}    
    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    # number of sites
    N = model_geometry.lattice.N

    # bonds of the lattice
    bonds = model_geometry.bond
   
    # initialize Hamiltonian and operator matrices
    H_vpars = Vector{Matrix{AbstractFloat}}()
    V       = Vector{Matrix{AbstractFloat}}()

    if pht == true
        # add s-wave pairing 
        add_pairing_symmetry!(
            "s", 
            determinantal_parameters, 
            optimize, 
            H_vpars, 
            V, 
            dims,
            N,
            pht
        )

        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        adding s-wave pairing matrix =>
        initial Δ_0 = $(determinantal_parameters.det_pars.Δ_0)
        """

        if optimize.Δ_0
            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            optimize Δ_0 = true
            """
        else
            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            optimize Δ_0 = false
            """
        end

        if dims > 1
            # add sPDW pairing
            add_pairing_symmetry!(
                "sPDW",
                determinantal_parameters,
                optimize,
                H_vpars,
                bonds,
                V,
                dims,
                N,
                pht
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding site-dependent s-wave pairing matrix =>
            initial Δ_spd = $(determinantal_parameters.det_pars.Δ_spd)
            """

            if optimize.Δ_spd
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_spd = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_spd = false
                """
            end

            # add d-wave pairing 
            add_pairing_symmetry!(
                "d", 
                determinantal_parameters, 
                optimize, 
                H_vpars, 
                V, 
                bonds,
                dims,
                N,
                pht
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding d-wave pairing matrix =>
            initial Δ_d = $(determinantal_parameters.det_pars.Δ_d)
            """

            if optimize.Δ_d
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_d = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_d = false
                """
            end

            # add dPDW pairing
            add_pairing_symmetry!(
                "dPDW",
                determinantal_parameters,
                optimize,
                H_vpars,
                V,
                bonds,
                dims,
                N,
                pht
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding site-dependent d-wave pairing matrix =>
            initial Δ_dpd = $(determinantal_parameters.det_pars.Δ_dpd)
            """

            if optimize.Δ_dpd
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_dpd = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_dpd = false
                """
            end
        end
    end

    # add in-plane magnetization term
    add_spin_order!(
        "spin-x",
        determinantal_parameters,
        optimize,
        H_vpars,
        V,
        dims,
        N,
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding spin-x matrix =>
    initial Δ_sx = $(determinantal_parameters.det_pars.Δ_sx)
    """

    if optimize.Δ_sx
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize Δ_sx = true
        """
    else
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize Δ_sx = false
        """
    end

    # add antiferromagnetic (Neél) term
    add_spin_order!(
        "spin-z", 
        determinantal_parameters, 
        optimize, 
        H_vpars, 
        V, 
        dims,
        N, 
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding spin-z matrix =>
    initial Δ_sz = $(determinantal_parameters.det_pars.Δ_sz)
    """

    if optimize.Δ_sz
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize Δ_sz = true
        """
    else
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize Δ_sz = false
        """
    end

    # add site-dependent spin term
    if dims > 1
        add_spin_order!(
            "site-dependent", 
            determinantal_parameters, 
            optimize, 
            H_vpars, 
            V, 
            dims,
            N, 
            pht
        )

        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        adding site-dependent spin matrix =>
        initial Δ_ssd = $(determinantal_parameters.det_pars.Δ_ssd)
        """

        if optimize.Δ_ssd
            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            optimize Δ_ssd = true
            """
        else
            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            optimize Δ_ssd = false
            """
        end
    end

    # add chemical potential term
    add_chemical_potential!(
        determinantal_parameters, 
        optimize, 
        H_vpars, 
        V, 
        N,
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding chemical potential matrix =>
    initial μ = $(determinantal_parameters.det_pars.μ)
    """

    if optimize.μ
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize μ = true
        """
    else
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize μ = false
        """
    end

    # add charge-density-wave term
    add_charge_order!(
        "density wave", 
        determinantal_parameters, 
        optimize, 
        H_vpars, 
        V, 
        dims,
        N,
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding charge density wave matrix =>
    initial Δ_cdw = $(determinantal_parameters.det_pars.Δ_cdw)
    """

    if optimize.Δ_cdw
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize Δ_cdw = true
        """
    else
        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        optimize Δ_cdw = false
        """
    end

    # add site-dependent charge term
    if dims > 1
        add_charge_order!(
            "site-dependent", 
            determinantal_parameters, 
            optimize, 
            H_vpars,
            V, 
            dims,
            N, 
            pht
        )

        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        adding site-dependent charge matrix =>
        initial Δ_csd = $(determinantal_parameters.det_pars.Δ_csd)
        """

        if optimize.Δ_csd
            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            optimize Δ_csd = true
            """
        else
            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            optimize Δ_csd = false
            """
        end
    end

    # @assert length(H_vpars) == determinantal_parameters.num_det_pars
    # @assert length(V) == determinantal_parameters.num_det_opts

    return sum(H_vpars), V
end


@doc raw"""

    add_pairing_symmetry!( symmetry::S, 
                           determinantal_parameters::DeterminantalParameters{I}, 
                           optimize::NamedTuple, 
                           H_vpars::Vector{Matrix{T}}, 
                           V::Vector{Matrix{T}}, 
                           dims::I,
                           N::I, 
                           pht::Bool ) where {S<:AbstractString, I<:Integer, T<:Number}

Adds a pairing term to the auxiliary Hamiltonian. 

- `symmetry::S`: type of pairing symmetry: "s", "d", "sPDW", or "dPDW".
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Matrix{T}}`: vector of variational Hamiltonian matrices.
- `V::Vector{Matrix{T}}`: vector of variational operators.
- `dims::I`: dimensions of the lattice. 
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model is particle-hole transformed.

"""
function add_pairing_symmetry!(
    symmetry::S,
    determinantal_parameters::DeterminantalParameters{I},
    optimize::NamedTuple,
    H_vpars::Vector{Matrix{T}},
    V::Vector{Matrix{T}},
    bonds::Vector{Vector{Bond{I}}},
    dims::I,
    N::I,
    pht::Bool
) where {S<:AbstractString, I<:Integer, T<:Number}
    twoN = 2 * N
    unit_cell = model_geometry.unit_cell
    lattice = model_geometry.lattice

    # helper caches used by multiple branches
    # cache positions and spin-indices for sites 1..N (used by PDW branches)
    positions = Vector{NTuple{dims, Int}}(undef, N)
    for i in 1:N
        loc = site_to_loc(i, unit_cell, lattice)[1]
        # convert to NTuple{dims,Int}
        if dims == 1
            positions[i] = (loc[1],)
        else
            positions[i] = (loc[1], loc[2])
        end
    end

    # cache spindices mapping for index -> (up,down)
    sp_up = Vector{Int}(undef, N)
    sp_dn = Vector{Int}(undef, N)
    for i in 1:N
        up, dn = get_spindices_from_index(i, model_geometry)
        sp_up[i] = up
        sp_dn[i] = dn
    end

    # cache linked indices for the s-wave construction (keeps original zero-based call)
    linked = Vector{Int}(undef, twoN)
    for i0 in 0:(twoN - 1)
        linked[i0 + 1] = get_linked_spindex(i0, N) + 1
    end

    if symmetry == "s"
        @assert pht == true

        Δ_0 = determinantal_parameters.det_pars.Δ_0
        # Use element type T for matrices (keeps consistency with H_vpars/V types)
        V_s = zeros(T, twoN, twoN)
        @inbounds for i0 in 0:(twoN - 1)
            r = i0 + 1
            c = linked[r]
            V_s[r, c] = one(T)   # 1.0 but typed
        end

        # push scaled matrix (dense as before)
        Hs = Δ_0 * V_s
        push!(H_vpars, Hs)
        if optimize.Δ_0
            push!(V, V_s)
        end

    elseif symmetry == "d"
        @assert pht == true
        @assert dims > 1

        Δ_d = determinantal_parameters.det_pars.Δ_d

        # build neighbor table once
        nbr_table = build_neighbor_table(bonds[1], unit_cell, lattice)
        nbrs = map_neighbor_table(nbr_table)

        # small map: displacement (as NTuple) -> sign
        disp_sign_map = Dict{NTuple{2,Int}, Int}((1,0)=>1, (0,1)=>-1, (-1,0)=>1, (0,-1)=>-1)

        V_dwave = zeros(T, twoN, twoN)
        # loop sites
        @inbounds for i in 1:N
            nn = nbrs[i][2]           # neighbor indices (assumed as in original)
            idn_idx = sp_dn[i]       # get spin-down index for i (cached)
            iup_idx = sp_up[i]
            for j in nn
                # get displacement (returns something indexable like [dx,dy])
                disp_arr = sites_to_displacement(i, j, unit_cell, lattice)
                # convert to NTuple{2,Int} for dictionary lookup
                disp = (disp_arr[1], disp_arr[2])
                dsgn = get(disp_sign_map, disp, 0)
                if dsgn != 0
                    # get other site's spin indices (cache access)
                    jup = sp_up[j]
                    jdn = sp_dn[j]
                    # add symmetric pairing elements (typed to T)
                    V_dwave[iup_idx, jdn] += T(dsgn)
                    V_dwave[jup, idn_idx] += T(dsgn)
                    V_dwave[jdn, iup_idx] += T(dsgn)
                    V_dwave[idn_idx, jup] += T(dsgn)
                end
            end
        end

        push!(H_vpars, Δ_d * V_dwave)
        if optimize.Δ_d
            push!(V, V_dwave)
        end

    elseif symmetry == "sPDW"
        @assert pht == true
        @assert dims > 1

        # neighbor table once
        nbr_table = build_neighbor_table(bonds[1], unit_cell, lattice)

        q_p = determinantal_parameters.det_pars.q_p
        Δ_spd = determinantal_parameters.det_pars.Δ_spd

        # preallocate V_spd as a Vector of matrices (no push!)
        V_spd = [zeros(T, twoN, twoN) for _ in 1:N]

        # fill V_spd from neighbor columns:
        # The original code used two ranges for x and y neighbors; respecting the same layout:
        # columns 1:N are x-neighbors and N+1:2*N are y-neighbors (as original)
        @inbounds for col in 1:(2*N)
            i = nbr_table[1, col]
            j = nbr_table[2, col]
            # set elements (use cached spin indices)
            ip = i
            jp = j
            # up indices are sp_up[i], etc., but original code used i and j for up-sector and j+N etc.
            # We'll use mapped spin indices:
            up_i = sp_up[ip]
            dn_i = sp_dn[ip]
            up_j = sp_up[jp]
            dn_j = sp_dn[jp]
            # set pairing entries (match original structure: V_spd[i][i, j + N] etc.)
            V_spd[ip][up_i, dn_j] = one(T)
            V_spd[jp][up_j, dn_i] = one(T)
            V_spd[ip][dn_j, up_i] = one(T)   # hermitian conjugate entry
            V_spd[jp][dn_i, up_j] = one(T)
        end

        # Preallocate H matrices typed to T
        H_sff = zeros(T, twoN, twoN)
        H_slo = zeros(T, twoN, twoN)

        # assemble H_sff and H_slo in a single loop; push V entries if optimizing
        @inbounds for i in 1:N
            pos = positions[i]
            # convert position tuple to vector for dot; dot with q_p may need promotion to real/complex
            # use a small helper dotprod that matches original sites_to_loc semantics
            # here we assume q_p is indexable and positions are integer tuples
            qp_dot = zero(eltype(q_p))
            for k in 1:dims
                qp_dot += q_p[k] * pos[k]
            end
            ff_phase = exp(im * qp_dot)
            lo_phase = ff_phase + exp(-im * qp_dot)
            # Δ_spd[i] times phases; Δ_spd may be complex/real (assume compatible)
            H_sff .+= (Δ_spd[i] * ff_phase) .* V_spd[i]
            H_slo .+= (Δ_spd[i] * lo_phase) .* V_spd[i]

            if optimize.Δ_spd
                push!(V, V_spd[i])
            end
        end

        push!(H_vpars, H_sff)
        push!(H_vpars, H_slo)

    elseif symmetry == "dPDW"
        @assert pht == true
        @assert dims > 1

        # neighbor table
        nbr_table = build_neighbor_table(bonds[1], unit_cell, lattice)

        q_p = determinantal_parameters.det_pars.q_p
        Δ_dpd = determinantal_parameters.det_pars.Δ_dpd

        # preallocate per-site matrices
        V_dpdx = [zeros(T, twoN, twoN) for _ in 1:N]
        V_dpdy = [zeros(T, twoN, twoN) for _ in 1:N]
        V_dpd  = [zeros(T, twoN, twoN) for _ in 1:N]

        # x-direction neighbors in columns 1:N
        @inbounds for col in 1:N
            i = nbr_table[1, col]
            j = nbr_table[2, col]
            up_i = sp_up[i]; dn_j = sp_dn[j]
            up_j = sp_up[j]; dn_i = sp_dn[i]
            V_dpdx[i][up_i, dn_j] = one(T)
            V_dpdx[j][up_j, dn_i] = one(T)
            V_dpdx[i][dn_j, up_i] = one(T)
            V_dpdx[j][dn_i, up_j] = one(T)
        end

        # y-direction neighbors in columns N+1:2N with negative sign
        @inbounds for col in (N+1):(2*N)
            i = nbr_table[1, col]
            j = nbr_table[2, col]
            up_i = sp_up[i]; dn_j = sp_dn[j]
            up_j = sp_up[j]; dn_i = sp_dn[i]
            V_dpdy[i][up_i, dn_j] = -one(T)
            V_dpdy[j][up_j, dn_i] = -one(T)
            V_dpdy[i][dn_j, up_i] = -one(T)
            V_dpdy[j][dn_i, up_j] = -one(T)
        end

        # combine x and y components into V_dpd
        @inbounds for i in 1:N
            V_dpd[i] .+= V_dpdx[i]
            V_dpd[i] .+= V_dpdy[i]
        end

        # allocate H accumulators
        H_dff = zeros(T, twoN, twoN)
        H_dlo = zeros(T, twoN, twoN)

        # assemble H matrices and optionally collect V_dpd for optimization
        @inbounds for i in 1:N
            pos = positions[i]
            qp_dot = zero(eltype(q_p))
            for k in 1:dims
                qp_dot += q_p[k] * pos[k]
            end
            ff_phase = exp(im * qp_dot)
            lo_phase = ff_phase + exp(-im * qp_dot)
            H_dff .+= (Δ_dpd[i] * ff_phase) .* V_dpd[i]
            H_dlo .+= (Δ_dpd[i] * lo_phase) .* V_dpd[i]

            if optimize.Δ_spd   # original had this likely intended for dPDW as well
                push!(V, V_dpd[i])
            end
        end

        push!(H_vpars, H_dff)
        push!(H_vpars, H_dlo)
    end

    return nothing
end



@doc raw"""

    add_spin_order!( order::S, 
                     determinantal_parameters::DeterminantalParameters{I}, 
                     optimize::NamedTuple, 
                     H_vpars::Vector{Matrix{T}}, 
                     V::Vector{Matrix{T}}, 
                     dims::I,
                     N::I, 
                     pht::Bool )::Nothing

Adds a spin ordering term to the auxiliary Hamiltonian. 

- `order::String`: type of spin order: "spin-x", "spin-z", or "site-dependent"
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Any}`: vector of variational Hamiltonian matrices.
- `V::Vector{Any}`: vector of variational operators.
- `dims::I`: dimensions of the lattice. 
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_spin_order!(
    order::S,
    determinantal_parameters::DeterminantalParameters{I},
    optimize::NamedTuple,
    H_vpars::Vector{Matrix{T}},
    V::Vector{Matrix{T}},
    dims::I,
    N::I,
    pht::Bool
) where {S<:AbstractString, I<:Integer, T<:Number}
    twoN = 2 * N

    # Precompute mapping from spindex -> site index and site coordinates (one call each)
    idxs = Vector{Int}(undef, twoN)
    # store coordinates as a small tuple (ix, iy) or (ix,)
    if dims == 1
        locs = Vector{Int}(undef, twoN)  # store ix only
        for s in 1:twoN
            idx = get_index_from_spindex(s, model_geometry)
            idxs[s] = idx
            loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
            # site_to_loc(...)[1][1] in your original code — capture single integer
            locs[s] = loc[1][1]
        end
    else
        locs = Vector{NTuple{2,Int}}(undef, twoN)
        for s in 1:twoN
            idx = get_index_from_spindex(s, model_geometry)
            idxs[s] = idx
            loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
            # capture (ix, iy)
            locs[s] = (loc[1][1], loc[1][2])
        end
    end

    if order == "spin-z"
        # Use scalars / element types from the parameter to allocate minimal objects
        Δ_sz = determinantal_parameters.det_pars.Δ_sz
        # mask vector as small int type to reduce memory (will be promoted when multiplied by Δ_sz)
        afm_vec = ones(Int8, twoN)

        if pht
            # stagger using precomputed locs; compute parity with isodd for speed
            if dims == 1
                for s in 1:twoN
                    ix = locs[s]  # Int
                    afm_vec[s] = isodd(ix) ? -1 : 1
                end
            else
                for s in 1:twoN
                    (ix, iy) = locs[s]
                    afm_vec[s] = isodd(ix + iy) ? -1 : 1
                end
            end

            # push Δ * Diagonal(afm_vec) directly (no intermediate dense zero matrix)
            V_afm = LinearAlgebra.Diagonal(afm_vec)           # this is cheap (stores diag only)
            push!(H_vpars, Δ_sz * V_afm)                      # scaled Diagonal pushed
            if optimize.Δ_sz
                push!(V, V_afm)
            end
        else
            # create center-symmetric vec for spin sectors: flip sign in second half
            afm_vec_neg = similar(afm_vec)
            # fill first half with 1/-1 based on parity then second half = - first half
            if dims == 1
                for s in 1:N
                    ix = locs[s]
                    v = isodd(ix) ? -1 : 1
                    afm_vec_neg[s] = v
                    afm_vec_neg[s + N] = -v  # flip sign for spin-down sector
                end
            else
                for s in 1:N
                    (ix, iy) = locs[s]
                    v = isodd(ix + iy) ? -1 : 1
                    afm_vec_neg[s] = v
                    afm_vec_neg[s + N] = -v
                end
            end

            V_afm_neg = LinearAlgebra.Diagonal(afm_vec_neg)
            push!(H_vpars, Δ_sz * V_afm_neg)
            if optimize.Δ_sz
                push!(V, V_afm_neg)
            end
        end

    elseif order == "spin-x"
        Δ_sx = determinantal_parameters.det_pars.Δ_sx
        # choose matrix element type according to Δ_sx
        Mty = typeof(Δ_sx)
        H_sx = zeros(Mty, twoN, twoN)  # dense as before, but typed to Mty

        # fill off-diagonals compactly
        half = N
        half_off = Δ_sx / 2
        for s in 1:half
            up_idx = s
            dn_idx = s + half
            H_sx[up_idx, dn_idx] += half_off
            H_sx[dn_idx, up_idx] += half_off
        end

        push!(H_vpars, H_sx)
        if optimize.Δ_sx
            push!(V, H_sx)
        end

    elseif order == "site-dependent"
        L = model_geometry.lattice.L
        Δ_ssd = determinantal_parameters.det_pars.Δ_ssd
        # number of shifts along unit cell direction 1
        nshifts = L[1]

        # process each shift immediately (avoid storing all ssd_vectors)
        for shift in 0:(nshifts - 1)
            # small-int mask vector
            ssd_vec = zeros(Int8, twoN)

            # set every L[1]th element according to pattern
            # loop through up-to 2*L[1] positions as original did, but stop at twoN
            # compute base indices and set 1 where valid
            for i in 1:(2 * L[1])
                idx = (i - 1) * L[1] + 1 + shift
                if idx > twoN
                    break
                end
                ssd_vec[idx] = 1
            end

            if pht
                # flip sign on spin-down half (in-place)
                @inbounds for j in (N+1):twoN
                    ssd_vec[j] = -ssd_vec[j]
                end
            end

            # compute staggering factor from locations and apply (in-place)
            if dims == 1
                for s in 1:twoN
                    ix = locs[s]
                    ssd_vec[s] = ssd_vec[s] * (isodd(ix) ? -1 : 1)
                end
            else
                for s in 1:twoN
                    (ix, iy) = locs[s]
                    ssd_vec[s] = ssd_vec[s] * (isodd(ix + iy) ? -1 : 1)
                end
            end

            # create diagonal and push scaled H_vpar
            V_ssd = LinearAlgebra.Diagonal(ssd_vec)
            push!(H_vpars, Δ_ssd[shift + 1] * V_ssd)  # Δ_ssd is 1-based indexed
            if optimize.Δ_ssd
                push!(V, V_ssd)
            end
        end
    end

    return nothing
end


@doc raw"""

    add_charge_order!( order::S, 
                       determinantal_parameters::DeterminantalParameters{I}, 
                       optimize::NamedTuple, 
                       H_vpars::Vector{Matrix{T}}, 
                       V::Vector{Matrix{T}}, 
                       dims::I,
                       N::I, 
                       pht::Bool ) where {S<:AbstractString, I<:Integer, T<:Number}

Adds a charge ordering term to the auxiliary Hamiltonian.

- `order::S`: type of spin order: "density wave" or "site-dependent"
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Matrix{T}}`: vector of variational Hamiltonian matrices.
- `V::Vector{Matrix{T}}`: vector of variational operators.
- `dims::I`: dimensions of the lattice. 
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_charge_order!(
    order::S,
    determinantal_parameters::DeterminantalParameters{I},
    optimize::NamedTuple,
    H_vpars::Vector{Matrix{T}},
    V::Vector{Matrix{T}},
    dims::I,
    N::I,
    pht::Bool
) where {S<:AbstractString, I<:Integer, T<:Number}
    twoN = 2 * N

    # Precompute mapping from spindex -> site index and site coordinates (one call each)
    idxs = Vector{Int}(undef, twoN)
    if dims == 1
        locs = Vector{Int}(undef, twoN)
        for s in 1:twoN
            idx = get_index_from_spindex(s, model_geometry)
            idxs[s] = idx
            loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
            locs[s] = loc[1][1]
        end
    else
        locs = Vector{NTuple{2,Int}}(undef, twoN)
        for s in 1:twoN
            idx = get_index_from_spindex(s, model_geometry)
            idxs[s] = idx
            loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
            locs[s] = (loc[1][1], loc[1][2])
        end
    end

    if order == "density wave"
        Δ_cdw = determinantal_parameters.det_pars.Δ_cdw
        # small-int mask
        cdw_vec = ones(Int8, twoN)

        if pht
            # prepare negated second half
            cdw_vec_neg = similar(cdw_vec)
            # set first half and second half with sign flip in-place
            if dims == 1
                for s in 1:N
                    ix = locs[s]
                    v = isodd(ix) ? -1 : 1
                    cdw_vec_neg[s] = v
                    cdw_vec_neg[s + N] = -v
                end
            else
                for s in 1:N
                    (ix, iy) = locs[s]
                    v = isodd(ix + iy) ? -1 : 1
                    cdw_vec_neg[s] = v
                    cdw_vec_neg[s + N] = -v
                end
            end

            V_cdw = LinearAlgebra.Diagonal(cdw_vec_neg)
            push!(H_vpars, Δ_cdw * V_cdw)
            if optimize.Δ_cdw
                push!(V, V_cdw)
            end
        else
            # non-pht: parity only
            if dims == 1
                for s in 1:twoN
                    ix = locs[s]
                    cdw_vec[s] = isodd(ix) ? -1 : 1
                end
            else
                for s in 1:twoN
                    (ix, iy) = locs[s]
                    cdw_vec[s] = isodd(ix + iy) ? -1 : 1
                end
            end

            V_cdw = LinearAlgebra.Diagonal(cdw_vec)
            push!(H_vpars, Δ_cdw * V_cdw)
            if optimize.Δ_cdw
                push!(V, V_cdw)
            end
        end

    elseif order == "site-dependent"
        L = model_geometry.lattice.L
        Δ_csd = determinantal_parameters.det_pars.Δ_csd
        nshifts = L[1]

        for shift in 0:(nshifts - 1)
            # small-int mask vector
            csd_vec = zeros(Int8, twoN)

            # set every L[1]th element according to pattern
            for i in 1:(2 * L[1])
                idx = (i - 1) * L[1] + 1 + shift
                if idx > twoN
                    break
                end
                csd_vec[idx] = 1
            end

            if pht
                @inbounds for j in (N + 1):twoN
                    csd_vec[j] = -csd_vec[j]
                end
            end

            # apply parity/stagger
            if dims == 1
                for s in 1:twoN
                    ix = locs[s]
                    csd_vec[s] = csd_vec[s] * (isodd(ix) ? -1 : 1)
                end
            else
                for s in 1:twoN
                    (ix, iy) = locs[s]
                    csd_vec[s] = csd_vec[s] * (isodd(ix + iy) ? -1 : 1)
                end
            end

            V_csd = LinearAlgebra.Diagonal(csd_vec)
            push!(H_vpars, Δ_csd[shift + 1] * V_csd)
            if optimize.Δ_csd
                push!(V, V_csd)
            end
        end
    end

    return nothing
end


@doc raw"""

    add_chemical_potential!( determinantal_parameters::DeterminantalParameters{I}, 
                             optimize::NamedTuple, 
                             H_vpars::Vector{Matrix{T}}, 
                             V::Vector{Matrix{T}}, 
                             N::I, 
                             pht::Bool )::Nothing

Adds a chemical potential term to the auxiliary Hamiltonian.

- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Matrix{T}}`: vector of variational Hamiltonian matrices.
- `V::Vector{Matrix{T}}`: vector of variational operators.
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_chemical_potential!(
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    H_vpars::Vector{Matrix{T}}, 
    V::Vector{Matrix{T}}, 
    N::I, 
    pht::Bool
) where {I<:Integer, T<:Number}
    # initial chemical potential value
    μ = determinantal_parameters.det_pars.μ

    μ_vec = fill(-1,2*N)
    if pht
        # account for minus sign 
        μ_vec_neg = copy(μ_vec)
        μ_vec_neg[N+1:2*N] .= -μ_vec_neg[N+1:2*N]

        # add chemical potential matrix
        H_μ = zeros(Number, 2*N, 2*N)
        V_μ_neg = LinearAlgebra.Diagonal(μ_vec_neg)
        H_μ += μ * V_μ_neg
        push!(H_vpars, H_μ)

        # if μ is being optimized, save Vμ matrix
        if optimize.μ
            push!(V, V_μ_neg)
        end
    else
        # add chemical potential matrix
        H_μ = zeros(Number, 2*N, 2*N)
        V_μ = LinearAlgebra.Diagonal(μ_vec)
        H_μ += μ * V_μ
        push!(H_vpars,H_μ)

        # if μ is being optimized, save Vμ matrix
        if optimize.μ
            push!(V, V_μ)
        end
    end

    return nothing
end


@doc raw"""

    diagonalize!( H::Matrix{T} ) where {T<:Number}

Returns all eigenenergies ``\varepsilon`` and all eigenstates ``\phi`` of the 
auxiliary Hamiltonian matrix. All eigenstates are stored in the columns of a 
unitary matrix `U_aux`. 

- `H::Matrix{T}`: auxiliary Hamiltonian matrix.

"""
function diagonalize!(
    H::Matrix{T}
) where {T<:Number}
    # check if Hamiltonian is Hermitian
    @assert ishermitian(H) 

    if T <: Real
        # in-place path exists
        F = eigen!(H)
    else
        # Complex Hermitian must use non-mutating eigen
        F = eigen(H)
    end

    return F.values, F.vectors  
end


@doc raw"""

    is_openshell( ε::Vector{E},  
                  Np::I ) where {E<:AbstractFloat, I<:Integer}

Checks whether a mean-field energy configuration is open shell.

- `ε::Vector{E}`: vector of mean-field energies.
- `Np::I`: total number of particles in the system.

"""
function is_openshell(
    ε::Vector{E}, 
    Np::I
) where {E<:AbstractFloat, I<:Integer}
    return ε[Np + 1] - ε[Np] < 1e-4 #1e-12
end


@doc raw"""

    get_variational_matrices( V::Vector{Matrix{T1}}, 
                              U_aux::Matrix{T2}, 
                              ε::Vector{E}, 
                              Np::I
                              model_geometry::ModelGeometry ) where {T1<:Number, T2<:Number, E<:AbstractFloat, I<:Integer}
    
Returns a set of variational parameter matrices ``A_k`` constructed from each variational operator ``V_k`` by computing 
```math
Q_k = \frac{(U^{\dagger}V_{k}U)_{\eta\nu}}{(\varepsilon_{\eta} - \varepsilon_{\nu})},
```
where ``\eta > N_p`` and ``\nu \leq N_p``.

- `V::Vector{Matrix{T1}`: vector of variational operators. 
- `U_aux::Matrix{T2}`: matrix which diagonalizes the auxiliary Hamiltonian.
- `ε::Vector{E}`: initial energies.
- `Np::I`: number of particles in the system. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_variational_matrices(
    V::Vector{Matrix{T1}}, 
    U_aux::Matrix{T2}, 
    ε::Vector{E}, 
    Np::I,
    model_geometry::ModelGeometry
) where {T1<:Number, T2<:Number, E<:AbstractFloat, I<:Integer}
    # number of lattice sites
    N = model_geometry.unit_cell.n * model_geometry.lattice.N

    # define perturbation mask
    ptmask = zeros(AbstractFloat, 2*N, 2*N)
        
    for η in 1:2*N
        for ν in 1:2*N
            if η > Np && ν <= Np
                ptmask[η, ν] = 1.0 / (ε[ν] - ε[η])
            end
        end
    end

    int_A = Vector{Matrix{<:Number}}()
    for v in V
        push!(int_A, U_aux * ((U_aux' * v * U_aux) .* ptmask) * U_aux')
    end

    return int_A
end


@doc raw"""

    get_tb_chem_pot( Ne::I, 
                     tight_binding_model::TightBindingModel{E}, 
                     model_geometry::ModelGeometry ) where {I<:Integer, E<:AbstractFloat}

Computes the appropriate chemical potential ``\mu`` for a non-interacting tight-binding model
on a finite lattice.

- `Ne::I`: total number of electrons.
- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_tb_chem_pot(
    Ne::I, 
    tight_binding_model::TightBindingModel{E}, 
    model_geometry::ModelGeometry
) where {I<:Integer, E<:AbstractFloat}
    # number of lattice sites
    N = model_geometry.lattice.N

    # bonds of the lattice
    bonds = model_geometry.bond
    
    # preallocate matrices
    H_t₀ = zeros(Complex, 2*N, 2*N)
    H_t₁ = zeros(Complex, 2*N, 2*N)
    H_t₂ = zeros(Complex, 2*N, 2*N)

    # hopping amplitudes
    t₀ = tight_binding_model.t₀
    t₁ = tight_binding_model.t₁
    t₂ = tight_binding_model.t₂

    # add nearest neighbor hopping
    nbr0 = build_neighbor_table(
        bonds[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # add nearest neighbor hopping
    for (i,j) in eachcol(nbr0)
        H_t₀[i,j] += -t₀
        if model_geometry.lattice.N > 2
            H_t₀[j,i] += -t₀
        end
    end
    for (i,j) in eachcol(nbr0 .+ N)    
        H_t₀[i,j] += -t₀
        if model_geometry.lattice.N > 2
            H_t₀[j,i] += -t₀
        end
    end

    if t₁ != 0.0
        # build next-nearest neighbor table
        nbr1 = build_neighbor_table(
            bonds[2],
            model_geometry.unit_cell,
            model_geometry.lattice
        )

        # add next-nearest neighbor hopping
        for (i,j) in eachcol(nbr1)
            H_t₁[i,j] += t₁
            H_t₁[j,i] += t₁
        end
        for (i,j) in eachcol(nbr1 .+ N)    
            H_t₁[i,j] += t₁
            H_t₁[j,i] += t₁
        end
    end

    if t₂ != 0.0
        # build third nearest neighbor table
        nbr2 = build_neighbor_table(
            bonds[3],
            model_geometry.unit_cell,
            model_geometry.lattice
        )

        # add third nearest neighbor hopping
        for (i,j) in eachcol(nbr2)
            H_t₂[i,j] += t₂
            H_t₂[j,i] += t₂
        end
        for (i,j) in eachcol(nbr2 .+ N)    
            H_t₂[i,j] += t₂
            H_t₂[j,i] += t₂
        end
    end

    # full tight-binding Hamiltonian
    H_tb = H_t₀ + H_t₁ + H_t₂

    # solve for eigenvalues
    ε_F, _ = diagonalize!(H_tb)

    # tight-binding chemical potential
    μ = 0.5 * (ε_F[Ne + 1] + ε_F[Ne])

    @debug """
    Hamiltonian::get_tb_chem_pot() :
    tight-binding chemical potential =>
    μ = $(μ)
    """

    return μ
end


















