@doc raw"""

    build_auxiliary_hamiltonian( tight_binding_model::TightBindingModel{E}, 
                                 determinantal_parameters::DeterminantalParameters{I}, 
                                 optimize::NamedTuple, 
                                 model_geometry::ModelGeometry, 
                                 pht::Bool;
                                 q_p::AbstractVector{T} = [0.0, 0.0] ) where {E<:AbstractFloat, I<:Integer}

Constructs an auxiliary Hamiltonian matrix ``H_{\mathrm{aux}}`` by combining the non-interacting
hopping matrix ``H_t`` with matrices of variational terms ``H_{\mathrm{var}}``.

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hol transformed.
- `q_p::AbstractVector{T} = [0.0, 0.0]`: pairing momentum for density wave pairing. 

"""
function build_auxiliary_hamiltonian(
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    pht::Bool;
    q_p = [0.0, 0.0]
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
        pht;
        q_p = q_p
    )

    return H_tb + H_var, V
end


@doc raw"""

    build_auxiliary_hamiltonian( tight_binding_model::TightBindingModel{E}, 
                                 determinantal_parameters::DeterminantalParameters{I}, 
                                 optimize::NamedTuple, 
                                 model_geometry::ModelGeometry, 
                                 twist_angles::AbstractRange{E},
                                 pht::Bool;
                                 q_p::AbstractVector{T} = [0.0, 0.0] ) where {E<:AbstractFloat, I<:Integer}

Constructs an auxiliary Hamiltonian matrix ``H_{\mathrm{aux}}^{\theta}`` for ``N_\theta`` twist angles
by combining the non-interacting hopping matrix ``H_t`` with matrices of variational terms ``H_{\mathrm{var}}``.

- `tight_binding_model::TightBindingModel{E}`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `twist_angles::AbstractRange{E}`: set of twist angles.
- `pht::Bool`: whether model is particle-hol transformed.
- `q_p::AbstractVector{T} = [0.0, 0.0]`: pairing momentum for density wave pairing. 

"""
function build_auxiliary_hamiltonian(
    tight_binding_model::TightBindingModel{E}, 
    determinantal_parameters::DeterminantalParameters{I}, 
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    twist_angles::AbstractRange{E},
    pht::Bool;
    q_p = [0.0, 0.0]
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
        pht;
        q_p = q_p
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
                                   pht::Bool;
                                   q_p::AbstractVector{T} = [0.0, 0.0] ) where {I<:Integer}

Constructs a set of ``2N`` by ``2N`` matrices for each variational parameter, where `N` is the number of 
lattice sites. Returns a total variational Hamiltonian matrix ``H_{\mathrm{var}}`` as well has a vector 
of operators ``V``.

- `determinantal_parameters::DeterminantalParameters{I}`: set of determinantal variational parameters.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `optimize::NamedTuple`: field of optimization flags.
- `pht::Bool`: whether model is particle-hole transformed. 
- `q_p::AbstractVector{T} = [0.0, 0.0]`: pairing momentum for density wave pairing. 

"""
function build_variational_hamiltonian(
    determinantal_parameters::DeterminantalParameters{I}, 
    model_geometry::ModelGeometry,
    optimize::NamedTuple, 
    pht::Bool;
    q_p = [0.0, 0.0]
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

    # get all determinantal parameters
    pars = determinantal_parameters.det_pars

    # Precompute mapping from spindex -> site index and site coordinates (one call each)
    idxs = Vector{Int}(undef, 2*N)
    if dims == 1
        locs = Vector{Int}(undef, 2*N)
        @inbounds for s in 1:2*N
            idx = get_index_from_spindex(s, model_geometry)
            idxs[s] = idx
            loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
            locs[s] = loc[1][1]
        end
    else
        locs = Vector{NTuple{2,Int}}(undef, 2*N)
        @inbounds for s in 1:2*N
            idx = get_index_from_spindex(s, model_geometry)
            idxs[s] = idx
            loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
            locs[s] = (loc[1][1], loc[1][2])
        end
    end 

    if pht == true
        # add s-wave pairing 
        add_pairing_symmetry!(
            "s", 
            optimize, 
            H_vpars, 
            V, 
            model_geometry,
            dims,
            bonds,
            locs,
            N,
            pht;
            q_p = q_p,
        )

        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        adding s-wave pairing matrix =>
        initial Δ_0 = $(pars.Δ_0)
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
            # add sLO pairing
            add_pairing_symmetry!(
                "sLO",
                optimize,
                H_vpars,
                V,
                model_geometry,
                dims,
                bonds,
                locs,
                N,
                pht;
                q_p = q_p
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding site-dependent s-wave pairing (Larkin-Ovchinnikov-type) matrix =>
            initial Δ_slo = $(pars.Δ_slo)
            """

            if optimize.Δ_slo
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_slo = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_slo = false
                """
            end

            # add sFF pairing
            add_pairing_symmetry!(
                "sFF",
                optimize,
                H_vpars,
                V,
                model_geometry,
                dims,
                bonds,
                locs,
                N,
                pht;
                q_p = q_p
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding site-dependent s-wave pairing (Fulde-Ferrell-type) matrix =>
            initial Δ_sff = $(pars.Δ_sff)
            """

            if optimize.Δ_sff
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_sff = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_sff = false
                """
            end

            # add d-wave pairing 
            add_pairing_symmetry!(
                "d", 
                optimize, 
                H_vpars, 
                V, 
                model_geometry,
                dims,
                bonds,
                locs,
                N,
                pht;
                q_p = q_p
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding d-wave pairing matrix =>
            initial Δ_d = $(pars.Δ_d)
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

            # add dLO pairing
            add_pairing_symmetry!(
                "dLO",
                optimize,
                H_vpars,
                V,
                model_geometry,
                dims,
                bonds,
                locs,
                N,
                pht;
                q_p = q_p
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding site-dependent d-wave (Larkin-Ovchinnikov-type) pairing matrix =>
            initial Δ_dlo = $(pars.Δ_dlo)
            """

            if optimize.Δ_dlo
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_dlo = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_dlo = false
                """
            end

            # add dFF pairing
            add_pairing_symmetry!(
                "dFF",
                optimize,
                H_vpars,
                V,
                model_geometry,
                dims,
                bonds,
                locs,
                N,
                pht;
                q_p = q_p
            )

            @debug """
            Hamiltonian::build_variational_hamiltonian() :
            adding site-dependent d-wave (Fulde-Ferrell-type) pairing matrix =>
            initial Δ_dff = $(pars.Δ_dff)
            """

            if optimize.Δ_dff
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_dff = true
                """
            else
                @debug """
                Hamiltonian::build_variational_hamiltonian() :
                optimize Δ_dff = false
                """
            end
        end
    end

    # add in-plane magnetization term
    add_spin_order!(
        "spin-x",
        optimize,
        H_vpars,
        V,
        model_geometry,
        locs,
        dims,
        N,
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding spin-x matrix =>
    initial Δ_sx = $(pars.Δ_sx)
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
        optimize, 
        H_vpars, 
        V, 
        model_geometry,
        locs,
        dims,
        N, 
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding spin-z matrix =>
    initial Δ_sz = $(pars.Δ_sz)
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
            optimize, 
            H_vpars, 
            V, 
            model_geometry,
            locs,
            dims,
            N, 
            pht
        )

        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        adding site-dependent spin matrix =>
        initial Δ_ssd = $(pars.Δ_ssd)
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
        optimize, 
        H_vpars, 
        V, 
        N,
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding chemical potential matrix =>
    initial μ = $(pars.μ)
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
        optimize, 
        H_vpars, 
        V, 
        model_geometry,
        locs,
        dims,
        N,
        pht
    )

    @debug """
    Hamiltonian::build_variational_hamiltonian() :
    adding charge density wave matrix =>
    initial Δ_cdw = $(pars.Δ_cdw)
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
            optimize, 
            H_vpars,
            V, 
            model_geometry,
            locs,
            dims,
            N, 
            pht
        )

        @debug """
        Hamiltonian::build_variational_hamiltonian() :
        adding site-dependent charge matrix =>
        initial Δ_csd = $(pars.Δ_csd)
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

    # allocate final Hamiltonian
    H_par = zero(H_vpars[1])

    # apply each variational parameter 
    hidx = 1
    @inbounds for p in pars
        if p isa AbstractVector
            for α in p
                H_par += α * H_vpars[hidx]
                hidx += 1
            end
        else
            H_par += p * H_vpars[hidx]
            hidx += 1
        end
    end 

    return H_par, V
end


@doc raw"""

    add_pairing_symmetry!( symmetry::S, 
                           optimize::NamedTuple, 
                           H_vpars::Vector{Matrix{T}}, 
                           V::Vector{Matrix{T}}, 
                           model_geometry::ModelGeometry,
                           bonds,
                           dims::I,
                           N::I, 
                           pht::Bool;
                           q_p::AbstractVector{T} = [0.0, 0.0],
                           locs::Vector{Tuple{I, I}} = Tuple{I,I}[] ) where {S<:AbstractString, I<:Integer, T<:Number}

Adds a pairing term to the auxiliary Hamiltonian along with its perturbative operator. 

- `symmetry::S`: type of pairing symmetry: "s", "d", "sLO", "sFF", "dLO", or "dFF".
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Matrix{T}}`: vector of variational Hamiltonian matrices.
- `V::Vector{Matrix{T}}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains lattice and unit cell quantities.
- `dims::I`: dimensions of the lattice. 
- `bonds`: lattice bonds.
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model is particle-hole transformed.
- `q_p::AbstractVector{T} = [0.0, 0.0]`: pairing momentum for density wave pairing.
- `locs::Vector{Tuple{I, I}} = Tuple{I,I}[]`: contains all positions (or locations) of each lattice site.

"""
function add_pairing_symmetry!(
    symmetry::S,
    optimize::NamedTuple,
    H_vpars::Vector{Matrix{T}},
    V::Vector{Matrix{T}},
    model_geometry::ModelGeometry,
    dims::I,
    bonds,
    locs,
    N::I,
    pht::Bool;
    q_p = [0.0, 0.0]
) where {S<:AbstractString, I<:Integer, T<:Number} 
    @assert pht == true

    twoN = 2 * N

    if symmetry == "s"
        V_s = zeros(T, twoN, twoN)

        @inbounds for r in 1:twoN
            c = get_linked_spindex(r - 1, N) + 1
            V_s[r, c] = one(T)
        end

        push!(H_vpars, V_s)
        if optimize.Δ_0
            push!(V, V_s)
        end
    elseif symmetry == "d"
        @assert dims > 1

        unit_cell = model_geometry.unit_cell
        lattice = model_geometry.lattice

        # cache spin indices
        sp_up = Vector{Int}(undef, N)
        sp_dn = Vector{Int}(undef, N)
        @inbounds for i in 1:N
            up, dn = get_spindices_from_index(i, model_geometry)
            sp_up[i] = up
            sp_dn[i] = dn
        end

        # build neighbor table
        nbr_table = build_neighbor_table(bonds[1], unit_cell, lattice)
        ncols = size(nbr_table, 2)
        
        V_dwave = zeros(T, twoN, twoN)

        @inbounds for col in 1:ncols
            i = nbr_table[1, col]
            j = nbr_table[2, col]

            # convert sites to displacement
            disp = sites_to_displacement(i, j, unit_cell, lattice)
            dx = disp[1]
            dy = disp[2]

            # d-wave phase: +1 for x-bonds, -1 for y-bonds
            dsgn = (dx != 0) - (dy != 0)
            dsgn == 0 && continue

            si_up = sp_up[i]
            si_dn = sp_dn[i]
            sj_up = sp_up[j]
            sj_dn = sp_dn[j]

            V_dwave[si_up, sj_dn] += T(dsgn)
            V_dwave[sj_up, si_dn] += T(dsgn)
            V_dwave[sj_dn, si_up] += T(dsgn)
            V_dwave[si_dn, sj_up] += T(dsgn)
        end

        push!(H_vpars, V_dwave)
        if optimize.Δ_d
            push!(V, V_dwave)
        end
    elseif symmetry == "sLO" || symmetry == "sFF"
        @assert dims > 1

        unit_cell = model_geometry.unit_cell
        lattice = model_geometry.lattice

        # build neighbor table
        nbr_table = build_neighbor_table(bonds[1], unit_cell, lattice)

        # cache spin indices
        sp_up = Vector{Int}(undef, N)
        sp_dn = Vector{Int}(undef, N)
        @inbounds for i in 1:N
            up, dn = get_spindices_from_index(i, model_geometry)
            sp_up[i] = up
            sp_dn[i] = dn
        end

        Vspdw = [zeros(T, twoN, twoN) for _ in 1:N]

        @inbounds for col in axes(nbr_table, 2)
            i = nbr_table[1, col]
            j = nbr_table[2, col]

            up_i = sp_up[i]
            dn_i = sp_dn[i]
            up_j = sp_up[j]
            dn_j = sp_dn[j]

            Vspdw[i][up_i, dn_j] += one(T)
            Vspdw[i][up_j, dn_i] += one(T)
            Vspdw[i][dn_j, up_i] += one(T)
            Vspdw[i][dn_i, up_j] += one(T)
        end

        # apply the appropriate phase for either the LO or FF state
        V_phase = Vector{Matrix{eltype(Vspdw[1])}}(undef, length(Vspdw))
        V_sign  = Vector{Matrix{eltype(Vspdw[1])}}(undef, length(Vspdw))

        for i in eachindex(Vspdw)
            V_phase[i] = similar(Vspdw[i])
            V_sign[i]  = similar(Vspdw[i])
        end

        @inbounds for i in 1:N
            r = locs[i]
            θ = q_p[1]*r[1] + q_p[2]*r[2]

            phase = if symmetry == "sLO"
                cis(θ) + cis(-θ)
            elseif symmetry == "sFF"
                cis(θ)
            end

            amp = abs(phase)
            sgn = iszero(amp) ? zero(phase) : phase / amp

            A  = Vspdw[i]
            Bp = V_phase[i]
            Bs = V_sign[i]

            @inbounds @simd for k in eachindex(A)
                a = A[k]
                Bp[k] = a * phase   # full phase for H_vpars
                Bs[k] = a * sgn     # unit phase for V
            end
        end

        append!(H_vpars, V_phase)

        if optimize.Δ_slo && symmetry == "sLO"
            append!(V, V_sign)
        end
        if optimize.Δ_sff && symmetry == "sFF"
            append!(V, V_sign)
        end
    elseif symmetry == "dLO" || symmetry == "dFF"
        @assert dims > 1

        unit_cell = model_geometry.unit_cell
        lattice = model_geometry.lattice

        # build neighbor table
        nbr_table = build_neighbor_table(bonds[1], unit_cell, lattice)

        # cache spin indices
        sp_up = Vector{Int}(undef, N)
        sp_dn = Vector{Int}(undef, N)
        @inbounds for i in 1:N
            up, dn = get_spindices_from_index(i, model_geometry)
            sp_up[i] = up
            sp_dn[i] = dn
        end

        Vdpdw = [zeros(T, twoN, twoN) for _ in 1:N]

        @inbounds for col in 1:(2N)
            i = nbr_table[1, col]
            j = nbr_table[2, col]

            sgn = col ≤ N ? one(T) : -one(T)

            up_i = sp_up[i]; dn_i = sp_dn[i]
            up_j = sp_up[j]; dn_j = sp_dn[j]

            Vi = Vdpdw[i]
            Vj = Vdpdw[j]

            Vi[up_i, dn_j] += sgn
            Vj[up_j, dn_i] += sgn
            Vi[dn_j, up_i] += sgn
            Vj[dn_i, up_j] += sgn
        end

        # apply the appropriate phase for either the LO or FF state
        V_phase = Vector{Matrix{eltype(Vdpdw[1])}}(undef, length(Vdpdw))
        V_sign  = Vector{Matrix{eltype(Vdpdw[1])}}(undef, length(Vdpdw))

        for i in eachindex(Vdpdw)
            V_phase[i] = similar(Vdpdw[i])
            V_sign[i]  = similar(Vdpdw[i])
        end

        @inbounds for i in 1:N
            r = locs[i]
            θ = q_p[1]*r[1] + q_p[2]*r[2]

            phase = if symmetry == "dLO"
                cis(θ) + cis(-θ)
            elseif symmetry == "dFF"
                cis(θ)
            end

            amp = abs(phase)
            sgn = iszero(amp) ? zero(phase) : phase / amp

            A  = Vdpdw[i]
            Bp = V_phase[i]
            Bs = V_sign[i]

            @inbounds @simd for k in eachindex(A)
                a = A[k]
                Bp[k] = a * phase   # full phase for H_vpars
                Bs[k] = a * sgn     # unit phase for V
            end
        end

        append!(H_vpars, V_phase)

        if optimize.Δ_dlo && symmetry == "dLO"
            append!(V, V_sign)
        end
        if optimize.Δ_dff && symmetry == "dFF"
            append!(V, V_sign)
        end    
    end

    return nothing
end


@doc raw"""

    add_spin_order!( order::S, 
                     optimize::NamedTuple, 
                     H_vpars::Vector{Matrix{T}}, 
                     V::Vector{Matrix{T}}, 
                     model_geometry::ModelGeometry,
                     locs::,
                     dims::I,
                     N::I, 
                     pht::Bool )::Nothing

Adds a spin ordering term to the auxiliary Hamiltonian along with its perturbative operator. 

- `order::String`: type of spin order: "spin-x", "spin-z", or "site-dependent"
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Any}`: vector of variational Hamiltonian matrices.
- `V::Vector{Any}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `locs::`: preallocated 
- `dims::I`: dimensions of the lattice. 
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_spin_order!(
    order::S,
    optimize::NamedTuple,
    H_vpars::Vector{Matrix{T}},
    V::Vector{Matrix{T}},
    model_geometry::ModelGeometry,
    locs,
    dims::I,
    N::I,
    pht::Bool
) where {S<:AbstractString, I<:Integer, T<:Number}
    twoN = 2 * N

    if order == "spin-z"
        afm_vec = Vector{Int8}(undef, twoN)

        if dims == 1
            @inbounds for s in 1:N
                ix = locs[s]
                v = isodd(ix) ? Int8(-1) : Int8(1)

                afm_vec[s]     = v
                afm_vec[s + N] = pht ? -v : v
            end
        else
            @inbounds for s in 1:N
                (ix, iy) = locs[s]
                v = isodd(ix + iy) ? Int8(-1) : Int8(1)

                afm_vec[s]     = v
                afm_vec[s + N] = pht ? -v : v
            end
        end

        V_afm = LinearAlgebra.Diagonal(afm_vec)
        push!(H_vpars, V_afm)

        if optimize.Δ_sz
            push!(V, V_afm)
        end
    elseif order == "spin-x"
        row = Vector{Int}(undef, 2N)
        col = Vector{Int}(undef, 2N)
        val = Vector{Float64}(undef, 2N)

        @inbounds for s in 1:N
            row[2s-1] = s
            col[2s-1] = s + N
            val[2s-1] = 1.0 # apply the 0.5 when the parameters are applied

            row[2s]   = s + N
            col[2s]   = s
            val[2s]   = 1.0
        end
        H_sx =  SparseArrays.sparse(row, col, val, twoN, twoN)

        push!(H_vpars, H_sx)                   
        if optimize.Δ_sx
            push!(V, H_sx)
        end
    elseif order == "site-dependent"
        L = model_geometry.lattice.L
        nshifts = L[1]

        # --- precompute staggering signs ---
        stag = Vector{Int8}(undef, twoN)
        if dims == 1
            @inbounds for s in 1:twoN
                ix = locs[s]
                stag[s] = isodd(ix) ? Int8(-1) : Int8(1)
            end
        else
            @inbounds for s in 1:twoN
                (ix, iy) = locs[s]
                stag[s] = isodd(ix + iy) ? Int8(-1) : Int8(1)
            end
        end

        # reusable buffer
        ssd_vec = Vector{Int8}(undef, twoN)
        fill!(ssd_vec, 0)

        for shift in 0:(nshifts - 1)

            # fill pattern + signs in one pass
            @inbounds begin
                k = 0
                while true
                    idx = shift + 1 + k * L[1]
                    idx > twoN && break

                    v = stag[idx]
                    if pht && idx > N
                        v = -v
                    end

                    ssd_vec[idx] = v
                    k += 1
                end
            end

            V_ssd = LinearAlgebra.Diagonal(ssd_vec)
            push!(H_vpars, V_ssd)
            if optimize.Δ_ssd
                push!(V, V_ssd)
            end

            # clear only touched entries
            @inbounds begin
                k = 0
                while true
                    idx = shift + 1 + k * L[1]
                    idx > twoN && break
                    ssd_vec[idx] = 0
                    k += 1
                end
            end
        end
    end

    return nothing
end


@doc raw"""

    add_charge_order!( order::S, 
                       optimize::NamedTuple, 
                       H_vpars::Vector{Matrix{T}}, 
                       V::Vector{Matrix{T}}, 
                       model_geometry::ModelGeometry,
                       locs::,
                       dims::I,
                       N::I, 
                       pht::Bool ) where {S<:AbstractString, I<:Integer, T<:Number}

Adds a charge ordering term to the auxiliary Hamiltonian along with its perturbative operator.

- `order::S`: type of spin order: "density wave" or "site-dependent"
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Matrix{T}}`: vector of variational Hamiltonian matrices.
- `V::Vector{Matrix{T}}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `locs::`:
- `dims::I`: dimensions of the lattice. 
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_charge_order!(
    order::S,
    optimize::NamedTuple,
    H_vpars::Vector{Matrix{T}},
    V::Vector{Matrix{T}},
    model_geometry::ModelGeometry,
    locs,
    dims::I,
    N::I,
    pht::Bool
) where {S<:AbstractString, I<:Integer, T<:Number}
    twoN = 2 * N

    if order == "density wave"
        cdw_vec = Vector{Int8}(undef, twoN)

        if dims == 1
            if pht
                @inbounds for s in 1:twoN
                    ix = locs[s]
                    cdw_vec[s] = isodd(ix) ? Int8(-1) : Int8(1)
                end
            else
                @inbounds for s in 1:N
                    ix = locs[s]
                    v = isodd(ix) ? Int8(-1) : Int8(1)
                    cdw_vec[s]     = v
                    cdw_vec[s + N] = -v
                end
            end
        else
            if pht
                @inbounds for s in 1:twoN
                    (ix, iy) = locs[s]
                    cdw_vec[s] = isodd(ix + iy) ? Int8(-1) : Int8(1)
                end
            else
                @inbounds for s in 1:N
                    (ix, iy) = locs[s]
                    v = isodd(ix + iy) ? Int8(-1) : Int8(1)
                    cdw_vec[s]     = v
                    cdw_vec[s + N] = -v
                end
            end
        end

        V_cdw = LinearAlgebra.Diagonal(cdw_vec)
        push!(H_vpars, V_cdw)
        if optimize.Δ_cdw
            push!(V, V_cdw)
        end
    elseif order == "site-dependent"
        L = model_geometry.lattice.L
        nshifts = L[1]

        for shift in 0:(nshifts - 1)

            csd_vec = Vector{Int8}(undef, twoN)
            fill!(csd_vec, 0)

            # stride-based index fill
            @inbounds for idx in (1 + shift):L[1]:twoN
                csd_vec[idx] = 1
            end

            if dims == 1
                if pht
                    @inbounds for s in 1:twoN
                        ix = locs[s]
                        stag = isodd(ix) ? Int8(-1) : Int8(1)
                        csd_vec[s] *= (s > N ? -stag : stag)
                    end
                else
                    @inbounds for s in 1:twoN
                        ix = locs[s]
                        csd_vec[s] *= isodd(ix) ? Int8(-1) : Int8(1)
                    end
                end
            else
                if pht
                    @inbounds for s in 1:twoN
                        (ix, iy) = locs[s]
                        stag = isodd(ix + iy) ? Int8(-1) : Int8(1)
                        csd_vec[s] *= (s > N ? -stag : stag)
                    end
                else
                    @inbounds for s in 1:twoN
                        (ix, iy) = locs[s]
                        csd_vec[s] *= isodd(ix + iy) ? Int8(-1) : Int8(1)
                    end
                end
            end

            V_csd = LinearAlgebra.Diagonal(csd_vec)
            push!(H_vpars, V_csd)
            if optimize.Δ_csd
                push!(V, V_csd)
            end
        end
    end

    return nothing
end


@doc raw"""

    add_chemical_potential!( optimize::NamedTuple, 
                             H_vpars::Vector{Matrix{T}}, 
                             V::Vector{Matrix{T}}, 
                             N::I, 
                             pht::Bool ) where {T<:Number, I<:Integer}

Adds a chemical potential term to the auxiliary Hamiltonian.

- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Matrix{T}}`: vector of variational Hamiltonian matrices.
- `V::Vector{Matrix{T}}`: vector of variational operators.
- `N::I`: number of sites in the lattice. 
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_chemical_potential!(
    optimize::NamedTuple, 
    H_vpars::Vector{Matrix{T}}, 
    V::Vector{Matrix{T}}, 
    N::I, 
    pht::Bool
) where {T<:Number, I<:Integer}
    μ_vec = fill(-1,2*N)
    if pht
        # account for minus sign 
        μ_vec_neg = copy(μ_vec)
        μ_vec_neg[N+1:2*N] .= -μ_vec_neg[N+1:2*N]

        # add chemical potential matrix
        H_μ = zeros(Number, 2*N, 2*N)
        V_μ_neg = LinearAlgebra.Diagonal(μ_vec_neg)
        H_μ += V_μ_neg
        push!(H_vpars, H_μ)

        # if μ is being optimized, save Vμ matrix
        if optimize.μ
            push!(V, V_μ_neg)
        end
    else
        # add chemical potential matrix
        H_μ = zeros(Number, 2*N, 2*N)
        V_μ = LinearAlgebra.Diagonal(μ_vec)
        H_μ += V_μ
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
                  Np::I,
                  pht::Bool ) where {E<:AbstractFloat, I<:Integer}

Checks whether a mean-field energy configuration is open shell.

- `ε::Vector{E}`: vector of mean-field energies.
- `Np::I`: total number of particles in the system.
- `pht::Bool`: whether the model is particle-hole transformed.

"""
function is_openshell(
    ε::Vector{E}, 
    Np::I,
    pht::Bool
) where {E<:AbstractFloat, I<:Integer}
    if pht
        return ε[Np] - ε[Np - 1] < 1e-4 
    else
        return ε[Np + 1] - ε[Np] < 1e-4 
    end 
end


@doc raw"""

    get_variational_matrices( V::Vector{Matrix{T1}}, 
                              U_aux::Matrix{T2}, 
                              ε::Vector{E}, 
                              Np::I
                              N::I ) where {T1<:Number, T2<:Number, E<:AbstractFloat, I<:Integer}
    
Returns a set of variational parameter matrices ``A_k`` constructed from each variational operator ``V_k`` by computing 
```math
Q_k = \frac{(U^{\dagger}V_{k}U)_{\eta\nu}}{(\varepsilon_{\eta} - \varepsilon_{\nu})},
```
where ``\eta > N_p`` and ``\nu \leq N_p``.

- `V::Vector{Matrix{T1}`: vector of variational operators. 
- `U_aux::Matrix{T2}`: matrix which diagonalizes the auxiliary Hamiltonian.
- `ε::Vector{E}`: initial energies.
- `Np::I`: number of particles in the system. 
- `N::I`: number of lattice sites.

"""
function get_variational_matrices(
    ptmask::Matrix{E},
    V::Vector{Matrix{T1}}, 
    U_aux::Matrix{T2}, 
    ε::Vector{E}, 
    Np::I,
    N::I
) where {T1<:Number, T2<:Number, E<:AbstractFloat, I<:Integer}
    # ### OLD ALGORTHIM ###
    # # populate perturbation mask
    # for η in 1:2*N
    #     for ν in 1:2*N
    #         if η > Np && ν <= Np
    #             ptmask[η, ν] = 1.0 / (ε[ν] - ε[η])
    #         end
    #     end
    # end
    #
    # # compute A matrices
    # As = Vector{Matrix{<:Number}}()
    # for v in V
    #     A = U_aux * ((U_aux' * v * U_aux) .* ptmask) * U_aux'
    #     push!(As, A)
    # end

    ### NEW ALGORITHM ###
    # populate perturbation theory mask
    @inbounds for ν in 1:Np
        εν = ε[ν]
        for η in (Np+1):(2*N)
            ptmask[η, ν] = inv(εν - ε[η])
        end
    end

    Nt  = size(U_aux, 1)
    Np1 = Np + 1

    Uocc  = @view U_aux[:, 1:Np]
    Uvirt = @view U_aux[:, Np1:Nt]

    As = Vector{Matrix{eltype(U_aux)}}(undef, length(V))

    T   = Matrix{eltype(U_aux)}(undef, Nt, Np)
    Q   = Matrix{eltype(U_aux)}(undef, Nt-Np, Np)
    Tmp = Matrix{eltype(U_aux)}(undef, Nt, Np)
    A   = Matrix{eltype(U_aux)}(undef, Nt, Nt)

    # # version 1
    # for (k, v) in pairs(V)
    #     mul!(T, v, Uocc)          # T = v * Uocc
    #     mul!(Q, Uvirt', T)        # Q = Uvirt' * T

    #     @inbounds for i in axes(Q,1), j in axes(Q,2)
    #         Q[i,j] /= (ε[j] - ε[i+Np])
    #     end

    #     mul!(Tmp, Uvirt, Q)       # Tmp = Uvirt * Q
    #     mul!(A, Tmp, Uocc')       # A = Tmp * Uocc'

    #     As[k] = copy(A)
    # end

    for (k, v) in pairs(V)
        mul!(T, v, Uocc)              # T = v * Uocc
        mul!(Q, Uvirt', T)            # Q = Uvirt' * T

        @inbounds for i in axes(Q,1), j in axes(Q,2)
            Q[i,j] /= (ε[j] - ε[i+Np])
        end

        mul!(Tmp, Uvirt, Q)           # Tmp = Uvirt * Q
        mul!(A, Tmp, Uocc')           # A = Tmp * Uocc'

        @inbounds @simd for i in eachindex(A)
            A[i] += conj(A[i])        # in-place symmetrization
        end

        As[k] = copy(A)               
    end

    return As
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


















