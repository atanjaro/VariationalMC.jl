@doc raw"""

    build_auxiliary_hamiltonian( tight_binding_model::TightBindingModel, 
                                 determinantal_parameters::DeterminantalParameters, 
                                 pht::Bool ) 

Constructs a complete Hamiltonian matrix by combining the non-interacting matrix with
matrices of variational terms.

- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model.
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hol transformed.

"""
function build_auxiliary_hamiltonian(
    tight_binding_model::TightBindingModel, 
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    model_geometry::ModelGeometry, 
    pht::Bool
)
    # hopping matrix
    H_tb = build_tight_binding_hamiltonian(
        tight_binding_model, 
        model_geometry, 
        pht
    )

    # variational matrices and operators
    H_var, V = build_variational_hamiltonian(
        determinantal_parameters, 
        optimize, 
        pht
    )

    return H_tb + H_var, V
end


@doc raw"""

    build_tight_binding_hamiltonian( tight_binding_model::TightBindingModel,
                                     model_geometry::ModelGeometry,
                                     pht::Bool ) 

Constructs a 2N by 2N Hamiltonian matrix where N is the number of lattice sites, 
given tight binding parameters t, and t'.

- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model. 
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function build_tight_binding_hamiltonian(
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry, 
    pht::Bool
)
    # number of sites
    N = model_geometry.unit_cell.n*model_geometry.lattice.N 

    # generate neighbor table
    nbr0 = build_neighbor_table(
        bonds[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )

    # initialize matrices
    H_t₀ = zeros(Complex, 2*N, 2*N)
    H_t₁ = zeros(Complex, 2*N, 2*N)

    # hopping parameters
    t₀ = tight_binding_model.t₀ 
    t₁ = tight_binding_model.t₁ 

    debug && println("Hamiltonian::build_tight_binding_hamiltonian() : ")
    debug && println("building tight-binding hopping matrix")
    debug && println("hopping : t₀ = ", t₀)
    debug && println("hopping : t₁ = ", t₁)
    debug && println("particle-hole transformation : ", pht)

    if pht == true
        # add nearest-neighbor hopping
        if Lx == 2 && Ly == 2 
            for (i,j) in eachcol(nbr0)
                H_t₀[i,j] += -t₀
            end
            for (i,j) in eachcol(nbr0 .+ N)    
                H_t₀[i,j] += t₀
            end
        # special case for Lx = 2 
        elseif Lx == 2 && Ly > Lx
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Ly)])
                H_t₀[i,j] += -t₀
                H_t₀[j,i] += -t₀
            end
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Ly)] .+ N)
                H_t₀[i,j] += t₀
                H_t₀[j,i] += t₀
            end 
        # special case for Ly = 2 
        elseif Ly == 2 && Lx > Ly
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Lx)])
                H_t₀[i,j] += -t₀
                H_t₀[j,i] += -t₀
            end
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Lx)] .+ N)
                H_t₀[i,j] += t₀
                H_t₀[j,i] += t₀
            end 
        else
            for (i,j) in eachcol(nbr0)
                H_t₀[i,j] += -t₀
                if model_geometry.lattice.N > 2
                    H_t₀[j,i] += -t₀
                else
                end
            end
            for (i,j) in eachcol(nbr0 .+ N)    
                H_t₀[i,j] += t₀
                if model_geometry.lattice.N > 2
                    H_t₀[j,i] += t₀
                else
                end
            end
        end

        # add next-nearest neighbor hopping
        nbr1 = build_neighbor_table(
            bonds[2],
            model_geometry.unit_cell,
            model_geometry.lattice
        )
        if Lx == 2 && Ly == 2
            for (i,j) in eachcol(nbr1)
                H_t₁[i,j] += t₁/2
            end
            for (i,j) in eachcol(nbr1 .+ N)    
                H_t₁[i,j] += -t₁/2
            end
        else
            for (i,j) in eachcol(nbr1)
                H_t₁[i,j] += t₁
                H_t₁[j,i] += t₁
            end
            for (i,j) in eachcol(nbr1 .+ N)    
                H_t₁[i,j] += -t₁
                H_t₁[j,i] += -t₁
            end
        end
    else
        # nearest neighbor hopping
        if Lx == 2 && Ly == 2 
            for (i,j) in eachcol(nbr0)
                H_t₀[i,j] += -t₀
            end
            for (i,j) in eachcol(nbr0 .+ N)    
                H_t₀[i,j] += -t₀
            end
        # special case for Lx = 2 
        elseif Lx == 2 && Ly > Lx
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Ly)])
                H_t₀[i,j] += -t₀
                H_t₀[j,i] += -t₀
            end
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Ly)] .+ N)
                H_t₀[i,j] += -t₀
                H_t₀[j,i] += -t₀
            end 
        # special case for Ly = 2 
        elseif Ly == 2 && Lx > Ly
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Lx)])
                H_t₀[i,j] += -t₀
                H_t₀[j,i] += -t₀
            end
            for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Lx)] .+ N)
                H_t₀[i,j] += -t₀
                H_t₀[j,i] += -t₀
            end  
        else
            for (i,j) in eachcol(nbr0)
                H_t₀[i,j] += -t₀
                if model_geometry.lattice.N > 2
                    H_t₀[j,i] += -t₀
                else
                end
            end
            for (i,j) in eachcol(nbr0 .+ N)    
                H_t₀[i,j] += -t₀
                if model_geometry.lattice.N > 2
                    H_t₀[j,i] += -t₀
                else
                end
            end
        end

        # add next-nearest neighbor hopping 
        nbr1 = build_neighbor_table(
            bonds[2],
            model_geometry.unit_cell,
            model_geometry.lattice
        )
        if Lx == 2 && Ly ==2 
            for (i,j) in eachcol(nbr1)
                H_t₁[i,j] += t₁/2
            end
            for (i,j) in eachcol(nbr1 .+ N)    
                H_t₁[i,j] += t₁/2
            end
        else
            for (i,j) in eachcol(nbr1)
                H_t₁[i,j] += t₁
                H_t₁[j,i] += t₁
            end
            for (i,j) in eachcol(nbr1 .+ N)    
                H_t₁[i,j] += t₁
                H_t₁[j,i] += t₁
            end
        end
    end

    return H_t₀ + H_t₁
end


@doc raw"""

    build_variational_hamiltonian( determinantal_parameters::DeterminantalParameters, 
                                   optimize::NamedTuple, 
                                   pht::Bool ) 

Constructs 2N by 2N matrices to be added to the non-interacting tight binding Hamiltonian for each variational parameter, 
where N is the number of lattice sites. Returns a vector of the sum of matrices and a vector of individual matrix terms.

- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `pht::Bool`: whether model is particle-hole transformed. 

"""
function build_variational_hamiltonian(
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    pht::Bool
)
    # dimensions
    dims = size(model_geometry.lattice.L)[1]
   
    # initialize Hamiltonian and operator matrices
    H_vpars = []
    V = []
    
    # add chemical potential term
    add_chemical_potential!(
        determinantal_parameters, 
        optimize, 
        H_vpars, 
        V, 
        model_geometry, 
        pht
    )

    debug && println("Hamiltonian::build_variational_hamiltonian() : ")
    debug && println("adding chemical potential matrix")
    debug && println("initial μ = ", determinantal_parameters.det_pars.μ)
    if optimize.μ
        debug && println("optimize = true")
    else
        debug && println("optimize = false")
    end

    # add s-wave term
    if pht == true
        add_pairing_symmetry!(
            "s", 
            determinantal_parameters, 
            optimize, 
            H_vpars, 
            V, 
            model_geometry, 
            pht
        )

        debug && println("Hamiltonian::build_variational_hamiltonian() : ")
        debug && println("adding s-wave pairing matrix")
        debug && println("initial Δ_0 = ", determinantal_parameters.det_pars.Δ_0)
        if optimize.Δ_0
            debug && println("optimize = true")
        else
            debug && println("optimize = false")
        end

        # add d-wave pairing 
        if dims > 1
            add_pairing_symmetry!(
                "d", 
                determinantal_parameters, 
                optimize, 
                H_vpars, 
                V, 
                model_geometry, 
                pht
            )

            debug && println("Hamiltonian::build_variational_hamiltonian() : ")
            debug && println("adding d-wave pairing matrix")
            debug && println("initial Δ_d = ", determinantal_parameters.det_pars.Δ_d)
            if optimize.Δ_d
                debug && println("optimize = true")
            else
                debug && println("optimize = false")
            end
        end
    end

    # add in-plane magnetization term
    # add_spin_order!(
    #     "spin-x"
    #     determinantal_parameters,
    #     optimize,
    #     H_vpars,
    #     V,
    #     model_geometry,
    #     pht
    # )

    # add antiferromagnetic (Neél) term
    add_spin_order!(
        "spin-z", 
        determinantal_parameters, 
        optimize, 
        H_vpars, 
        V, 
        model_geometry, 
        pht
    )

    debug && println("Hamiltonian::build_variational_hamiltonian() : ")
    debug && println("adding spin-z matrix")
    debug && println("initial Δ_afm = ", determinantal_parameters.det_pars.Δ_afm)
    if optimize.Δ_afm
        debug && println("optimize = true")
    else
        debug && println("optimize = false")
    end

    # add charge-density-wave term
    add_charge_order!(
        "density wave", 
        determinantal_parameters, 
        optimize, 
        H_vpars, 
        V, 
        model_geometry, 
        pht
    )

    debug && println("Hamiltonian::build_variational_hamiltonian() : ")
    debug && println("adding charge density wave matrix")
    debug && println("initial Δ_cdw = ", determinantal_parameters.det_pars.Δ_cdw)
    if optimize.Δ_cdw
        debug && println("optimize = true")
    else
        debug && println("optimize = false")
    end

    # add site-dependent charge term
    if dims > 1
        add_charge_order!(
            "site-dependent", 
            determinantal_parameters, 
            optimize, 
            H_vpars,
            V, 
            model_geometry, 
            pht
        )

        debug && println("Hamiltonian::build_variational_hamiltonian() : ")
        debug && println("adding site-dependent charge matrix")
        debug && println("initial Δ_sdc = ", determinantal_parameters.det_pars.Δ_sdc)
        if optimize.Δ_sdc
            debug && println("optimize = true")
        else
            debug && println("optimize = false")
        end
    end

    # add site-dependent spin term
    if dims > 1
        add_spin_order!(
            "site-dependent", 
            determinantal_parameters, 
            optimize, 
            H_vpars, 
            V, 
            model_geometry, 
            pht
        )

        debug && println("Hamiltonian::build_variational_hamiltonian() : ")
        debug && println("adding site-dependent spin matrix")
        debug && println("initial Δ_sds = ", determinantal_parameters.det_pars.Δ_sds)
        if optimize.Δ_sds
            debug && println("optimize = true")
        else
            debug && println("optimize = false")
        end
    end

    @assert length(H_vpars) == determinantal_parameters.num_det_pars
    @assert length(V) == determinantal_parameters.num_det_opts

    return sum(H_vpars), V
end


@doc raw"""

    add_pairing_symmetry!( symmetry::String, 
                           determinantal_parameters::DeterminantalParameters, 
                           optimize::NamedTuple, 
                           H_vpars::AbstractMatrix{<:Complex}, 
                           V, 
                           model_geometry::ModelGeometry, 
                           pht::Bool )::Nothing

Adds a pairing symmetry term to the auxiliary Hamiltonian. 

- `symmetry::String`: type of pairing symmetry: "s" or "d"
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Any}`: vector of variational Hamiltonian matrices.
- `V::Vector{Any}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities. 
- `pht::Bool`: whether model is particle-hole transformed.

"""
function add_pairing_symmetry!(
    symmetry::String, 
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    H_vpars::Vector{Any}, 
    V::Vector{Any}, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Nothing
    # lattice dimensions
    dims = size(model_geometry.lattice.L)[1]

    # lattice sites
    N = model_geometry.lattice.N

    # add s-wave pairing 
    if symmetry == "s"
        @assert pht == true

        # s-wave parameter
        Δ_0 = determinantal_parameters.det_pars.Δ_0

        # add s-wave symmetry
        H_s = zeros(Complex, 2*N, 2*N) 
        V_s = copy(H_s)
        for i in 0:(2 * N - 1)
            V_s[i + 1, get_linked_spindex(i, N) + 1] = 1.0
        end

        # add s-wave matrix
        H_s += Δ_0 * V_s
        push!(H_vpars,H_s)

        if optimize.Δ_0
            push!(V, V_s)
        end
    elseif symmetry == "d"
        @assert pht == true
        @assert dims > 1

        # d-wave parameter
        Δ_d = determinantal_parameters.det_pars.Δ_d

        # create neighbor table
        nbr_table = build_neighbor_table(
            bonds[1],
            model_geometry.unit_cell,
            model_geometry.lattice
        )

        # maps neighbor table to dictionary of bonds and neighbors                                
        nbrs = map_neighbor_table(nbr_table)

        # Predefine the displacement-to-sign map
        disp_sign_map = Dict([1,0] => 1, [0,1] => -1, [-1,0] => 1, [0,-1] => -1)

        # initial variational operator matrix
        H_d = zeros(Complex, 2*N, 2*N)
        V_dwave = zeros(AbstractFloat, 2*N, 2*N)

        # add d-wave symmetry
        for i in 1:N
            # get all neighbors of site i
            nn = nbrs[i][2]

            # loop over neighbors
            for j in nn
                # Find the displacement between site i and one of its neighbors j
                disp = sites_to_displacement(i, j, unit_cell, lattice)

                # Lookup sign of Δd
                dsgn = get(disp_sign_map, disp, 0)  # Default to 0 if no match

                if dsgn != 0
                    # println(dsgn > 0 ? "+1" : "-1")

                    # Store spin-down indices
                    idn_idx = get_spindices_from_index(i, model_geometry)[2]
                    jdn_idx = get_spindices_from_index(j, model_geometry)[2]

                    # Add elements to variational operator
                    V_dwave[i, jdn_idx] = dsgn
                    V_dwave[j, idn_idx] = dsgn
                    V_dwave[jdn_idx, i] = dsgn
                    V_dwave[idn_idx, j] = dsgn
                end
            end
        end

        # add d-wave matrix
        H_d += Δ_d * V_dwave
        push!(H_vpars,H_d)

        # if Δ_d is being optimized, store Vdwave matrix
        if optimize.Δ_d
            push!(V, V_dwave)
        end
    end

    return nothing
end


@doc raw"""

    add_spin_order!( order::String, 
                     determinantal_parameters::DeterminantalParameters, 
                     optimize::NamedTuple, 
                     H_vpars::Vector{Any}, 
                     V::Vector{Any}, 
                     model_geometry::ModelGeometry, 
                     pht::Bool )::Nothing

Adds a spin order term to the auxiliary Hamiltonian. 

- `order::String`: type of spin order: "spin-z" or "site-dependent"
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Any}`: vector of variational Hamiltonian matrices.
- `V::Vector{Any}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_spin_order!(
    order::String, 
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    H_vpars::Vector{Any}, 
    V::Vector{Any}, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Nothing
    # lattice sites
    N = model_geometry.lattice.N

    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    if order == "spin-z"
        afm_vec = fill(1,2*N)

        # antiferromagnetic parameter
        Δ_afm = determinantal_parameters.det_pars.Δ_afm

        if pht
            # stagger
            for s in 1:2*N
                # get proper site index
                idx = get_index_from_spindex(s, model_geometry)

                # 1D
                if dims == 1
                    # get site coordinates
                    ix = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1][1]

                    # apply phase
                    afm_vec[s] *= (-1)^(ix)
                # 2D
                elseif dims == 2
                    # get site coordinates
                    (ix, iy) = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1]

                    # apply phase
                    afm_vec[s] *= (-1)^(ix + iy)
                end
            end

            # add afm matrix
            H_afm = zeros(Complex, 2*N, 2*N)
            V_afm = LinearAlgebra.Diagonal(afm_vec)
            H_afm += Δ_afm * V_afm
            push!(H_vpars, H_afm)

            # if Δ_afm is being optimized, store Vafm matrix
            if optimize.Δ_afm
                push!(V, V_afm)
            end
        else
            # additional sign flip in spin down sector
            afm_vec_neg = copy(afm_vec)
            afm_vec_neg[N+1:2*N] .= -afm_vec_neg[N+1:2*N]

            # stagger
            for s in 1:2*N
                # get proper site index
                idx = get_index_from_spindex(s, model_geometry)

                # 1D
                if dims == 1
                    # get site coordinates
                    ix = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1][1]

                    # apply phase
                    afm_vec_neg[s] *= (-1)^(ix)
                # 2D
                elseif dims == 2
                    # get site coordinates
                    (ix, iy) = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1]

                    # apply phase
                    afm_vec_neg[s] *= (-1)^(ix + iy)
                end
            end

            # add afm matrix
            H_afm = zeros(Complex, 2*N, 2*N)
            V_afm_neg = LinearAlgebra.Diagonal(afm_vec_neg)
            H_afm += Δ_afm * V_afm_neg
            push!(H_vpars, H_afm)

            # if Δ_afm is being optimized, store Vafm matrix
            if optimize.Δ_afm
                push!(V, V_afm_neg)
            end
        end
    # elseif order == "spin-x"
    #     # spin-x parameter
    #     Δ_x = determinantal_parameters.det_pars.Δ_x

    #     # spin-x matrix (off-diagonal in spin)
    #     H_sx = zeros(Complex, 2*N, 2*N)

    #     for s in 1:N
    #         up_idx = s
    #         dn_idx = s + N

    #         # Add spin-flip terms: S^x = (1/2)(|↑⟩⟨↓| + |↓⟩⟨↑|)
    #         H_sx[up_idx, dn_idx] += Δ_x / 2
    #         H_sx[dn_idx, up_idx] += Δ_x / 2
    #     end

    #     # Add to variational Hamiltonian
    #     push!(H_vpars, H_sx)

    #     # If Δ_x is being optimized, store the matrix
    #     if optimize.Δ_x
    #         push!(V, H_sx)
    #     end
    elseif order == "site-dependent"
        # lattice dimensions
        L = model_geometry.lattice.L

        Δ_sds = determinantal_parameters.det_pars.Δ_sds
        sds_vectors = []
        if pht
            for shift in 0:(L[1]-1)
                vec = zeros(Int, 2 * N)
                for i in 1:2*L[1]
                    idx = (i - 1) * L[1] + 1 + shift
                    if idx <= 2 * N
                        vec[idx] = 1
                    end
                end
                # Flip the sign of the last N elements
                vec[N+1:end] .*= -1
                push!(sds_vectors, vec)
            end
        else
            for shift in 0:(L[1]-1)
                vec = zeros(Int, 2 * N)
                for i in 1:2*L[1]
                    idx = (i-1)*L[1] + 1 + shift
                    if idx <= 2 * N  
                        vec[idx] = 1
                    end
                end
                push!(sds_vectors, vec)  
            end
        end

        for (i, sds_vec) in enumerate(sds_vectors)
            H_sds = zeros(AbstractFloat, 2*N, 2*N)
            sds_vec_neg = copy(sds_vec)
            sds_vec_neg[N+1:2*N] .= -sds_vec_neg[N+1:2*N]
        
            for s in 1:2*N
                idx = get_index_from_spindex(s, model_geometry)
                loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
                sds_vec_neg[s] *= (-1)^(loc[1][1] + loc[1][2])
            end
        
            V_sds_neg = LinearAlgebra.Diagonal(sds_vec_neg)
            H_sds += Δ_sds[i] * V_sds_neg
            push!(H_vpars, H_sds)
        
            if optimize.Δ_sds
                push!(V, V_sds_neg)
            end
        end
    end

    return nothing
end


@doc raw"""

    add_charge_order!( order::String, 
                       determinantal_parameters::DeterminantalParameters, 
                       optimize::NamedTuple, 
                       H_vpars::Vector{Any}, 
                       V::Vector{Any}, 
                       model_geometry::ModelGeometry, 
                       pht::Bool )::Nothing

Adds a charge order term to the auxiliary Hamiltonian.

- `order::String`: type of spin order: "density wave" or "site-dependent"
- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Any}`: vector of variational Hamiltonian matrices.
- `V::Vector{Any}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_charge_order!(
    order::String, 
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    H_vpars::Vector{Any}, 
    V::Vector{Any}, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Nothing
    # lattice sites
    N = model_geometry.lattice.N

    # dimensions
    dims = size(model_geometry.lattice.L)[1]

    if order == "density wave"
        # charge density wave parameter
        Δ_cdw = determinantal_parameters.det_pars.Δ_cdw

        # diagonal vector
        cdw_vec = fill(1,2*N)

        # account for particle-hole transformation
        if pht
            # sign flip in the spin-down sector
            cdw_vec_neg = copy(cdw_vec)
            cdw_vec_neg[N+1:2*N] .= -cdw_vec_neg[N+1:2*N]

            # stagger
            for s in 1:2*N
                # get proper site index
                idx = get_index_from_spindex(s, model_geometry)

                # 1D
                if dims == 1
                    # get site coordinates
                    ix = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1][1]

                    # apply phase
                    cdw_vec_neg[s] *= (-1)^(ix)
                # 2D
                elseif dims == 2
                    # get site coordinates
                    (ix, iy) = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1]

                    # apply phase
                    cdw_vec_neg[s] *= (-1)^(ix + iy)
                end
            end

            # add cdw matrix
            H_cdw = zeros(Complex, 2*N, 2*N) 
            V_cdw = LinearAlgebra.Diagonal(cdw_vec_neg)
            H_cdw += Δ_cdw * V_cdw
            push!(H_vpars, H_cdw)

            # if Δ_cdw is being optimized, save Vcdw matrix
            if optimize.Δ_cdw
                push!(V, V_cdw)
            end
        else
            # stagger
            for s in 1:2*N
                # get proper site index
                idx = get_index_from_spindex(s, model_geometry)

                # 1D
                if dims == 1
                    # get site coordinates
                    ix = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1][1]

                    # apply phase
                    cdw_vec[s] *= (-1)^(ix)
                # 2D
                elseif dims == 2
                    # get site coordinates
                    (ix, iy) = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)[1]

                    # apply phase
                    cdw_vec[s] *= (-1)^(ix + iy)
                end
            end

            # add cdw matrix
            H_cdw = zeros(Complex, 2*N, 2*N) 
            V_cdw = LinearAlgebra.Diagonal(cdw_vec)
            H_cdw += Δ_cdw * V_cdw
            push!(H_vpars, H_cdw)

            # if Δ_cdw is being optimized, save Vcdw matrix
            if optimize.Δ_cdw
                push!(V, V_cdw)
            end
        end
    elseif order == "site-dependent"
        # lattice dimensions
        L = model_geometry.lattice.L
        Δ_sdc = determinantal_parameters.det_pars.Δ_sdc
        sdc_vectors = []

        if pht
            for shift in 0:(L[1]-1)
                vec = zeros(Int, 2 * N)
                for i in 1:2*L[1]
                    idx = (i - 1) * L[1] + 1 + shift
                    if idx <= 2 * N
                        vec[idx] = 1
                    end
                end
                vec[N+1:end] .*= -1  # Flip sign of last N elements
                push!(sdc_vectors, vec)
            end
        else
            for shift in 0:(L[1]-1)
                vec = zeros(Int, 2 * N)  
                for i in 1:2*L[1]
                    idx = (i-1)*L[1] + 1 + shift
                    if idx <= 2 * N  
                        vec[idx] = 1
                    end
                end
                push!(sdc_vectors, vec)  
            end
        end

        for (i, sdc_vec) in enumerate(sdc_vectors)
            H_sdc = zeros(AbstractFloat, 2*N, 2*N)
            V_sdc = LinearAlgebra.Diagonal(sdc_vec)
            H_sdc += Δ_sdc[i] .* V_sdc
            push!(H_vpars, H_sdc)

            # if Δ_sdc is being optimized, store Vsdc matrix
            if optimize.Δ_sdc
                push!(V, V_sdc)
            end
        end
    end

    return nothing
end


@doc raw"""

    add_chemical_potential!( determinantal_parameters::DeterminantalParameters, 
                             optimize::NamedTuple, 
                             H_vpars::Vector{Any}, 
                             V::Vector{Any}, 
                             model_geometry::ModelGeometry, 
                             pht::Bool )::Nothing

Adds a chemical potential term to the auxiliary Hamiltonian.

- `determinantal_parameters::DeterminantalParameters`: set of determinantal variational parameters.
- `optimize::NamedTuple`: field of optimization flags.
- `H_vpars::Vector{Any}`: vector of variational Hamiltonian matrices.
- `V::Vector{Any}`: vector of variational operators.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model if particle-hole transformed.

"""
function add_chemical_potential!(
    determinantal_parameters::DeterminantalParameters, 
    optimize::NamedTuple, 
    H_vpars::Vector{Any}, 
    V::Vector{Any}, 
    model_geometry::ModelGeometry, 
    pht::Bool
)::Nothing
    # lattice sites
    N = model_geometry.lattice.N

    # initial chemical potential value
    μ = determinantal_parameters.det_pars.μ

    μ_vec = fill(-1,2*N)
    if pht
        # account for minus sign 
        μ_vec_neg = copy(μ_vec)
        μ_vec_neg[N+1:2*N] .= -μ_vec_neg[N+1:2*N]

        # add chemical potential matrix
        H_μ = zeros(Complex, 2*N, 2*N)
        V_μ_neg = LinearAlgebra.Diagonal(μ_vec_neg)
        H_μ += μ * V_μ_neg
        push!(H_vpars, H_μ)

        # if μ is being optimized, save Vμ matrix
        if optimize.μ
            push!(V, V_μ_neg)
        end
    else
        # add chemical potential matrix
        H_μ = zeros(Complex, 2*N, 2*N)
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

    diagonalize( H::Matrix{ComplexF64} )::Tuple{Vector{Float64}, Matrix{ComplexF64}}

Returns all eigenenergies and all eigenstates of the mean-field Hamiltonian, 
the latter being stored in the columns of a matrix. Convention: H(diag) = U⁺HU.

- `H::AbstractMatrix{<:Complex}`: auxiliary Hamiltonian.

"""
function diagonalize(
    H::AbstractMatrix{<:Complex}
)
    # check if Hamiltonian is Hermitian
    @assert ishermitian(H) == true
    F = eigen(H)    
    return F.values, F.vectors  
end


@doc raw"""

    is_openshell( ε::Vector{Float64},  
                  Np::Int )::Bool

Checks whether a energy configuration is open shell.

"""
function is_openshell(
    ε::Vector{Float64}, 
    Np::Int
)::Bool
    return ε[Np + 1] - ε[Np] < 0.0001
end


@doc raw"""

    get_variational_matrices( V::Vector{Any}, 
                              U_int::Matrix{ComplexF64}, 
                              ε::Vector{Float64}, 
                              model_geometry::ModelGeometry )::Vector{Any} 
    
Returns variational parameter matrices Aₖ from the corresponding Vₖ. Computes 
Qₖ = (U⁺VₖU)_(ην) / (ε_η - ε_ν), for η > Nₚ and ν ≤ Nₚ and is 0 otherwise (η and ν run from 1 to 2N).

- `V::Vector{Any}`: vector of variational operators. 
- `U_int::Matrix{ComplexF64}`: matrix which diagonalizes the auxiliary Hamiltonian.
- `ε::Vector{Float64}`: initial energies.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_variational_matrices(
    V::Vector{Any}, 
    U_int::Matrix{ComplexF64}, 
    ε::Vector{Float64}, 
    model_geometry::ModelGeometry
)::Vector{Any}
    # number of lattice sites
    N = model_geometry.unit_cell.n * model_geometry.lattice.N

    # define perturbation mask
    ptmask = zeros(Float64, 2*N, 2*N)
        
    for η in 1:2*N
        for ν in 1:2*N
            if η > Ne && ν <= Ne
                ptmask[η, ν] = 1.0 / (ε[ν] - ε[η])
            end
        end
    end

    int_A = []
    for v in V
        push!(int_A, U_int * ((U_int' * v * U_int) .* ptmask) * U_int')
    end

    return int_A
end


@doc raw"""

    get_tb_chem_pot( Ne::Int64, 
                     tight_binding_model::TightBindingModel, 
                     model_geometry::ModelGeometry )::Float64

For a tight-binding model that has not been particle-hole transformed, returns the  
chemical potential.

- `Ne::Int64`: total number of electrons.
- `tight_binding_model::TightBindingModel`: parameters for a non-interacting tight-binding model.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.

"""
function get_tb_chem_pot(
    Ne::Int64, 
    tight_binding_model::TightBindingModel, 
    model_geometry::ModelGeometry
)::Float64
    @assert pht == false

    # number of lattice sites
    N = model_geometry.lattice.N
    
    # preallocate matrices
    H_t₀ = zeros(Complex, 2*N, 2*N)
    H_t₁ = zeros(Complex, 2*N, 2*N)

    # hopping amplitudes
    t₀ = tight_binding_model.t₀
    t₁ = tight_binding_model.t₁

    # add nearest neighbor hopping
    nbr0 = build_neighbor_table(
        bonds[1],
        model_geometry.unit_cell,
        model_geometry.lattice
    )
    # special case for Lx, Ly = 2
    if Lx == 2 && Ly == 2 
        for (i,j) in eachcol(nbr0)
            H_t₀[i,j] += -t₀
        end
        for (i,j) in eachcol(nbr0 .+ N)    
            H_t₀[i,j] += -t₀
        end
    # special case for Lx = 2 
    elseif Lx == 2 && Ly > Lx
        for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Ly)])
            H_t₀[i,j] += -t₀
            H_t₀[j,i] += -t₀
        end
        for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Ly)] .+ N)
            H_t₀[i,j] += -t₀
            H_t₀[j,i] += -t₀
        end 
    # special case for Ly = 2 
    elseif Ly == 2 && Lx > Ly
        for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Lx)])
            H_t₀[i,j] += -t₀
            H_t₀[j,i] += -t₀
        end
        for (i,j) in eachcol(nbr0[:,1:(size(nbr0,2) - Lx)] .+ N)
            H_t₀[i,j] += -t₀
            H_t₀[j,i] += -t₀
        end  
    else
        for (i,j) in eachcol(nbr0)
            H_t₀[i,j] += -t₀
            if model_geometry.lattice.N > 2
                H_t₀[j,i] += -t₀
            else
            end
        end
        for (i,j) in eachcol(nbr0 .+ N)    
            H_t₀[i,j] += -t₀
            if model_geometry.lattice.N > 2
                H_t₀[j,i] += -t₀
            else
            end
        end
    end
    # add next nearest neighbor hopping
    nbr1 = build_neighbor_table(
        bonds[2],
        model_geometry.unit_cell,
        model_geometry.lattice
    )
    if Lx == 2 && Ly ==2 
        for (i,j) in eachcol(nbr1)
            H_t₁[i,j] += 0.5 * t₁
        end
        for (i,j) in eachcol(nbr1 .+ N)    
            H_t₁[i,j] += 0.5 * t₁
        end
    else
        for (i,j) in eachcol(nbr1)
            H_t₁[i,j] += t₁
            H_t₁[j,i] += t₁
        end
        for (i,j) in eachcol(nbr1 .+ N)    
            H_t₁[i,j] += t₁
            H_t₁[j,i] += t₁
        end
    end

    # full tight-binding Hamiltonian
    H_tb = H_t₀ + H_t₁

    # solve for eigenvalues
    ε_F, Uₑ = diagonalize(H_tb)

    # tight-binding chemical potential
    μ = 0.5 * (ε_F[Ne + 1] + ε_F[Ne])

    debug && println("Hamiltonian::get_tb_chem_pot() : ")
    debug && println("tight-binding chemical potential")
    debug && println("μ = ", μ)

    return μ
end


















