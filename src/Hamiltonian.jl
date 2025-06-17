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
    N = model_geometry.unit_cell.n * model_geometry.lattice.N 
    Lx = model_geometry.lattice.L[1]
    Ly = model_geometry.lattice.L[2]

    # generate nearest neighbor table
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
        @assert Lx > 2 && Ly > 2

        # add nearest-neighbor hopping
        for (i,j) in eachcol(nbr0)
            H_t₀[i,j] += -t₀
            if N > 2
                H_t₀[j,i] += -t₀
            end
        end
        for (i,j) in eachcol(nbr0 .+ N)    
            H_t₀[i,j] += t₀
            if N > 2
                H_t₀[j,i] += t₀
            end
        end

        if t₁ != 0
            # generate next-nearest neighbor table
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
                H_t₁[i,j] += -t₁
                H_t₁[j,i] += -t₁
            end
        end
    else
        @assert Lx > 2 && Ly > 2

        # add nearest neighbor hopping
        for (i,j) in eachcol(nbr0)
            H_t₀[i,j] += -t₀
            if N > 2
                H_t₀[j,i] += -t₀
            end
        end
        for (i,j) in eachcol(nbr0 .+ N)    
            H_t₀[i,j] += -t₀
            if N > 2
                H_t₀[j,i] += -t₀
            end
        end

        if t₁ != 0
            # generate next-nearest neighbor table
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

    if pht == true
        # add s-wave pairing 
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

        if dims > 1
            # add sPDW pairing
            add_pairing_symmetry!(
                "sPDW",
                determinantal_parameters,
                optimize,
                H_vpars,
                V,
                model_geometry,
                pht
            )

            debug && println("Hamiltonian::build_variational_hamiltonian() : ")
            debug && println("adding site-dependent s-wave pairing matrix")
            debug && println("initial Δ_spd = ", determinantal_parameters.det_pars.Δ_spd)
            if optimize.Δ_spd
                debug && println("optimize = true")
            else
                debug && println("optimize = false")
            end

            # add d-wave pairing 
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

            # add dPDW pairing
            add_pairing_symmetry!(
                "dPDW",
                determinantal_parameters,
                optimize,
                H_vpars,
                V,
                model_geometry,
                pht
            )

            debug && println("Hamiltonian::build_variational_hamiltonian() : ")
            debug && println("adding site-dependent d-wave pairing matrix")
            debug && println("initial Δ_dpd = ", determinantal_parameters.det_pars.Δ_dpd)
            if optimize.Δ_dpd
                debug && println("optimize = true")
            else
                debug && println("optimize = false")
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
        model_geometry,
        pht
    )

    debug && println("Hamiltonian::build_variational_hamiltonian() : ")
    debug && println("adding spin-x matrix")
    debug && println("initial Δ_sx = ", determinantal_parameters.det_pars.Δ_sx)
    if optimize.Δ_sz
        debug && println("optimize = true")
    else
        debug && println("optimize = false")
    end

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
    debug && println("initial Δ_sz = ", determinantal_parameters.det_pars.Δ_sz)
    if optimize.Δ_sz
        debug && println("optimize = true")
    else
        debug && println("optimize = false")
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
        debug && println("initial Δ_ssd = ", determinantal_parameters.det_pars.Δ_ssd)
        if optimize.Δ_ssd
            debug && println("optimize = true")
        else
            debug && println("optimize = false")
        end
    end

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
        debug && println("initial Δ_csd = ", determinantal_parameters.det_pars.Δ_csd)
        if optimize.Δ_csd
            debug && println("optimize = true")
        else
            debug && println("optimize = false")
        end
    end

    # @assert length(H_vpars) == determinantal_parameters.num_det_pars
    # @assert length(V) == determinantal_parameters.num_det_opts

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

- `symmetry::String`: type of pairing symmetry: "s", "d", "sPDW", or "dPDW
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

        # if Δ_d is being optimized, store V_dwave matrix
        if optimize.Δ_d
            push!(V, V_dwave)
        end
    elseif symmetry == "sPDW"
        @assert pht == true
        @assert dims > 1
        
        # create neighbor table
        nbr_table = build_neighbor_table(
            bonds[1],
            model_geometry.unit_cell,
            model_geometry.lattice
        )

        # pairing momentum
        q_p = determinantal_parameters.det_pars.q_p

        # s-wave pair fields
        Δ_spd = determinantal_parameters.det_pars.Δ_spd

        # initialize sPDW matrices
        H_sff = zeros(Complex, 2*N, 2*N)     
        H_slo = zeros(Complex, 2*N, 2*N)     
        V_spd = Matrix{Complex}[]
        for i in 1:N
            push!(V_spd, zeros(Complex, 2*N, 2*N))
        end

        # loop over neighbors in the x-direction 
        for col in 1:N
            i, j = nbr_table[:, col]

            # loop over x-direction 
            V_spd[i][i, j + N] = 1.0 
            V_spd[j][j, i + N] = 1.0 

            # Hermitian conjugate
            V_spd[i][j + N, i] = 1.0 
            V_spd[j][i + N, j] = 1.0 
        end

        # loop over neighbors in the y-direction 
        for col in N+1:2*N
            i, j = nbr_table[:, col]

            V_spd[i][i, j + N] = 1.0 
            V_spd[j][j, i + N] = 1.0 

            # Hermitian conjugate
            V_spd[i][j + N, i] = 1.0 
            V_spd[j][i + N, j] = 1.0 
        end

        # get all lattice site positions
        positions = []
        for i in 1:N
            pos = site_to_loc(i, unit_cell, lattice)[1]
            push!(positions, pos)
        end

        for i in 1:N
            # construct Fulde-Ferrell term
            ff_phase = exp(im*dot(q_p, positions[i]))
            H_sff += Δ_spd[i] * ff_phase * V_spd[i]      

            # construct Larkin-Ovchinnikov term
            lo_phase = exp(im*dot(q_p, positions[i])) + exp(-im*dot(q_p, positions[i]))
            H_slo += Δ_spd[i] * lo_phase * V_spd[i]                

            # if Δ_spd is being optimized, store each V_spd matrix
            if optimize.Δ_spd
                push!(V, V_spd[i])
            end
        end

        # add Fulde-Ferrell sPDW matrix
        push!(H_vpars, H_sff)

        # add Larkin-Ovchinnikov sPDW matrix
        push!(H_vpars, H_slo)

        # # add complete sPDW matrix
        # H_spdw = H_sff + Hslo
        # push!(H_vpars, H_spdw)
    elseif symmetry == "dPDW"
        @assert pht == true
        @assert dims > 1

        # create neighbor table
        nbr_table = build_neighbor_table(
            bonds[1],
            model_geometry.unit_cell,
            model_geometry.lattice
        )

        # pairing momentum
        q_p = determinantal_parameters.det_pars.q_p

        # d-wave pair fields
        Δ_dpd = determinantal_parameters.det_pars.Δ_dpd

        # initialize dPDW matrices
        H_dff = zeros(Float32, 2*N, 2*N)     
        H_dlo = zeros(Float32, 2*N, 2*N)   
        V_dpd = Matrix{Float32}[]
        V_dpdx = Matrix{Float32}[]
        V_dpdy = Matrix{Float32}[]
        for i in 1:N
            push!(V_dpd, zeros(Float32, 2*N, 2*N))
            push!(V_dpdx, zeros(Float32, 2*N, 2*N))
            push!(V_dpdy, zeros(Float32, 2*N, 2*N))
        end

        # loop over neighbors in the x-direction 
        for col in 1:N
            i, j = nbr_table[:, col]

            # loop over x-direction 
            V_dpdx[i][i, j + N] = 1.0
            V_dpdx[j][j, i + N] = 1.0

            # Hermitian conjugate
            V_dpdx[i][j + N, i] = 1.0
            V_dpdx[j][i + N, j] = 1.0
        end

        # loop over neighbors in the y-direction 
        for col in N+1:2*N
            i, j = nbr_table[:, col]

            V_dpdy[i][i, j + N] = -1.0
            V_dpdy[j][j, i + N] = -1.0

            # Hermitian conjugate
            V_dpdy[i][j + N, i] = -1.0
            V_dpdy[j][i + N, j] = -1.0
        end

        # complete V_dpd matrix
        for i in 1:N
            V_dpd[i] += V_dpdx[i]
            V_dpd[i] += V_dpdy[i] 
        end

        # get all lattice site positions
        positions = []
        for i in 1:N
            pos = site_to_loc(i, unit_cell, lattice)[1]
            push!(positions, pos)
        end

        for i in 1:N
            # construct Fulde-Ferrell term
            ff_phase = exp(im*dot(q_p, positions[i]))
            H_dff += Δ_dpd[i] * ff_phase * V_dpd[i]       

            # construct Larkin-Ovchinnikov term
            lo_phase = exp(im*dot(q_p, positions[i])) + exp(-im*dot(q_p, positions[i]))
            H_dlo += Δ_dpd[i] * lo_phase * V_dpd[i]                

            # if Δ_dpd is being optimized, store each V_dpd matrix
            if optimize.Δ_spd
                push!(V, V_dpd[i])
            end
        end

        # add Fulde-Ferrell dPDW matrix
        push!(H_vpars, H_dff)

        # add Larkin-Ovchinnikov dPDW matrix
        push!(H_vpars, H_dlo)

        # # add complete dPDW matrix
        # H_dpdw = H_dff + Hdlo
        # push!(H_vpars, H_dpdw)
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

- `order::String`: type of spin order: "spin-x", "spin-z", or "site-dependent"
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
        Δ_sz = determinantal_parameters.det_pars.Δ_sz

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
            H_afm += Δ_sz * V_afm
            push!(H_vpars, H_afm)

            # if Δ_sz is being optimized, store Vafm matrix
            if optimize.Δ_sz
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
            H_afm += Δ_sz * V_afm_neg
            push!(H_vpars, H_afm)

            # if Δ_sz is being optimized, store Vafm matrix
            if optimize.Δ_sz
                push!(V, V_afm_neg)
            end
        end
    elseif order == "spin-x"
        # spin-x parameter
        Δ_sx = determinantal_parameters.det_pars.Δ_sx

        # spin-x matrix (off-diagonal in spin)
        H_sx = zeros(Complex, 2*N, 2*N)

        for s in 1:N
            up_idx = s
            dn_idx = s + N

            # Add spin-flip terms: S^x = (1/2)(|↑⟩⟨↓| + |↓⟩⟨↑|)
            H_sx[up_idx, dn_idx] += Δ_sx / 2
            H_sx[dn_idx, up_idx] += Δ_sx / 2
        end

        # Add to variational Hamiltonian
        push!(H_vpars, H_sx)

        # If Δ_x is being optimized, store the matrix
        if optimize.Δ_sx
            push!(V, H_sx)
        end
    elseif order == "site-dependent"
        # lattice dimensions
        L = model_geometry.lattice.L

        Δ_ssd = determinantal_parameters.det_pars.Δ_ssd
        ssd_vectors = []
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
                push!(ssd_vectors, vec)
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
                push!(ssd_vectors, vec)  
            end
        end

        for (i, ssd_vec) in enumerate(ssd_vectors)
            H_ssd = zeros(AbstractFloat, 2*N, 2*N)
            ssd_vec_neg = copy(ssd_vec)
            ssd_vec_neg[N+1:2*N] .= -ssd_vec_neg[N+1:2*N]
        
            for s in 1:2*N
                idx = get_index_from_spindex(s, model_geometry)
                loc = site_to_loc(idx, model_geometry.unit_cell, model_geometry.lattice)
                ssd_vec_neg[s] *= (-1)^(loc[1][1] + loc[1][2])
            end
        
            V_ssd_neg = LinearAlgebra.Diagonal(ssd_vec_neg)
            H_ssd += Δ_ssd[i] * V_ssd_neg
            push!(H_vpars, H_ssd)
        
            if optimize.Δ_ssd
                push!(V, V_ssd_neg)
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
        Δ_csd = determinantal_parameters.det_pars.Δ_csd
        csd_vectors = []

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
                push!(csd_vectors, vec)
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
                push!(csd_vectors, vec)  
            end
        end

        for (i, csd_vec) in enumerate(csd_vectors)
            H_csd = zeros(AbstractFloat, 2*N, 2*N)
            V_csd = LinearAlgebra.Diagonal(csd_vec)
            H_csd += Δ_csd[i] .* V_csd
            push!(H_vpars, H_csd)

            # if Δ_csd is being optimized, store Vcsd matrix
            if optimize.Δ_csd
                push!(V, V_csd)
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
                  Ne::Int )::Bool

Checks whether a energy configuration is open shell.

"""
function is_openshell(
    ε::Vector{Float64}, 
    Ne::Int
)::Bool
    return ε[Ne + 1] - ε[Ne] < 0.0001
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


















