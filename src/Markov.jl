@doc raw"""

    local_fermion_update!(  detwf::DeterminantalWavefunction{T, V, E1, I},
                            model_geometry::ModelGeometry,
                            Np::I,  
                            δW::E2,
                            n_stab_W::I,
                            rng::AbstractRNG    ) where {T<:Number, V, E1<:Number, I<:Integer, E2<:AbstractFloat}

Performs a local update to the fermionic sector by attempting to move a particle ``\beta`` at a randomly 
selected lattice site ``k`` to another random site ``l``, with the move accepted or rejected
by the Meteropolis algorithm.

- `detwf::DeterminantalWavefunction{T, V, E, I}`: current determinantal variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `δW::E`: error threshold for the Green's function.
- `n_stab_W::I`: frequency of Green's function stabilization steps.
- `rng::AbstractRNG`: random number generator.

"""
function local_fermion_update!(
    detwf::DeterminantalWavefunction{T, V, E1, I},
    model_geometry::ModelGeometry,
    Np::I,  
    δW::E2,
    n_stab_W::I,
    rng::AbstractRNG
) where {T<:Number, V, E1<:Number, I<:Integer, E2<:AbstractFloat}
    acceptances = 0.0
    rejections = 0.0

    @debug """
    Markov:local_fermion_update!() :
    Starting new Monte Carlo step!
    """

    # attempt configuration update
    acceptance = metropolis_step(
        detwf, 
        model_geometry,
        Np, 
        δW, 
        n_stab_W, 
        rng
    )

    if acceptance == "accepted"
        acceptances += 1.0
    elseif acceptance == "rejected"
        rejections += 1.0
    end

    # calculate local acceptance rate
    rate_calc = acceptances / (acceptances + rejections)
    if isnan(rate_calc)
        acceptance_rate = 0.0
    else
        acceptance_rate = rate_calc
    end

    return acceptance_rate, detwf
end


@doc raw"""

    local_fermion_update!(  detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                            jastrow_factor::JastrowFactor{E2}, 
                            model_geometry::ModelGeometry, 
                            jastrow_parameters::JastrowParameters{S, K, V, I}, 
                            Np::I, 
                            δW::E2, 
                            δT::E2, 
                            n_stab_W::I,
                            n_stab_T::I,
                            rng::AbstractRNG,
                            pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Performs a local update to the fermionic sector by attempting to move a particle ``\beta`` at a randomly 
selected lattice site ``k`` to another random site ``l``, with the move accepted or rejected
by the Meteropolis algorithm.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current determinantal variational wavefunction.
- `jastrow_facor::JastrowFactor{E}`: current Jastrow factor.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow variational parameters.
- `Np::I`: total number of particles in the system.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `pht::Bool`: whether model is particle-hole transformed.
- `n_stab_W::I`: frequency of Green's function stabilization steps.
- `n_stab_T::I`: frequency of T vector stabilization steps.
- `δW::E`: error threshold for the Green's function.
- `δT::E`: error treshold for the Jastrow T vector.
- `rng::AbstractRNG`: random number generator.

"""
function local_fermion_update!(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2}, 
    model_geometry::ModelGeometry, 
    jastrow_parameters::JastrowParameters{S, K, V, I}, 
    Np::I, 
    δW::E2, 
    δT::E2, 
    n_stab_W::I,
    n_stab_T::I,
    rng::AbstractRNG,
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    acceptances = 0.0
    rejections = 0.0

    @debug """
    Markov:local_fermion_update!() :
    Starting new Monte Carlo step!
    """

    # attempt configuration update
    acceptance = metropolis_step(
        detwf, 
        jastrow_factor, 
        model_geometry,
        jastrow_parameters, 
        Np, 
        δW, 
        δT,
        n_stab_W, 
        n_stab_T,
        rng,
        pht
    )

    if acceptance == "accepted"
        acceptances += 1.0
    elseif acceptance == "rejected"
        rejections += 1.0
    end

    # calculate local acceptance rate
    rate_calc = acceptances / (acceptances + rejections)
    if isnan(rate_calc)
        acceptance_rate = 0.0
    else
        acceptance_rate = rate_calc
    end

    return acceptance_rate, detwf, jastrow_factor
end


@doc raw"""

    local_fermion_update!(  detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                            jastrow_factor_1::JastrowFactor{E2}, 
                            jastrow_factor_2::JastrowFactor{E2},
                            model_geometry::ModelGeometry, 
                            jastrow_parameters_1::JastrowParameters{S, K, V, I},  
                            jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
                            Np::I, 
                            δW::E2, 
                            δT::E2,
                            n_stab_W::I,
                            n_stab_T::I,
                            rng::AbstractRNG,
                            pht::Bool   ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Performs a local update to the fermionic sector by attempting to move a particle ``\beta`` at a randomly 
selected lattice site ``k`` to another random site ``l``, with the move accepted or rejected
by the Meteropolis algorithm.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current determinantal variational wavefunction.
- `jastrow_factor_1::Jastrow{E}`: first Jastrow factor.
- `jastrow_factor_2::Jastrow{E}`: second Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow variational parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow variational parameters.
- `Np::I`: total number of particles in the system.
- `δW::E`: error threshold for the Green's function.
- `δT::E`: error treshold for the Jastrow T vector.
- `n_stab_W::I`: frequency of Green's function stabilization steps.
- `n_stab_T::I`: frequency of T vector stabilization steps.
- `rng::AbstractRNG`: random number generator.
- `pht::Bool`: whether model is particle-hole tranformed.

"""
function local_fermion_update!(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor_1::JastrowFactor{E2}, 
    jastrow_factor_2::JastrowFactor{E2},
    model_geometry::ModelGeometry, 
    jastrow_parameters_1::JastrowParameters{S, K, V, I},  
    jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
    Np::I, 
    δW::E2, 
    δT::E2,
    n_stab_W::I,
    n_stab_T::I,
    rng::AbstractRNG,
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    acceptances = 0.0
    rejections = 0.0

    @debug """
    Markov:local_fermion_update!() :
    Starting new Monte Carlo step!
    """

    # attempt configuration update
    acceptance = metropolis_step(
        detwf, 
        jastrow_factor_1,  
        jastrow_factor_2, 
        model_geometry,
        jastrow_parameters_1,
        jastrow_parameters_2, 
        Np, 
        δW, 
        δT,
        n_stab_W, 
        n_stab_T,  
        rng,
        pht
    )

    if acceptance == "accepted"
        acceptances += 1.0
    elseif acceptance == "rejected"
        rejections += 1.0
    end

    # calculate local acceptance rate
    rate_calc = acceptances / (acceptances + rejections)
    if isnan(rate_calc)
        acceptance_rate = 0.0
    else
        acceptance_rate = rate_calc
    end

    return acceptance_rate, detwf, jastrow_factor_1, jastrow_factor_2
end


@doc raw"""

    metropolis_step(    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        model_geometry::ModelGeometry,
                        Np::I, 
                        δW::E2,
                        n_stab_W::I, 
                        rng::AbstractRNG    ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}

Determines the acceptance criterion for a particle ``\beta`` at lattice site ``k`` moving to a randomly selected site ``l``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current determinantal variational wavefunction.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `Np::I`: total number of particles in the system.
- `δW::E`: error threshold for the Green's function.
- `n_stab_W::I`: frequency of Green's function stabilization steps.
- `rng::AbstractRNG`: random number generator.

"""
function metropolis_step(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    model_geometry::ModelGeometry,
    Np::I, 
    δW::E2,
    n_stab_W::I, 
    rng::AbstractRNG
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat}
    # propose a random move
    markov_move = propose_random_move(
        Np, 
        detwf.pconfig, 
        model_geometry, 
        rng
    ) 

    if markov_move.possible == false
        @debug """
        Markov::metropolis_step() : 
        hop impossible!
        """

        return "impossible"
    else
        # get wavefunction ratio (element of the equal-time Green's function)
        Rₛ = real(detwf.W[markov_move.l, markov_move.particle])

        # acceptance probability
        acceptance_prob = Rₛ^2   # TODO: include phase

        @debug """
        Markov::metropolis_step() :
        hop possible =>
        R_s = $(Rₛ)
        acceptance_prob = $(acceptance_prob)
        """

        if acceptance_prob > rand(rng)
            @debug """
            Markov::metropolis_step() : 
            hop accepted!
            """

            # hop the particle
            hop!(
                markov_move, 
                detwf.pconfig
            )

            # perform rank-1 update to W matrix
            update_equal_time_greens!(
                markov_move, 
                detwf, 
                model_geometry, 
                Np, 
                n_stab_W, 
                δW
            ) 

            return "accepted"
        else
            @debug """
            Markov::metropolis_step() : 
            hop rejected!
            """

            return "rejected"
        end
    end
end


@doc raw"""

    metropolis_step(    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        jastrow_factor::JastrowFactor{E2},  
                        model_geometry::ModelGeometry,
                        jastrow_parameters::JastrowParameters{S, K, V, I},
                        Np::I, 
                        δW::E2, 
                        δT::E2,
                        n_stab_W::I, 
                        n_stab_T::I,  
                        rng::AbstractRNG,
                        pht::Bool ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Determines the acceptance criterion for a particle ``\beta`` at lattice site ``k`` moving to a randomly selected site ``l``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `jastrow_factor::JastrowFactor{E}`: current Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `jastrow_parameters::JastrowParameters{S, K, V, I}`: current set of Jastrow parameters.
- `Np::I`: total number of particles in the system.
- `δW::E`: error threshold for the Green's function.
- `δT::E`: error threshold for the T vector.
- `n_stab_W::I`: frequency of Green's function stabilization steps.
- `n_stab_T::I`: frequency of T vector stabilization steps.
- `rng::AbstractRNG`: random number generator.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function metropolis_step(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor::JastrowFactor{E2},  
    model_geometry::ModelGeometry,
    jastrow_parameters::JastrowParameters{S, K, V, I},
    Np::I, 
    δW::E2, 
    δT::E2,
    n_stab_W::I, 
    n_stab_T::I,  
    rng::AbstractRNG,
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    # propose a random move
    markov_move = propose_random_move(
        Np, 
        detwf.pconfig, 
        model_geometry, 
        rng
    ) 

    if markov_move.possible == false
        @debug """
        Markov::metropolis_step() : 
        hop impossible!
        """

        return "impossible"
    else
        # check particle spin
        spin = get_spindex_type(
            markov_move.k, 
            model_geometry
        )

        # get element of T vector
        @assert jastrow_parameters.jastrow_type == "e-den-den" || jastrow_parameters.jastrow_type == "e-spn-spn"
        Rⱼ = get_fermionic_jastrow_ratio(
            markov_move, 
            jastrow_parameters, 
            jastrow_factor, 
            pht, 
            spin, 
            model_geometry
        ) 

        # get wavefunction ratio (element of the equal-time Green's function)
        Rₛ = real(detwf.W[markov_move.l, markov_move.particle])

        # acceptance probability
        acceptance_prob = Rⱼ^2 * Rₛ^2      

        @debug """
        Markov::metropolis_step() :
        hop possible =>
        R_j = $(Rⱼ)
        R_s = $(Rₛ)
        acceptance_prob = $(acceptance_prob)
        """

        if acceptance_prob > rand(rng)
            @debug """
            Markov::metropolis_step() : 
            hop accepted!
            """

            # hop the particle
            hop!(markov_move, detwf.pconfig) 

            # perform rank-1 update to W matrix
            update_equal_time_greens!(
                markov_move, 
                detwf, model_geometry, 
                Np, 
                n_stab_W, 
                δW
            ) 

            # update T vector
            update_fermionic_Tvec!(
                markov_move, 
                spin, 
                jastrow_parameters, 
                jastrow_factor, 
                detwf,
                model_geometry, 
                n_stab_T, 
                δT, 
                pht
            )

            return "accepted"
        else
            @debug """
            Markov::metropolis_step() : 
            hop rejected!
            """

            return "rejected"
        end 
    end
end


@doc raw"""

    metropolis_step(    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
                        jastrow_factor_1::JastrowFactor{E2},  
                        jastrow_factor_2::JastrowFactor{E2}, 
                        model_geometry::ModelGeometry,
                        jastrow_parameters_1::JastrowParameters{S, K, V, I},
                        jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
                        Np::I, 
                        δW::E2, 
                        δT::E2,
                        n_stab_W::I, 
                        n_stab_T::I, 
                        rng::AbstractRNG,
                        pht::Bool   ) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}

Determines the acceptance criterion for a particle ``\beta`` at lattice site ``k`` moving to a randomly selected site ``l``.

- `detwf::DeterminantalWavefunction{T, Q, E, I}`: current variational wavefunction.
- `jastrow_factor_1::JastrowFactor{E}`: first Jastrow factor.
- `jastrow_factor_2::JastrowFactor{E}`: second Jastrow factor.
- `model_geometry::ModelGeometry`: contains unit cell and lattice quantities.
- `jastrow_parameters_1::JastrowParameters{S, K, V, I}`: first set of Jastrow parameters.
- `jastrow_parameters_2::JastrowParameters{S, K, V, I}`: second set of Jastrow parameters.
- `δW::E`: error threshold for the Green's function.
- `δT::E`: error threshold for the T vector.
- `Np::I`: total number of particles in the system.
- `n_stab_W::I`: frequency of Green's function stabilization steps.
- `n_stab_T::I`: frequency of T vector stabilization steps.
- `rng::AbstractRNG`: random number generator.
- `pht::Bool`: whether model is particle-hole transformed.

"""
function metropolis_step(
    detwf::DeterminantalWavefunction{T, Q, E1, I}, 
    jastrow_factor_1::JastrowFactor{E2},  
    jastrow_factor_2::JastrowFactor{E2}, 
    model_geometry::ModelGeometry,
    jastrow_parameters_1::JastrowParameters{S, K, V, I},
    jastrow_parameters_2::JastrowParameters{S, K, V, I}, 
    Np::I, 
    δW::E2, 
    δT::E2,
    n_stab_W::I, 
    n_stab_T::I, 
    rng::AbstractRNG,
    pht::Bool
) where {T<:Number, Q, E1<:Number, I<:Integer, E2<:AbstractFloat, S<:AbstractString, K, V}
    # propose a random move
    markov_move = propose_random_move(
        Np, 
        detwf.pconfig, 
        model_geometry, 
        rng
    ) 

    if markov_move.possible == false
        @debug """
        Markov::metropolis_step() : 
        hop impossible!
        """

        return "impossible"
    else
        # check particle spin
        spin = get_spindex_type(
            markov_move.k, 
            model_geometry
        )

        @assert jastrow_parameters_1.jastrow_type == "e-den-den" || jastrow_parameters_1.jastrow_type == "e-spn-spn"
        @assert jastrow_parameters_2.jastrow_type == "e-den-den" || jastrow_parameters_2.jastrow_type == "e-spn-spn"

        # get element of the first T vector
        Rⱼ₁ = get_fermionic_jastrow_ratio(
            markov_move, 
            jastrow_parameters_1, 
            jastrow_factor_1, 
            pht, 
            spin, 
            model_geometry
        )

        # get element of the second T vector
        Rⱼ₂ = get_fermionic_jastrow_ratio(
            markov_move, 
            jastrow_parameters_2, 
            jastrow_factor_2, 
            pht, 
            spin, 
            model_geometry
        )

        # get wavefunction ratio (element of the equal-time Green's function)
        Rₛ = real(detwf.W[markov_move.l, markov_move.particle])

        # acceptance probability
        acceptance_prob = Rⱼ₁^2 * Rⱼ₂^2 * Rₛ^2       # TODO: include phase

        @debug """
        Markov::metropolis_step() :
        hop possible =>
        R_j1 = $(Rⱼ₁)
        R_j2 = $(Rⱼ₂)
        R_s = $(Rₛ)
        acceptance_prob = $(acceptance_prob)
        """

        if acceptance_prob > rand(rng)
            @debug """
            Markov::metropolis_step() : 
            hop accepted!
            """

            # hop the particle
            hop!(markov_move, detwf.pconfig)

            # perform rank-1 update to W matrix
            update_equal_time_greens!(
                markov_move, 
                detwf, model_geometry, 
                Np, 
                n_stab_W, 
                δW
            ) 

            # update T vectors
            update_fermionic_Tvec!(
                markov_move, 
                spin, 
                jastrow_parameters_1, 
                jastrow_factor_1, 
                detwf,
                model_geometry, 
                n_stab_T, 
                δT, 
                pht
            )

            update_fermionic_Tvec!(
                markov_move, 
                spin, 
                jastrow_parameters_2, 
                jastrow_factor_2, 
                detwf,
                model_geometry, 
                n_stab_T, 
                δT, 
                pht
            )

            return "accepted"
        else
            @debug """
            Markov::metropolis_step() : 
            hop rejected!
            """

            return "rejected"
        end
    end
end


