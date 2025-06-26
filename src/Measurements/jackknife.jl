@doc raw"""

    jackknife( g::Function,
               samples...;
               # KEYWORD ARGUMENTS
               bias_corrected = true,
               jackknife_samples = similar.(samples),
               jackknife_g = similar(samples[1]) )

Propagates errors through the evaluation of a function `g` given the binned `samples`,
returning both the mean and error. If the keyword argument `bias = true`, then the 
``\mathcal{O}(1/N)`` bias is corrected. The keyword arguments `jackknife_samples` and 
`jackknife_g` can be passed to avoid temporary memory allocations.

"""
function jackknife(
    g::Function,
    samples...;
    # KEYWORD ARGUMENTS
    bias_corrected = true,
    jackknife_samples = similar.(samples),
    jackknife_g = similar(samples[1])
)

    # get sample size
    N = length(jackknife_g)

    # iterate over input variables
    for i in eachindex(samples)

        # calculate mean current input variable
        x̄ = mean(samples[i])

        # iterate over samples
        for j in eachindex(samples[i])

            # calculate the mean of the j'th jackknife sample by updating the mean to
            # reflect removing the j'th sample
            jackknife_samples[i][j] = (N*x̄ - samples[i][j])/(N-1)
        end
    end

    # evaluate the input function using the jackknife sample means
    @. jackknife_g = g(jackknife_samples...)

    # calculate jackknife mean
    ḡ = mean(jackknife_g)

    # calculate jackkife error
    Δg = sqrt( (N-1) * varm(jackknife_g, ḡ, corrected=false) )

    # correct O(1/N) bias, usually doesn't matter as error scales as O(1/sqrt(N))
    # and is typically much larger than the bias
    if bias_corrected
        Ḡ = g(map(mean, samples)...)
        ḡ = N * Ḡ - (N-1) * ḡ
    end

    return ḡ, Δg
end


# @doc raw"""

#     jackknife_stats( bin_sample::Vector{Float64} )

# Performs jacknife resampling to obtain the standard deviation per bin.

# """
# function jackknife_stats(
#     bin_samples::Vector{Float64}
# )
#     N = length(bin_samples)
#     if N ≤ 1
#         return (mean=NaN, std=NaN)
#     end

#     jk_means = [mean(deleteat!(copy(bin_samples), i)) for i in 1:N]
#     jk_mean = mean(jk_means)
#     jk_std = sqrt((N - 1) / N * sum((jk_means .- jk_mean).^2))

#     return (mean=jk_mean, std=jk_std)
# end