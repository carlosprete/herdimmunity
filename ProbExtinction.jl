# Probability of extinction (no large outbreak) for a Negative Binomial distribution with mean R0 and dispersion parameter k
using Roots

"""
    ProbExtinction(R0, k)
Finds the probability of extinction (i.e., of a large outbreak taking place) assuming a negative binomial distribution for the number of secondary cases for each infected with mean `R0` and dispersion `k`.
"""
function ProbExtinction(R0,k)
    ## Probability generation function
    g(s) = (1 + R0*(1-s)/k) ^(-k) - s
    ## Find probability of extinction
    if R0 < 1
        p = 1.0
    else
        p = find_zero(g, (0.0, 1 - 1e-6))
    end
    return p
end
