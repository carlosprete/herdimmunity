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

# Final attack rate
"""
    attackrate(R0, k)
Returns the final attack rate for an epidemic with reproduction number `R0` and dispersion factor `k`.  An homogeneous population corresponds to `k == Inf`.
"""
function attackrate(R0, k)
    if k == Inf
        AR = find_zero(ar -> ar + log(1 - ar)/R0, (1e-9, 1.0))
    else
        G1 = x -> (1 + (R0/k) * (1 - x))^(-k)
        if k == 1
            G0 = (x, p0) -> p0 + (1-p0) * (1 - log(1+R0*(1-x))/log(1+R0))
            G0′unnorm = R0 / log(1 + R0)


        else
            G0 = (x, p0) -> p0 + (1-p0) * (1 - (1 - R0*x / (R0 + k))^(1-k))/(1 - (k / (R0 + k))^(1-k))
            G0′unnorm = ((1 - k)/k) * R0 / ((k/(R0 + k))^(k-1) - 1)
            println("G0' $(G0′unnorm)")
        end
        if G0′unnorm > R0
            p0 = 1 - R0 / G0′unnorm
        else
            p0 = 0.0
        end
        println("p0 $p0")
        uk = find_zero(u -> G1(u) - u, (0, 1-1e-9))
        println("uk $uk")
        AR = 1 - G0(uk, p0)
    end
    return AR
end

function G1negbin(x, R0, k)
    return (1 + (R0/k) * (1 - x))^(-k)
end
