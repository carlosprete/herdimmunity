using Distributions, Statistics

"""
    function discretestochasticSEIRSeroRevNewVariant!(upost, upre, p, t)
Discrete-time stochastic solver with seroreversion - computes the state at next time instant `t+Δt`, `upost`, given current state `upre`.  `upost = model(upre, t)`.  `p` is a parameter vector.
"""
function discretestochasticSEIRSeroRevNewVariant!(upost, upre, p, t)
    sind = p[:S]
    eind = p[:E]
    iind = p[:I]
    rspind = p[:Rsp]     # fraction of R that does not serorevert
    rlsind = p[:Rl]     # fraction of R that loses immunity against first variant
    rvarind = p[:Rvar] # fraction of R that loses immunity only to new variant
    srind = p[:Sr]       # fraction of S that was infected by one of the variants, but lost immunity to both
    svarind = p[:Svar] # fraction of S that was infected by variant 0, but is not immune to variant 1
    evarind = p[:Evar]
    ivarind = p[:Ivar]
    ATi = p[:ContactMatrix]        # Mixing matrix
    σ = p[:IncubationRate]      # incubation rate (ρ in the paper)
    μ = p[:RecoveryRate]        # recovery rate (γ in the paper)
    p0 = p[:LossImmProb]        # probability of loss of immunity (to both variants)
    ar = p[:LossImmRate]        # rate of seroreversion (loss of immunity - to both variants)
    p1 = p[:LossImmProbVar]       # probability of loss of immunity only to second variant
    ar1 = p[:LossImmRateVar]    # rate of immunity loss only against second variant
    N = p[:Population]         # total population
    Δt = p[:StepSize]     # stepsize (in days)
    Δtnorm = -Δt / N
    TimeIntroduction = p[:TimeIntroduction] # time new variant is introduced
#    α = p[:NPI]   # function to compute social distancing index
    q = get(caseopt, :q, 1.0)
    s = view(upre, sind)
    sr = view(upre, srind)
    svar = view(upre, svarind)
    e = view(upre, eind)
    i = view(upre, iind)
    rsp = view(upre, rspind)
    rls = view(upre, rlsind)
    rvar = view(upre, rvarind)

    Ai = ATi(t, i.^q)
    F = 1 .- exp.(Δtnorm.*(Ai))
    B = rand.(Binomial.(s, F))
    Bsr = rand.(Binomial.(sr, F))
    C = rand.(Binomial.(e, 1 - exp(-σ*Δt)))
    D = rand.(Binomial.(i, 1 - exp(-μ*Δt)))
    F .= rand.(Binomial.(rls, 1 - exp(-ar*Δt)))
    G = rand.(Binomial.(rvar, 1- exp(-ar1*Δt)))

    upost[sind] = s .- B
    upost[srind] = sr .- Bsr .+ F # Before introduction of new variant, Sr decreases only due to first variant
    upost[svarind] = svar .+ G # Before introduction of new variant, Svar only increases.
    upost[eind] = e .+ B .+ Bsr .- C
    upost[iind] = i .+ C .- D
    Dsp = zeros(length(D), 3)
    for ind in eachindex(D)
        Dsp[ind,:] = rand(Multinomial(D[ind], [p0, p1, 1-p0-p1]))
        # Dsp[:,1]: part of D that will eventually loose immunity against variants 1 (and 2)
        # Dsp[:,2]: part of D that will eventually loose immunity against variant 2 only
        # Dsp[:,3]: part of D that will not loose immunity.
    end
    upost[rspind] = rsp .+ Dsp[:, 3]
    upost[rlsind] = rls .+ Dsp[:, 1] .- F
    upost[rvarind] = rvar .+ Dsp[:, 2] .- G

    # New variant appears at t == TimeIntroduction
    if t == TimeIntroduction
        N0var = p[:N0var] # Number of cases for new variant
        N0vareff = sum(upost[sind]) >= N0var ? N0var : sum(upost[sind])
        πv = p[:PopDistribution]
        paciente0 = rand(DiscreteNonParametric(p[:Evar], πv), N0vareff)
        for pac in paciente0
           upost[pac] += 1
        end
        upost[sind] .-= upost[evarind]
    elseif t > TimeIntroduction
        R0fac = p[:R0fac] # (R0 for new variant)/R0

        evar = view(upre, evarind)
        ivar = view(upre, ivarind)
        Ai = (R0fac*Δtnorm) .* ATi(t, ivar.^q)
        F .= 1 .- exp.(Ai)
        B .= rand.(Binomial.(sr, F))
        C .= rand.(Binomial.(svar, F))

        upost[evarind] .= rand.(Binomial.(upost[sind], F))
        upost[ivarind] .= rand.(Binomial.(evar, 1 - exp(-σ*Δt)))
        D .= rand.(Binomial.(ivar, 1 - exp(-μ*Δt)))

        upost[sind] .-= upost[evarind]
        upost[srind] .-= B
        upost[svarind] .-= C
        upost[evarind] .+= evar .+ B .+ C .- upost[ivarind]
        upost[ivarind] .+= ivar .- D
        upost[rspind] .+= D # We are assuming here that recovery from the second variant gives full immunity

    end




    for k = 1:length(upost)
        if upost[k] < 0
            println("Negative intermediate value!")
            upost[k] = 0
        end
    end

    return upost
end

"""
    function discretestochasticSEIRNewVariant!(upost, upre, p, t)
Discrete-time stochastic solver - computes the state at next time instant `t+Δt`, `upost`, given current state `upre`.  `upost = model(upre, t)`.  `p` is a parameter vector.
"""
function discretestochasticSEIRNewVariant!(upost, upre, p, t)
    sind = p[:S]
    eind = p[:E]
    iind = p[:I]
    rind = p[:R]
    evarind = p[:Evar]
    ivarind = p[:Ivar]
    rvarind = p[:Rvar]
    ATi = p[:ContactMatrix]        # Mixing matrix
    σ = p[:IncubationRate]      # incubation rate (ρ in the paper)
    μ = p[:RecoveryRate]        # recovery rate (γ in the paper)
    N = p[:Population]         # total population
    Δt = p[:StepSize]          # stepsize (in days)
#    α = p[:NPI]   # function to compute social distancing index
    q = get(caseopt, :q, 1.0)
    TimeIntroduction = p[:TimeIntroduction] # time new variant is introduced

    s = view(upre, sind)
    e = view(upre, eind)
    i = view(upre, iind)
    r = view(upre, rind)
    Δtnorm = -Δt / N
    A = Δtnorm .* ATi(t, i.^q)
    upost[eind] .= rand.(Binomial.(s, 1 .- exp.(A)))
    upost[iind] .= rand.(Binomial.(e, 1 - exp(-σ*Δt)))
    upost[rind] .= rand.(Binomial.(i, 1 - exp(-μ*Δt)))

    upost[sind] .= s .- upost[eind]
    upost[eind] .+= e .- upost[iind]
    upost[iind] .+= i .- upost[rind]
    upost[rind] .+= r
    if t == TimeIntroduction
        N0var = p[:N0var] # Number of cases for new variant
        N0vareff = sum(upost[sind]) >= N0var ? N0var : sum(upost[sind])
        πv = p[:PopDistribution]
        paciente0 = rand(DiscreteNonParametric(p[:Evar], πv), N0vareff)
        for pac in paciente0
           upost[pac] += 1
        end
        upost[sind] .-= upost[evarind]
    end



    if t > TimeIntroduction
        R0fac = p[:R0fac] # (R0 for new variant)/R0

        evar = view(upre, evarind)
        ivar = view(upre, ivarind)
        rvar = view(upre, rvarind)
        A = (R0fac*Δtnorm) .* ATi(t, ivar.^q)
        upost[evarind] .= rand.(Binomial.(upost[sind], 1 .- exp.(A)))
        upost[ivarind] .= rand.(Binomial.(evar, 1 - exp(-σ*Δt)))
        upost[rvarind] .= rand.(Binomial.(ivar, 1 - exp(-μ*Δt)))

        upost[sind] .-= upost[evarind]
        upost[evarind] .+= evar .- upost[ivarind]
        upost[ivarind] .+= ivar .- upost[rvarind]
        upost[rvarind] .+= rvar

    end


    for k = 1:length(upost)
        if upost[k] < 0
            upost[k] = 0
        end
    end

    return upost
end


"""
    function discretestochasticSEIR!(upost, upre, p, t)
Discrete-time stochastic solver - computes the state at next time instant `t+Δt`, `upost`, given current state `upre`.  `upost = model(upre, t)`.  `p` is a parameter vector.
"""
function discretestochasticSEIR!(upost, upre, p, t)
    sind = p[:S]
    eind = p[:E]
    iind = p[:I]
    rind = p[:R]
    ATi = p[:ContactMatrix]        # Mixing matrix
    σ = p[:IncubationRate]      # incubation rate (ρ in the paper)
    μ = p[:RecoveryRate]        # recovery rate (γ in the paper)
    N = p[:Population]         # total population
    Δt = p[:StepSize]          # stepsize (in days)
#    α = p[:NPI]   # function to compute social distancing index
    q = get(caseopt, :q, 1.0)

    s = view(upre, sind)
    e = view(upre, eind)
    i = view(upre, iind)
    r = view(upre, rind)
    Δtnorm = -Δt / N
    A = Δtnorm .* ATi(t, i.^q)
    upost[eind] .= rand.(Binomial.(s, 1 .- exp.(A)))
    upost[iind] .= rand.(Binomial.(e, 1 - exp(-σ*Δt)))
    upost[rind] .= rand.(Binomial.(i, 1 - exp(-μ*Δt)))

    upost[sind] .= s .- upost[eind]
    upost[eind] .+= e .- upost[iind]
    upost[iind] .+= i .- upost[rind]
    upost[rind] .+= r

    for k = 1:length(upost)
        if upost[k] < 0
            upost[k] = 0
        end
    end

    return upost
end

"""
    function discretestochasticSEIRSeroRev!(upost, upre, p, t)
Discrete-time stochastic solver with seroreversion - computes the state at next time instant `t+Δt`, `upost`, given current state `upre`.  `upost = model(upre, t)`.  `p` is a parameter vector.
"""
function discretestochasticSEIRSeroRev!(upost, upre, p, t)
    sind = p[:S]
    eind = p[:E]
    iind = p[:I]
    rspind = p[:Rsp]     # fraction of R that does not serorevert
    rsnind = p[:Rsn]     # fraction of R that seroreverts, but remains immune
    rlsind = p[:Rl]     # fraction of R that loses immunity
    srind = p[:Sr]       # fraction of S that was infected at least once
    ATi = p[:ContactMatrix]        # Mixing matrix
    σ = p[:IncubationRate]      # incubation rate (ρ in the paper)
    μ = p[:RecoveryRate]        # recovery rate (γ in the paper)
    ϵ0 = p[:SeroRevProb]        # probability of seroreversion
    p0 = p[:LossImmProb]        # probability of loss of immunity
    ar = p[:LossImmRate]        # rate of seroreversion (loss of immunity)
    N = p[:Population]         # total population
    Δt = p[:StepSize]          # stepsize (in days)
#    α = p[:NPI]   # function to compute social distancing index
    q = get(caseopt, :q, 1.0)
    s = view(upre, sind)
    sr = view(upre, srind)
    e = view(upre, eind)
    i = view(upre, iind)
    rsp = view(upre, rspind)
    rsn = view(upre, rsnind)
    rls = view(upre, rlsind)
    Ai = ATi(t, i.^q)
    F = 1 .- exp.(-(Δt/N).*(Ai))

    B = rand.(Binomial.(s, F))
    Bsr = rand.(Binomial.(sr, F))
    C = rand.(Binomial.(e, 1 - exp(-σ*Δt)))
    D = rand.(Binomial.(i, 1 - exp(-μ*Δt)))
    F = rand.(Binomial.(rls, 1 - exp(-ar*Δt)))

    upost[sind] = s .- B
    upost[srind] = sr .- Bsr .+ F
    upost[eind] = e .+ B .+ Bsr .- C
    upost[iind] = i .+ C .- D
    Dsp = rand.(Binomial.(D, 1-ϵ0))
    Dls = rand.(Binomial.(D-Dsp, p0))
    upost[rspind] = rsp .+ Dsp
    upost[rsnind] = rsn .+ D .- Dsp .- Dls
    upost[rlsind] = rls .+ Dls .- F

    for k = 1:length(upost)
        if upost[k] < 0
            upost[k] = 0
        end
    end

    return upost
end
"""
    function discretesolve(u0, md, p, tspan)
Discrete-time solver - computes state evolution in tspan[1]...tspan[2], with initial condition `u0` and model `md` (e.g., `md=discretestochasticSEIR!`).
`p` is a parameter vector.
"""
function discretesolve(u0, md, p, tspan)
    t = tspan[1]:p[:StepSize]:tspan[2]
    u = [zeros(Int64, length(u0)) for row in 1:length(t)]
    u[1] = copy(u0)
    for n in 2:length(t)
        md(u[n], u[n-1], p, t[n-1])
    end
    return u,t
end
