function stochasticsimulation(caseopt, tspan, Nsim)

    # Note that the text in Britton's paper is misleading, the definition of A was taken out from their software (ref [14] in the paper).


    # This computes the distribution in the total number of infected pacients in Nsim realizations
    Nsteps = floor(Int, (tspan[2] - tspan[1]) / discretestep(caseopt[:StepSize])) + 1
    Ninfected = zeros(Int, Nsim)
    d0 = firstday(caseopt[:FirstDay])
    u = zeros(Int, Nsteps)
    t = zeros(Nsteps)
    s = zeros(Int, Nsteps, Nsim)
    e = zeros(Int, Nsteps, Nsim)
    i = zeros(Int, Nsteps, Nsim)
    Pop = totalpopulation(caseopt[:Population])
    p = Dict() # Attention: in some models, matrix A is random - p[:ContactMatrix] from the simulation will have only the last A

    Rt = Vector{Array{Float64}}(undef,Nsim)

    OutputTotals = Vector{Dict}(undef, Nsim)
    @showprogress "Computing..." for sim = 1:Nsim
        md, πv, A, p, u0 = model(caseopt)
        u, t = discretesolve(u0, md, p, tspan)
        # u is the state, with the S, E, I, R compartments.  u[p[:S]] are the states related to the S compartment, etc.  Note that there may be several S compartments, one for each age and activity group (and similarly for E, I, R). In the models with seroreversion, there may be more compartments --- see SEIRmodel.jl
        OutputTotals[sim] = Dict(compartment => sumind(u, p[compartment]) for compartment in p[:IndexNames])
        if ComputeRt
            Rt[sim] = ReplacementNumber(t, u, p, caseopt[:Model])
        end

        s[:, sim] = sumind(u, p[:S])
        e[:, sim] = sumind(u, p[:E])
        i[:, sim] = sumind(u, p[:I])
        # This is the percentage of people that get infected
        Ninfected[sim] = totalpopulation(caseopt[:Population]) .- s[end, sim]
        if PLOTtxt
            println("")
            println(lineplot(t, (Pop .- s[:,sim])*100/Pop, title="sim = $(sim), Attack rate = $(round((Pop-s[end,sim])*100/Pop,digits=2))"))
        end
    end

    d = DateTime(d0) .+ Second.(round.(Int, t * 24 * 3600))
    # number of infected at each time instant (only the last simulation)
    if ComputeRt
        return s, e, i, t, d, Ninfected, p, OutputTotals, Rt
    else
        return s, e, i, t, d, Ninfected, p, OutputTotals
    end
end
