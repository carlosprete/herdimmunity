## Simulation about heterogeneity and herd immunity
using OrdinaryDiffEq, PyPlot, Match, DSP, JLD
include("SEIRmodel.jl")
include("ContactMatrices.jl")
include("plotresults.jl")
## Init
SAVE = true
SAVEPLOT = true

Simulation = :SP_NoAge_EstimatedRt
#dinit = Date("2020-01-01")
caseopt = caseoptions(Simulation)
dinit = firstday(caseopt[:FirstDay])
N0 = 5000 * round(totalpopulation(caseopt[:Population])/totalpopulation(:Manaus), digits = 2)
caseopt[:N0] = N0

dir = "Results/"*string(Simulation)*"_$(dinit)_N0=$(get(caseopt,:N0,1))__$(caseopt[:NPI]== :None ? "" : "NPI_")$((get(caseopt,:NPIprediction,"")))_$(today())/"
if SAVE || SAVEPLOT
    try
        mkdir(dir)
    catch e
    end
end
# Note that the text in the paper is misleading, the definition of A was taken out from their software (ref [14] in the paper).
#Δdate = Day(17)
#dinit = dinit + Δdate # For Manaus only
days = Day.([0])#, 5, 10, 15, 20, 25, 30]) #, 35, 40, 45, 50, 55, 60])
Pop = totalpopulation(caseopt[:Population])
tspan = (0.0, 365.0)

for dd in days
    global p, dint
    d0 = dinit + dd
    caseopt[:FirstDay] = d0
    md,πv,A,p,u0 = model(caseopt)


    prob = ODEProblem(md, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol = 1e-15, abstol = 1e-15)


    i = sumind(sol.u, p[:I])
    s = sumind(sol.u, p[:S])
    e = sumind(sol.u, p[:E])
    r = sumind(sol.u, p[:R])

    # This is the percentage of high-activity people that get infected


    d = DateTime(d0) .+ Second.(round.(Int,sol.t * 24 * 3600))

    figure(1),clf()

    PyPlot.semilogy(d, Pop*i, linewidth = 2)
    xticks(rotation=30)
    infectedtitle = string(caseopt[:Population])*(Pop == 1 ? "Fraction" : "No.") * " of infected people"
    title(infectedtitle)
    if Pop == 1
        ylabel(L"i(t)")
    else
        ylabel(L"I(t)")
    end
    xlabel("Date")
    grid()
    show()

    fig = figure(2)
    clf()
    ax1 = plt.axes()

    color = "tab:red"
    ax1.set_xlabel("Date")
    ax1.set_ylabel(L"$i(t)$", color = color)
    ax1.plot(d, i, lw = 2, "-", color = color, label = "Fraction of infected people")
    ax1.tick_params(axis = "y", labelcolor = color)
    HerdImmunityPoint = findmax(i)[2]
    HerdImmunityLevel = 1.0 - 1.0 / Rzero(caseopt[:R0])
    #ax1.plot([d[HerdImmunityPoint], d[HerdImmunityPoint]],[0,maximum(i)],":k")
    xticks(rotation=30)

    ax1.legend(loc="center right")
    ax2 = ax1.twinx()


    color = "tab:blue"
    ax2.set_ylabel(L"1-s(t)", color=color)  # we already handled the x-label with ax1
    ax2.plot(d, 1 .- s, lw = 2, label = "Cummulative infected", color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    xlabel("Date")
    # title("Fraction of infected population (left); Cumulative infected (right)")
    title("Evolution of epidemic for "*string(caseopt[:Population])*", N0 = $(N0) $(N0 == 1 ? "case" : "cases") at d0 = $(d0)")



    #ax2.plot(d[[1, end]], [1.0, 1.0]*HerdImmunityLevel, lw = 2, ":", label = "Homogeneous herd immunity level")
    #HerdImmunityPointHomog = findall((1 .- s) .>= HerdImmunityLevel)
    #if length(HerdImmunityPointHomog) > 1
    #    HerdImmunityPointHomog = HerdImmunityPointHomog[1]
    #    ax2.plot(d[HerdImmunityPointHomog],(1 - s[HerdImmunityPointHomog]), "ok", label = "Cummulative value at peak")
    #end
    # Maxi = findall([i[j+1]-i[j] < 0 && i[j] - i[j-1]>0 for j in 2:length(i)-1])
    # if length(Maxi) > 1
    #     ki = findall(Maxi .> HerdImmunityPoint)
    #     if length(ki) > 0
    #         HIP2 = Maxi[ki[end]]
    #         HIP2 = findmax(i[HIP2-10:min(HIP2+10,length(d))])[2]+HIP2-11
    #         ax2.plot([d[HIP2], d[HIP2]],[0,1 - s[HIP2]],":k", label="Peak attained")
    #     end
    # end
    ax2.legend()
    grid()
    show();
    if SAVEPLOT
        savefig(dir*string(caseopt[:Population])*"Evolution_N0_$(N0)_d0_$(d0).svg")
    end

    figure(3), clf()
    dint = d0 .+ Day.(tspan[1]:tspan[2])
    PyPlot.plot(dint, p[:NPI].(tspan[1]:tspan[2]), lw=3)
    #PyPlot.plot(dint, 1 .- p[:NPI].(tspan[1]:tspan[2]), lw=3)
    hweek = PolynomialRatio(fill(1/7, 7),[1])
    #socialindexfiltered = filt(hweek, 1 .- p[:NPI].(tspan[1]:tspan[2]))
    socialindexfiltered = filt(hweek, p[:NPI].(tspan[1]:tspan[2]))

    PyPlot.plot(dint[1:end-3], socialindexfiltered[4:end])
    xticks(rotation=30)
    xlabel("Date")
    #title("Social distancing index, "*string(caseopt[:Population]))
    title("Rt, "*string(caseopt[:Population]))
    savefig(dir*string(caseopt[:Population])*"Rt_N0_$(N0)_$(d0).svg")

    if SAVE
        # The macro @save has problems with caseopt sometimes.
        JLD.save(dir*"Data$(d0).jld", "caseopt", caseopt, "s", s, "e", e, "i", i, "r", r, "dint", dint, "dinit", dinit, "d0", d0, "days", days, "t", sol.t)
    end
end
