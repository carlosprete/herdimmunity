using PyPlot
@isdefined(PLOTtxt) ? true : PLOTtxt = false

if PLOTtxt
    using UnicodePlots
end
function plotresults(dir, Simulation, caseopt, s, i, d, Ninfected, PlotSocialDist = false; αd = d -> 1.0)
    Nsteps,Nsim = size(s)
    colors = ("b","g","r","c","m","y","tab:blue","tab:orange","tab:green","tab:red",
    "tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan")
    Ncolors = length(colors)
    figure(1)
    clf()
    Pop = totalpopulation(caseopt[:Population])

    for sim in 1:Nsim
        PyPlot.semilogy(d, i[:,sim], linewidth = 2)
    end

    infectedtitle = "No. of infected people"
    PyPlot.title(infectedtitle)
    PyPlot.ylabel(L"I(t)")
    PyPlot.xlabel("Date")
    xticks(rotation=30)
    grid()
    if SAVEPLOT
        savefig(dir*string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$(string(get(caseopt,:NPIprediction,"")))_LossImmunityProb_$(string(caseopt[:SeroRevProb]))_LossImmuntyRate_$(string(caseopt[:LossImmRate]))_SymmetricSQRT_$(today())Infected.svg")
    end
    fig = figure(2)
    clf()
    outbreaks = findall(s[end,:] .< 0.95Pop) # finds the simulations with outbreaks to plot
    ax1 = plt.axes()

    red = "tab:red"
    ax1.set_xlabel("Date")
    ax1.set_ylabel(L"$I(t)$", color = red)
    ax1.tick_params(axis = "y", labelcolor = red)
    xticks(rotation=30)

    #legend(loc="center right")

    ax2 = ax1.twinx()
    blue = "tab:blue"
    ax2.set_ylabel("Cumulative infected", color=blue)  # we already handled the x-label with ax1
    ax2.tick_params(axis="y", labelcolor=blue)
    PyPlot.xlabel("Date")
    PyPlot.title("Infected population (left); Cumulative infected (right)")


    # This is the homogeneous, theoretical herd immunity level.
    HerdImmunityHomogeneousLevel = (1.0 - 1.0 / Rzero(caseopt[:R0]))*100
    ax2.plot(d[[1, end]], [1.0, 1.0]*HerdImmunityHomogeneousLevel, lw = 2, ":", label = "Homogeneous herd immunity level")

    color = 1
    for sim in outbreaks

        ax1.plot(d, i[:,sim], lw = 2, color = colors[mod1(color,Ncolors)], "-", label = "Infected people")

        ax2.plot(d, (Pop .- s[:,sim])*100/Pop, color = colors[mod1(color,Ncolors)], lw = 2, label = "Cummulative infected") #, color=color)

        color += 1




    #
        HerdImmunityPointHomog = findall((Pop .- s[:,sim]) .>= HerdImmunityHomogeneousLevel*Pop/100)
        if length(HerdImmunityPointHomog) > 1
            HerdImmunityPointHomog = HerdImmunityPointHomog[1]
            ax2.plot(d[HerdImmunityPointHomog],(Pop - s[HerdImmunityPointHomog,sim])*100/Pop, "ok", label = "Homogeneous herd immunity point")
        end

    # Maxi = findall([i[j+1]-i[j] < 0 && i[j] - i[j-1]>0 for j in 2:length(i)-1])
    # if length(Maxi) > 1
    #     ki = findall(Maxi .> HerdImmunityPoint)
    #     if length(ki) > 0
    #         HIP2 = Maxi[ki[end]]
    #         HIP2 = findmax(i[HIP2-10:min(HIP2+10,length(d))])[2]+HIP2-11
    #         ax2.plot([d[HIP2], d[HIP2]],[0,1 - s[HIP2]],":k", label="Herd immunity attained")
    #     end
    # end
    end
    grid(true)
    if SAVEPLOT
        savefig(dir*string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$(string(get(caseopt,:NPIprediction,"")))_LossImmunityProb_$(string(caseopt[:SeroRevProb]))_LossImmuntyRate_$(string(caseopt[:LossImmRate]))_SymmetricSQRT_$(today())Cummulative.svg")
    end
    #legend(loc = "best")
    show();
    if PlotSocialDist
        figure(3), clf()
        d0 = firstday(caseopt[:FirstDay])
        dint = d0 .+ Day.(tspan[1]:tspan[2])
        PyPlot.plot(dint, 1 .- αd.(tspan[1]:tspan[2]), lw = 3)
        hweek = PolynomialRatio(fill(1 / 7, 7), [1])
        socialindexfiltered = filt(hweek, 1 .-  p[:NPI].(tspan[1]:tspan[2]))
        PyPlot.plot(dint[1:end-3], socialindexfiltered[4:end])

        PyPlot.xlabel("Date")
        PyPlot.title("Social distancing index")
        xticks(rotation=30)
        if SAVEPLOT
            savefig(dir*string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$(string(get(caseopt,:NPIprediction,"")))_LossImmunityProb_$(string(caseopt[:SeroRevProb]))_LossImmuntyRate_$(string(caseopt[:LossImmRate]))_SymmetricSQRT_$(today())SocialDistancing.svg")
        end
    end

    figure(4);clf()
    PyPlot.hist(Ninfected/totalpopulation(caseopt[:Population]), density=false, bins=-0.005:0.01:1.005, rwidth = 0.5, align="mid")
    PyPlot.xlabel("Percentage of total population infected")
    PyPlot.title("Number of cases in $Nsim simulations")
    show()
    if SAVEPLOT
        savefig(dir*string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$(string(get(caseopt,:NPIprediction,"")))_LossImmunityProb_$(string(caseopt[:SeroRevProb]))_LossImmuntyRate_$(string(caseopt[:LossImmRate]))_SymmetricSQRT_$(today())Distribution.svg")
    end
end
if PLOTtxt
    function plottext(caseopt, s, i, t)
        Pop = totalpopulation(caseopt[:Population])
        outbreaks = findall(s[end,:] .< 0.95Pop) # finds the simulations with outbreaks to plot
        unplt=lineplot(t, (Pop .- s[:,1])*100/Pop)
        println(unplt)
        for sim = 2:length(outbreaks)
            println(lineplot!(unplt, t, (Pop .- s[:,outbreaks[sim]])*100/Pop))
        end
    end
end
