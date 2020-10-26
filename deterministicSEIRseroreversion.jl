## Simulation about heterogeneity and herd immunity
using OrdinaryDiffEq, PyPlot, Match
include("SEIRmodel.jl")
include("ContactMatrices.jl")
include("plotresults.jl")
## Init

Simulation = :SPhomogeneous_nonNPI_Continuous


# Note that the text in the paper is misleading, the definition of A was taken out from their software (ref [14] in the paper).
caseopt = caseoptions(Simulation)
md,πv,A,p,u0 = model(caseopt)


tspan = (0.0, 184.0)
prob = ODEProblem(md, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol = 1e-15, abstol = 1e-15)


i = sumind(sol.u, p[:I])
s = sumind(sol.u, p[:S])

# This is the percentage of high-activity people that get infected


d0 = firstday(caseopt[:FirstDay])
d = DateTime(d0) .+ Second.(round.(Int,sol.t * 24 * 3600))

figure(1),clf()
Pop = totalpopulation(caseopt[:Population])

PyPlot.semilogy(d, Pop*i, linewidth = 2)

infectedtitle = (Pop == 1 ? "Fraction" : "No.") * " of infected people"
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
ax1.plot([d[HerdImmunityPoint], d[HerdImmunityPoint]],[0,maximum(i)],":k")


legend(loc="center right")
ax2 = ax1.twinx()


color = "tab:blue"
ax2.set_ylabel(L"1-s(t)", color=color)  # we already handled the x-label with ax1
ax2.plot(d, 1 .- s, lw = 2, label = "Cummulative infected", color=color)
ax2.tick_params(axis="y", labelcolor=color)

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
xlabel("Date")
title("Fraction of infected population (left); Cumulative infected (right)")



ax2.plot(d[[1, end]], [1.0, 1.0]*HerdImmunityLevel, lw = 2, ":", label = "Homogeneous herd immunity level")
HerdImmunityPointHomog = findall((1 .- s) .>= HerdImmunityLevel)
if length(HerdImmunityPointHomog) > 1
    HerdImmunityPointHomog = HerdImmunityPointHomog[1]
    ax2.plot(d[HerdImmunityPointHomog],(1 - s[HerdImmunityPointHomog]), "ok", label = "Homogeneous herd immunity point")
end
Maxi = findall([i[j+1]-i[j] < 0 && i[j] - i[j-1]>0 for j in 2:length(i)-1])
if length(Maxi) > 1
    ki = findall(Maxi .> HerdImmunityPoint)
    if length(ki) > 0
        HIP2 = Maxi[ki[end]]
        HIP2 = findmax(i[HIP2-10:min(HIP2+10,length(d))])[2]+HIP2-11
        ax2.plot([d[HIP2], d[HIP2]],[0,1 - s[HIP2]],":k", label="Herd immunity attained")
    end
end
legend(loc = "best")
grid()
show();


figure(3), clf()
dint = d0 .+ Day.(tspan[1]:tspan[2])
PyPlot.plot(dint, 1 .- p[:NPI].(tspan[1]:tspan[2]), lw=3)
using DSP
hweek = PolynomialRatio(fill(1/7, 7),[1])
socialindexfiltered = filt(hweek, 1 .- p[:NPI].(tspan[1]:tspan[2]))
PyPlot.plot(dint[1:end-3], socialindexfiltered[4:end])

xlabel("Date")
title("Social distancing index")
