## Finds probability of an epidemic taking place for a given situation
using PyPlot, ProgressMeter, JLD, Dates, DSP, DataStructures, Distributions, StatsBase
# Plot results?
PLOT = true
SAVEPLOT = true
TextPlot = false
PLOTtxt = false
if PLOTtxt
    using UnicodePlots
end
include("SEIRmodel.jl")   # Continuous-time model and options
include("DiscreteSEIR.jl") # Discrete-time model
include("plotresults.jl")
include("stochasticsimulations.jl")
include("ProbExtinction.jl")
##
# Configuration of the simulation -- choose one of the sets of parameter choices for your simulation.  Create new codes for new simulations, this way it'll be easy to reproduce simulations for a paper.
Simulation = :SP_NoAge_Dispersion_NoNPI_Discrete #:Manaus_1_5_Dispersion_NoNPI_Discrete #:SPHomog_InLoco_Discrete
caseopt = caseoptions(Simulation)
N0 = 1
caseopt[:N0] = N0
caseopt[:q] = 1.0
# Time span for the simulation (in days)
tspan = (0, 365)

# Number of different runs (to estimate the probability of an outbreak, or the interval of possible final values)
Nsim = 10_000
println()
## Main program - first, create variable with chosen set of parameters
k = [0.04, 0.07, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
theorecticalProbOutbreak = 1 .- ProbExtinction.(Rzero(caseopt[:R0]), k)
ProbOutbreak = zeros(length(k))
AttackRate = zeros(length(k))
@time @showprogress "Outer Loop..." for n in 1:length(k)
    caseopt[:Dispersion] = k[n]
    caseopt[:N0] = N0
    s,e,i,t,d,Ninfected,p = stochasticsimulation(caseopt, tspan, Nsim)
    ProbOutbreak[n] = count(s[end,:] .< 0.95 * totalpopulation(caseopt[:Population])) / Nsim
    AttackRate[n] = (totalpopulation(caseopt[:Population]) - mean(s[end,findall(s[end,:] .< 0.95 * totalpopulation(caseopt[:Population]))])) * 100 / totalpopulation(caseopt[:Population])
    println(ProbOutbreak[n],"\t",AttackRate[n])
end

if gethostname() == "tesla"
    ioff()
end
dir = "Results/OutbreakProb_$(Simulation)_Nsim_$(Nsim)_$(today())/"
mkdir(dir)
figure(1);clf()
stem(log10.(k), ProbOutbreak, label = "Simulation")
stem(log10.(k), theorecticalProbOutbreak, label = "Theorectical", "r", markerfmt = "r^")
xlabel(L"Dispersion factor $k$")
ylabel(L"$P(\{$Outbreak$\})$")
title("Outbreak probability vs. dispersion factor")
legend()
savefig(dir*"OutbreakProb_$(Simulation)_Nsim_$(Nsim).svg")
figure(2);clf()
stem(k, AttackRate)
xlabel(L"Dispersion factor $k$")
ylabel("Final infected")
title("Attack rate vs. dispersion factor")
if SAVEPLOT
    savefig(dir*"AttackRate_$(Simulation)_Nsim_$(Nsim).svg")
end
JLD.save(dir*"OutbreakProb_$(Simulation)_Nsim_$(Nsim).jld", "ProbOutbreak", ProbOutbreak, "AttackRate", AttackRate, "Simulation", Simulation, "N0", N0, "Nsim", Nsim, "caseopt", caseopt, "k", k)
