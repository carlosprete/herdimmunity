## Simulation about heterogeneity and herd immunity
using PyPlot, ProgressMeter, JLD, Dates, DSP, DataStructures
# Plot results?
PLOT = true
SAVEPLOT = false
TextPlot = false
PLOTtxt = false
if PLOTtxt
    using UnicodePlots
end
include("SEIRmodel.jl")   # Continuous-time model and options
include("DiscreteSEIR.jl") # Discrete-time model
include("plotresults.jl")
include("stochasticsimulations.jl")
##
# Configuration of the simulation -- choose one of the sets of parameter choices for your simulation.  Create new codes for new simulations, this way it'll be easy to reproduce simulations for a paper.
Simulation = :Manaus_NoAge_Dispersion_1_5_NoNPI_Fractal_Discrete #:Manaus_1_5_Dispersion_NoNPI_Discrete #:SPHomog_InLoco_Discrete
caseopt = caseoptions(Simulation)
caseopt[:N0] = 100
caseopt[:q] = 1.0
caseopt[:Dispersion] = 1.5
#caseopt[:ActivityVector] = :None
#caseopt[:ActivityStructure] = :None
# Time span for the simulation (in days)
tspan = (0, 365)

# Number of different runs (to estimate the probability of an outbreak, or the interval of possible final values)
Nsim = 10

## Main program - first, create variable with chosen set of parameters



@time s,e,i,t,d,Ninfected,p = stochasticsimulation(caseopt, tspan, Nsim)
# The macro @save has problems with caseopt sometimes.
JLD.save( string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$((get(caseopt,:NPIprediction,"")))_LossImmunityProb_$(string(caseopt[:SeroRevProb]))_LossImmuntyRate_$(string(caseopt[:LossImmRate]))_SymmetricSQRT_$(today()).jld", "caseopt", caseopt, "s", s, "e", e, "i", i, "Ninfected", Ninfected, "d", d, "t", t)

if PLOT
    plotresults(Simulation, caseopt, s, i, d, Ninfected, true; αd = p[:NPI])
end

if TextPlot
    plottext(caseopt, s, i, t)
end
