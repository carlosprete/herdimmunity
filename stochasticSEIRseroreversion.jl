## Simulation about heterogeneity and herd immunity
using PyPlot, ProgressMeter, JLD, Dates, DSP, DataStructures
# Plot results?
PLOT = true
SAVE = false
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
Simulation = :Manaus_Quarantine_Dispersion_SchoolClosure_LossImm_Discrete
#:SP_Quarantine_LowDispersion_GovSP_Discrete
#:Manaus_Age_NoNPI_Discrete
#:ManausHomog_NoNPI_Discrete #:SP_Quarantine_LowDispersion_GovSP_Slow_SeroRev_Discrete
#:SP_Quarantine_LowDispersion_GovSP_Discrete
        #:SPDispersion_NoNPI_Discrete
        #:Manaus_1_5_Dispersion_NoNPI_Discrete #:SP_Quarantine_LowDispersion_GovSP_Slow_SeroRev_Discrete #:SP_NoAge_Dispersion_EstimatedRt_Discrete  #:SPHomog_InLoco_Discrete
caseopt = caseoptions(Simulation)
#caseopt[:N0] = 1
caseopt[:Dispersion] = 5.0
#caseopt[:NPI] = :SPGov_double
#caseopt[:LossImmRate] = 1 / 180
#caseopt[:SeroRevProb] = 0.3
# caseopt[:LossImmProb] = 1.0
#caseopt[:Dispersion] = 20.0
#caseopt[:R0] = 3.5
#caseopt[:Dispersion] = 1.5
#caseopt[:ActivityVector] = :None
#caseopt[:ActivityStructure] = :None
# Time span for the simulation (in days)
tspan = (0, 1.5*365)

# Number of different runs (to estimate the probability of an outbreak, or the interval of possible final values)
Nsim = 10

## Main program - first, create variable with chosen set of parameters



@time s,e,i,t,d,Ninfected,p = stochasticsimulation(caseopt, tspan, Nsim)

dir = "Results/"*string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$((get(caseopt,:NPIprediction,"")))_LossImmunityProb_$(string(caseopt[:SeroRevProb]))_LossImmunityRate_$(string(caseopt[:LossImmRate]))_SymmetricSQRT_$(today()))/"
if SAVE
    mkdir(dir)
    # The macro @save has problems with caseopt sometimes.
    JLD.save(dir*"Data.jld", "caseopt", caseopt, "s", s, "e", e, "i", i, "Ninfected", Ninfected, "d", d, "t", t)
end
if PLOT
    plotresults(dir, Simulation, caseopt, s, i, d, Ninfected, true; αd = p[:NPI])
end

if TextPlot
    plottext(caseopt, s, i, t)
end
