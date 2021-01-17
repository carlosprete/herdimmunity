## Simulation about heterogeneity and herd immunity
using PyPlot, ProgressMeter, JLD, Dates, DSP, DataStructures
# Plot results?
PLOT = true
SAVE = true
SAVEPLOT = true
TextPlot = false
PLOTtxt = false
ComputeRt = true

if PLOTtxt
    using UnicodePlots
end
include("SEIRmodel.jl")   # Continuous-time model and options
include("DiscreteSEIR.jl") # Discrete-time model
include("plotresults.jl")
include("stochasticsimulations.jl")
##
# Configuration of the simulation -- choose one of the sets of parameter choices for your simulation.  Create new codes for new simulations, this way it'll be easy to reproduce simulations for a paper.
Simulation = :Manaus_Dispersion_NPI_SeroRev_Discrete
#:Manaus_NoDispersion_NoAge_NPI_SeroRev_Discrete
#:Manaus_NoDispersion_NPI_SeroRev_Discrete
#:SP_Quarantine_LowDispersion_GovSP_Discrete
#:Manaus_Age_NoNPI_Discrete
#:ManausHomog_NoNPI_Discrete #:SP_Quarantine_LowDispersion_GovSP_Slow_SeroRev_Discrete
#:SP_Quarantine_LowDispersion_GovSP_Discrete
        #:SPDispersion_NoNPI_Discrete
        #:Manaus_1_5_Dispersion_NoNPI_Discrete #:SP_Quarantine_LowDispersion_GovSP_Slow_SeroRev_Discrete #:SP_NoAge_Dispersion_EstimatedRt_Discrete  #:SPHomog_InLoco_Discrete
caseopt = caseoptions(Simulation)
#caseopt[:N0] = 1
#caseopt[:Dispersion] = 0.9
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
Nsim = 100

## Main program - first, create variable with chosen set of parameters



@time s,e,i,t,d,Ninfected,p,OT, Rt = stochasticsimulation(caseopt, tspan, Nsim)

NewVar = (occursin("NewVariant", String(caseopt[:Model])) ? "_R0var=" * string(caseopt[:R0var]) : "")

LossImmTC = caseopt[:LossImmRate] != :None ? string(round(1 / caseopt[:LossImmRate], digits = 2)) : ""

dir = "Results/"*string(Simulation)*"_Nsim=$(Nsim)_N0=$(get(caseopt,:N0,1))_Dispersion_$(dispersion_factor(caseopt[:Dispersion]))_$(string(get(caseopt,:Quarantine,"")))_$(caseopt[:NPI]== :None ? "" : "NPI_")$((get(caseopt,:NPIprediction,"")))$(NewVar)_LossImmunityProb_$(string(caseopt[:LossImmProb]))_LossImmunityTimeConst_$(LossImmTC)_SymmetricSQRT_$(now()))/"
if SAVE
    mkdir(dir)
    # The macro @save has problems with caseopt sometimes.
    JLD.save(dir*"Data.jld", "caseopt", caseopt, "s", s, "e", e, "i", i, "Ninfected", Ninfected, "d", d, "t", t, "Rt", Rt, "OT", OT)
end
if :Ivar in p[:IndexNames] # If there are 2 variants, add the infected for each one.
    for sim in 1:Nsim
        i[:, sim] .+= OT[sim][:Ivar]
    end
end
if PLOT
    plotresults(dir, Simulation, caseopt, s, i, d, Rt, Ninfected, true; αd = p[:NPI])
end

if TextPlot
    plottext(caseopt, s, i, t)
end
