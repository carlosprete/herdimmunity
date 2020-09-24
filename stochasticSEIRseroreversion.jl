## Simulation about heterogeneity and herd immunity
using PyPlot, ProgressMeter
include("SEIRmodel.jl")   # Continuous-time model and options
include("DiscreteSEIR.jl") # Discrete-time model

##
# Configuration of the simulation -- choose one of the sets of parameter choices for your simulation.  Create new codes for new simulations, this way it'll be easy to reproduce simulations for a paper.
Simulation = :SPHeter_SPGov_Discrete #:ManausHomog_NoNPI_Discrete #:SPHomog_InLoco_Discrete
##
# Chooses the appropriate set of parameters for each simulation
function caseoptions(Simulation)
    caseopt = @match Simulation begin

        :SPHeter_SPGov_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP,
            :MixingMatrix => :BrittonScience2020,
            :ActivityVector => :Masks,
            :ActivityStructure => :Masks,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :ManausHomog_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :Manaus,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :MixingMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
            :ManausHomog_InLoco_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :LossImmRate => :None,
                :SeroRevProb => :None,
                :LossImmProb => :None,
                :Population => :Manaus,
                :FirstDay => :Manaus,
                :R0 => :AM,
                :NPI => :ManausInLoco,
                :RecoveryRate => :BrittonScience2020,
                :IncubationRate => :BrittonScience2020,
                :AgeStructure => :None,
                :MixingMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end


            :SPHomog_InLoco_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :LossImmRate => :None,
                :SeroRevProb => :None,
                :LossImmProb => :None,
                :Population => :SP,
                :FirstDay => :SP,
                :R0 => :SP,
                :NPI => :SPInLoco,
                :RecoveryRate => :BrittonScience2020,
                :IncubationRate => :BrittonScience2020,
                :AgeStructure => :None,
                :MixingMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end

            :SPHomog_NPI_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :LossImmRate => :None,
                :SeroRevProb => :None,
                :LossImmProb => :None,
                :Population => :SP,
                :FirstDay => :SP,
                :R0 => :SP,
                :NPI => :SPGov,
                :RecoveryRate => :BrittonScience2020,
                :IncubationRate => :BrittonScience2020,
                :AgeStructure => :None,
                :MixingMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end
            :SPHomog_noNPI_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :LossImmRate => :None,
                :SeroRevProb => :None,
                :LossImmProb => :None,
                :Population => :SP,
                :FirstDay => :SP,
                :R0 => :SP,
                :NPI => :None,
                :RecoveryRate => :BrittonScience2020,
                :IncubationRate => :BrittonScience2020,
                :AgeStructure => :None,
                :MixingMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end

    end
    return caseopt
end

## Main program - first, create variable with chosen set of parameters
caseopt = caseoptions(Simulation)

# Note that the text in Britton's paper is misleading, the definition of A was taken out from their software (ref [14] in the paper).
md,πv,A,p,u0 = model(caseopt)

# Time span for the simulation (in days)
tspan = (0, 1.5*365)

# Number of different runs (to estimate the probability of an outbreak, or the interval of possible final values)
Nsim = 1
# This computes the distribution in the total number of infected pacients in Nsim realizations
Nsteps = floor(Int, (tspan[2]-tspan[1]) / discretestep(caseopt[:StepSize]))
Ninfected = zeros(Int, Nsim)
d0 = firstday(caseopt[:FirstDay])
u = zeros(Int, Nsteps)
t = zeros(Nsteps)
s = zeros(Int, Nsteps)
@showprogress "Computing..." for sim in 1:Nsim
    global u, s, t, Ninfected
    u,t = discretesolve!(u0, md, p, tspan)
    # u is the state, with the S, E, I, R compartments.  u[p[:S]] are the states related to the S compartment, etc.  Note that there may be several S compartments, one for each age and activity group (and similarly for E, I, R). In the models with seroreversion, there may be more compartments --- see SEIRmodel.jl
    s = sumind(u, p[:S])
    # This is the percentage of people that get infected
    Ninfected[sim] = totalpopulation(caseopt[:Population]) .- s[end]
end


d = DateTime(d0) .+ Second.(round.(Int, t * 24 * 3600))
i = sumind(u, p[:I])  # number of infected at each time instant (only the last simulation)

figure(1)
#clf()
Pop = totalpopulation(caseopt[:Population])

PyPlot.semilogy(d, i, linewidth = 2)

infectedtitle = "No. of infected people"
title(infectedtitle)
ylabel(L"I(t)")
xlabel("Date")
grid()

fig = figure(2)
clf()
ax1 = plt.axes()

color = "tab:red"
ax1.set_xlabel("Date")
ax1.set_ylabel(L"$I(t)$", color = color)
ax1.plot(d, i, lw = 2, "-", color = color, label = "Infected people")
ax1.tick_params(axis = "y", labelcolor = color)
# This is the homogeneous, theoretical herd immunity level.
HerdImmunityHomogeneousLevel = (1.0 - 1.0 / Rzero(caseopt[:R0]))*100


legend(loc="center right")
ax2 = ax1.twinx()


color = "tab:blue"
ax2.set_ylabel("Cumulative infected", color=color)  # we already handled the x-label with ax1
ax2.plot(d, (Pop .- s)*100/Pop, lw = 2, label = "Cummulative infected", color=color)
ax2.tick_params(axis="y", labelcolor=color)

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
xlabel("Date")
title("Infected population (left); Cumulative infected (right)")


#
ax2.plot(d[[1, end]], [1.0, 1.0]*HerdImmunityHomogeneousLevel, lw = 2, ":", label = "Homogeneous herd immunity level")
HerdImmunityPointHomog = findall((Pop .- s) .>= HerdImmunityHomogeneousLevel*Pop/100)
if length(HerdImmunityPointHomog) > 1
    HerdImmunityPointHomog = HerdImmunityPointHomog[1]
    ax2.plot(d[HerdImmunityPointHomog],(Pop - s[HerdImmunityPointHomog])*100/Pop, "ok", label = "Homogeneous herd immunity point")
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

figure(4);clf()
PyPlot.hist(Ninfected/totalpopulation(caseopt[:Population]), density=false, bins=-0.005:0.01:1.005, rwidth = 0.5, align="mid")
xlabel("Percentage of total population infected")
