# SEIR model for COVID-19 simulations

These files implement a few models for infectious diseases based on the classical SEIR model.  The program is based on different modules that allow the user to define and save the parameters necessary to reproduce a particular simulation in the form of a dictionary recording each option.  All models are variations on the two basic models described below.

## Basic models

1. Continuous-time deterministic SEIR

   The base of all the implemented continuous-time, deterministic models is the SEIR model defined by the equations
   $$
   \frac{d S}{dt} = -\beta \frac{IS}{N},
   $$

   $$
   \frac{dE}{dt} = \beta\frac{IS}{N} - \sigma E
   $$

   $$
   \frac{d I}{dt} = \sigma E - \mu I,
   $$

   $$
   \frac{dR}{dt} = \mu I,
   $$

   where $S(t)$ is the number of susceptible people at time $t$, $E$ is the number of exposed people (that have contracted the disease, but are not yet contagious), $I$ is the number of infected people (that are contagious), and $R$ is the number of recovered or dead.  Parameter $\beta$ is the contact rate, that is, the rate at which an infected person infects susceptibles, $\sigma$ is the incubation rate ($\tau_\sigma = 1/\sigma$ is the average incubation time), $\mu$ is the recovery rate ($\tau_\mu = 1/\mu$ is the average recovery time).

2. Discrete-time stochastic SEIR[^Based on Lekone and Finkenstädt, "Statistical inference in a stochastic epidemic SEIR model with control intervention: Ebola as a case study", _Biometrics_ vol. 62, Dec. 2006]

   $$
   S(t+\Delta t) = S(t)-\Delta S,\\
   E(t+\Delta t) = E(t) + \Delta S - \Delta E,\\
   I(t+\Delta t) = I(t) + \Delta E - \Delta I,\\
   R(t + \Delta t) = R(t) + \Delta I,
   $$
   where $\Delta t$ is the step size for the simulation, and the increments $\Delta S$, $\Delta E$ and $\Delta I$ follow binomial distributions  $\Delta S\sim \text{Bin}(S(t), p_{\text{exp}}$, $\Delta E \sim \text{Bin}(E(t), p_{\text{inc}})$, $\Delta I(t) \sim \text{Bin}(I(t), p_{\text{rec}})$ with exposure, incubation and recovery probabilities given by
   $$
   p_{\text{exp}} = 1 - e^{-\frac{\beta(t)}{N}\Delta t I(t)},\\
   p_{\text{inc}} = 1 - e^{-\sigma\Delta t},\\
   p_{\text{rec}} = 1 - e^{-\mu \Delta t}.
   $$
   

## Augmented models

Several features were added to these basic models:

1. A function $\alpha(t)$ that multiplies $\beta$, to form a simple model for the effect of non-pharmaceutical interventions (NPIs), such as social distancing, use of masks, school closures, etc.

2. The division of each compartment, $S$, $E$, $I$, $R$, into several subcompartments  each for a different age group (so that $S$, $E$, $I$ and $R$ are now vectors with dimensions equal to the number of age groups).  The interactions between age groups is described by a contact matrix $A$, such that
   $$
   \frac{dS}{dt} = - A^T I \frac{S}{N}
   $$
   for the continuous model.  For the discrete model, the exporsure probability now depends on the age compartments, and becomes vector such that
   $$
   p_{\text{exp}} = 1 - \exp.\left(-A^T \frac{I}{N} \Delta t\right),
   $$
   where $\exp.(x)$ should be understood as a broadcast in Julia, that is, a element-by-element operation.

   The contact matrix can depend on time.  In particular, if the contact matrices from [Prem 2020] are used, which provides separate contact matrices for home, school, work, and others, NPI measures such as closing of schools can be taken into account in the simulations.

3. The division of each compartment into _activity level_ groups to allow the introduction of dispersion[^See Lloyd-Smith et al., Superspreading and the effect of individual variation on disease emergence, Nature, 2005, and Britton et al, A mathematical model reveals the influence of population heterogeneity on herd immunity to SARS-CoV-2, Science, 2020].  Each activity level subgroup corresponds to a different number of contacts for that subgroup.  This difference in number of contacts corresponds to the individual reproductive number $\nu$ described in [Lloyd-Smith, 2005], and can be defined either in an ad-hoc manner, as described in [Britton, 2020], or using a gamma distribution, adapting the proposal in [Lloyd-Smith, 2005].

   Note that age and activity groups can be used simultaneously.

4. A fractal model adapting [Abbasi 2020] was also included, changing the update equations to
   $$
   \frac{dS}{dt} = - A^T I(0)^(1-q)I(t)^q \frac{S}{N},
   $$
   where $q$ is a parameter associated to Tsallis statistics.  The discrete-time model was similarly modified.

The options for the simulations are chosen through a dictionary, with is used by several functions to define the parameters of a particular simulation.  Saving the values of this option dictionary provides a simple way of reproducing the same simulation again.

The structure of the dictionary is (note that some options can be given as a Julia `Symbol`, as `:Population => :Manaus`, to be converted to values in special functions to be explained later, or directly as numerical values, as `:Population => 2_219_580` ).

```{julia}
:Simulation_Name => begin
  	Dict(:Model => :SEIRDiscrete,   # Solver - continuous or discrete-time
            :Population => :Manaus, # Total population for this simulation
            :FirstDay => :Manaus,   # Day corresponding to `t=0` in the simulation 
            :R0 => :SP,						  # The value of the base (average) R0 
            :NPI => :None,				  # Will NPI information be used?
            :RecoveryRate => :BrittonScience2020, # Recovery rate μ
            :IncubationRate => :BrittonScience2020, # Incubation rate σ
            :AgeStructure => :None,  # Should age structure of the population be used?
            :ContactMatrix => :None, # Choice of contact matrix between age groups
            :ActivityVector => :Superspreaders, # Individual contact number distribution ν among people of same age group
            :ActivityStructure => :Superspreaders, # Frequency of each contact number group
            :Dispersion => :None, # Choice of dispersion factor k for the simulation
            :q => 0.13,           # Fractal index
            :N0 => 100,           # Number of initial infected patients in the population
            :InitCond => :Discrete, # Choice of initial condition for the simulation
            :StepSize => :QuarterDay # Step size
            :Quarantine => :None,  # Uses time-varying contact matrix considering NPIs
            :LossImmRate => :None, # Rate of loss of immunity
            :SeroRevProb => :None, # Probability of seroreversion
            :LossImmProb => :None  # Probability of loss of immunity given that patient seroreverted
  			)
end
```

We shall give details of all these options in the following.





