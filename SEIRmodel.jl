using LinearAlgebra, Match, Dates, Pandas, Distributions, StatsBase
include("ContactMatrices.jl")


⊗ = (x, y) -> kron(x, y)


##
# Chooses the appropriate set of parameters for each simulation
"""
    function caseoptions(Simulation)
Returns a dictionary with the options chosen for a particular simulation.
"""
function caseoptions(Simulation)
    caseopt = @match Simulation begin
        :Manaus_NoAge_Continuous => begin
            Dict(:Model => :SEIRAgeContinuous,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => 3.0,
            :NPI => :None,
            :NPIprediction => :KeepLast,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :Dispersion => :None,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Continuous,
            :StepSize => :QuarterDay)
        end
        :SP_NoAge_EstimatedRt => begin
            Dict(:Model => :SEIRAgeContinuous,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => 1.0,
            :NPI => :RtSP_29out2020,
            :NPIprediction => :KeepLast, # For this simulation, this means no social distancing
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :Dispersion => :None,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Continuous,
            :StepSize => :QuarterDay)
        end

        :Manaus_NoAge_EstimatedRt => begin
            Dict(:Model => :SEIRAgeContinuous,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :N0 => 1,
            :R0 => 1.0,
            :NPI => :RtManaus_29out2020,
            :NPIprediction => :KeepLast,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :Dispersion => :None,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Continuous,
            :StepSize => :QuarterDay)
        end

        :SP_NoAge_EstimatedRt_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :Jan1st,
            :R0 => 1.0,
            :NPI => :RtSP_29out2020,
            :NPIprediction => :KeepLast,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :Dispersion => :None,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :Manaus_NoAge_EstimatedRt_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Jan18,
            :R0 => 1.0,
            :NPI => :RtManaus_29out2020,
            :NPIprediction => :KeepLast,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :Dispersion => :None,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :Manaus_Quarantine_1_5_Dispersion_InLoco_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :ManausSchoolClosure,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :N0 => 1,
            :R0 => :AM,
            :NPI => :ManausInLoco,
            :NPIprediction => :Baseline,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :Manaus75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :Small_k,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :SP_NoAge_Dispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
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
            :ContactMatrix => :None,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => 0.2,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :Manaus_NoAge_Dispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :AM,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => 1.5,
            :q => 1.0,
            :N0 => 1,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :SP_NoAge_Dispersion_1_5_NoNPI_Fractal_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
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
            :ContactMatrix => :None,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => 1.5,
            :q => 0.13,
            :N0 => 100,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_Quarantine_LowDispersion_GovSP_Fractal_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :q => 0.13,
            :N0 => 10,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :SP_Quarantine_LowDispersion_GovSP_LossProb_0_005_SeroRev_Discrete => begin
            Dict(:Model => :SEIRSeroRevDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => 1 / (4*30), # four months
            :SeroRevProb => 0.005,
            :LossImmProb => 1.0,
            :N0 => 100,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :NPIprediction => :Baseline,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_Quarantine_LowDispersion_GovSP_LossProb_0_1_SeroRev_Discrete => begin
            Dict(:Model => :SEIRSeroRevDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => 1 / (6*30), # one year
            :SeroRevProb => 1 / 10,
            :LossImmProb => 1.0,
            :N0 => 100,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :NPIprediction => :Baseline,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :SP_Quarantine_LowDispersion_GovSP_High_SeroRev_Discrete => begin
            Dict(:Model => :SEIRSeroRevDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => 1 / 365, # one year
            :SeroRevProb => 1 / 2,
            :LossImmProb => 1.0,
            :N0 => 100,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_Quarantine_LowDispersion_GovSP_Slow_SeroRev_Discrete => begin
            Dict(:Model => :SEIRSeroRevDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => 1 / (8*30), # Eight months
            :SeroRevProb => 1 / 500,
            :LossImmProb => 1.0,
            :N0 => 100,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_Quarantine_LowDispersion_GovSP_SeroRev_Discrete => begin
            Dict(:Model => :SEIRSeroRevDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => 1 / (4*30), # Four months
            :SeroRevProb => 1 / 10_000,
            :LossImmProb => 1.0,
            :N0 => 100,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_Quarantine_LowDispersion_GovSP_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :SP_SchoolClosureFullReturn,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :N0 => 10,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilSeparate,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :SP_LowDispersion_GovSP_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPGov,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :SP_Age_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :Dispersion => :None,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_LowDispersion_InLoco_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPInLoco,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :Manaus_1_5_Dispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :N0 => 1,
            :R0 => :AM,
            :NPI => :None,
            :NPIprediction => :Baseline,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :Manaus75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :Small_k,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :Manaus_1_5_Dispersion_InLoco_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :AM,
            :NPI => :ManausInLoco,
            :NPIprediction => :Baseline,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :Manaus75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :Small_k,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SP_LowDispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SPDispersion_InLoco_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :SPInLoco,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SPDispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :SP,
            :FirstDay => :SP,
            :R0 => :SP,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :SP75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end



        :ManausDispersion_InLoco_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :AM,
            :NPI => :ManausInLoco,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :Manaus75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :Manaus_LowDispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :AM,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :Manaus75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020_max,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
        :ManausDispersion_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :AM,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :Manaus75,
            :ContactMatrix => :BrazilFull,
            :ActivityVector => :Superspreaders,
            :ActivityStructure => :Superspreaders,
            :Dispersion => :COVID_Endo2020,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :SPHeter_SPGov_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
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
            :ContactMatrix => :BrittonScience2020,
            :ActivityVector => :Masks,
            :ActivityStructure => :Masks,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end

        :ManausHomog_NoNPI_Discrete => begin
            Dict(:Model => :SEIRDiscrete,
            :Quarantine => :None,
            :LossImmRate => :None,
            :SeroRevProb => :None,
            :LossImmProb => :None,
            :Population => :Manaus,
            :FirstDay => :Manaus,
            :R0 => :AM,
            :NPI => :None,
            :RecoveryRate => :BrittonScience2020,
            :IncubationRate => :BrittonScience2020,
            :AgeStructure => :None,
            :ContactMatrix => :None,
            :ActivityVector => :None,
            :ActivityStructure => :None,
            :InitCond => :Discrete,
            :StepSize => :QuarterDay)
        end
            :ManausHomog_InLoco_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :Quarantine => :None,
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
                :ContactMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end


            :SPHomog_InLoco_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :Quarantine => :None,
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
                :ContactMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end

            :SPHomog_NPI_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :Quarantine => :None,
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
                :ContactMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end
            :SPHomog_noNPI_Discrete => begin
                Dict(:Model => :SEIRDiscrete,
                :Quarantine => :None,
                :LossImmRate => :None,
                :SeroRevProb => :None,
                :LossImmProb => :None,
                :Population => :SP,
                :FirstDay => :SP,
                :R0 => :SP,
                :NPI => :None,
                :Dispersion => :None,
                :RecoveryRate => :BrittonScience2020,
                :IncubationRate => :BrittonScience2020,
                :AgeStructure => :None,
                :ContactMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end
            _ => begin
                println("No option set defined for $(Simulation). Using SP no NPI, no age, no dispersion.")
                Dict(:Model => :SEIRDiscrete,
                :Quarantine => :None,
                :LossImmRate => :None,
                :SeroRevProb => :None,
                :LossImmProb => :None,
                :Population => :SP,
                :FirstDay => :SP,
                :Dispersion => :None,
                :R0 => :SP,
                :NPI => :None,
                :RecoveryRate => :BrittonScience2020,
                :IncubationRate => :BrittonScience2020,
                :AgeStructure => :None,
                :ContactMatrix => :None,
                :ActivityVector => :None,
                :ActivityStructure => :None,
                :InitCond => :Discrete,
                :StepSize => :QuarterDay)
            end
    end
    return caseopt
end


"""
    function secondary80(R0, k)
Returns the fraction `p80` of infected responsible for 80% of secondary transmissions, given values of the reproduction number `R0` and the dispersion factor `k`, following [Endo 2020].
"""
function secondary80(R0, k)
    X = 0
    NB = NegativeBinomial(k, k/(R0+k))
    S = pdf(NB, X)
    Me = 0*S / R0
    while Me < 0.2
        X += 1
        prob = pdf(NB, X)
        Me += X * prob / R0
        S += prob
    end
    return 1 - S
end

"""
    function searchdate(t, d0)
Returns the date corresponding to a given time `t` in the simulation.  `d0` corresponds to `t = 0`.
"""
function searchdate(t, d0)
    return d0 + Day(floor(Int,t))
end
"""
    function InstantaneousContactMatrix(caseopt, A, α)
Returns a function that computes `A(t)'*i`, whene `A(t)' is the instantaneous contact matrix for a given time `t`, and `i` is the vector of infected at time `t`.  Assumes `A` is a vector of matrices.  If `length(A) > 1`, assumes the order of matrices is:
    `A[1] = Ahome`, `A[2] = Aschool`, `A[3] = Awork`, `A[4] = Aothers`. Depending on the case, the last matrices may be added together, for example, if one does not wish to consider the contact matrices for work and others separately, the third entry of `A` could be `Awork + Aothers`, whithout a fourth entry.
Applies NPI factor using function `α` (see  function `socialdistancing`).
"""
function InstantaneousContactMatrix(caseopt, A, α)
    d0 = firstday(caseopt[:FirstDay])
    @match caseopt[:Quarantine] begin
        :ManausSchoolClosure => begin
        A[2] .= A[1] .+ A[2]
        q = get(caseopt, :q, 1.0)
        i0q = get(caseopt,:N0, 1.0)^(1-q)
        (t, i) ->
        begin
            d = searchdate(t, d0)
            i0qα = i0q * α(t)
            if d < Date("2020-03-16") || d >= Date("2020-08-10")
                return (i0qα .*  (A[3] * i)) .+ i0q .* (A[2] * i)
            else
                return (i0qα .* (A[3] * i)) .+ i0q .* (A[1] * i)
            end
        end


    end
        :SP_SchoolClosureFullReturn => begin
                                A[2] .= A[1] .+ A[2]
                                q = get(caseopt, :q, 1.0)
                                i0q = get(caseopt,:N0, 1.0)^(1-q)
                                (t, i) ->
                                begin
                                    d = searchdate(t, d0)
                                    if d < Date("2020-03-24") || d >= Date("2020-10-05")
                                        return ((i0q * α(t)) .*  (A[3] * i)) .+ i0q .* (A[2] * i)
                                    else
                                        return ((i0q * α(t)) .* (A[3] * i)) .+ i0q .* (A[1] * i)
                                    end
                                end
                            end
        :None => (t, i) -> α(t) .* (sum(A) * i)
         _ =>  (t, i) -> α(t) .* (sum(A) * i) # Constant contact matrix
    end
end


"""
    DistributePopulation(Pop, πv)
Returns a vector with the total population `Pop` distributed according to the relative frequencies given in `πv`, making sure that the sum of the entries is indeed equal to `Pop`.
"""
function DistributePopulation(Pop, πv)
    x = zeros(Int, length(πv))
    x = round.(Int,Pop*πv) # segfault if .*πv
    while sum(x) != Pop
       srand = rand(1:length(πv))
       x[srand] += sum(x) < Pop ? 1 : -1
       if x[srand] < 0
           x[srand] = 0
       end
    end
    return x
end

"""
    dispersion_factor(case)
Returns the dispersion paramater `k` for superspreader models.
"""
function dispersion_factor(case)
    @match case begin
        n::Real => n
        :SARS_Lloyd2005 => 0.15 # J. O. Lloyd-Smith, S. J. Schreiber, P. E. Kopp, and W. M. Getz, “Superspreading and the effect of individual variation on disease emergence,” Nature, vol. 438, no. 7066, Art. no. 7066, Nov. 2005, doi: 10.1038/nature04153.

        :COVID_Endo2020 => 0.10 # A. Endo, Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group, S. Abbott, A. J. Kucharski, and S. Funk, “Estimating the overdispersion in COVID-19 transmission using outbreak sizes outside China,” Wellcome Open Res, vol. 5, p. 67, Jul. 2020, doi: 10.12688/wellcomeopenres.15842.3.
        :COVID_Endo2020_max => 0.2 # Upper limit of CrI
        :COVID_Endo2020_min => 0.04 # Lower limit of CrI
        :Small_k => 1.5
        :None => Inf
        _ => Inf # No dispersion
    end
end

"""
    function discretestep(case)
Returns the stepsize (in days) for the discrete-time simulations.
"""
function discretestep(case)
    @match case begin
        Δt::Number => Δt
        :QuarterDay => 0.25
        :HalfDay => 0.5
        :Day => 1.0
        _ => begin
                println("Choice for :StepSize not found")
                0.25
            end
    end
end

"""
    function lossofimmunityrate(case)
Returns the rate at which patients lose immunity (simple seroreversion model).
"""
function lossofimmunityrate(case)
    @match case begin
        n::Real => n
        :None  => 0.0
        :PrevalenceSP => (0.8865, 0.65)
        :PrevalenceManaus => (0.7886, 0.85)
        :PrevalenceMin => 0.8795
        _ => begin
                println("Choice for :LossImmRate not found")
                0.0
            end
    end
end
"""
    function seroreversionprobability(case)
Returns the probability of seroreversion.
"""
function seroreversionprobability(case)
    @match case begin
        n::Real, if 0 ≤ n ≤ 1 end => n
        :None => 0.0
        :PrevalencePaper => 1.0
        _ => begin
                println("Choice for :SeroRevProb not found")
                0.0
            end
    end
end

"""
    function lossofimmunityprobability(case)
Returns the probability of loss of immunity.
"""
function lossofimmunityprobability(case)
    @match case begin
        n::Real, if 0 ≤ n ≤ 1 end => n
        :None => 0.0
        :TenPercent => 0.1
        :TwentyPercent => 0.2
        _ => begin
                println("Choice for :LossImmProb not found")
                0.0
            end
    end
end

"""
    function totalpopulation(case)
Returns the total population for each simulation.
"""
function totalpopulation(case)
    @match case begin
        n::Integer => n
        :BrittonScience2020   =>  1 # Normalized continuous-time model
        :SP =>  12_325_232 # https://agenciadenoticias.ibge.gov.br/agencia-sala-de-imprensa/2013-agencia-de-noticias/releases/28668-ibge-divulga-estimativa-da-populacao-dos-municipios-para-2020
        :Manaus => 2_219_580 # Como acima
        _ => begin
                println("Choice for :Population not found")
                7_594_000_000
            end
    end
end

"""
    function firstday(case)
Returns the day of first case for a given case.
"""
function firstday(case)
    @match case begin
        d::Date => d
        :Jan1st => Date("2020-01-01")
        :Jan15 => Date("2020-01-015")
        :Jan18 => Date("2020-01-018")
        :Feb1st => Date("2020-02-01")
        :Feb18 => Date("2020-02-18")
        :BrittonScience2020 => Date("2019-12-01")
        :Brazil       ||
        :SP    =>  Date("2020-02-25")
        :Manaus => Date("2020-03-13")
        _       => begin
                        println("Choice for :FirstDay not found")
                        Date("2019-12-01")
                    end
    end
end


"""
    function searchsocialdistancing(t, d, d0, αd)
Returns the value of the social distancing index (essentially, the reduction in R0) for a given day.  `t` is a continuous variable, `d` is a vector of dates for which social distancing indices are available, `d0` is the first day, `αd` is a vector of same length as `d` with the social distancing indices.
"""
function searchsocialdistancing(
    t,
    d,
    d0,
    αd,
    endoption = :Baseline,
    type = :SocialDistancing,
)
    @match type begin
        :Rt => begin
            if d[1] > d0 + Day(floor(Int, t)) # before first day
                return 3.0
            elseif d[end] < d0 + Day(floor(Int, t)) # There are many choices of what to do after the last value of social distancing index.
                @match endoption begin
                    :Baseline => return 3.0 # No social distancing after last data
                    :KeepLast => return αd[end] # keep last value
                    :Lockdown => return 0.0
                    n::Real,
                    if 0 ≤ endoption
                    end => αd[end] - (t - (d[end] - d0).value) * endoption # gradual return to baseline.
                end
            else
                n0 = findfirst(x -> x >= d0 + Day(floor(Int, t)), d)
                return αd[n0]
            end
        end
        _ || :SocialDistancing => begin
            if d[1] > d0 + Day(floor(Int, t)) # before first day
                return 1.0
            elseif d[end] < d0 + Day(floor(Int, t)) # There are many choices of what to do after the last value of social distancing index.
                @match endoption begin
                    :Baseline => return 1.0 # No social distancing after last data
                    :KeepLast => return min(1.0 - αd[end], 1.0) # keep last value
                    :Lockdown => return 0.0
                    n::Real,
                    if 0 ≤ endoption
                    end =>
                        min(1.0 - αd[end] + (t - (d[end] - d0).value) * endoption, 1.0) # gradual return to baseline.
                end
            else
                n0 = findfirst(x -> x >= d0 + Day(floor(Int, t)), d)
                return 1.0 - αd[n0]
            end
        end
    end
end
"""
    function Rzero(case)
Returns the value of R0 used in each simulation.  The values here for the states are from William's paper.
"""
function Rzero(case)
    @match case begin
        n::Number        => n
        :BrittonScience2020 || :LewisPaper2020   =>  2.5
        :Brazil          =>  3.1
        :SP              =>  2.9
        :RJ              =>  2.9
        :AM              =>  2.6
        :CE              =>  1.9
        _                =>  begin # Average estimate in the world.
                                println("Choice for :R0 not found")
                                2.5
                            end
    end
end

"
    function socialdistancing(case)
Returns a function to compute the social index index.  Some options read data obtained from different sources.
`case = :None` -> α = 1.0
"
function socialdistancing(caseopt)
    @match caseopt[:NPI] begin
        :None   =>  (t -> 1.0) # constant equal to 1.0
        :SPGovRaw   => begin # From https://www.saopaulo.sp.gov.br/coronavirus/isolamento
        # Remember to remove first blank row, if necessary.
            df = read_excel("NPI_Data/IsolamentoGovSP15out2020.xlsx")
            dfnames = columns(df)
            d = [Date(dfnames[i]*"20", "dd/mm/yyyy") for i in 5:size(df)[2]]
            #d = [Date(dfnames[i][end-4:end]*"/2020", "dd/mm/yyyy") for i in 4:size(df)[2]]
            αd = Array(query(df, :(Município1 == "SÃO PAULO")))[5:end]
            # basal level -> minimum of period immediately before NPIs implemented.
            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), αd, get(caseopt, :NPIprediction, :Baseline))
        end
        :SPGov   => begin # From https://www.saopaulo.sp.gov.br/coronavirus/isolamento
        # Remember to remove first blank row, if necessary.
            df = read_excel("NPI_Data/IsolamentoGovSP15out2020.xlsx")
            dfnames = columns(df)
            d = [Date(dfnames[i]*"20", "dd/mm/yyyy") for i in 5:size(df)[2]]
            #d = [Date(dfnames[i][end-4:end]*"/2020", "dd/mm/yyyy") for i in 4:size(df)[2]]
            αd = Array(query(df, :(Município1 == "SÃO PAULO")))[5:end]
            # basal level -> minimum of period immediately before NPIs implemented.
            basal_level = minimum(αd[1:4])
            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), (αd .- basal_level)./(1-basal_level))
        end
        :SPInLocoRaw => begin

                    df = read_csv("NPI_Data/InLoco-isolation_Manaus_SaoPaulo.csv")
                    dfnames = columns(df)
                    d = Date.(Array(query(df, :(state_name == "São Paulo" && city_name == "São Paulo"))[[:dt]]))
                    αd = Array(query(df, :(state_name == "São Paulo" && city_name == "São Paulo"))[[:isolated]])
                    # Minimum index for the period immediately before isolation begins

                    return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), αd, get(caseopt, :NPIprediction, :Baseline))
                end
                :ManausInLocoRaw => begin
                    df = read_csv("NPI_Data/InLoco-isolation_Manaus_SaoPaulo.csv")
                    dfnames = columns(df)
                    d = Date.(Array(query(df, :(state_name == "Amazonas" && city_name == "Manaus"))[[:dt]]))
                    αd = Array(query(df, :(state_name == "Amazonas" && city_name == "Manaus"))[[:isolated]])

                    return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), αd, get(caseopt, :NPIprediction, :Baseline))
                end
        :SPInLoco => begin

            df = read_csv("NPI_Data/InLoco-isolation_Manaus_SaoPaulo.csv")
            dfnames = columns(df)
            d = Date.(Array(query(df, :(state_name == "São Paulo" && city_name == "São Paulo"))[[:dt]]))
            αd = Array(query(df, :(state_name == "São Paulo" && city_name == "São Paulo"))[[:isolated]])
            # Minimum index for the period immediately before isolation begins
            basal_level = minimum(αd[1:40])

            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), (αd .- basal_level)./(1-basal_level))
        end
        :ManausInLoco => begin
            df = read_csv("NPI_Data/InLoco-isolation_Manaus_SaoPaulo.csv")
            dfnames = columns(df)
            d = Date.(Array(query(df, :(state_name == "Amazonas" && city_name == "Manaus"))[[:dt]]))
            αd = Array(query(df, :(state_name == "Amazonas" && city_name == "Manaus"))[[:isolated]])
            basal_level = minimum(αd[1:40])
            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), (αd .- basal_level) ./ (1-basal_level))
        end
        :RtManaus_29out2020 => begin
            df = read_csv("Rt/Rt ManausTeste.csv")
            dfnames = columns(df)
            d = Date.(Array(df["dates"]))
            αd = Array(df["Rt"])
            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), #fill(1.0-3.0*0.25,length(d)),
                αd,
                         get(caseopt,:NPIprediction,:KeepLast), :Rt)
        end
        :RtSP_29out2020 => begin
            df = read_csv("Rt/Rt São PauloTeste.csv")
            dfnames = columns(df)
            d = Date.(Array(df["dates"]))
            αd = Array(df["Rt"])
            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), #fill(1.0-3.0*0.25,length(d)),
                αd,
                         get(caseopt,:NPIprediction,:KeepLast), :Rt)
        end
        _ => begin
                println("Choice for :NPI not found.")
                t -> 1.0
            end
    end
end

"
    function recoveryrate(case)
Returns the inverse of the recovery time constant, depending on the value of `case`.
"
function recoveryrate(case)
    @match case begin
        :BrittonScience2020 =>  1/4
        :DaviesLancet2020   => 1/5 # Adding clinical and subclinical average times
        _ => begin
                    println("Choice for :RecoveryRate not found")
                    1/5
                end
    end
end

"
    function incubationrate(case)
Returns the inverse of the incubation time constant, depending on the value of `case`.
"
function incubationrate(case)
    @match case begin
        :BrittonScience2020 => 1/3
        :DaviesLancet2020 => 1/4
        _ => begin
                    println("Choice for :IncubationRate not found")
                    1/4
                end
    end
end
"""
    function age_structure(case)
Returns a vector with the probabilities for each age range.  If using Britton's contact matrix, the ranges are [0-5, 6-12, 13-19, 20-39, 40-59, 60-].
Returns 1.0 if age structure is not included in the simulation.
"""
function age_structure(case)
    πv = @match case begin
        :BrittonScience2020   =>   [0.0725, 0.0866, 0.1124, 0.3323, 0.2267, 0.1695] #[0-5, 6-12, 13-19, 20-39, 40-59, 60-]
        :Brazil       =>   [0.087758779
                            0.116819264
                            0.125311918
                            0.335384766
                            0.226792041
                            0.107933238] # From https://sidra.ibge.gov.br/tabela/2020",0, whole country
        :SP60 || :SP   =>  [0.076912061
                             0.100387104
                             0.105932843
                             0.351962947
                             0.24768619
                             0.117119] # Only SP city
        :RMSP75        => [102897 # 0 a 4
                           105038 # 5 a 9
                           120339 # 10 a 14
                           76399+52161 # 15 a 19
                           149348 # 20 a 24
                           141394 # 25 a 29
                           120163 # 30 a 34
                           101333 # 35 a 39
                           88757 # 40 a 44
                           72071 # 45 a 49
                           58520 # 50 a 54
                           44413 # 55 a 59
                           55413 # 60 a 64
                           23519 # 65 a 69
                           17552 # 70 a 74
                           11671+11108+2356+248] # 75+
        :SP75          => [1307653 # 0
                           1398122 # 5
                           1611871 # 10
                           938835+615310 # 15
                           1758995 # 20
                           1865192 # 25
                           1752736 # 30
                           1557670 # 35
                           1431684 # 40
                           1289237 # 45
                           1138259 # 50
                           920705 # 55
                           1180237 # 60
                           484906 # 65
                           368478 # 70
                           256250+253752+37488+1501] # 75+
        :Manaus75        =>  # Censo 2010 - https://sidra.ibge.gov.br/tabela/1378
                                [196298 # 0 a 4
                                       202326 # 5 a 9
                                       218668 # 10 a 14
                                       128714+79586 # 15 a 19
                                       208301 # 20 a 24
                                       210450 # 25 a 29
                                       190581 # 30 a 34
                                       157532 # 35 a 39
                                       130470 # 40 a 44
                                       105759 # 45 a 49
                                       85616 # 50 a 54
                                       63085 # 55 a 59
                                       75632 # 60 a 64
                                       30993 # 65 a 69
                                       22299 # 70 a 74
                                       14626+13606+2540+233] # 75+

        :BrazilFull     =>  readpopulationstructure("Brazil")
        :None           => 1.0
        _               => begin
                                println("Choice for :AgeStructure not found")
                                1.0
                            end
     end
     return(πv ./ sum(πv))
 end

"""
    function ContactMatrix(case)
Returns the contact matrix between ages for each case, or 1.0 if no age structure.  Note that this matrix will be further normalized to the desired R0.
"""
function ContactMatrix(case)
     @match case begin
        :BrittonScience2020  => begin
                            A = [2.2257 0.4136 0.2342 0.4539 0.2085 0.1506;
                            0.4139 3.6140 0.4251 0.4587 0.2712 0.1514;
                            0.2342 0.4257 2.9514 0.6682 0.4936 0.1972;
                            0.4539 0.4592 0.6676 0.9958 0.6510 0.3300;
                            0.2088 0.2706 0.4942 0.6508 0.8066 0.4341;
                            0.1507 0.1520 0.1968 0.3303 0.4344 0.7136]
                            ((A+A')./2,)
                        end
        :BrazilFull => begin
        # From
        #K. Prem et al., “Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era,” medRxiv, p. 2020.07.22.20159772, Jul. 2020, doi: 10.1101/2020.07.22.20159772.
        # Age groups: 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75+
                            Afull, Ahome, Aschool, Awork, Aothers = readcontactmatrices("BRA")
                            return (((Afull+Afull'))./2,)
                        end
        :BrazilSeparate => begin
                            Afull, Ahome, Aschool, Awork, Aothers = readcontactmatrices("BRA")
                            return (Ahome+Ahome')/2, (Aschool+Aschool')./2,((Awork+Aothers+(Awork+Aothers)'))./2
                        end

        :None => begin
                    A = Array{Float64,2}(undef, 1, 1)
                    A[1,1] = 1.0
                    return (A,)
                end
        _ => begin
                    println("Choice for :ContactMatrix not found")
                    A = Array{Float64,2}(undef, 1, 1)
                    A[1,1] = 1.0
                    return (A,)
                end
    end
end
"""
    function activity_vector(caseopt)
Returns the vector with the increase factors in contact (activity level) in heterogeneous models.  For example, a vector `[4.0, 1.0]` means that each age structure and class will have individuals with 4 times more contacts than the others.  The distribution of people in each set is given by `activity_percentages(case)`.
"""
function activity_vector(caseopt)
    @match caseopt[:ActivityVector] begin
        # High, average, low
        :Masks               =>  begin
                                    act = [4.0, 1.0, 0.25]
                                    #act = act / (act⋅activity_percentages(caseopt[:Ac:ActivityStructure]))
                        end
        :BrittonScience2020  =>  begin
                                act = [2.0, 1.0, 0.5]
                        end
        :BrittonNormalized => begin # average is one
                                act = [2.0, 1.0, 0.5]
                                act = act / (act⋅activity_percentages(caseopt[:ActivityStructure]))
                        end
        :SPInLoco       =>  begin # Using Inloco mobility data
                                act = [1.15595, 1.06564, 1.0]
                                act = act / (act⋅activity_percentages(caseopt[:ActivityStructure]))
                        end
        :Superspreaders => begin
                            k = dispersion_factor(caseopt[:Dispersion])

                            Γ = Gamma(k,Rzero(caseopt[:R0])/k) #Gamma(α, θ) # Gamma distribution
                            Pop = totalpopulation(caseopt[:Population])
                            agefreqs = age_structure(caseopt[:AgeStructure])
                            Nages = length(agefreqs)
                            x = DistributePopulation(Pop, agefreqs)
                            R0values = Vector{Vector{Int64}}(undef,Nages)
                            R0freqs = Vector{Vector{Float64}}(undef,Nages)
                            Nreal = 0
                            for j in 1:Nages
                                realizations = round.(Int, rand(Γ, x[j])) # R0 values, rounded to closest integer to reduce size of matrices
                                propmap = proportionmap(realizations)

                                perm = sortperm(collect(keys(propmap)))
                                R0values[j] = collect(keys(propmap))[perm] # R0 values that were sampled
                                R0freqs[j] = collect(values(propmap))[perm] # corresponding frequencies
                                Nreal += length(R0values[j]) # number of entries
                            end
                            return R0values,R0freqs,Nreal

                        end



        :None         => [1.0]
        _             => begin
                            println("Choice for :ActivityVector not found")
                            [1.0]
                        end
    end
end

"""
    function activity_percentages(case)
Returns the vector with the percentages of each different contact (activity) level group in contact in heterogeneous models.  See also `activity_vector`.
This function is not used if `case == :Superspreaders`.
"""
function activity_percentages(case)
    @match case begin
            :BrittonScience2020       =>  [0.25, 0.5, 0.25] # High, average, low
            :SPInloco     =>  [0.7, 0.2, 0.1] # From Inloco mobility data
            :Masks        =>  [0.3, 0.2, 0.5]
            :None         =>  [1.0] # homogeneous
            _             => begin
                                println("Choice for :ActivityStructure not found")
                                [1.0]
                            end

    end
end


"""
    function πvector(caseopt, actprobs)
Computes the complete vector of age+contact level groups (a Kronecker product of the age and contact level structures) when `:ActivityVector = :Superspreders`.
`actprobs` is the list of probabilities for each activity level, returned by `Amatrix`.
"""
function πvector(caseopt, actprobs)
    @match caseopt[:ActivityVector] begin
                :Superspreaders => begin
                        agestr = age_structure(caseopt[:AgeStructure])
                        reduce(vcat, [actprobs[i] * agestr[i] for i in 1:length(agestr)]) # Joins result into a single vector
                    end
    end
end

"""
    function πvector(caseopt)
Computes the complete vector of age+contact level groups (a Kronecker product of the age and contact level structures).
"""
function πvector(caseopt)
    @match caseopt[:ActivityVector] begin
                :Superspreaders => begin
                    println("If the `:Superspreaders` model is being used, funcion `πvector` should be called with two arguments: `caseopt` and `actprobs`.  `actprobs` is the list of probabilities for each activity level, returned by `Amatrix`.")
                    [1.0]
                end

                _ => kron(age_structure(caseopt[:AgeStructure]),
                activity_percentages(caseopt[:ActivityStructure]))
            end
end

"""
    function Amatrix(caseopt)
Computes the contact matrix considering age and contact level groups.
In the case of models considering superspreaders, returns also a vector of vectors actprobs, such that `actprobs[agegroup]` is the list of relative frequencies for the activity levels sampled for `agegroup`.
"""
function Amatrix(caseopt)

    Avec = ContactMatrix(caseopt[:ContactMatrix])
    Nmatrices = length(Avec)
    A = Vector{Array{Float64,2}}(undef,Nmatrices)

    @match caseopt[:ActivityVector] begin
        :Superspreaders => begin
                            Na = size(Avec[1])[1]
                            actvalues,actprobs,Nact = activity_vector(caseopt)
                            Nk = maximum(length.(actvalues))
                            actprobsext = Vector{Vector{Float64}}(undef, Na)
                            actvaluessqrt = Vector{Vector{Float64}}(undef, Na)
                            for ia in 1:Nmatrices
                                A[ia] = zeros(Na*Nk, Na*Nk)
                                Atemp = Avec[ia]
                                for r in 1:Na
                                    Nkr = length(actvalues[r])
                                    actprobsext[r] = vcat(actprobs[r],zeros(Nk-Nkr))
                                    # This option assumes superspreaders are not supersusceptible
                                    #A[1+(r-1)*Nk:r*Nk, :] .= Atemp[r,:]' ⊗ (vcat(actvalues[r],zeros(Nk-Nkr)) * ones(1, Nk))
                                    # This options assumes superspreaders are also supersusceptible
                                    actvaluessqrt[r] = vcat(sqrt.(actvalues[r]),zeros(Nk-Nkr))
                                end
                                for r in 1:Na
                                    for l in 1:Na
                                        A[ia][1+(r-1)*Nk:r*Nk, 1+(l-1)*Nk:l*Nk] .= Atemp[r,l] .* (actvaluessqrt[r] * actvaluessqrt[l]')
                                    end
                                end
                            end
                            return A,actprobsext
                        end
                        _ =>    begin
                            for ia in 1:length(Avec)
                                A[ia] = Avec[ia] ⊗ (activity_vector(caseopt)*activity_vector(caseopt)')
                            end
                            return A,activity_vector(caseopt)
                        end
    end
end
"""
    function NormalizationFactor(caseopt, A, πvec)
Computes the normalization factor `κ` such that the maximum eigenvalue of `κ*A*πvec` is equal to `Rzero(caseopt[:R0])`.
"""
function NormalizationFactor(caseopt, A, πvec)
    μ = recoveryrate(caseopt[:RecoveryRate])
    sπ = diagm(sqrt.(πvec))  # Recall that λ(A*D) = λ(D^(1/2)*A*D^(1/2)) if D is diagonal.
    eig = eigmax(Symmetric(sπ*A*sπ)) ./ μ # Only works if A is symmmetric, otherwise replace by maximum(abs.(eigvals(A*diagm(πvec))))
    κ = Rzero(caseopt[:R0]) / maximum(real.(eig))
    return κ
end

"""
    function NormalizedContactMatrix(caseopt)
Returns the normalized contact matrix used in the SEIR model (corresponding to the `β` parameter in the standard, unstructured model).  The matrix is normalized so that the corresponding `R0` is the one returned by `Rzero(caseopt[:R0])`.
"""
function NormalizedContactMatrix(caseopt)
    A,πvec = @match caseopt[:ActivityVector] begin
        :Superspreaders => begin
            A,actprobs = Amatrix(caseopt)
            πvec = πvector(caseopt, actprobs)
            κ = NormalizationFactor(caseopt, sum(A), πvec)
            κ .* A,πvec
        end
        _ => begin
            πvec = πvector(caseopt)
            A,temp = Amatrix(caseopt)
            κ = NormalizationFactor(caseopt, A[1], πvec)
            κ .* A,πvec
        end
    end
    return A,πvec
end

"""
    function indices(caseopt, πvec)
Returns the indices for the susceptible (`:S`), exposed (`:E`), infected (`:I`) and recovered/dead (`:R`) entries in the state vector.  Models with more compartments will include other variables (e.g., seroreversion, loss of immunity).
`caseopt` defines the set of parameters to be used.
`πvec` contains the proportions for the population in each `S` compartment.
"""
function indices(caseopt, πvec)
    # Number of compartments in the model
    compartments = 4 # S-E-I-R
    compartments += @match caseopt[:SeroRevProb] begin # Check if Seroreversion is possible
                :None => begin
                            compnames = (:S, :E, :I, :R)
                            0 # No seroreversion, no loss of immunity - S-E-I-R
                        end
                _ => begin @match caseopt[:LossImmProb] begin
                                :None => begin
                                        compnames = (:S, :Sr, :E, :I, :Rsp, :Rsn)
                                        1 # S-E-I-Rseropositive-Rseroneg (Rseroneg = fraction of R that seroreverts)
                                    end
                                _ => begin
                                        compnames = (:S, :Sr, :E, :I, :Rsp, :Rsn, :Rl)
                                        3 # S-Srecovered-E-I-Rseropositive-Rseroneg-Rloss (Srecovered - susceptible that was infected at least once, Rloss - Recovered that will lose immunity)
                                    end
                            end
                    end
    end
    Nlevels = length(πvec)
    return Dict(compnames[i] => (i-1)*Nlevels+1:i*Nlevels for i in 1:compartments)
end

"""
    function initcond(caseopt,πv,sind,eind)
Returns a vector with the initial conditions for each case.
"""
function initcond(caseopt,πv,sind,eind)
    @match caseopt[:InitCond] begin
        :Continuous        => begin
                                ϵ = caseopt[:N0] / totalpopulation(caseopt[:Population])
                                x = zeros(length(indices(caseopt,πv))*length(sind))
                                x[eind] .= ϵ .* πv
                                x[sind] .= (1-ϵ) .* πv
                                return x
                            end
        :BrittonScience2020  =>  begin
                                ϵ = 1e-4
                                x = zeros(length(indices(caseopt,πv))*length(sind))
                                x[eind] .= ϵ .* πv
                                x[sind] .= (1-ϵ) .* πv
                                return x
                            end
        :Discrete => begin
                            N0 = get(caseopt, :N0, 1) # No. of seed cases
                            Pop = totalpopulation(caseopt[:Population])
                            x = zeros(Int, length(indices(caseopt,πv))*length(sind))
                            x[sind] .= DistributePopulation( Pop, πv)
                            pacient0 = rand(DiscreteNonParametric(eind,πv), N0)
                            for pac = 1:N0
                               x[pacient0[pac]] += 1
                            end
                            x[sind] .-= x[eind]
                            for k = 1:length(x)
                               if x[k] < 0
                                   x[k] = 0
                               end
                            end

                            return x
                        end
        _ => begin
                println("Choice for :InitCond not found")
                zeros(x = zeros(Int, length(indices(caseopt,πv))*length(sind)))
            end

    end
end


"""
    function seirseroreversion!(du, u, p, t)
Continuous,time model with seroreversion.
 `du` - derivative for state `u` at time `t`. `p` -> parameter vector:
`p[:S]: sind`
`p[:E]: eind`
`p[:I]: iind`
`p[:Rsp]: rspind` - fraction of R that does not serorevert
`p[:Rsn]: rslind` - fraction of R that seroreverts, but remains immune
`p[:Rsl]: rslind` - fraction of R that seroreverts and becomes susceptible again
`p[:Sr]: srind`   - fraction of S that has been infected at least once
"""
function seirseroreversion!(du, u, p, t)

    sind = p[:S]
    eind = p[:E]
    iind = p[:I]
    rspind = p[:Rsp]      # fraction of R that does not serorevert
    rsnind = p[:Rsn]     # fraction of R that seroreverts, but remains immune
    rlsind = p[:Rl]     # fraction of R that loses immunity
    srind = p[:Sr] # R that becomes again S
    A = p[:ContactMatrix]
    σ = p[:IncubationRate]
    μ = p[:RecoveryRate]
    ϵ0 = p[:SeroRevProb]
    p0 = p[:LossImmProb]
    ar = p[:LossImmRate]
    #α = p[:NPI]

    s = view(u, sind)
    sr = view(u, srind)
    e = view(u, eind)
    i = view(u, iind)
    rsp = view(u, rspind)
    rsn = view(u, rsnind)
    rls = view(u, rlsind)
    tempAi = A(t,i)
    tempse = tempAi .* s
    tempsre = tempAi .* sr
    tempei = σ .* e
    tempir = μ .* i
    temprs = ar .* rls

    du[sind] .=  -tempse
    du[srind] .= temprs .- tempsre
    du[eind] .=  tempse .+ tempsre .- tempei
    du[iind] .= tempei .- tempir
    du[rspind] .= (1-ϵ0)*tempir
    du[rsnind] .= ϵ0*(1-p0)*tempir
    du[rslind] .= ϵ0*p0*tempir .- temprs

end

"""
    function seir!(du, u, p, t)
Continuous-time model without seroreversion.
`du` - derivative for state `u` at time `t`. `p` -> parameter vector.
`p[:S]: sind`
`p[:E]: eind`
`p[:I]: iind`
`p[:R]: rind`
"""
function seir!(du, u, p, t)

    sind = p[:S]
    eind = p[:E]
    iind = p[:I]
    rind = p[:R]
    A = p[:ContactMatrix]
    σ = p[:IncubationRate]
    μ = p[:RecoveryRate]
    #α = p[:NPI]

    s = view(u, sind)
    e = view(u, eind)
    i = view(u, iind)
    r = view(u, rind)
    tempse = (A(t,i) .* s)
    tempei = σ .* e
    tempir = μ .* i

    du[sind] .=  -tempse
    du[eind] .=  tempse .- tempei
    du[iind] .= tempei .- tempir
    du[rind] .= tempir

end
"""
    function model(caseopt)
Returns the model for a particular simulation (both continuous or discrete-time).
"""
function model(caseopt)
    #println("Model: ", caseopt[:Model])

    @match caseopt[:Model] begin
        :SEIRAgeContinuous    =>  begin
                μ = recoveryrate(caseopt[:RecoveryRate])
                σ = incubationrate(caseopt[:IncubationRate])
                α = socialdistancing(caseopt)

                A,πv = NormalizedContactMatrix(caseopt)
                funcA = InstantaneousContactMatrix(caseopt, A, α)
                Indices = indices(caseopt,πv)
                u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
                p = merge(Indices, Dict(:ContactMatrix => funcA, :IncubationRate => σ, :RecoveryRate => μ, :NPI => α))
                return seir!,πv,A,p,u0
        end
        :SeroRevContinuous => begin
            μ = recoveryrate(caseopt[:RecoveryRate])
            σ = incubationrate(caseopt[:IncubationRate])
            α = socialdistancing(caseopt)
            ar = lossofimmunityrate(caseopt[:LossImmRate])
            ϵ0 = seroreversionprobability(caseopt[:SeroRevProb])
            p0 = lossofimmunityprobability(caseopt[:LossImmProb])
            A,πv = NormalizedContactMatrix(caseopt)
            funcA = InstantaneousContactMatrix(caseopt, A, α)

            #println(maximum(real.(eig))*κ)
            Indices = indices(caseopt,πv)
            u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
            p = merge(Indices, Dict(:ContactMatrix => funcA, :IncubationRate => σ, :RecoveryRate => μ, :SeroRevProb => ϵ0, :LossImmProb => p0, :LossImmRate => ar, :NPI => α))

            return seirseroreversion!,πv,A,p,u0
        end
        :SEIRSeroRevDiscrete  => begin
            #println("Model: Discrete SEIR with Seroreversion")
            μ = recoveryrate(caseopt[:RecoveryRate])
            σ = incubationrate(caseopt[:IncubationRate])
            α = socialdistancing(caseopt)
            ar = lossofimmunityrate(caseopt[:LossImmRate])
            ϵ0 = seroreversionprobability(caseopt[:SeroRevProb])
            p0 = lossofimmunityprobability(caseopt[:LossImmProb])
            Δt = discretestep(caseopt[:StepSize])
            pop = totalpopulation(caseopt[:Population])

            A,πv = NormalizedContactMatrix(caseopt)
            funcA = InstantaneousContactMatrix(caseopt, A, α)

            Indices = indices(caseopt,πv)
            u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
            p = merge(Indices, Dict(:ContactMatrix => funcA, :IncubationRate => σ, :RecoveryRate => μ, :SeroRevProb => ϵ0, :LossImmProb => p0, :LossImmRate => ar, :Population => pop, :StepSize => Δt, :NPI => α))

            return discretestochasticSEIRSeroRev!,πv,A,p,u0
        end
        :SEIRDiscrete  => begin
            #println("Model: Discrete SEIR")
            μ = recoveryrate(caseopt[:RecoveryRate])
            σ = incubationrate(caseopt[:IncubationRate])
            α = socialdistancing(caseopt)
            Δt = discretestep(caseopt[:StepSize])
            pop = totalpopulation(caseopt[:Population])

            A,πv = NormalizedContactMatrix(caseopt)
            funcA = InstantaneousContactMatrix(caseopt, A, α)

            Indices = indices(caseopt,πv)
            u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
            p = merge(Indices, Dict(:ContactMatrix => funcA, :IncubationRate => σ, :RecoveryRate => μ, :Population => pop, :StepSize => Δt, :NPI => α))

            return discretestochasticSEIR!,πv,A,p,u0
        end
        _ => println("Model ", caseopt[:Model], " not found.")
    end
end

function sumind(u, ind) # u is Vector{Vector{Union{Int64,Float64}}}
    s = [sum(u[i][ind]) for i in 1:length(u)]
end
