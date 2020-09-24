using LinearAlgebra, Match, Dates, Pandas

"""
    function discretestep(case)
Returns the stepsize (in days) for the discrete-time simulations.
"""
function discretestep(case)
    @match case begin
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
function searchsocialdistancing(t, d, d0, αd)
    if d[1] > d0 + Day(floor(Int,t)) # before first day
        return 1.0
    elseif d[end] < d0 + Day(floor(Int,t)) # There are many choices of what to do after the last value of social distancing index.  We should turn this into an option.
        #return 1.0 - αd[1] # Returns to first measured social distancing index
        #return 1.0 # No social distancing
        #return min(1.0 - αd[end], 1.0) # Last value of social distancing
        return min(1.0-αd[end] + (t-(d[end]-d0).value) * 0.0008, 1.0) # gradual decay
    else
        n0 = findfirst(x -> x >= d0 + Day(floor(Int,t)), d)
        return 1.0-αd[n0]
    end
end
"""
    function Rzero(case)
Returns the value of R0 used in each simulation.  The values here for the states are from William's paper.
"""
function Rzero(case)
    @match case begin
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
        :SPGov   => begin # From https://www.saopaulo.sp.gov.br/coronavirus/isolamento
            df = read_excel("NPI_Data/IsolamentoGovSP01set2020.xlsx")
            dfnames = columns(df)
            d = [Date(dfnames[i]*"20", "dd/mm/yyyy") for i in 4:size(df)[2]]
            #d = [Date(dfnames[i][end-4:end]*"/2020", "dd/mm/yyyy") for i in 4:size(df)[2]]
            αd = Array(query(df, :(Município1 == "SÃO PAULO")))[4:end]
            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), αd)
        end
        :SPInLoco => begin
            df = read_csv("NPI_Data/InLoco-isolation_Manaus_SaoPaulo.csv")
            dfnames = columns(df)
            d = Date.(Array(query(df, :(state_name == "São Paulo" && city_name == "São Paulo"))[[:dt]]))
            αd = Array(query(df, :(state_name == "São Paulo" && city_name == "São Paulo"))[[:isolated]])

            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), αd)
        end
        :ManausInLoco => begin
            df = read_csv("NPI_Data/InLoco-isolation_Manaus_SaoPaulo.csv")
            dfnames = columns(df)
            d = Date.(Array(query(df, :(state_name == "Amazonas" && city_name == "Manaus"))[[:dt]]))
            αd = Array(query(df, :(state_name == "Amazonas" && city_name == "Manaus"))[[:isolated]])

            return t -> searchsocialdistancing(t, d, firstday(caseopt[:FirstDay]), αd)
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
        :SP             =>  [0.076912061
                             0.100387104
                             0.105932843
                             0.351962947
                             0.24768619
                             0.117119] # Only SP city
        :None           => 1.0
        _               => begin
                                println("Choice for :AgeStructure not found")
                                1.0
                            end
     end
     return(πv ./ sum(πv))
 end

"""
    function MixingMatrix(case)
Returns the contact matrix between ages for each case, or 1.0 if no age structure.  Note that this matrix will be further normalized to the desired R0.
"""
function MixingMatrix(case)
     @match case begin
         :BrittonScience2020  => [2.2257 0.4136 0.2342 0.4539 0.2085 0.1506;
                            0.4139 3.6140 0.4251 0.4587 0.2712 0.1514;
                            0.2342 0.4257 2.9514 0.6682 0.4936 0.1972;
                            0.4539 0.4592 0.6676 0.9958 0.6510 0.3300;
                            0.2088 0.2706 0.4942 0.6508 0.8066 0.4341;
                            0.1507 0.1520 0.1968 0.3303 0.4344 0.7136]
        :None => [1.0]
        _ => begin
                    println("Choice for :MixingMatrix not found")
                    [1.0]
                end
    end
end
"""
    function activity_vector(case)
Returns the vector with the increase factors in contact (activity level) in heterogeneous models.  For example, a vector `[4.0, 1.0]` means that each age structure and class will have individuals with 4 times more contacts than the others.  The distribution of people in each set is given by `activity_percentages(case)`.
"""
function activity_vector(case)
    @match case begin
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
        :SPInloco       =>  begin # Using Inloco mobility data
                            act = [1.15595, 1.06564, 1.0]
                            act = act / (act⋅activity_percentages(caseopt[:ActivityStructure]))
                        end
        :None         => [1.0]
        _             => begin
                            println("Choice for :AtivityVector not found")
                            [1.0]
                        end
    end
end

"""
    function activity_percentages(case)
Returns the vector with the percentages of each different contact (activity) level group in contact in heterogeneous models.  See also `activity_vector(case)`.
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
    function πvector(caseopt)
Computes the complete vector of age+contact level groups (a Kronecker product of the age and contact level structures).
"""
function πvector(caseopt)
    kron(age_structure(caseopt[:AgeStructure]),
            activity_percentages(caseopt[:ActivityStructure]))
end

"""
    function Amatrix(caseopt)
Computes the contact matrix considering age and contact level groups.
"""
function Amatrix(caseopt)
    #println(activity_vector(case))
    kron(MixingMatrix(caseopt[:MixingMatrix]),
                activity_vector(caseopt[:ActivityVector])*activity_vector(caseopt[:ActivityVector])')
end


"""
    function Āmatrix(caseopt)
This is the matrix actually used in the SEIR model (equivalent to the β parameter in some texts.)
"""
function Āmatrix(caseopt)
    πvec = πvector(caseopt)
    Ā = Amatrix(caseopt) * diagm(πvec)
    return Ā,πvec
end

"""
    function indices(caseopt)
Returns the indices for the susceptible (`:S`), exposed (`:E`), infected (`:I`) and recovered/dead (`:R`) entries in the state vector.  Models with more compartments will include other variables (e.g., seroreversion, loss of immunity).
`caseopt` defines the set of parameters to be used.
"""
function indices(caseopt)
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
    Nlevels = length(πvector(caseopt))
    return Dict(compnames[i] => (i-1)*Nlevels+1:i*Nlevels for i in 1:compartments)
end

"""
    function initcond(caseopt,πv,sind,eind)
Returns a vector with the initial conditions for each case.
"""
function initcond(caseopt,πv,sind,eind)
    @match caseopt[:InitCond] begin
        :BrittonScience2020  =>  begin
                                ϵ = 1e-4
                                x = zeros(length(indices(caseopt))*length(sind))
                                x[eind] .= ϵ .* πv
                                x[sind] .= (1-ϵ) .* πv
                                return x
                            end
        :Discrete => begin
                            N0 = 1
                            Pop = totalpopulation(caseopt[:Population])
                            x = zeros(Int, length(indices(caseopt))*length(sind))
                            x[sind]=round.(Int,Pop*πv) # segfault if .*πv
                            while sum(x[sind]) != Pop
                               srand = rand(sind)
                               x[srand] += sum(x[sind]) < Pop ? 1 : -1
                               if x[srand] < 0
                                   x[srand] = 0
                               end
                            end
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
                zeros(x = zeros(Int, length(indices(caseopt))*length(sind)))
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
    A = p[:MixingMatrix]
    σ = p[:IncubationRate]
    μ = p[:RecoveryRate]
    ϵ0 = p[:SeroRevProb]
    p0 = p[:LossImmProb]
    ar = p[:LossImmRate]
    α = p[:NPI]

    s = view(u, sind)
    sr = view(u, srind)
    e = view(u, eind)
    i = view(u, iind)
    rsp = view(u, rspind)
    rsn = view(u, rsnind)
    rls = view(u, rlsind)
    tempse = α(t) .* (A'*i) .* s
    tempsre = α(t) .* (A'*i) .* sr
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
    A = p[:MixingMatrix]
    σ = p[:IncubationRate]
    μ = p[:RecoveryRate]
    α = p[:NPI]

    s = view(u, sind)
    e = view(u, eind)
    i = view(u, iind)
    r = view(u, rind)
    tempse = α(t) .* ((A'*i) .* s)
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
    println("Model: ", caseopt[:Model])

    @match caseopt[:Model] begin
        :SEIRAgeContinuous    =>  begin
                A = Amatrix(caseopt)
                μ = recoveryrate(caseopt[:RecoveryRate])
                σ = incubationrate(caseopt[:IncubationRate])
                α = socialdistancing(caseopt)
                Ā,πv = Āmatrix(caseopt)
                eig = eigen(Ā).values ./ μ
                κ = Rzero(caseopt[:R0]) / maximum(real.(eig))
                A .= κ .* A
                #println(maximum(real.(eig))*κ)
                Indices = indices(caseopt)
                u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
                p = merge(Indices, Dict(:MixingMatrix => A, :IncubationRate => σ, :RecoveryRate => μ, :NPI => α))
                return seir!,πv,A,p,u0
        end
        :SeroRevContinuous => begin
            A = Amatrix(caseopt)
            μ = recoveryrate(caseopt[:RecoveryRate])
            σ = incubationrate(caseopt[:IncubationRate])
            α = socialdistancing(caseopt)
            ar = lossofimmunityrate(caseopt[:LossImmRate])
            ϵ0 = seroreversionprobability(caseopt[:SeroRevProb])
            p0 = lossofimmunityprobability(caseopt[:LossImmProb])
            Ā,πv = Āmatrix(caseopt)
            eig = eigen(Ā).values ./ μ
            κ = Rzero(caseopt[:R0]) / maximum(real.(eig))
            A .= κ .* A
            #println(maximum(real.(eig))*κ)
            Indices = indices(caseopt)
            u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
            p = merge(Indices, Dict(:MixingMatrix => A, :IncubationRate => σ, :RecoveryRate => μ, :SeroRevProb => ϵ0, :LossImmProb => p0, :LossImmRate => ar, :NPI => α))

            return seirseroreversion!,πv,A,p,u0
        end
        :SEIRSeroRevDiscrete  => begin
            println("Model: Discrete SEIR with Seroreversion")
            A = Amatrix(caseopt)
            μ = recoveryrate(caseopt[:RecoveryRate])
            σ = incubationrate(caseopt[:IncubationRate])
            α = socialdistancing(caseopt)
            ar = lossofimmunityrate(caseopt[:LossImmRate])
            ϵ0 = seroreversionprobability(caseopt[:SeroRevProb])
            p0 = lossofimmunityprobability(caseopt[:LossImmProb])
            Δt = discretestep(caseopt[:StepSize])
            pop = totalpopulation(caseopt[:Population])

            Ā,πv = Āmatrix(caseopt)
            eig = eigen(Ā).values ./ μ
            κ = Rzero(caseopt[:R0]) / maximum(real.(eig))
            A .= κ .* A
            Indices = indices(caseopt)
            u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
            p = merge(Indices, Dict(:MixingMatrix => A, :IncubationRate => σ, :RecoveryRate => μ, :SeroRevProb => ϵ0, :LossImmProb => p0, :LossImmRate => ar, :Population => pop, :StepSize => Δt, :NPI => α))

            return discretestochasticSEIRSeroRev!,πv,A,p,u0
        end
        :SEIRDiscrete  => begin
            println("Model: Discrete SEIR")
            A = Amatrix(caseopt)
            μ = recoveryrate(caseopt[:RecoveryRate])
            σ = incubationrate(caseopt[:IncubationRate])
            α = socialdistancing(caseopt)
            Δt = discretestep(caseopt[:StepSize])
            pop = totalpopulation(caseopt[:Population])

            Ā,πv = Āmatrix(caseopt)
            eig = eigen(Ā).values ./ μ
            κ = Rzero(caseopt[:R0]) / maximum(real.(eig))
            A .= κ .* A
            Indices = indices(caseopt)
            u0 = initcond(caseopt, πv, Indices[:S], Indices[:E])
            p = merge(Indices, Dict(:MixingMatrix => A, :IncubationRate => σ, :RecoveryRate => μ, :Population => pop, :StepSize => Δt, :NPI => α))

            return discretestochasticSEIR!,πv,A,p,u0
        end
        _ => println("Model ", caseopt[:Model], " not found.")
    end
end

function sumind(u, ind) # u is Vector{Vector{Union{Int64,Float64}}}
    s = [sum(u[i][ind]) for i in 1:length(u)]
end
