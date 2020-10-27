# Read synthetic contact matrices
using RData, DataFramesMeta

"""
    function readcontactmatrices(region)
Reads the contact matrices from `region` (3-letter code) using data from
K. Prem et al., “Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era,” medRxiv, p. 2020.07.22.20159772, Jul. 2020.

Returns the four contact matrices: full, home, school, work, and others.
"""
function readcontactmatrices(region)
    #dir = "PremSyntheticContactMatrices/generate_synthetic_matrices/output/syntheticcontactmatrices2020/"
    dir = "ContactMatrices/syntheticcontactmatrices2020/"
    vars = load(dir * "overall/contact_all.rdata")
    Afull = vars["contact_all"][region]

    vars_school = load(dir * "overall/contact_school.rdata")
    Aschool = vars_school["contact_school"][region]

    vars_home = load(dir * "overall/contact_home.rdata")
    Ahome = vars_home["contact_home"][region]

    vars_others = load(dir * "overall/contact_others.rdata")
    Aothers = vars_others["contact_others"][region]

    vars_work = load(dir * "overall/contact_work.rdata")
    Awork = vars_work["contact_work"][region]
    return Afull, Ahome, Aschool, Awork, Aothers
end
"""
    function readpopulationstructure(region)
Reads the population structure (by age) `region` (3-letter code) using data from
K. Prem et al., “Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era,” medRxiv, p. 2020.07.22.20159772, Jul. 2020.

Returns a vector with population percentages per age group.  Age groups are
0-4, 5-9, ..., 70-74, 75+.
"""
function readpopulationstructure(countryname)
    dir = "PremSyntheticContactMatrices/generate_synthetic_matrices/input/pop/"
    vars = load(dir * "popratio.rdata")
    popratiofull = Array(@where(vars["popratio"], :countryname .== countryname)[1,4:end-1]) # last column is total population, discard
    popratio75 = [popratiofull[1:15];sum(popratiofull[16:end])] # population ratio includes groups up to 100 yo., but contact matrices agregate all 75+ year olds in a single group.
    return popratio75
end
