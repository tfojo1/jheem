
##--------------------------------------------------------##
##-- THE FUNCTION TO ACTUALLY GET A POPULATION SKELETON --##
##--          (formulated as a private helper)          --##
##--------------------------------------------------------##

#if value=='sequential'
# value is set to 1:prod(dims)
do.get.population.skeleton <- function(dim.names, value)
{
    dims = sapply(dim.names, length)
    if (length(value)==1 && !is.na(value) && value=='sequential')
        value = 1:prod(dims)

    array(value, dim=dims, dimnames=dim.names)
}

##--------------------------------------------------##
##-- THE GENERAL GET POPULATION SKELETON FUNCTION --##
##--------------------------------------------------##

#'@title Get a skeleton for a population with indicated dimensions
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the shared elements of both HIV-positive and HIV-negative populations (age, race, subpopulation, sex*risk) for the specified JHEEM model
#'
#'@family functions to build population array skeletons
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param value The value to fill the array with. If value=='sequential', the array will be population with numbers 1:<size of the array>
#'@param age,race,subpopulation,sex,risk,non.hiv.subset,continuum,cd4,hiv.subset,age.from,race.from,subpopulation.from,sex.from,risk.from,non.hiv.subset.from,continuum.from,cd4.from,hiv.subset.from,age.to,race.to,subpopulation.to,sex.to,risk.to,non.hiv.subset.to,continuum.to,cd4.to,hiv.subset.to Indicators for whether these strata should be incorporated into the returned skeleton
#'
#'@seealso \code{\link{get.population.skeleton}}, \code{\link{get.hiv.negative.population.skeleton}}, \code{\link{get.hiv.positive.population.skeleton}}
#'
#'@export
get.population.skeleton <- function(jheem,
                                    value=0,
                                    age=F,
                                    race=F,
                                    subpopulation=F,
                                    sex=F,
                                    risk=F,
                                    non.hiv.subset=F,
                                    continuum=F,
                                    cd4=F,
                                    hiv.subset=F,
                                    age.from=F,
                                    race.from=F,
                                    subpopulation.from=F,
                                    sex.from=F,
                                    risk.from=F,
                                    non.hiv.subset.from=F,
                                    continuum.from=F,
                                    cd4.from=F,
                                    hiv.subset.from=F,
                                    age.to=F,
                                    race.to=F,
                                    subpopulation.to=F,
                                    sex.to=F,
                                    risk.to=F,
                                    non.hiv.subset.to=F,
                                    continuum.to=F,
                                    cd4.to=F,
                                    hiv.subset.to=F)
{
    dim.names = get.dimnames(jheem,
                             age=age, race=race, subpopulation=subpopulation,
                             sex=sex, risk=risk, non.hiv.subset=non.hiv.subset,
                             continuum.of.care=continuum, cd4=cd4, hiv.subset=hiv.subset,
                             age.from=age.from, race.from=race.from, subpopulation.from=subpopulation.from,
                             sex.from=sex.from, risk.from=risk.from, non.hiv.subset.from=non.hiv.subset.from,
                             continuum.of.care.from=continuum.from, cd4.from=cd4.from, hiv.subset.from=hiv.subset.from,
                             age.to=age.to, race.to=race.to, subpopulation.to=subpopulation.to,
                             sex.to=sex.to, risk.to=risk.to, non.hiv.subset.to=non.hiv.subset.to,
                             continuum.of.care.to=continuum.to, cd4.to=cd4.to, hiv.subset.to=hiv.subset.to)

    do.get.population.skeleton(dim.names, value)
}

##----------------------------------------------------------------##
##-- GET POPULATION SKELETON FUNCTIONS FOR SPECIFIC POPULATIONS --##
##----------------------------------------------------------------##

#'@title Get a skeleton HIV-negative population
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the HIV-negative population for the specified JHEEM model
#'
#'@family functions to build population array skeletons
#'@inheritParams get.population.skeleton
#'
#'@export
get.hiv.negative.population.skeleton <- function(jheem, value=0)
{
    get.population.skeleton(jheem, value=value,
                            age = T, race = T, subpopulation = T, sex = T, risk = T,
                            non.hiv.subset = T, continuum = F, cd4 = F, hiv.subset = F)
}

#'@title Get a skeleton HIV-positive population
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the HIV-positive population for the specified JHEEM model
#'
#'@family functions to build population array skeletons
#'@inheritParams get.population.skeleton
#'
#'@seealso \code{\link{get.population.skeleton}}, \code{\link{get.general.population.skeleton}}, \code{\link{get.hiv.negative.population.skeleton}}
#'
#'@export
get.hiv.positive.population.skeleton <- function(jheem, value=0)
{
    get.population.skeleton(jheem, value=value,
                            age = T, race = T, subpopulation = T, sex = T, risk = T,
                            non.hiv.subset = F, continuum = T, cd4 = T, hiv.subset = T)
}

#'@title Get a skeleton for the population characteristics shared by both HIV-negative and HIV-positive
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the shared elements of both HIV-positive and HIV-negative populations (age, race, subpopulation, sex*risk) for the specified JHEEM model
#'
#'@family functions to build population array skeletons
#'@inheritParams get.population.skeleton
#'
#'@seealso \code{\link{get.population.skeleton}}, \code{\link{get.hiv.negative.population.skeleton}}, \code{\link{get.hiv.positive.population.skeleton}}
#'
#'@export
get.general.population.skeleton <- function(jheem, value=0)
{
    get.population.skeleton(jheem, value=value,
                            age = T, race = T, subpopulation = T, sex = T, risk = T,
                            non.hiv.subset = F, continuum = F, cd4 = F, hiv.subset = F)
}


##---------------------------------------------------------------------##
##-- GET POPULATION SKELETON FUNCTIONS FOR SPECIFIC PARAMETER ARRAYS --##
##---------------------------------------------------------------------##

#'@title Get a skeleton contact array
#'
#'@family functions for transmission contact arrays
#'@family functions to build population array skeletons
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the contact transmission array: [age from, race from, subpopulation from, sex from, risk from, age to, race to, subpopulation to, sex to, risk to, non-hiv subset to]
#'
#'@inheritParams get.population.skeleton
#'@param all.dimensions If true, returns a contact array with all dimensions (age, race, subpopulation, sex, risk - both to and from) in it. If false, returns only the selected dimensions
#'
#'@seealso \code{\link{set.transmission.contact.array}}
#'
#'@export
get.contact.array.skeleton <- function(jheem, value=0,
                                       age=F,
                                       race=F,
                                       subpopulation=F,
                                       sex=F,
                                       risk=F,
                                       non.hiv.subset=F,
                                       age.from=age,
                                       race.from=race,
                                       subpopulation.from=subpopulation,
                                       sex.from=sex,
                                       risk.from=risk,
                                       age.to=age.from,
                                       race.to=race.from,
                                       subpopulation.to=subpopulation.from,
                                       sex.to=sex.from,
                                       risk.to=risk.from,
                                       non.hiv.subset.to=non.hiv.subset,
                                       all.dimensions=!age.from && !race.from && !subpopulation.from && !sex.from && !risk.from && !age.to && !race.to && !subpopulation.to && !sex.to && !risk.to && !non.hiv.subset.to
)
{
    get.population.skeleton(jheem, value=value,
                            age.from = age.from || all.dimensions,
                            race.from = race.from || all.dimensions,
                            subpopulation.from = subpopulation.from || all.dimensions,
                            sex.from = sex.from || all.dimensions,
                            risk.from = risk.from || all.dimensions,
                            age.to = age.to || all.dimensions,
                            race.to = race.to || all.dimensions,
                            subpopulation.to = subpopulation.to || all.dimensions,
                            sex.to = sex.to || all.dimensions,
                            risk.to = risk.to || all.dimensions,
                            non.hiv.subset.to = non.hiv.subset.to || all.dimensions)
}
