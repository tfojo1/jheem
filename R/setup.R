
##---------------------------------------##
##-- THE CONSTRUCTOR FOR JHEEM OBJECTS --##
##---------------------------------------##

NUM_HIV_STRATA = 8
NUM_NONHIV_STRATA = 6
NUM_GENERAL_STRATA = 5

#'@title Initializes a new controller for JHEEM
#'
#'@param version An identified of the version, however that is specified by the user. Recorded here for convenience of access later
#'@param age.cutoffs The cut points for age stratifications in the model, in years. Will result in length(age.cutoffs)-1 strata: age.cutoffs[1] (inclusive) to age.cutoffs[2] (exclusive), age.cutoffs[2] to age.cutoffs[3], ..., age.cutoffs[n-1] to age.cutoffs[n]
#'
#'
#'@export
initialize.jheem <- function(version,
                             age.cutoffs=NULL,
                             race.strata=NULL,
                             locations=NULL,
                             subpopulations=NULL,
                             sex.strata=NULL,
                             risk.strata=NULL,
                             continuum.of.care.states=NULL,
                             cd4.strata=NULL,
                             hiv.subsets=NULL,
                             nonhiv.subsets=NULL,
                             first.diagnosed.hiv.continuum.states=NULL,
                             all.diagnosed.hiv.continuum.states=first.diagnosed.hiv.continuum.states,
                             transmission.route.names='sexual',
                             verbose=T,
                             new.diagnoses.keep.dimensions=c('age','race','subpopulation','sex','risk')
)
{
    jheem = list(version=version)

    ##-- FILL IN MISSING STRATA --##
    if (is.null(age.cutoffs))
        age.cutoffs=c(0,Inf)
    if (is.null(race.strata))
        race.strata='all_races'
    if (is.null(locations))
        locations = 'all_locations'
    if (is.null(subpopulations))
        subpopulations = 'all_subpopulations'
    if (is.null(sex.strata))
        sex.strata = 'all_sexes'
    if (is.null(risk.strata))
        risk.strata = 'all_risk'
    if (is.null(continuum.of.care.states))
        continuum.of.care.states = 'all_continuum'
    if (is.null(cd4.strata))
        cd4.strata = 'all_cd4'
    if (is.null(hiv.subsets))
        hiv.subsets = 'all_hiv_positive'
    if (is.null(nonhiv.subsets))
        nonhiv.subsets = 'all_hiv_negative'

    if (class(race.strata) != 'character' || any(is.na(race.strata)) || length(race.strata)==0)
        stop("The race.strata argument must be a character vector with at least one element and no NA values")
    if (class(locations) != 'character' || any(is.na(locations)) || length(locations)==0)
        stop("The locations argument must be a character vector with at least one element and no NA values")
    if (class(subpopulations) != 'character' || any(is.na(subpopulations)) || length(subpopulations)==0)
        stop("The subpopulations argument must be a character vector with at least one element and no NA values")
    if (class(sex.strata) != 'character' || any(is.na(sex.strata)) || length(sex.strata)==0)
        stop("The sex.strata argument must be a character vector with at least one element and no NA values")
    if (class(risk.strata) != 'character' || any(is.na(risk.strata)) || length(risk.strata)==0)
        stop("The risk.strata argument must be a character vector with at least one element and no NA values")
    if (class(continuum.of.care.states) != 'character' || any(is.na(continuum.of.care.states)) || length(continuum.of.care.states)==0)
        stop("The continuum.of.care.states argument must be a character vector with at least one element and no NA values")
    if (class(cd4.strata) != 'character' || any(is.na(cd4.strata)) || length(cd4.strata)==0)
        stop("The cd4.strata argument must be a character vector with at least one element and no NA values")
    if (class(hiv.subsets) != 'character' || any(is.na(hiv.subsets)) || length(hiv.subsets)==0)
        stop("The hiv.subsets argument must be a character vector with at least one element and no NA values")
    if (class(nonhiv.subsets) != 'character' || any(is.na(nonhiv.subsets)) || length(nonhiv.subsets)==0)
        stop("The nonhiv.subsets argument must be a character vector with at least one element and no NA values")

    if ((class(age.cutoffs) != 'integer' && class(age.cutoffs) != 'numeric') || any(is.na(age.cutoffs)) || length(age.cutoffs)==0)
        stop("The age.cutoffs argument must be a numeric or integer vector with at least one element and no NA values")

    ##-- TRANSMISSION ROUTES --##
    if (is.null(transmission.route.names) || length(transmission.route.names)==0)
        stop("There must be at least one route of transmission specified in transition.route.names")
    if (class(transmission.route.names) != 'character' || any(is.na(transmission.route.names)))
        stop("The transmission.route.names argument must be a character vector with no NA values")
    if (any(table(transmission.route.names)>2))
        stop("The transmission.route.names argument can not have repeated names")

    jheem$transmission.routes = transmission.route.names

    ##-- AGES --##
    jheem$age = make.age.strata(age.cutoffs)

    ##-- OTHER STRATA --##
    jheem$race = race.strata
    jhem$locations = locations
    jheem$subpopulations = subpopulations
    jheem$sex = sex.strata
    jheem$risk.strata = risk.strata
    jheem$continuum.of.care = continuum.of.care.states
    jheem$cd4.strata = cd4.strata
    jheem$hiv.subsets = hiv.subsets
    jheem$nonhiv.subsets = nonhiv.subsets


    ##-- CHECK TO MAKE SURE NO DIM NAMES ARE REPEATED --##
    all.dimnames = get.dimnames.all(jheem)
    all.dimnames.flat = unlist(all.dimnames)
    all.dimnames.names = unlist(lapply(names(all.dimnames), function(name){rep(name, length(all.dimnames[[name]]))}))
    tabled.dimnames = table(all.dimnames.flat)
    if (any(tabled.dimnames>1))
    {
        error = paste0("A dimension name can only appear once (across all dimensions) for the JHEEM. The following name",
                       ifelse(sum(tabled.dimnames)>1, 's appear', ' appears'), " multiple times: ")
        for (name in names(tabled.dimnames)[tabled.dimnames>1])
        {
            error = paste0(error, "'", name, "' (in dimensions for ",
                           paste0(all.dimnames.names[all.dimnames.flat==name], collapse=", "),
                           ")")
        }
        stop(error)
    }


    ##-- CHECK TO MAKE SURE NO DIM NAMES CONTAIN PROTECTED WORDS --##
#    contain.from = grepl('.from$', all.dimnames.flat)
#    if (any(contain.from))
#        stop(paste0("No dimension names can end with '.from'. The following name",
#                    ifelse(sum(contain.from)>1, 's do', ' does'), ": ",
#                    paste0(all.dimnames.flat[contain.from], " (from the dimension for ", all.dimname.names[contain.from], ")", collapse=", ")))

#    contain.to = grepl('.to$', all.dimnames.flat)
#    if (any(contain.from))
#        stop(paste0("No dimension names can end with '.to'. The following name",
#                    ifelse(sum(contain.to)>1, 's do', ' does'), ": ",
#                    paste0(all.dimnames.flat[contain.to], " (from the dimension for ", all.dimname.names[contain.to], ")", collapse=", ")))


    ##-- SET UP THE PARAMETERS STUMP --##
    jheem$parameters = list(NUM_AGE_STRATA = length(jheem$age$labels),
                            NUM_RACE_STRATA = length(jheem$race),
                            NUM_SUBPOPULATIONS = length(jheem$subpopulations),
                            NUM_SEX_STRATA = length(jheem$sex),
                            NUM_RISK_STRATA = length(jheem$risk.strata),
                            NUM_CONTINUUM_STATES = length(jheem$continuum.of.care),
                            NUM_CD4_STRATA = length(jheem$cd4.strata),
                            NUM_HIV_SUBSETS = length(jheem$hiv.subsets),
                            NUM_NONHIV_SUBSETS = length(jheem$nonhiv.subsets),
                            NUM_TRANSMISSION_ROUTES = length(jheem$transmission.routes))

    jheem$parameters$NUM_GENERAL_STATES = jheem$parameters$NUM_AGE_STRATA * jheem$parameters$NUM_RACE_STRATA *
        jheem$parameters$NUM_SUBPOPULATIONS * jheem$parameters$NUM_SEX_STRATA * jheem$parameters$NUM_RISK_STRATA
    jheem$parameters$NUM_HIV_STATES = jheem$parameters$NUM_GENERAL_STATES *
        jheem$parameters$NUM_CONTINUUM_STATES * jheem$parameters$NUM_CD4_STRATA * jheem$parameters$NUM_HIV_SUBSETS
    jheem$parameters$NUM_NONHIV_STATES = jheem$parameters$NUM_GENERAL_STATES * jheem$parameters$NUM_NONHIV_SUBSETS

    jheem$parameters$NUM_HIV_STATES_PER_AGE_AND_RACE = jheem$parameters$NUM_SUBPOPULATIONS *
        jheem$parameters$NUM_SEX_STRATA * jheem$parameters$NUM_RISK_STRATA * jheem$parameters$NUM_CONTINUUM_STATES *
        jheem$parameters$NUM_CD4_STRATA * jheem$parameters$NUM_HIV_SUBSETS
    jheem$parameters$NUM_NONHIV_STATES_PER_AGE_AND_RACE = jheem$parameters$NUM_SUBPOPULATIONS *
        jheem$parameters$NUM_SEX_STRATA * jheem$parameters$NUM_RISK_STRATA * jheem$parameters$NUM_NONHIV_SUBSETS


    ##-- SET-UP FOR TRACKED TRANSITIONS --##
    jheem = init.transition.tracking(jheem)
    jheem$DIAGNOSED_CONTINUUM_STATES = jheem$UNDIAGNOSED_CONTINUUM_STATES =
        jheem$FIRST_DIAGNOSED_CONTINUUM_STATES = jheem$ALL_DIAGNOSED_HIV_CONTINUUM_STATES = NULL



    if (is.null(first.diagnosed.hiv.continuum.states))
    {
        print("Unable to track new diagnoses when 'first.diagnosed.hiv.continuum.states' is NULL")
    }
    else
    {
        first.diagnosed.hiv.continuum.states = unique(first.diagnosed.hiv.continuum.states)
        invalid.continuum.states = first.diagnosed.hiv.continuum.states[sapply(first.diagnosed.hiv.continuum.states, function(state){
            all(jheem$continuum.of.care != state)
        })]

        if (length(invalid.continuum.states)>0)
            stop(paste0("The following 'first.diagnosed.continuum.states' were not actually specified as part of the HIV continuum of care: ",
                        paste0("'", invalid.continuum.states, "'", collapse=', ')))

        jheem$FIRST_DIAGNOSED_CONTINUUM_STATES = first.diagnosed.hiv.continuum.states
    }

    all.diagnosed.hiv.continuum.states = union(first.diagnosed.hiv.continuum.states, all.diagnosed.hiv.continuum.states)
    if (is.null(all.diagnosed.hiv.continuum.states))
    {
        print("Unable to track new diagnoses when 'all.diagnosed.hiv.continuum.states' is NULL")
    }
    else
    {
        all.diagnosed.hiv.continuum.states = unique(all.diagnosed.hiv.continuum.states)
        invalid.continuum.states = all.diagnosed.hiv.continuum.states[sapply(all.diagnosed.hiv.continuum.states, function(state){
            all(jheem$continuum.of.care != state)
        })]

        if (length(invalid.continuum.states)>0)
            stop(paste0("The following 'all.diagnosed.continuum.states' were not actually specified as part of the HIV continuum of care: ",
                        paste0("'", invalid.continuum.states, "'", collapse=', ')))

        jheem$DIAGNOSED_CONTINUUM_STATES = all.diagnosed.hiv.continuum.states
    }

    if (!is.null(jheem$FIRST_DIAGNOSED_CONTINUUM_STATES) && !is.null(jheem$DIAGNOSED_CONTINUUM_STATES))
    {

        if (length(jheem$DIAGNOSED_CONTINUUM_STATES)==jheem$parameters$NUM_CONTINUUM_STATES)
        {
            print("All HIV continuum of care states were designated as first diagnosed states. Unable to track new diagnoses")
            jheem$UNDIAGNOSED_CONTINUUM_STATES = NULL
        }
        else
            jheem$UNDIAGNOSED_CONTINUUM_STATES = setdiff(jheem$continuum.of.care, jheem$DIAGNOSED_CONTINUUM_STATES)

        jheem = set.tracked.transition.hiv.positive(jheem, name='new_diagnoses',
                                                    transition.dimension = 'continuum',
                                                    from=jheem$UNDIAGNOSED_CONTINUUM_STATES,
                                                    to=jheem$FIRST_DIAGNOSED_CONTINUUM_STATES,
                                                    keep.dimensions=new.diagnoses.keep.dimensions)
    }

    ##-- SET UP FOR CONTACT MATRICES --##
    jheem$parameters$CONTACT_TYPE_FROM_DIMENSIONS = lapply(jheem$transmission.routes, function(route){NULL})

    ##-- SET UP FOR BIRTH PROPORTIONS --##
    jheem = init.for.birth.proportions(jheem)

    ##-- SET UP AGING RATES --##
    jheem$parameters$AGING_RATES = 1/jheem$age$spans

    ##-- SET UP DEFAULT TRACKING DIMENSIONS (ALL) --#
    jheem = set.track.incidence.dimensions(jheem)
    jheem = set.track.mortality(jheem)

    ##-- RETURN THE OBJECT --##

    class(jheem) = 'jheem'
    jheem
}

#'@title Make Age Strata based on cutoffs
#'
#'@description makes length(age.cutoffs)-1 age strata
#'
#'@param age.cutoffs A numeric vector of ages, where age.cutoffs[1] is the lower bound of the first age bracket, age.cutoffs[length(age.cutoffs)] is the upper bound of the last age bracket
#'
#'@return a list with the following components: $labels, $lowers, $uppers, and $spans
#'
#'@export
make.age.strata <- function(age.cutoffs)
{
    lowers = age.cutoffs[-length(age.cutoffs)]
    uppers = age.cutoffs[-1]
    labels = sapply(1:length(lowers), function(i){
        if (uppers[i]==Inf)
            paste0(lowers[i], '+ years')
        else
            paste0(lowers[i], '-', uppers[i]-1, ' years')
    })

    names(labels) = names(lowers) = names(uppers) = labels

    list(labels=labels,
         lowers=lowers,
         uppers=uppers,
         spans=uppers-lowers)
}




check.transmission.route.names <- function(jheem, route.names)
{
    invalid.mask = sapply(route.names, function(name){
        all(name != jheem$transmission.routes)
    })

    if (any(invalid.mask))
    {
        invalid.names = route.names[invalid.mask]

        stop(paste0(get.character.list(paste0("'", invalid.names, "'")),
                    ifelse(length(invalid.names)==1, ' is not a', ' are not'),
                    " pre-specified transmission route",
                    ifelse(length(invalid.names)==1, '', 's'),
                    " for the JHEEM object. The valid transmission route",
                    ifelse(length(jheem$transmission.routes)==1, ' is', 's are'),
                    get.character.list(paste0("'", jheem$transmission.routes, "'"))))
    }
}

##--------------------------##
##-- POPULATION SKELETONS --##
##--------------------------##





##--------------------------------------------------------------##
##-- STRATIFYING and OTHERWISE MANIPULATING POPULATION ARRAYS --##
##--------------------------------------------------------------##

#'@title Substratify a population
#'
#'@description Returns an array such that each cell in 'population' is multiplied by the array given by 'stratification' to yield an array of dim=c(dim(population), dim(stratification))
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param value The value to fill the array with
#'
#'@seealso \code{\link{get.population.skeleton}}, \code{\link{get.hiv.negative.population.skeleton}}, \code{\link{get.hiv.positive.population.skeleton}}
#'
#'@export
stratify.population <- function(population,
                                stratification)
{
    if (is.null(dim(stratification)))
        stratification = array(stratification, dim=length(stratification), dimnames=list(names(stratification)))

    array(as.numeric(population) * rep(as.numeric(stratification), each=prod(dim(population))),
          dim=c(dim(population), dim(stratification)),
          dimnames=c(dimnames(population), dimnames(stratification)))
}



#'@title Takes a population and splits males into heterosexual males and MSM
#'
#'@description Returns an array with same dimensions and dimnames as population.to.stratify, except that the males have been split into heterosexual males and MSM
#'
#'@param population.to.stratify A population where dimnames for one dimension includes male.name
#'@param msm.multiplier,heterosexual.multiplier The numbers by which to multiply the males in population.to.stratify to yield the corresponding populations in the stratified population
#'@param male.name The name given to males in a dimension of population.to.stratify
#'@param heterosexual.male.name,msm.name The names to be given to heterosexual males and MSM in the return value
#'
#'@seealso \code{\link{expand.population.to.hiv.negative}}
#'
#'@export
stratify.males.by.orientation <- function(population.to.stratify,
                                          msm.multiplier=1,
                                          heterosexual.multiplier=1,
                                          male.name='male',
                                          heterosexual.male.name='heterosexual_male',
                                          msm.name='msm')
{
    #Figure out which current dimension contains males
    orig.dim = dim(population.to.stratify)
    dim.names = dimnames(population.to.stratify)
    sex.dim.mask = sapply(dim.names, function(dnames){
        any(dnames==male.name)
    })
    if (!any(sex.dim.mask))
        stop(paste0("None of the dimnames for population.to.stratify contain the field '", male.name, "'"))
    sex.dim = (1:length(dim.names))[sex.dim.mask]
    male.index.within.dim = index.of(male.name, dim.names[[sex.dim]])

    #Set up new dimensions to split males
    new.dim.names = dim.names
    non.male.sex.mappings = 1:orig.dim[sex.dim]
    non.male.sex.mappings[non.male.sex.mappings>male.index.within.dim]  = non.male.sex.mappings[non.male.sex.mappings>male.index.within.dim] + 1
    non.male.sex.mappings = non.male.sex.mappings[-male.index.within.dim]
    new.dim.names[[sex.dim]] = c(dim.names[[sex.dim]][(1:orig.dim[sex.dim])<male.index.within.dim],
                                 heterosexual.male.name, msm.name,
                                 dim.names[[sex.dim]][(1:orig.dim[sex.dim])>male.index.within.dim])
    new.dims = dim(population.to.stratify)
    new.dims[sex.dim] = new.dims[sex.dim] + 1

    #Fold up the population.to.stratify into three dimensions, where the middle is the sex dimension
    dims.before = prod(orig.dim[1:length(orig.dim)<sex.dim])
    dims.after = prod(orig.dim[1:length(orig.dim)>sex.dim])
    dim(population.to.stratify) = c(dims.before, orig.dim[sex.dim], dims.after)

    #Map to a stratified population
    rv = array(NA, dim=c(dims.before, orig.dim[sex.dim]+1, dims.after))
    #non-males
    rv[,non.male.sex.mappings,] = population.to.stratify[,-male.index.within.dim,]
    #heterosexual males
    rv[,male.index.within.dim,] = population.to.stratify[,male.index.within.dim,] * heterosexual.multiplier
    #MSM
    rv[,male.index.within.dim+1,] = population.to.stratify[,male.index.within.dim,] * msm.multiplier

    #Hydrate the rv into an array with desired dimensions and return
    dim(rv) = new.dims
    dimnames(rv) = new.dim.names
    rv
}





##-------------------##
##-- INITIAL STATE --##
##-------------------##

#'@title Set the initial HIV-positive population
#'
#'@family functions to set model births
#'
#'@description A convenience function that bundles \code{\link{set.initial.population.hiv.positive}} and \code{\link{set.initial.population.hiv.negastive}} into one function call
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param init.hiv.negative,init.hiv.positive Arrays representing the initial hiv-positive and hiv-negative populations
#'
#'@export
set.initial.populations <- function(jheem,
                                    init.hiv.negative,
                                    init.hiv.positive)
{
    jheem = set.initial.population.hiv.negative(jheem, init.hiv.negative)
    jheem = set.initial.population.hiv.positive(jheem, init.hiv.positive)

    jheem
}

#'@title Set the initial HIV-positive population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param init An array representing the initial hiv-positive population
#'
#'@seealso \code{\link{set.initial.state.hiv.negative}}, \code{\link{get.hiv.positive.population.skeleton}}
#'
#'@export
set.initial.population.hiv.positive <- function(jheem,
                                           init)
{
    needed.dim.names = get.dimnames.hiv(jheem)
    if (!array.matches.dimnames(init, needed.dim.names))
        stop("The initial state must have ", length(needed.dim.names), " dimensions: ",
             paste0(paste0(names(needed.dim.names), " (", sapply(needed.dim.names, length), " state", c('','s')[1+(sapply(needed.dim.names, length)>1)] ,")"), collapse=', '))

    jheem$initial.state.hiv.positive = init

    jheem
}

#'@title Set the initial HIV-negative population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param init An array representing the initial hiv-negative population
#'
#'@seealso \code{\link{set.initial.state.hiv.positive}}, \code{\link{get.hiv.negative.population.skeleton}}
#'
#'@export
set.initial.population.hiv.negative <- function(jheem,
                                           init)
{
    needed.dim.names = get.dimnames.nonhiv(jheem)
    if (!array.matches.dimnames(init, needed.dim.names))
        stop("The initial state must have ", length(needed.dim.names), " dimensions: ",
             paste0(paste0(names(needed.dim.names), " (", sapply(needed.dim.names, length), " state", c('','s')[1+(sapply(needed.dim.names, length)>1)] ,")"), collapse=', '))

    jheem$initial.state.hiv.negative = init

    jheem
}


#'@title Get a pair of seeded initial states
#'
#'@description Seeds cases into an HIV-positive population, and subtracts those cases from the corresponding strata in the HIV-negative population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param cases.per.subpopulation Either (1) a numeric vector with length the number of subpopulations, where each value represents the number of cases to be seeded into each subpopulation, or (2) a single number, which will be the number of cases seeded into each subpopulation
#'@param seed.to.races The races into which cases will be seeded
#'@param seed.to.ages The age strata into which cases will be seeded
#'@param seed.to.risk.strata The sex*risk strata into which cases will be seeded
#'
#'@seealso \code{\link{set.initial.state.hiv}}, \code{\link{get.hiv.negative.population.skeleton}}
#'
#'@export
get.seeded.initial.populations <- function(jheem,
                                      nonhiv.init,
                                      cases.per=1,
                                      seed.to.ages=jheem$age$labels[ceiling(length(jheem$age$labels)/2)],
                                      seed.to.races=jheem$race[1],
                                      seed.to.subpopulations=jheem$subpopulations[1],
                                      seed.to.sexes=jheem$sex[1],
                                      seed.to.risk.strata=jheem$risk.strata[1],
                                      seed.to.continuum.of.care=jheem$continuum.of.care[1],
                                      seed.to.cd4=jheem$cd4.strata[1],
                                      seed.to.hiv.subset=jheem$hiv.subsets[1],
                                      seed.from.non.hiv.subset=jheem$nonhiv.subsets[1])
{
    nonhiv.dimnames = get.dimnames.nonhiv(jheem)
    dim(nonhiv.init) = sapply(nonhiv.dimnames, length)
    dimnames(nonhiv.init) = nonhiv.dimnames

    hiv.init = get.hiv.positive.population.skeleton(jheem)

    for (age in seed.to.ages)
    {
        for (race in seed.to.races)
        {
            for (subpopulation in seed.to.subpopulations)
            {
                for (sex in seed.to.sexes)
                {
                    for (risk in seed.to.risk.strata)
                    {
                        for (continuum in seed.to.continuum.of.care)
                        {
                            for (cd4 in seed.to.cd4)
                            {
                                for (hiv.subset in seed.to.hiv.subset)
                                {
                                    hiv.init[age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset] = cases.per
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cases.from = rowSums(hiv.init, dims=5)
    from.proportions = rep(0, jheem$parameters$NUM_NONHIV_SUBSETS)
    names(from.proportions) = jheem$nonhiv.subsets
    from.proportions[seed.from.non.hiv.subset] = 1/length(seed.from.non.hiv.subset)
    cases.from = array(as.numeric(cases.from) * rep(from.proportions, each=prod(dim(cases.from))),
                       dim=dim(nonhiv.init), dimnames=dimnames(nonhiv.init))

    if (any(cases.from > nonhiv.init))
        stop('More cases are requested for some strata than are present in the hiv-negative population')

    list(hiv.positive=hiv.init,
         hiv.negative=nonhiv.init - cases.from)
}



##----------------------##
##-- FIXED POPULATION --##
##----------------------##




##---------------##
##-- FERTILITY --##
##---------------##

##-----------------------------##
##-- TIME-VARYING PARAMETERS --##
##-----------------------------##

set.time.varying.parameter.value <- function(jheem,
                                             parameter.name,
                                             parameter.value,
                                             times,
                                             stepwise=F)
{
    jheem = create.time.varying.parameter.slot(jheem, parameter.name, stepwise = stepwise)

    for (one.time in times)
    {
        jheem$parameters$time.varying.parameters[[parameter.name]]$values = c(jheem$parameters$time.varying.parameters[[parameter.name]]$value, list(parameter.value))
        jheem$parameters$time.varying.parameters[[parameter.name]]$times = c(jheem$parameters$time.varying.parameters[[parameter.name]]$times, one.time)
    }

    o = order(jheem$parameters$time.varying.parameters[[parameter.name]]$times)
    jheem$parameters$time.varying.parameters[[parameter.name]]$values = jheem$parameters$time.varying.parameters[[parameter.name]]$values[o]
    jheem$parameters$time.varying.parameters[[parameter.name]]$times = jheem$parameters$time.varying.parameters[[parameter.name]]$times[o]

    jheem
}

set.time.varying.parameter.stepwise <- function(jheem,
                                                parameter.name,
                                                stepwise)
{
    jheem = create.time.varying.parameter.slot(jheem, parameter.name)
    jheem$parameters$time.varying.parameters[[parameter.name]]$stepwise=stepwise

    jheem
}

create.time.varying.parameter.slot <- function(jheem, parameter.name, stepwise=F)
{
    if (is.null(jheem$parameters$time.varying.parameters[[parameter.name]]))
        jheem$parameters$time.varying.parameters[[parameter.name]] = list(values=NULL,
                                                                          times=numeric(),
                                                                          stepwise=stepwise)
    else
        jheem$parameters$time.varying.parameters[[parameter.name]]$stepwise = stepwise

    jheem
}

#Parameters we are going to need:
# NUM_AGE_STRATA
# FERTILITY_RATES (and anchor times)
# MODEL_BIRTHS
#   a boolean indicating whether we explicitly model births based on birth rates
# BIRTH_PROPORTIONS (and anchor times)
#   indexed[race, subpopulation, risk_stratum, nonhiv_subset]
#   BIRTH_PROPORTIONS[r,s,,] sums to 1, and represents the fractional
#   distribution of new births in to the race=r, subpopulation=s stratum across
#   levels of risk_strata and nonhiv_subset
# GENERAL_MORTALITY_FOR_HIV_NEGATIVE (and anchor times)
#   Indexed same as HIV- population
# GENERAL_MORTALITY_FOR_HIV_POSITIVE (and anchor times)
#   Indexed same as HIV+ population
# HIV_SPECIFIC_MORTALITY (and anchor times)
#   Indexed same as HIV+ population
# TRANSMISSIBILITY (and anchor times)
#   Indexed same as HIV+ population
# SUSCEPTIBILITY (and anchor times)
#   Indexed same as HIV- population
# CONTACT_MATRIX (and anchor times)
#   a NUM_HIV_STATES x NUM_NONHIV_STATES matrix
# NEW_INFECTION_PROPORTIONS (and anchor times)
#   indexed [age,race,subpopulation,risk*sex,nonhiv_subset,continuum,cd4,hiv_subset]
#   NEW_INFECTION_PROPORTIONS[a,r,s,rf,nhs,,,] sums to 1,
#   and represents the distribution of new infections from the
#   age=a, race=r, subpopulation=s, risk*sex=rf, nonhiv_subset=nhs stratum of uninfected
#   across the levels of hiv+
# HIV_POSITIVE_TRANSITION_MATRICES (and anchor times)
#   for each a,r
#       HIV_POSITIVE_TRANSITION_MATRICES[a,r,,] is a 1 x NUM_HIV_STATES_PER_AGE_AND_RACE x NUM_HIV_STATES_PER_AGE_AND_RACE matrix, where
#       HIV_POSITIVE_TRANSITION_MATRICES[a,r,i,j] represents the rate at which individuals move from i to j
#       flattened(hiv.positive[a,r,,,,,]) %*% HIV_POSITIVE_TRANSITION_MATRICES[a,r,,] is the movement to the flattened population
#       flattened(hiv.positive[a,r,,,,,]) * rowSums(HIV_POSITIVE_TRANSITION_MATRICES[a,r,,]) is the movement from the flattened population
# HIV_NEGATIVE_TRANSITION_MATRICES
#   as for HIV_POSITIVE_TRANSITION_MATRICES (and anchor times)
# FIX_STRATA_SIZES
#   a boolean indicating whether we use fixed strata sizes
# TARGET_STRATUM_SIZE (and anchor times)
#   a numeric of length NUM_FIXED_STRATA indicating the size of each
# FIXED_ENTRY_PROPORTIONS (and anchor times)
#   a numeric of length (NUM_HIV_STATES + NUM_NONHIV_STATES) indicating the proportion of fixed entries into the stratum
#   to which new entries into that stratum should be put into
# FIXED_STRATUM
#   a numeric of length (NUM_HIV_STATES + NUM_NONHIV_STATES) indicating which fixed-size stratum the flattened compartment goes into
# IS_STRATUM_FIXED_PROPORTION
#   a boolean of length NUM_FIXED_STRATA indicating whether the stratum is fixed-proportion (as opposed to fixed-entry)
# NUM_FIXED_STRATA
# NUM_HIV_STATES
# NUM_NONHIV_STATES
# NUM_AGE_STRATA, NUM_RACE_STRATA, NUM_SUBPOPULATIONS, NUM_RISK_STRATA
# NUM_CONTINUUM_STATES, NUM_CD4_STRATA, NUM_HIV_SUBSETS
# NUM_NONHIV_SUBSETS
# FEMALE_INDICES
#   The indices in the sex*risk stratum that represent females (who can give birth)
# AGING.RATES
# TRACKED_HIV_TRANSITIONS_FROM, TRACKED_HIV_TRANSITIONS_TO
# TRACKED_NONHIV_TRANSITIONS_FROM, TRACKED_NONHIV_TRANSITIONS_TO
# NUM_HIV_STATES_PER_AGE_AND_RACE, NUM_NONHIV_STATES_PER_AGE_AND_RACE
