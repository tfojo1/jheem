
##--------------------------------------------------------------------------------##
##--                          POPULATION_GROWTH.R                               --##
##--                                                                            --##
##-- This file contains functions to govern the way the population grows
##--  either by new births into the population or fixing population strata size --##                                              --##
##-- This comprises: fertility, birth proportions, and fixing population size   --##
##--                                                                            --##
##--------------------------------------------------------------------------------##


##---------------##
##-- FERTILITY --##
##---------------##


#'@title Sets fertility rates for both the HIV-negative and HIV-positive population
#'
#'@family functions to specify births in the population
#'@family functions to set fertility rates
#'
#'@description *calling with a single scalar value for fertility would be equivalent to setting a birth rate (rather than a fertility rate) for the entire population
#'
#'@param jheem A JHEEM object
#'@param fertility An array of fertility rates, whose dimensions are a subset of [age, race, subpopulation, sex, risk]. ***Rates should be per person per year
#'@param time The time at which these fertility rates apply, if time-varying
#'
#'@export
set.fertility <- function(jheem,
                          fertility,
                          time=-Inf)
{
    jheem = set.fertility.hiv.negative(jheem, fertility, time)
    jheem = set.fertility.hiv.positive(jheem, fertility, time)

    jheem
}

#'@title Sets fertility rates for the HIV-negative population
#'
#'@family functions to specify births in the population
#'@family functions to set fertility rates
#'
#'@description *calling with a single scalar value for fertility would be equivalent to setting a birth rate (rather than a fertility rate) for the entire population
#'
#'@param jheem A JHEEM object
#'@param fertility An array of fertility rates, whose dimensions are a subset of [age, race, subpopulation, sex, risk, non.hiv.subset]. ***Rates should be per person per year
#'@param time The time at which these fertility rates apply, if time-varying
#'
#'@export
set.fertility.hiv.negative <- function(jheem,
                                       fertility,
                                       time=-Inf)
{
    fertility = expand.population.to.hiv.negative(jheem, fertility)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'FERTILITY_RATES_NONHIV',
                                     parameter.value = fertility,
                                     time = time)
}

#'@title Sets fertility rates for the HIV-positive population
#'
#'@family functions to specify births in the population
#'@family functions to set fertility rates
#'
#'@description *calling with a single scalar value for fertility would be equivalent to setting a birth rate (rather than a fertility rate) for the entire population
#'
#'@param jheem A JHEEM object
#'@param fertility An array of fertility rates, whose dimensions are a subset of [age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset]. ***Rates should be per person per year
#'@param time The time at which these fertility rates apply, if time-varying
#'
#'@export
set.fertility.hiv.positive <- function(jheem,
                                       fertility,
                                       time=-Inf)
{
    fertility = expand.population.to.hiv.positive(jheem, fertility)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'FERTILITY_RATES_HIV',
                                     parameter.value = fertility,
                                     time = time)
}

##---------------------------##
##-- SET BIRTH PROPORTIONS --##
##---------------------------##

#'@title Sets the proportions according to which new births are distributed across the population assuming no maternal-fetal transmission of HIV
#'
#'@family functions to specify birth proportions
#'@family functions to specify births in the population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param proportions An array of birth proportions, typically created by calling \code{\link{create.birth.proportions.hiv.negative}}. Assumed to be the birth proportions for both HIV-negative and HIV-positive parents
#'@param time The time at which these proportions apply, if time-varying
#'
#'@export
set.birth.proportions.no.maternal.transmission <- function(jheem,
                                                           proportions,
                                                           time=-Inf)
{
    jheem = set.birth.proportions.hiv.negative(jheem, proportions, time)
    jheem = set.birth.proportions.hiv.positive(jheem,
                                               fraction.births.infected=0,
                                               proportions.for.uninfected.births = proportions,
                                               time=time)

    jheem
}

#'@title Sets the proportions according to which new births from HIV-negative parents are distributed across the HIV-negative population
#'
#'@family functions to specify birth proportions
#'@family functions to specify births in the population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param proportions An array of birth proportions, typically created by calling \code{\link{create.birth.proportions.hiv.negative}}
#'@param time The time at which these proportions apply, if time-varying
#'
#'@export
set.birth.proportions.hiv.negative <- function(jheem,
                                               proportions,
                                               time=-Inf)
{
    do.set.birth.proportions(jheem,
                             proportions=proportions,
                             suffix="NONHIV",
                             required.to.dim.names=get.to.hiv.negative.birth.proportion.dimnames(jheem),
                             allowed.from.dim.names=get.from.hiv.negative.birth.proportion.dimnames(jheem),
                             time=time)
}

#'@title Sets the proportions according to which new births from HIV-positive are distributed across the HIV-negative and HIV-positive populations
#'
#'@family functions to specify birth proportions
#'@family functions to specify births in the population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param proportions.for.uninfected.births,proportions.for.infected.births Arrays of birth proportions for uninfected and infected offpring, typically created by calling \code{\link{create.birth.proportions.hiv.positive.to.negative}} or \code{\link{create.birth.proportions.hiv.positive.to.positive}}
#'@param fraction.births.infected An array that matches the HIV-positive population, indicating what fraction of offspring are born with HIV. Dimensions should be a subset of [age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset]
#'@param time The time at which these proportions apply, if time-varying
#'
#'@export
set.birth.proportions.hiv.positive <- function(jheem,
                                               fraction.births.infected,
                                               proportions.for.uninfected.births,
                                               proportions.for.infected.births=NULL,
                                               time=-Inf)
{
    if (!is.null(fraction.births.infected))
        fraction.births.infected = expand.population.to.hiv.positive(jheem, fraction.births.infected)


    if (!is.null(fraction.births.infected))
        jheem = set.time.varying.parameter.value(jheem,
                                                 parameter.name = 'FRACTION_BIRTHS_INFECTED',
                                                 parameter.value = fraction.births.infected,
                                                 time = time)

    if (!is.null(proportions.for.uninfected.births))
        jheem = do.set.birth.proportions(jheem,
                                         proportions=proportions.for.uninfected.births,
                                         suffix="HIV_TO_NONHIV",
                                         required.to.dim.names=get.to.hiv.negative.birth.proportion.dimnames(jheem),
                                         allowed.from.dim.names=get.from.hiv.positive.birth.proportion.dimnames(jheem),
                                         time=time)

    if (!is.null(proportions.for.infected.births))
        jheem = do.set.birth.proportions(jheem,
                                         proportions=proportions.for.infected.births,
                                         suffix="HIV_TO_HIV",
                                         required.to.dim.names=get.to.hiv.positive.birth.proportion.dimnames(jheem),
                                         allowed.from.dim.names=get.from.hiv.positive.birth.proportion.dimnames(jheem),
                                         time=time)

    jheem
}

#internal helper to set up dimnames and add proportions parameter value
do.set.birth.proportions <- function(jheem,
                                     proportions,
                                     suffix,
                                     required.to.dim.names,
                                     allowed.from.dim.names,
                                     time)
{
    required.to.dim.name.names = names(required.to.dim.names)
    allowed.from.dim.name.names = names(allowed.from.dim.names)
    dim.name.names = names(dimnames(proportions))

    #Hydrate up to make sure we include race
    proportions = expand.population(proportions, get.dimnames.by.name(jheem, union('race.from', dim.name.names)))
    dim.name.names = names(dimnames(proportions))

    #Check to make sure all required are present
    missing.to.dim.name.names = setdiff(required.to.dim.name.names, dim.name.names)
    if (length(missing.to.dim.name.names)>0)
        stop("The following '.to' dimension",
             ifelse(length(missing.to.dim.name.names)>1, "s are", " is"),
             " missing from the given proportions array: ",
             paste0("'", missing.to.dim.name.names, "'", collapse=", "))

    #Check to make sure only allowed are present
    extra.dim.name.names = setdiff(dim.name.names, union(required.to.dim.name.names, allowed.from.dim.name.names))
    if (length(extra.dim.name.names)>1)
        stop("The following dimension",
             ifelse(length(extra.dim.name.names)>1, "s are", " is"),
             " not allowed in the given proportions array: ",
             paste0("'", extra.dim.name.names, "'", collapse=", "))

    #If we haven't previously specified this type of birth proportions, save the from dimensions
    # Otherwise, check to make sure these from dimensions match those from dimensions

    #From dims - the numeric indices into vector[age.from, race.from, supopulation.from, sex.from, risk.from, <other from indices>]
    from.dims = (1:length(allowed.from.dim.name.names))[sapply(allowed.from.dim.name.names, function(name){
        any(dim.name.names==name)
    })]
    #race dim - the index into from dims of the race.from field
    race.dim = index.of('race.from', allowed.from.dim.name.names[from.dims])

    #Normalize so all proportions sum to 1
    present.from.dim.name.names = intersect(dim.name.names, allowed.from.dim.name.names)
    denominators = apply(proportions, present.from.dim.name.names, sum)
    denominators = expand.population(denominators, dimnames(proportions))
    proportions = proportions / as.numeric(denominators)

    #Save array dims into the jheem
    if (is.null(jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS[[suffix]]))
    {
        jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS[[suffix]] = from.dims
    }
    else if (!all(from.dims==jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS[[suffix]]))
    {
        stop("The 'from' dimensions for the given birth proportions array (",
             paste0("'", allowed.from.dimensions[from.dims], "'", collapse=', '),
             ") do not match the 'from' dimensions previously specified for '", suffix, "' birth proportions arrays (",
             paste0("'", names(jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS[[suffix]]), "'", collapse=', '),
             ")")
    }

    #Set the time-varying parameter
    proportions.name = paste0('BIRTH_PROPORTIONS_', suffix)
    set.time.varying.parameter.value(jheem,
                                     parameter.name = proportions.name,
                                     parameter.value = proportions,
                                     time = time)
}

#sets up the lists needed for birth proportions
init.for.birth.proportions <- function(jheem)
{
    jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS = list()

    jheem
}

##---------------------------##
##-- SET BIRTH PROPORTIONS --##
##---------------------------##

#'@title Create a Birth Proportions array based on sex distribution by race from a given population
#'
#'@family functions to create birth proportions arrays
#'@family functions to specify births in the population
#'
#'@param jheem A JHEEM object
#'@param population The population based on which to determine sex distribution by race. Must have at least a dimension for sex
#'
#'@description Creates a birth proportions array with the marginals of sex by race drawn from the population. The returned value is valid as birth proportions for either HIV-negative or for uninfected births to HIV-positive
#'
#'@details Assumes that children are born into the same subpopulation as their parents, that all children are born into the first risk state (unless otherwise specified) and the first non.hiv.subset (unless otherwise specified), and that distributions of sex by race are determined by the marginals on the given population
#'
#'@export
create.birth.proportions.from.population <- function(jheem,
                                                     population,
                                                     risk.born.into.state=1,
                                                     non.hiv.subset.born.into.state=1)
{
    #Identify the dim names and flag race and sex
    dim.name.names = identify.dimnames.in.array(jheem, population, single.occurence.as='')
    dim.name.names[dim.name.names=='race'] = 'race.from'
    dim.name.names[dim.name.names=='sex'] = 'sex.to'
    population = set.dim.name.names(population, dim.name.names)

    #marginalize out all but race and sex
    if (any(dim.name.names=='race.from'))
        keep.dimensions = c('race.from', 'sex.to')
    else
        keep.dimensions = 'sex.to'

    population = apply(population, c('race.from','sex.to'), sum)

    #hydrate up to a birth proportions array
    create.birth.proportions.hiv.negative(jheem, population,
                                          subpopulation.same.as.parents = T,
                                          risk.born.into.state = risk.born.into.state,
                                          non.hiv.subset.born.into.state = non.hiv.subset.born.into.state)
}

#'@title Create a skeleton Birth Proportions for births into the HIV-negative population (from either HIV-negative or HIV-positive parents)
#'
#'@family functions to create birth proportions arrays
#'@family functions to specify births in the population
#'
#'@param jheem A JHEEM object
#'@param value The default value with which to populate the array
#'
#'@description Creates an empty birth proportions array for births into the HIV-negative population
#'
#'@export
get.birth.proportions.to.hiv.negative.skeleton <- function(jheem, value=0)
{
    get.population.skeleton(jheem,
                            race.from=T,
                            subpopulation.from=T,
                            subpopulation.to=T,
                            sex.to=T,
                            risk.to=T,
                            non.hiv.subset.to=T,
                            value=value)
}

#'@title Create a valid Birth Proportions array for births to HIV-negative parents
#'
#'@family functions to create birth proportions arrays
#'@family functions to specify births in the population
#'
#'@inheritParams do.create.birth.proportions
#'
#'@inherit do.create.birth.proportions details
#'
#'@export
create.birth.proportions.hiv.negative <- function(jheem, ...,
                                                  subpopulation.same.as.parents=F,
                                                  risk.same.as.parents=F,
                                                  non.hiv.subset.same.as.parents=F,
                                                  subpopulations.born.into.state=NULL,
                                                  sex.born.into.state=NULL,
                                                  risk.born.into.state=NULL,
                                                  non.hiv.subset.born.into.state=NULL)
{
    do.create.birth.proportions(jheem,
                                required.to.dimnames=get.to.hiv.negative.birth.proportion.dimnames(jheem),
                                ...,
                                subpopulation.same.as.parents=subpopulation.same.as.parents,
                                risk.same.as.parents=risk.same.as.parents,
                                non.hiv.subset.same.as.parents=non.hiv.subset.same.as.parents,
                                subpopulations.born.into.state=subpopulations.born.into.state,
                                sex.born.into.state=sex.born.into.state,
                                risk.born.into.state=risk.born.into.state,
                                non.hiv.subset.born.into.state=non.hiv.subset.born.into.state)
}

#'@title Create a valid Birth Proportions array for HIV-uninfected births to HIV-positive parents
#'
#'@family functions to create birth proportions arrays
#'@family functions to specify births in the population
#'
#'@inheritParams do.create.birth.proportions
#'
#'@inherit do.create.birth.proportions details
#'
#'@export
create.birth.proportions.hiv.positive.to.negative <- function(jheem, ...,
                                                                 subpopulation.same.as.parents=F,
                                                                 risk.same.as.parents=F,
                                                                 subpopulations.born.into.state=NULL,
                                                                 sex.born.into.state=NULL,
                                                                 risk.born.into.state=NULL,
                                                                 non.hiv.subset.born.into.state=NULL)
{
    do.create.birth.proportions(jheem,
                                required.to.dimnames=get.to.hiv.negative.birth.proportion.dimnames(jheem),
                                ...,
                                subpopulation.same.as.parents=subpopulation.same.as.parents,
                                risk.same.as.parents=risk.same.as.parents,
                                subpopulations.born.into.state=subpopulations.born.into.state,
                                sex.born.into.state=sex.born.into.state,
                                risk.born.into.state=risk.born.into.state,
                                non.hiv.subset.born.into.state=non.hiv.subset.born.into.state)
}

#'@title Create a valid Birth Proportions array for HIV-infected births to HIV-positive parents
#'
#'@family functions to create birth proportions arrays
#'@family functions to specify births in the population
#'
#'@inheritParams do.create.birth.proportions
#'
#'@inherit do.create.birth.proportions details
#'
#'@export
create.birth.proportions.hiv.positive.to.positive <- function(jheem, ...,
                                                                 subpopulation.same.as.parents=F,
                                                                 risk.same.as.parents=F,
                                                                 hiv.subset.same.as.parents=F,
                                                                 subpopulations.born.into.state=NULL,
                                                                 sex.born.into.state=NULL,
                                                                 risk.born.into.state=NULL,
                                                                 continuum.born.into.state=NULL,
                                                                 cd4.born.into.state=NULL,
                                                                 hiv.subset.born.into.state=NULL)
{
    do.create.birth.proportions(jheem,
                                required.to.dimnames=get.to.hiv.positive.birth.proportion.dimnames(jheem),
                                ...,
                                subpopulation.same.as.parents=subpopulation.same.as.parents,
                                risk.same.as.parents=risk.same.as.parents,
                                hiv.subset.same.as.parents=hiv.subset.same.as.parents,
                                subpopulations.born.into.state=subpopulations.born.into.state,
                                sex.born.into.state=sex.born.into.state,
                                risk.born.into.state=risk.born.into.state,
                                continuum.born.into.state=continuum.born.into.state,
                                cd4.born.into.state=cd4.born.into.state,
                                hiv.subset.born.into.state=hiv.subset.born.into.state)
}

#'
#'@param jheem A instance of a JHEEM
#'@param ...
#'@param subpopulation.same.as.parents,risk.same.as.parents,non.hiv.subset.same.as.parents,hiv.subset.same.as.parents Boolean values indicating whether births from each of these dimensions should be into the same state as the parents
#'@param subpopulations.born.into.state,sex.born.into.state,risk.born.into.state,non.hiv.subset.born.into.state,continuum.born.into.state,cd4.born.into.state,hiv.subset.born.into.state The name or index of the state within each dimension into which all new births should go. A value of NULL indicates that births into a specific state should not be used for the dimension
#'
#'@details A full birth proportions array will have a set of '.from' dimensions (representing the demographics of parents) and '.to' dimensions (representing the demographics of newborns). While the full set of '.to' dimensions must be specified for a valid birth proportions array, any or all of the '.from' dimensions may be omitted (in which case, the birth proportions will be assumed to be the same for the omitted dimensions). Specifically:
#'\itemize{
#'\item{\strong{Births to HIV-negative parents} must have a subset of [age.from, race.from, subpopulation.from, sex.from, risk.from, non.hiv.subset.from] and all of [subpopulation.to, sex.to, risk.to, non.hiv.subset.to]}
#'\item{\strong{Uninfected births to HIV-positive parents} must have a subset of [age.from, race.from, subpopulation.from, sex.from, risk.from, continuum.from, cd4.from, hiv.subset.from] and all of [subpopulation.to, sex.to, risk.to, non.hiv.subset.to]}
#'\item{\strong{Infected births to HIV-positive parents} must have a subset of [age.from, race.from, subpopulation.from, sex.from, risk.from, continuum.from, cd4.from, hiv.subset.from] and all of [subpopulation.to, sex.to, risk.to, continuum.to, cd4.to, hiv.subset.to]}
#'}
#'(There is no age.to because all births are into the first age stratum, and there is no race.to because it is assumed all newborns will be of the same race as their parents)
#'The create.birth.proportions functions are convenience functions used to easily construct valid birth proportions arrays from minimal information.
#'They require the user to specify how each of the '.to' dimensions will be constructed in the birth proportions array. The form of each '.to' dimension may be specified by one of:
#'\enumerate{
#'\item{An explicit array passed as one of the ... arguments. These arrays can have any combination of the requisite '.from' and '.to' dimensions (the functions will assume that if two dimensions corresponding to, say, sex, are present, the first is the .from and the second is the .to. If only one is present, it is assumed to be the .to)}
#'\item{An indication that newborns go into the same dimension as their parents. Eg, subpopulation.to.same.as.parents=T means that subpopulation.to wil equal subpopulation.from (unless otherwise specified by one of the explicitly passed arrays)}
#'\item{An indication of which state from the dimension all newborns go into. Eg, subpopulation.to.state=1 would indicate that all newborns go into the first subpopulation state}
#'\item{If a dimension only has one state (eg, just one subpopulation) then no specification is needed}
#'}
#'The arguments must specify using one of these means for each of the needed .to dimensions.
#'
#'@keywords internal
do.create.birth.proportions <- function(jheem,
                                        required.to.dimnames,
                                        ...,
                                        subpopulation.same.as.parents=F,
                                        risk.same.as.parents=F,
                                        non.hiv.subset.same.as.parents=F,
                                        hiv.subset.same.as.parents=F,
                                        subpopulations.born.into.state=NULL,
                                        sex.born.into.state=NULL,
                                        risk.born.into.state=NULL,
                                        non.hiv.subset.born.into.state=NULL,
                                        continuum.born.into.state=NULL,
                                        cd4.born.into.state=NULL,
                                        hiv.subset.born.into.state=NULL)
{
    #-- Figure out dimnames in arrays in ... --#
    array.marginals = lapply(list(...), function(one.marginal)
    {
        source.dim.name.names = identify.dimnames.in.array(jheem, one.marginal,
                                                           age.single.as='.from',
                                                           race.single.as='.from',
                                                           subpopulation.single.as='.to',
                                                           sex.single.as='.to',
                                                           risk.single.as='.to',
                                                           non.hiv.subset.single.as='.to',
                                                           continuum.single.as='.to',
                                                           cd4.single.as='.to',
                                                           hiv.subset.single.as='.to')
        one.marginal = set.dim.name.names(one.marginal, source.dim.name.names)
        one.marginal
    })

    #-- Organize the dim names from different components to be checked --#
    #Required dim names
    required.to.dim.name.names = names(required.to.dimnames)

    #Dims with only one state
    length.one.dim.name.names = required.to.dim.name.names[sapply(required.to.dimnames, length)==1]


    #Dims same from parents to children
    no.change.mask = c(subpopulation.to=subpopulation.same.as.parents, risk.to=risk.same.as.parents,
                       non.hiv.subset.to=non.hiv.subset.same.as.parents, hiv.subset.to=hiv.subset.same.as.parents)
    dim.name.names.no.change = names(no.change.mask)[no.change.mask]

    #Dims into a specified state
    specified.states = list(subpopulation.to=subpopulations.born.into.state,
                             sex.to=sex.born.into.state,
                             risk.to=risk.born.into.state,
                             non.hiv.subset.to=non.hiv.subset.born.into.state,
                             continuum.to=continuum.born.into.state,
                             cd4.to=cd4.born.into.state,
                             hiv.subset.to=hiv.subset.born.into.state)
    specified.states = specified.states[!sapply(specified.states, is.null)]
    dim.name.names.specified.state = names(specified.states)

    #Dims specified in arrays
    dim.name.names.present.in.arrays = unique(unlist(sapply(array.marginals, function(one.marginal){
        names(dimnames(one.marginal))
    })))


    #-- Check to make sure all required .to dim names have been specified one way or another --#
    unspecified.required = setdiff(required.to.dim.name.names,
                                   union(union(length.one.dim.name.names, dim.name.names.present.in.arrays),
                                         union(dim.name.names.no.change, dim.name.names.specified.state)))
    if (length(unspecified.required)>0)
        stop("Births into the following '.to' dimension",
             ifelse(length(unspecified.required)>1, "s have", " has"),
             " not been specified in any way: ",
             paste0("'", unspecified.required, "'", collapse=", "))


    #-- Check to make sure that .to dim names have not been specified in more than one way --#
    overlap.no.change.arrays = intersect(dim.name.names.present.in.arrays, dim.name.names.no.change)
    overlap.specified.arrays = intersect(dim.name.names.present.in.arrays, dim.name.names.specified.state)
    overlap.no.change.specified = intersect(dim.name.names.no.change, dim.name.names.specified.state)

    if (length(overlap.no.change.arrays)>0)
        stop("Births into the following '.to' dimension",
             ifelse(length(overlap.no.change.arrays)>1, "s have", " has"),
             " been specified both in an array and as no change from parents to offspring: ",
             paste0("'", overlap.no.change.arrays, "'", collapse=", "),
             ". You must specify in one way or the other, but not both.")

    if (length(overlap.specified.arrays)>0)
        stop("Births into the following '.to' dimension",
             ifelse(length(overlap.specified.arrays)>1, "s have", " has"),
             " been specified both in an array and as being born into particular states: ",
             paste0("'", overlap.specified.arrays, "'", collapse=", "),
             ". You must specify in one way or the other, but not both.")

    if (length(overlap.no.change.specified)>0)
        stop("Births into the following '.to' dimension",
             ifelse(length(overlap.no.change.specified)>1, "s have", " has"),
             " been specified both as no change from parents to offspring and as being born into particular states: ",
             paste0("'", overlap.no.change.specified, "'", collapse=", "),
             ". You must specify in one way or the other, but not both.")


    #-- Hydrate up the no changes and the specified states as arrays --#

    no.change.marginals = lapply(dim.name.names.no.change, function(to.name){
        create.no.change.transition.array(get.dimnames.by.name(jheem, gsub('.to', '', to.name)))
    })

    specified.state.marginals = lapply(dim.name.names.specified.state, function(to.name){
        rv = do.get.population.skeleton(get.dimnames.by.name(jheem, to.name), value=0)
        rv[specified.states[[to.name]]] = 1
        rv
    })


    #-- Pull dimensions in any of the component arrays and set up the rv --#

    all.marginals = c(array.marginals,
                      no.change.marginals,
                      specified.state.marginals)

    target.dim.name.names = unique(unlist(sapply(all.marginals, function(one.marginal){
        names(dimnames(one.marginal))
    })))
    target.dim.name.names = unique(c(target.dim.name.names, length.one.dim.name.names))

    target.dim.names = get.dimnames.by.name(jheem, target.dim.name.names)
    rv = do.get.population.skeleton(target.dim.names, value=1)


    #-- Multiply all component arrays together to form the rv --#
    for (one.marginal in all.marginals)
    {
        #hydrate the marginal array up to the full size
        one.marginal = expand.population(one.marginal, target.dim.names)
        #multiply into the rv
        rv = rv * one.marginal
    }

    #-- Normalize so that birth proportions across any stratum of .from sum to 1 --#
    if (length(target.dim.name.names) == length(required.to.dimnames))
        denominators = sum(rv)
    else
        denominators = rowSums(rv, dims=length(target.dim.name.names)-length(required.to.dimnames))

    denominators[denominators==0] = 1
    rv = rv / as.numeric(denominators)


    #-- Return it --#
    #rv
}

##-----------##
##-- AGING --##
##-----------##

#'@title Sets aging rates for both the HIV-negative and HIV-positive population
#'
#'@description Equivalent to calling set.aging.hiv.negative(jheem, aging.rates, time) and set.aging.hiv.positive(jheem, aging.rates, time). In general, aging rates can be approximated by the proportion of individuals in an age bracket who are in the last year of the age bracket
#'@family Functions to set aging rates
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param aging An array of aging rates, whose dimensions are a subset of the dimensions shared by both the HIV-negative and HIV-positive populations (age, race, subpopulation, sex, risk factor). ***Rates should be per person per year
#'@param time The time at which these mortality rates apply, if time-varying
#'
#'@export
set.aging <- function(jheem,
                      aging.rates,
                      time=-Inf)
{
    jheem = set.aging.hiv.positive(jheem, aging.rates=aging.rates, time=time)
    jheem = set.aging.hiv.negative(jheem, aging.rates=aging.rates, time=time)

    jheem
}

#'@title Sets aging rates for HIV-positive population
#'
#'@family Functions to set aging rates
#'
#'@inheritParams set.aging
#'
#'@export
set.aging.hiv.positive <- function(jheem,
                                   aging.rates,
                                   time=-Inf)
{
    aging.rates = expand.population.to.hiv.positive(jheem, aging.rates)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'AGING_RATES_HIV',
                                     parameter.value = aging.rates,
                                     time = time)
}

#'@title Sets aging rates for HIV-negative population
#'
#'@family Functions to set aging rates
#'
#'@inheritParams set.aging
#'
#'@export
set.aging.hiv.negative <- function(jheem,
                                   aging.rates,
                                   time=-Inf)
{
    aging.rates = expand.population.to.hiv.negative(jheem, aging.rates)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'AGING_RATES_NONHIV',
                                     parameter.value = aging.rates,
                                     time = time)
}

#'@title Get a default array of aging rates for the HIV-positive population
#'
#'@description Assumes that ages are distributed uniformly across the age bracket, such that the aging rate for an age bracket is 1 / the number of years in that age bracket
#'
#'@family Functions to set aging rates
#'
#'@inheritParams set.aging
#'
#'@export
get.aging.skeleton.hiv.positive <- function(jheem)
{
    get.hiv.positive.population.skeleton(jheem, value=1/jheem$age$spans)
}

#'@title Get a default array of aging rates for the HIV-negative population
#'
#'@description Assumes that ages are distributed uniformly across the age bracket, such that the aging rate for an age bracket is 1 / the number of years in that age bracket
#'
#'@family Functions to set aging rates
#'
#'@inheritParams set.aging
#'
#'@export
get.aging.skeleton.hiv.negative <- function(jheem)
{
    get.hiv.negative.population.skeleton(jheem, value=1/jheem$age$spans)
}

##---------------##
##-- MORTALITY --##
##---------------##

#'@title Sets mortality rates for both the HIV-negative and HIV-positive population
#'
#'@description Equivalent to calling set.general.mortality.hiv.negative(jheem, mortality, time) and set.general.mortality.hiv.positive(jheem, mortality, time)
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param mortality An array of mortality rates, whose dimensions are a subset of the dimensions shared by both the HIV-negative and HIV-positive populations (age, race, subpopulation, sex, risk factor). ***Rates should be per person per year
#'@param time The time at which these mortality rates apply, if time-varying
#'
#'@seealso \code{\link{set.general.mortality.hiv.negative}}, \code{\link{set.general.mortality.hiv.positive}}
#'
#'@export
set.general.mortality <- function(jheem,
                                  mortality,
                                  time=-Inf)
{
    jheem = set.general.mortality.hiv.negative(jheem, mortality, time)
    jheem = set.general.mortality.hiv.positive(jheem, mortality, time)

    jheem
}

#'@title Set general mortality rates for HIV-positive population (not including HIV-specific mortality)
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param mortality An array of mortality rates, whose dimensions are a subset of the dimensions of the HIV-negative population. ***Rates should be per person per year
#'@param time The time at which these mortality rates apply, if time-varying
#'
#'@seealso \code{\link{set.general.mortality}}, \code{\link{set.general.mortality.hiv.negative}}, \code{\link{expand.population.to.hiv.positive}}
#'
#'@export
set.general.mortality.hiv.positive <- function(jheem,
                                               mortality,
                                               time=-Inf)
{
    mortality = expand.population.to.hiv.positive(jheem, mortality)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'GENERAL_MORTALITY_FOR_HIV_POSITIVE',
                                     parameter.value = mortality,
                                     time = time)
}

#'@title Set mortality rates for HIV-negative population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param mortality An array of mortality rates, whose dimensions are a subset of the HIV-negative population. ***Rates should be per person per year
#'@param time The time at which these mortality rates apply, if time-varying
#'
#'@seealso \code{\link{set.hiv.specific.mortality}}, \code{\link{set.general.mortality}}, \code{\link{set.general.mortality.hiv.positive}}, \code{\link{expand.population.to.hiv.negative}}
#'
#'@export
set.general.mortality.hiv.negative <- function(jheem,
                                               mortality,
                                               time=-Inf)
{
    mortality = expand.population.to.hiv.negative(jheem, mortality)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'GENERAL_MORTALITY_FOR_HIV_NEGATIVE',
                                     parameter.value = mortality,
                                     time = time)
}


#'@title Set HIV-specific mortality rates for HIV-positive population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param mortality An array of mortality rates, whose dimensions are a subset of the dimensions of the HIV-negative population. ***Rates should be per person per year
#'@param time The time at which these mortality rates apply, if time-varying
#'
#'@seealso \code{\link{set.general.mortality}}, \code{\link{set.general.mortality.hiv.positive}}, \code{\link{expand.population.to.hiv.positive}}
#'
#'@export
set.hiv.specific.mortality <- function(jheem,
                                       mortality,
                                       time=-Inf)
{
    mortality = expand.population.to.hiv.positive(jheem, mortality)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'HIV_SPECIFIC_MORTALITY',
                                     parameter.value = mortality,
                                     time = time)
}

#'@title Set whether the JHEEM should track HIV-specific and overall mortality
#'
#'@family function to set up JHEEM mortality
#'
#'@description A convenience function that bundles \code{\link{set.track.hiv.specific.mortality}} and \code{\link{set.track.overall.mortality}} into one function call
#'
#'@inheritParams set.general.mortality
#'@param track.hiv.specific.mortality,track.overall.hiv.positive.mortality,track.overall.hiv.negative.mortality Booleans indicating whether or not to track each type for mortality
#'
#'@export
set.track.mortality <- function(jheem,
                                track.hiv.specific.mortality,
                                track.overall.hiv.positive.mortality,
                                track.overall.hiv.negative.mortality)
{
    jheem = set.track.hiv.specific.mortality(jheem, track.hiv.specific.mortality)
    jheem = set.track.overall.mortality(jheem,
                                        track.hiv.positive = track.overall.hiv.positive.mortality,
                                        track.hiv.negative = track.overall.hiv.negative.mortality)

    jheem
}



##--------------------------------##
##-- IMMIGRATION and EMIGRATION --##
##--------------------------------##

#'@title Sets immigration rates
#'
#'@family functions to specify immigration
#'
#'@details Immigration rates are interpreted similarly to birth rates: they are used to calculate the total number of immigrants, which are then distributed through either the HIV-negative or HIV-positive populations according to separately set proportions (set using \code{\link{set.immigration.proportions}})
#'
#'@param jheem A JHEEM object
#'@param immigration.rates An array of immigration rates, whose dimensions are a subset of [age, race, subpopulation, sex, risk].
#'@param is.relative If true, the rates are interpreted as per person per year. If false, they are interpreted as absolute numbers
#'@param time The time at which these immigration rates apply, if time-varying
#'
#'@export
set.immigration.rates <- function(jheem,
                                  immigration.rates,
                                  is.relative=T,
                                  time=-Inf)
{
    if ((is.null(dimnames(immigration.rates)) || is.null(names(dimnames(immigration.rates)))) &&
        !is.scalar(immigration.rates))
        stop("The 'immigration.rates' argument must have named dimnames or be a scalar value")

    if (is.scalar(immigration.rates))
        jheem = set.immigration.rate.dimensions(jheem)
    else
    {
        dim.name.names = names(dimnames(immigration.rates))
        jheem = set.immigration.rate.dimensions(jheem,
                                                age=any(dim.name.names=='age'),
                                                race=any(dim.name.names=='race'),
                                                subpopulation=any(dim.name.names=='subpopulation'),
                                                sex=any(dim.name.names=='sex'),
                                                risk=any(dim.name.names=='risk'))
    }

    jheem = set.time.varying.parameter.value(jheem,
                                             parameter.name = 'IMMIGRATION_RATES',
                                             parameter.value = immigration.rates,
                                             time = time)

    jheem = set.time.varying.parameter.value(jheem,
                                             parameter.name = 'USE_RELATIVE_IMMIGRATION',
                                             parameter.value = is.relative,
                                             time = time,
                                             stepwise = T)

    jheem
}

#'@export
set.immigration.rate.dimensions <- function(jheem,
                                            age=F,
                                            race=F,
                                            subpopulation=F,
                                            sex=F,
                                            risk=F)
{
    target.dim.name.names = get.dimnames(jheem, age=age, race=race, subpopulation=subpopulation, sex=sex, risk=risk)
    general.dim.name.names = names(get.dimnames.general(jheem))
    from.dims = (1:length(general.dim.name.names))[sapply(general.dim.name.names, function(name){
        any(name==target.dim.name.names)
    })]

    if (is.null(jheem$parameters$IMMIGRATION_FROM_DIMENSIONS))
        jheem$parameters$IMMIGRATION_FROM_DIMENSIONS = from.dims
    else if (!all(from.dims == parameters$IMMIGRATION_FROM_DIMENSIONS))
        stop(paste0("The dimensions from immigration (",
                    ifelse(length(from.dims)==0, "no dimensions",
                        paste0(general.dim.name.names[from.dims], collapse=", ")),
                    ") do not match the dimensions set by a previous call to set.immigration.rates or set.immigration.rate.dimensions (",
                    ifelse(length(parameters$IMMIGRATION_FROM_DIMENSIONS), "no dimensions",
                           paste0(general.dim.name.names[parameters$IMMIGRATION_FROM_DIMENSIONS])),
                    ")"))

    if (length(target.dim.name.names)==0)
    {
        jheem$parameters$IMMIGRATION_NONHIV_PROPORTIONS_INDICES = 1
        jheem$parameters$IMMIGRATION_HIV_PROPORTIONS_INDICES = 1
    }
    else
    {
        indexed.general = get.population.skeleton(jheem, value='sequential',
                                                  age=any(target.dim.name.names=='age'),
                                                  race=any(target.dim.name.names=='race'),
                                                  subpopulation=any(target.dim.name.names=='subpopulation'),
                                                  sex=any(target.dim.name.names=='sex'),
                                                  risk=any(target.dim.name.names=='risk'))
        jheem$parameters$IMMIGRATION_NONHIV_PROPORTIONS_INDICES = as.numeric(expand.population.to.hiv.negative(indexed.general))
        jheem$paramaters$IMMIGRATION_HIV_PROPORTIONS_INDICES = as.numeric(expand.population.to.hiv.positive(indexed.general))
    }

    jheem
}



#'@title Sets the proportions according to which new births are distributed across the population assuming no maternal-fetal transmission of HIV
#'
#'@family functions to specify birth proportions
#'@family functions to specify births in the population
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param proportions An array of birth proportions, typically created by calling \code{\link{create.birth.proportions.hiv.negative}}. Assumed to be the birth proportions for both HIV-negative and HIV-positive parents
#'@param time The time at which these proportions apply, if time-varying
#'
#'@export
set.immigration.proportions <- function(jheem,
                                        use.population.proportions,
                                        hiv.negative.proportions=NULL,
                                        hiv.positive.proportions=NULL,
                                        time=-Inf)
{
    if (is.null(jheem$parameters$IMMIGRATION_FROM_DIMENSIONS))
        stop("Prior to calling set.immigration.proportions, you must specify which dimensions are used for immigration with a call to set.immigration.rate.dimensions")

    #Check dimnames of given proportions
    if (!use.population.proportions)
    {
        nonhiv.is.zero = is.null(hiv.negative.proportions) || is.scalar(hiv.negative.proportions, 0)
        hiv.is.zero = is.null(hiv.positive.proportions) || is.scalar(hiv.positive.proportions, 0)

        if (nonhiv.is.zero && hiv.is.zero)
            stop("hiv.negative.proportions and hiv.positive.proportions cannot both be zero/NULL")

        if (nonhiv.is.zero)
            hiv.negative.proportions = expand.population.to.hiv.negative(jheem, 0)
        else if (!array.matches.dimnames(hiv.negative.proportions, get.dimnames.nonhiv(jheem), allow.missing.length.one.dims=T))
            stop("hiv.negative.proportions must have the full set of dimensions for the HIV-negative population: ",
                 paste0(names(get.dimnames.nonhiv(jheem)), collapse=', '))

        if (hiv.is.zero)
            hiv.positive.proportions = expand.population.to.hiv.negative(jheem, 0)
        else if (!array.matches.dimnames(hiv.positive.proportions, get.dimnames.hiv(jheem), allow.missing.length.one.dims=T))
            stop("hiv.positive.proportions must have the full set of dimensions for the HIV-positive population: ",
                 paste0(names(get.dimnames.hiv(jheem)), collapse=', '))
    }
    else
    {
        if (!is.null(hiv.negative.proportions))
            warning('You have passed a non-null value to hiv.negative.proportions, but since use.population.proportions is true, this value will be ignored')

        if (!is.null(hiv.negative.proportions))
            warning('You have passed a non-null value to hiv.negative.proportions, but since use.population.proportions is true, this value will be ignored')
    }

    #Normalize
    if (!use.population.proportions)
    {
        denominators = rowSums(hiv.negative.proportions, dims=5) + rowSums(hiv.positive.proportions, dims=5)
        if (any(denominators==0))
            stop("The given hiv.negative.proportions and hiv.positive.proportions sum to zero for some strata of ",
                 paste0(names(dimnames(get.dimnames.general(jheem))[jheem$parameters$IMMIGRATION_FROM_DIMENSIONS], collapse=' x ')),
                 ". The proportions must sum to a non-zero value for each stratum.")

        hiv.negative.proportions = hiv.negative.proportions / as.numeric(denominators)
        hiv.positive.proportions = hiv.positive.proportions / as.numeric(denominators)
    }

    #Set values

    jheem = set.time.varying.parameter.value(jheem,
                                             parameter.name='USE_SET_IMMGRATION_PROPORTIONS',
                                             parameter.value=!use.population.proportions,
                                             time=time)

    if (!use.population.proportions)
    {
        jheem = set.time.varying.parameter.value(jheem,
                                                 parameter.name='IMMIGRATION_PROPORTIONS_NONHIV',
                                                 parameter.value=hiv.negative.proportions,
                                                 time=time,
                                                 stepwise=T)

        jheem = set.time.varying.parameter.value(jheem,
                                                 parameter.name='IMMIGRATION_PROPORTIONS_HIV',
                                                 parameter.value=hiv.positive.proportions,
                                                 time=time,
                                                 stepwise=T)
    }

    jheem
}

#'@title Sets emigration rates for both the HIV-negative and HIV-positive population
#'
#'@description Equivalent to calling set.emigration.hiv.negative(jheem, mortality, time) and set.emigration.hiv.positive(jheem, mortality, time)
#'
#'@family functions to specify emigration
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param emigration An array of emigration rates, whose dimensions are a subset of the dimensions shared by both the HIV-negative and HIV-positive populations (age, race, subpopulation, sex, risk factor). ***Rates should be per person per year
#'@param time The time at which these emigration rates apply, if time-varying
#'
#'@export
set.emigration <- function(jheem,
                           emigration,
                           time=-Inf)
{
    jheem = set.emigration.hiv.negative(jheem, emigration, time)
    jheem = set.emigration.hiv.positive(jheem, emigration, time)

    jheem
}

#'@title Set emigration rates for HIV-positive population
#'
#'@family functions to specify emigration
#'
#'@param emigration An array of emigration rates, whose dimensions are a subset of the dimensions of the HIV-negative population. ***Rates should be per person per year
#'@inheritParams
#'
#'@export
set.emigration.hiv.positive <- function(jheem,
                                               emigration,
                                               time=-Inf)
{
    emigration = expand.population.to.hiv.positive(jheem, emigration)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'EMIGRATION_FOR_HIV_POSITIVE',
                                     parameter.value = emigration,
                                     time = time)
}

#'@title Set emigration rates for HIV-negative population
#'
#'@family functions to specify emigration
#'
#'@param mortality An array of mortality rates, whose dimensions are a subset of the HIV-negative population. ***Rates should be per person per year
#'@inheritParams
#'
#'@export
set.emigration.hiv.negative <- function(jheem,
                                        emigration,
                                        time=-Inf)
{
    emigration = expand.population.to.hiv.negative(jheem, emigration)

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'EMIGRATION_FOR_HIV_NEGATIVE',
                                     parameter.value = emigration,
                                     time = time)
}

##------------------------##
##-- FIXED BY MORTALITY --##
##------------------------##

#'@title Sets which dimensions define strata that will have sizes kept constant
#'
#'@family functions to fix population strata sizes
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param fix.age,fix.race,fix.subpopulation,fix.sex,fix.risk Boolean values indicating whether or not these dimensions define the strata whose sizes are fixed
#'
#'@details Fixing strata sizes allows the JHEEM to adjust rates of change such that the overall size of the population within strata stays constant. The strata may be defined by any subset of the five dimensions shared by both HIV-negative and HIV-positive subpopulations: <age, race, subpopulation, sex, risk>.
#'Note that two separate functions must be called to use this ability: \code{\link{set.fixed.size.strata}} to indicte which dimensions define fixed strata and \code{\link{set.keep.strata.sizes.constant}} to actually indicate at what time to fix strata this way
#'
#'
#'@export
set.fixed.size.strata <- function(jheem,
                                  fix.age=F,
                                  fix.race=F,
                                  fix.subpopulation=F,
                                  fix.sex=F,
                                  fix.risk=F)
{
    fixed.skeleton = get.population.skeleton(jheem, value='sequential',
                                             age=fix.age,
                                             race=fix.race,
                                             subpopulation=fix.subpopulation,
                                             sex=fix.sex,
                                             risk=fix.risk)
    fixed.dim.name.names = names(dimnames(fixed.skeleton))

    general.dim.name.names = names(get.dimnames.general(jheem))
    matched.dim.name.names = sapply(general.dim.name.names, function(name){any(name==fixed.dim.name.names)})

    jheem$parameters$CONSTANT_SIZE_DIMENSIONS = (1:length(general.dim.name.names))[matched.dim.name.names]

    jheem$parameters$HIV_NEGATIVE_INDICES_INTO_CONSTANT_SIZE_STRATA = as.numeric(expand.population.to.hiv.negative(jheem, fixed.skeleton))
    jheem$parameters$HIV_POSITIVE_INDICES_INTO_CONSTANT_SIZE_STRATA = as.numeric(expand.population.to.hiv.positive(jheem, fixed.skeleton))

    jheem
}

#'@title Sets whether to keep population strata sizes constant for a given time
#'
#'@family functions to fix population strata sizes
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param fix.strata.sizes A boolean indicated whether to use fixed strata sizes
#'@param time The time following which the setting applies
#'
#'@inherit set.fixed.size.strata details
#'
#'@export
set.keep.strata.sizes.constant <- function(jheem,
                                           fix.strata.sizes,
                                           time=-Inf)
{
    if (is.null(jheem$parameters$CONSTANT_SIZE_DIMENSIONS))
        warning("In order to fix strata sizes, you must also call 'set.fixed.size.strata' to indicate which dimensions define the strata to fix.")

    jheem = set.time.varying.parameter.value(jheem,
                                             parameter.name = 'KEEP_STRATA_SIZE_CONSTANT',
                                             parameter.value = fix.strata.sizes,
                                             time = time,
                                             stepwise=T)

    jheem
}

##----------------------------##
##-- FIXED POPULATION SIZES --##
##----------------------------##


#'@title Sets whether to use fixed population strata for a given time
#'
#'@family functions to fix population strata sizes
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param fix.population A boolean indicated whether to use fixed population strata
#'@param time The time following which the setting applies
#'
#'@seealso \code{\link{set.fixed.population.strata}}, \code{\link{set.fixed.population.size}}
#'
#'@export
set.use.fixed.population <- function(jheem,
                                     fix.population,
                                     time=-Inf)
{
    jheem = set.time.varying.parameter.value(jheem,
                                     parameter.name = 'FIX_STRATA_SIZES',
                                     parameter.value = fix.population,
                                     time = time,
                                     stepwise=T)

    #for now, we assume that if we are fixing population sizes, we are not modeling births
    jheem = set.time.varying.parameter.value(jheem,
                                     parameter.name = 'MODEL_BIRTHS',
                                     parameter.value = !fix.population,
                                     time = time,
                                     stepwise=T)

    jheem
}





#'@title Set which strata to use for fixed size populations
#'
#'@family functions to fix population strata sizes
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param fix.age,fix.race,fix.subpopulation,fix.sex,fix.risk Boolean values indicating whether or not to fix these dimensions
#'
#'@seealso \code{\link{set.use.fixed.population}}, \code{\link{set.fixed.population.size}}
#'
#'@export
set.fixed.population.strata <- function(jheem,
                                        fix.age=F,
                                        fix.race=F,
                                        fix.subpopulation=F,
                                        fix.sex=F,
                                        fix.risk=F
)
{
    #Make sure we haven't already set this
    if (!is.null(jheem$parameters$NUM_FIXED_STRATA) && !is.null(jheem$parameters$time.varying.parameters$TARGET_STRATUM_SIZE))
        stop("Fixed population sizes have already been set using the previously specified population strata. These cannot be reset")

    #Calculate the values
    dim.names = get.dimnames(jheem, age=fix.age, race=fix.race, subpopulation=fix.subpopulation, sex=fix.sex, risk=fix.risk)
    N = prod(sapply(dim.names, length))

    pop.map = array(1:N,
                    dim=sapply(dim.names, length),
                    dimnames=dim.names)

    pop.map.hiv.negative = expand.population.to.hiv.negative(jheem, pop.map)
    pop.map.hiv.positive = expand.population.to.hiv.positive(jheem, pop.map)


    # Set the parameters
    jheem$parameters$NUM_FIXED_STRATA = N
    jheem$parameters$FIXED_STRATUM_HIV = pop.map.hiv.positive
    jheem$parameters$FIXED_STRATUM_NONHIV = pop.map.hiv.negative

    jheem$parameters$FIXED_STRATA = names(dim.names)
#jheem$parameters$FIXED_STRATA = sort(index.of(names(dim.names), names(get.dimnames.general(jheem))))


    # Fixed entry proportions

    first.age.strata = get.hiv.negative.population.skeleton(jheem)
    first.age.strata[1,,,,,] = 1
    fixed.stratum.contains.first.age = as.logical(sapply(1:N, function(i){max(first.age.strata[pop.map.hiv.negative==i])}))

    jheem$parameters$IS_STRATUM_FIXED_PROPORTION = !fixed.stratum.contains.first.age

    # Return
    jheem
}

#'@title Set the population to which to fix sizes
#'
#'
#'@family functions to fix population strata sizes
#'
#'@description Must first call \code{\link{set.fixed.population.strata}}, which indicates which dimensions to fix population along. Calling set.fixed.population.size indicates what levels to fix the population to
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param hiv.negative.population,hiv.positive.population Arrays representing the HIV-negative and HIV-populations to use in calculating the sizes for fixed strata
#'@param time The time at which the fixed population is pegged
#'
#'@seealso \code{\link{set.fixed.population.strata}}, \code{\link{set.use.fixed.population}}
#'
#'@export
set.fixed.population.size <- function(jheem,
                                      hiv.negative.population=jheem$initial.state.hiv.negative,
                                      hiv.positive.population=jheem$initial.state.hiv.positive,
                                      fixed.birth.proportions=create.birth.proportions.from.population(jheem, hiv.negative.population),
                                      time=-Inf

)
{
    # Check to make sure we're set up
    if (is.null(jheem$parameters$NUM_FIXED_STRATA))
        stop("You must call 'set.fixed.population.strata' prior to calling 'set.fixed.population.size'")

    #Check to make sure the populations fit
    needed.dim.names = get.dimnames.nonhiv(jheem)
    if (is.null(hiv.negative.population) || !array.matches.dimnames(hiv.negative.population, needed.dim.names))
        stop("The hiv.negative.population state must have ", length(needed.dim.names), " dimensions: ",
             paste0(paste0(names(needed.dim.names), " (", sapply(needed.dim.names, length), " state", c('','s')[1+(sapply(needed.dim.names, length)>1)] ,")"), collapse=', '))

    needed.dim.names = get.dimnames.hiv(jheem)
    if (is.null(hiv.positive.population) || !array.matches.dimnames(hiv.positive.population, needed.dim.names))
        stop("The hiv.positive.population state must have ", length(needed.dim.names), " dimensions: ",
             paste0(paste0(names(needed.dim.names), " (", sapply(needed.dim.names, length), " state", c('','s')[1+(sapply(needed.dim.names, length)>1)] ,")"), collapse=', '))


    # Calculate the stratum sizes
    target.stratum.size = sapply(1:jheem$parameters$NUM_FIXED_STRATA, function(i){
        sum(hiv.negative.population[jheem$parameters$FIXED_STRATUM_NONHIV==i]) +
            sum(hiv.positive.population[jheem$parameters$FIXED_STRATUM_HIV==i])
    })

    #Calculate the entry proportions
        #only include the population from specified born into strata
#    pop.for.proportions = get.hiv.negative.population.skeleton(jheem)
#    access(pop.for.proportion, race=fixed.births.into.race, subpopulation=fixed.births.into.subpopulation, sex=fixed.births.into.sex, risk=fixed.births.into.risk, non.hiv.subset=fixed.births.into.non.hiv.subset) =
#        access(hiv.negative.population, race=fixed.births.into.race, subpopulation=fixed.births.into.subpopulation, sex=fixed.births.into.sex, risk=fixed.births.into.risk, non.hiv.subset=fixed.births.into.non.hiv.subset)
#    denominator.for.proportions.by.stratum = sapply(1:jheem$parameters$NUM_NONHIV_STATES, function(i){
#        if (jheem$parameters$IS_STRATUM_FIXED_PROPORTION[jheem$parameters$FIXED_STRATUM_NONHIV[i]])
#            1
#        else
#            sum(pop.for.proportions[jheem$parameters$FIXED_STRATUM_NONHIV==jheem$parameters$FIXED_STRATUM_NONHIV[i]])
#    })

#    if (any(denominator.for.proportions.by.stratum==0))
#        stop('You must specify the fixed.births.into.X arguments such that each fixed stratum of the initial hiv negative population is nonzero in size')

#    hiv.negative.proportions = sapply(1:jheem$parameters$NUM_NONHIV_STATES, function(i){
#        pop.for.proportions[i] / denominator.for.proportions.by.stratum[jheem$parameters$FIXED_STRATUM_NONHIV[i]]
#    })

    fixed.prop.dim.name.names = c('race.from','subpopulation.to','sex.to','risk.to','non.hiv.subset.to')
    fixed.birth.proportions = apply(fixed.birth.proportions, fixed.prop.dim.name.names, sum)
    fixed.birth.proportions = expand.population(fixed.birth.proportions, target.dim.names = get.dimnames.by.name(jheem, c('age',fixed.prop.dim.name.names)))
    denominator.for.proportions.by.stratum = sapply(1:jheem$parameters$NUM_NONHIV_STATES, function(i){
        if (jheem$parameters$IS_STRATUM_FIXED_PROPORTION[jheem$parameters$FIXED_STRATUM_NONHIV[i]])
            1
        else
            sum(fixed.birth.proportions[jheem$parameters$FIXED_STRATUM_NONHIV==jheem$parameters$FIXED_STRATUM_NONHIV[i]])
    })

    if (any(denominator.for.proportions.by.stratum==0))
            stop('You must specify the fixed.birth.proportions such that each fixed stratum within fixed.birth.proportions has at least one nonzero component')

    hiv.negative.proportions = sapply(1:jheem$parameters$NUM_NONHIV_STATES, function(i){
            fixed.birth.proportions[i] / denominator.for.proportions.by.stratum[i]
    })

    fixed.entry.proportions = c(rep(0,jheem$parameters$NUM_HIV_STATES),
                                hiv.negative.proportions)

    # Set the parameters
    jheem = set.time.varying.parameter.value(jheem,
                                             parameter.name = 'TARGET_STRATUM_SIZE',
                                             parameter.value = target.stratum.size,
                                             time = time,
                                             stepwise=F)

    jheem = set.time.varying.parameter.value(jheem,
                                             parameter.name = 'FIXED_ENTRY_PROPORTIONS',
                                             parameter.value = fixed.entry.proportions,
                                             time = time,
                                             stepwise=F)

    # Return
    jheem
}

##-----------------------------##
##-- HELPERS ON THE DIMNAMES --##
##-----------------------------##

get.from.hiv.negative.birth.proportion.dimnames <- function(jheem)
{
    get.dimnames(jheem, age.from=T, race.from=T, subpopulation.from=T, sex.from=T, risk.from=T, non.hiv.subset.from=T)
}

get.to.hiv.negative.birth.proportion.dimnames <- function(jheem)
{
    get.dimnames(jheem, subpopulation.to=T, sex.to=T, risk.to=T, non.hiv.subset.to=T)
}

get.from.hiv.positive.birth.proportion.dimnames <- function(jheem)
{
    get.dimnames(jheem,
                 age.from=T, race.from=T, subpopulation.from=T, sex.from=T, risk.from=T,
                 continuum.of.care.from=T, cd4.from=T, hiv.subset.from=T)
}

get.to.hiv.positive.birth.proportion.dimnames <- function(jheem)
{
    get.dimnames(jheem, subpopulation.to=T, sex.to=T, risk.to=T, continuum.of.care.to=T, cd4.to=T, hiv.subset.to=T)
}
