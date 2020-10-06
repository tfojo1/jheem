
##-----------------------##
##-----------------------##
##-- PARAMETER SETTERS --##
##-----------------------##
##-----------------------##

#'@title Sets the global HIV transmission rate
#'
#'@family functions to specify transmission
#'
#'@param jheem An initialized JHEEM object
#'@param transmission.route.names The names of the transmission routes to which the transmissibility should apply (a subset of the 'transmission.route.names' given to \code{\link{initialize.jheem}})
#'@param transmission.rate a single numeric value representing the global transmission rate, which can be used to tune the epidemic across all strata
#'@param time The time at which this transmission rate applies, if time-varying
#'
#'@details The transmission rate (for one route of transmission) from stratum a (hiv-positive) to stratum b (hiv-negative) is a product of five components (four of which are set as parameters, and one of which - HIV prevalence - is determined by running the model):
#'\itemize{
#'  \item {\strong{A global transmission rate} - applies to all strata equally (\code{\link{set.global.transmission.rate}})}
#'  \item {\strong{Transmissibility} - all factors which depend only on the characteristics of the HIV-positive stratum a (\code{\link{set.transmissibility}})}
#'  \item {\strong{Contact rate} - all factors which depend on both the general characteristics of the HIV-positive stratum a (ie age, race, subpopulation, sex, risk factor) and the characteristics of the HIV-negative stratum b (\code{\link{set.transmission.contact.array}})}
#'  \item {\strong{Susceptibility} - all factors which depend only on the characteristics of the HIV-negative stratum b (\code{\link{set.susceptibility}})}
#'  \item {\strong{Prevalence of HIV} within the general characteristics of the HIV-positive stratum a (ie age, race, subpopulation, sex, risk factor). This is not set, but is estiated at model run time}
#'}
#'
#'@export
set.global.transmission.rate <- function(jheem,
                                         transmission.rate,
                                         transmission.route.names=jheem$transmission.routes,
                                         time=-Inf)
{
    check.transmission.route.names(jheem, transmission.route.names)

    if ((!class(transmission.rate)=='numeric' && !class(transmission.rate)=='integer') ||
        length(transmission.rate) > 1 || is.na(transmission.rate))
        stop("transmission.rate must be a scalar, numeric, non-NA value")

    for (transmission.route in transmission.route.names)
    {
        jheem = set.time.varying.parameter.value(jheem,
                                                 parameter.name = paste0('GLOBAL_TRANSMISSION_RATE_', index.of(transmission.route, jheem$transmission.routes)),
                                                 parameter.value = transmission.rate,
                                                 time = time)
    }

    jheem
}

#'@title Sets transmissibility for the HIV-positive population
#'
#'@description Transmissibilty is a characteristic of the HIV-positive population, a multiplier indicating how likely an HIV-positive individual in each stratum is to transmit infection (per partnership) relative to a baseline. ie a value of 1 would not change the transmissibility from the global (baseline) rate
#'
#'@family functions to specify transmission
#'
#'@param transmissibility An array of transmissibility, whose dimensions match the HIV-positive population <age, race, subpopulation, sex, risk factor, continuum of care, cd4, hiv subset>.
#'@param time The time at which these transmissibilities apply, if time-varying
#'
#'@inheritParams set.global.transmission.rate
#'@inherit set.global.transmission.rate details
#'
#'@seealso \code{\link{get.hiv.positive.population.skeleton}}, \code{\link{create.transmissibility.array.from.marginals}}
#'
#'@export
set.transmissibility <- function(jheem,
                                 transmissibility,
                                 transmission.route.names=jheem$transmission.routes,
                                 time=-Inf)
{
    check.transmission.route.names(jheem, transmission.route.names)
    transmissibility = expand.population.to.hiv.positive(jheem, transmissibility)

    for (transmission.route in transmission.route.names)
    {
        jheem = set.time.varying.parameter.value(jheem,
                                     parameter.name = paste0('TRANSMISSIBILITY_', index.of(transmission.route, jheem$transmission.routes)),
                                     parameter.value = transmissibility,
                                     time = time)
    }

    jheem
}

#'@title Sets susceptibility to infection for the HIV-negative population
#'
#'@description Susceptibility is a characteristic of the HIV-negative population, a multiplier indicating how likely an HIV-negative individual in each stratum is to aquire infection (per partnership) relative to a baseline. ie a value of 1 would not change the susceptibility from the global (baseline) rate
#'
#'@family functions to specify transmission
#'
#'@param susceptibility An array of susceptibility rates, whose dimensions match the HIV-negative population [age, race, subpopulation, sex, risk, non.hiv.subset].
#'@param time The time at which these susceptibilities apply, if time-varying
#'
#'@inheritParams set.global.transmission.rate
#'@inherit set.global.transmission.rate details
#'
#'@seealso \code{\link{create.susceptibility.array.from.marginals}}, \code{\link{get.hiv.negative.population.skeleton}}
#'
#'@export
set.susceptibility <- function(jheem,
                               susceptibility,
                               transmission.route.names=jheem$transmission.routes,
                               time=-Inf)
{
    check.transmission.route.names(jheem, transmission.route.names)
    susceptibility = expand.population.to.hiv.negative(jheem, susceptibility)

    for (transmission.route in transmission.route.names)
    {
        jheem = set.time.varying.parameter.value(jheem,
                                     parameter.name = paste0('SUSCEPTIBILITY_', index.of(transmission.route, jheem$transmission.routes)),
                                     parameter.value = susceptibility,
                                     time = time)
    }

    jheem
}

#'@title Sets the transmission contact array between HIV-negative and HIV-positive populations
#'
#'@description A contact array is a function of both the hiv-infected population and hiv-negative population
#'
#'@family functions for transmission contact arrays
#'@family functions to specify transmission
#'
#'@param contact.array An contact array, with some subset of dimensions: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]. from->to follows the direction of HIV transmission
#'@param time The time at which these susceptibilities apply, if time-varying
#'
#'@inheritParams set.global.transmission.rate
#'@inherit set.global.transmission.rate details
#'
#'@seealso \code{\link{create.contact.array.from.marginals}}, \code{\link{get.contact.array.skeleton}}
#'
#'@export
set.transmission.contact.array <- function(jheem,
                                           contact.array,
                                           transmission.route.names=jheem$transmission.routes,
                                           time=-Inf)
{
    check.transmission.route.names(jheem, transmission.route.names)

    contact.array = expand.to.contact.array(jheem, contact.array)
    names.of.dims = names(dimnames(contact.array))

    for (transmission.route in transmission.route.names)
    {
        route.index = index.of(transmission.route, jheem$transmission.routes)
        contact.type = paste0('CONTACT_', route.index)

        contact.type.from.dims = get.contact.from.dim.indices(jheem, contact.array)

        if (length(jheem$parameters$CONTACT_TYPE_FROM_DIMENSION)==0 ||
            is.null(jheem$parameters$CONTACT_TYPE_FROM_DIMENSIONS[[route.index]]))
        {
            jheem$parameters$CONTACT_TYPE_FROM_DIMENSIONS[[route.index]] = contact.type.from.dims
        }
        else if (!all(contact.type.from.dims==jheem$parameters$CONTACT_TYPE_FROM_DIMENSIONS[[route.index]]))
        {
            stop("The 'from' dimensions for the given contact array (",
                 paste0("'", names(contact.type.from.dims), "'", collapse=', '),
                 ") do not match the 'from' dimensions previously specified for transmission route '",
                 transmission.route, "' (",
                 paste0("'", names(jheem$parameters$CONTACT_TYPE_FROM_DIMENSIONS[[route.index]]), "'", collapse=', '),
                 ")")
        }

        if (length(contact.type.from.dims)==0)
            dim(contact.array) = jheem$parameters$NUM_NONHIV_STATES
        else
            dim(contact.array) = c(prod(sapply(get.dimnames.general(jheem)[contact.type.from.dims], length)),
                                   jheem$parameters$NUM_NONHIV_STATES)

        jheem = set.time.varying.parameter.value(jheem,
                                                 parameter.name = contact.type,
                                                 parameter.value = contact.array,
                                                 time = time)
    }

    jheem
}

#'@title Sets the proportions into which new infections from strata of age x race x subpopulation x sex x risk are distributed
#'
#'@family functions handling new infection proportions
#'@family functions to specify transmission
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param proportions An contact array, with some subset of dimensions: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]. from->to follows the direction of HIV transmission
#'@param time The time at which these susceptibilities apply, if time-varying
#'
#'@details A new infection proportions array should have dimensions [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]. For each combination of age x race x subpopulation x sex x risk x non.hiv.subset, proportions[age,race,subpopulation,sex,risk,non.hiv.subset,,,] should sum to 1 and represent the proportions according to which new infections from the age x race x subpopulation x sex x risk x non.hiv.subset are distributed across strata of continuum x cd4 x hiv.subset
#'
#'@export
set.new.infection.proportions <- function(jheem,
                                          proportions,
                                          time=-Inf)
{
    proportions = expand.to.new.infection.proportions(jheem, proportions)

    #check to make sure they are valid (sum to one) proportions (within numerical error)
    sums = rowSums(proportions, dims=6)
    tolerance = 0.001

    if (max(abs(1-sums))>tolerance)
        stop('proportions[age,race,subpopulation,sex,risk,non.hiv.subset,,,] must sum to 1 for each combination of age, race, subpopulation, sex, risk, and non.hiv.subset')
    if (any(proportions < -tolerance))
        stop('proportions must be greater than or equal to zero')

    #normalize proportions to be >=0 and sum to one (in case there was numerical error)
    proportions = apply(proportions, 1:length(dim(proportions)), function(x){pmax(x,0)})
    proportions = proportions / rep(as.numeric(sums), prod(dim(proportions))/prod(dim(sums)))

    set.time.varying.parameter.value(jheem,
                                     parameter.name = 'NEW_INFECTION_PROPORTIONS',
                                     parameter.value = proportions,
                                     time = time)
}




##-------------------------------------------------##
##-------------------------------------------------##
##-- HELPER FUNCTIONS TO SET UP PARAMETER ARRAYS --##
##-------------------------------------------------##
##-------------------------------------------------##


##---------------------------------------------##
##-- TRANSMISSIBILITY ARRAY SET-UP FUNCTIONS --##
##---------------------------------------------##

#'@title Get a skeleton transmissibility array
#'
#'@family functions for building transmissibility array
#'@family generators of population skeletons
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the HIV-positive population: [age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset]
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param value The value to fill the array with
#'
#'@seealso \code{\link{set.transmissibility}}
#'
#'@export
get.transmissibility.array.skeleton <- function(jheem, value=0)
{
    get.hiv.positive.population.skeleton(jheem, value=value)
}

#'@title Makes a transmissibility array by multiplying multiple marginal arrays
#'
#'@inherit set.transmissibility description
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param ... One or more arrays whose dimensions are a subset of the dimensions of the HIV-positive population: [age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset]
#'
#'@seealso \code{\link{set.transmissibility}}, \code{\link{get.hiv.positive.population.skeleton}}
#'
#'@export
create.transmissibility.array.from.marginals <- function(jheem,
                                                         ...)
{
    create.product.array.from.marginals(get.dimnames.hiv(jheem), ...)
}


##-------------------------------------------##
##-- SUSCEPTIBILITY ARRAY SET-UP FUNCTIONS --##
##-------------------------------------------##


#'@title Get a skeleton susceptibility array
#'
#'@family functions for building susceptibility array
#'@family generators of population skeletons
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the HIV-negative population: [age, race, subpopulation, sex, risk, non.hiv.subset]
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param value The value to fill the array with
#'
#'@seealso \code{\link{set.transmissibility}}
#'
#'@export
get.susceptibility.array.skeleton <- function(jheem, value=0)
{
    get.hiv.negative.population.skeleton(jheem, value=value)
}

#'@title Makes a susceptibility array by multiplying multiple marginal arrays
#'
#'@inherit set.susceptibility description
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param ... One or more arrays whose dimensions are a subset of the dimensions of the HIV-negative population: [age, race, subpopulation, sex, risk, non.hiv.subset]
#'
#'@seealso \code{\link{set.susceptibility}}, \code{\link{get.hiv.negative.population.skeleton}}
#'
#'@export
create.susceptibility.array.from.marginals <- function(jheem,
                                                       ...)
{
    create.product.array.from.marginals(get.dimnames.nonhiv(jheem), ...)
}


##-------------------------------------------------##
##-- TRANSMISSION CONTACT ARRAY SET-UP FUNCTIONS --##
##-------------------------------------------------##


#'@title Create a transmission contact array from one or more marginal contact arrays
#'
#'@family functions for transmission contact arrays
#'@family functions to build arrays from one or more marginal arrays
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param ... One or more arrays whose dimensions are a subset of: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]
#'
#'@seealso \code{\link{set.transmission.contact.array}}, \code{\link{get.contact.array.skeleton}}
#'
#'@export
create.contact.array.from.marginals <- function(jheem,
                                                ...)
{
    marginals = list(...)

    #Find out which dimensions are present in any of the given marginals
    # Match dimensions by name, if names are present
    # Otherwise, match by the dimnames
    #  Assume that if a set of dimnames occurs twice in marginal, it is 'from' and 'to' in that order
    #  Assume that if a set of dimnames occurs only once in marginal, it is 'to'

    dim.names.present = list()

    for (i in 1:length(marginals))
    {
        dim.name.names = identify.dimnames.in.array(jheem, marginals[[i]], single.occurence.as='.to')
        marginals[[i]] = set.dim.name.names(marginals[[i]], dim.name.names)

        dim.names.present = union(dim.names.present, names(dimnames(marginals[[i]])))
    }

    #Build the return value skeleton
    rv = get.contact.array.skeleton(jheem, 1,
                                    age.from = any(dim.names.present=='age.from'),
                                    race.from = any(dim.names.present=='race.from'),
                                    subpopulation.from = any(dim.names.present=='subpopulation.from'),
                                    sex.from = any(dim.names.present=='sex.from'),
                                    risk.from = any(dim.names.present=='risk.from'),
                                    age.to = any(dim.names.present=='age.to'),
                                    race.to = any(dim.names.present=='race.to'),
                                    subpopulation.to = any(dim.names.present=='subpopulation.to'),
                                    sex.to = any(dim.names.present=='sex.to'),
                                    risk.to = any(dim.names.present=='risk.to'),
                                    non.hiv.subset.to = any(dim.names.present=='non.hiv.subset.to'))
    target.dim.names = dimnames(rv)

    #for each marginal, expand up the population, then multiply it into the return value
    for (marginal.array in marginals)
    {
        expanded.marginal.array = expand.population(marginal.array, target.dim.names)
        rv = rv * expanded.marginal.array
    }

    rv
}


#'@title Multiply a transmission contact array by rates for the 'from' population
#'
#'@family functions for transmission contact arrays
#'
#'@description Fundamentally, a contact array can be thought of as taking the outer product of a 'from array' with a 'to array'. This function equates to multiplying the from array portion prior to taking the outer product with the to array portion
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param contact.array A contact array, with dimensions: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]
#'@param from.multiplier An array with dimensions matching the general population: [age, race, subpopulation, sex, risk]
#'
#'@inherit set.new.infection.proportions details
#'
#'@export
multiply.contact.array.from <- function(jheem,
                                        contact.array,
                                        from.multiplier)
{
    stop('Not implemented')
    from.names = get.dimnames.general(jheem)[get.contact.from.dim.indices(jheem, contact.array)]
    from.multiplier = expand.population(from.multiplier, from.names)

    array(as.numeric(contact.array) *
              rep(as.numeric(from.multiplier), prod(dim(contact.array))/prod(dim(from.multiplier))),
          dim=dim(contact.array), dimnames=dimnames(contact.array))
}

#'@title Multiply a transmission contact array by rates for the 'to' population
#'
#'@family functions for transmission contact arrays
#'
#'@description Fundamentally, a contact array can be thought of as taking the outer product of a 'from array' with a 'to array'. This function equates to multiplying the to array portion prior to taking the outer product with the from array portion
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param contact.array A contact array, with dimensions: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]
#'@param from.multiplier An array with dimensions matching the HIV-negative population: [age, race, subpopulation, sex, risk, non.hiv.subset]
#'
#'@export
multiply.contact.array.to <- function(jheem,
                                      contact.array,
                                      to.multiplier)
{
    stop('Not implemented')

    to.multiplier = expand.population.to.hiv.negative(jheem, to.multiplier)

    array(as.numeric(contact.array) *
              rep(as.numeric(to.multiplier), each=prod(dim(contact.array))/prod(dim(to.multiplier))),
          dim=dim(contact.array), dimnames=dimnames(contact.array))
}


##-------------------------------------------------##
##-- TRANSMISSION CONTACT ARRAY SET-UP FUNCTIONS --##
##-------------------------------------------------##


#'@title Get a skeleton new infection proportions array
#'
#'@family functions handling new infection proportions
#'
#'@description Gets an array, with all values set to value, with dimensions and names corresponding to the new infection proportions array: [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param value The value to fill the array with.
#'
#'@inherit set.new.infection.proportions details
#'
#'@export
get.new.infection.proportions.skeleton <- function(jheem, value=0)
{
    get.population.skeleton(jheem, value=value,
                            age = T, race = T, subpopulation = T, sex = T, risk = T,
                            non.hiv.subset = T, continuum = T, cd4 = T, hiv.subset = T)
}

#'@title Get a new infection proportions array where new infections are seeded uniformly
#'
#'@family functions handling new infection proportions
#'
#'@description The return value will be a new infection proportions array such that the new infection proportions are the same across all strata of age x race x subpopulation x sex x risk x non.hiv.subset, with new infections being distributed uniformly across the specified continuum, cd4, and hiv.subset specified
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param initial.continuum,initial.cd4,initial.hiv.subset The index (indices) or name(s) of the dimensions into which new infections are to be distributed uniformly
#'
#'@inherit set.new.infection.proportions details
#'
#'@export
get.uniform.new.infection.proportions <- function(jheem,
                                                  initial.continuum=jheem$continuum.of.care[1],
                                                  initial.cd4=jheem$cd4.strata[1],
                                                  initial.hiv.subset=jheem$hiv.subsets[1])
{
    sub.array = get.population.skeleton(jheem, continuum = T, cd4 = T, hiv.subset = T)
    access(sub.array, continuum=initial.continuum, cd4=initial.cd4, hiv.subset=initial.hiv.subset) =
        1 / length(initial.continuum) / length(initial.cd4) / length(initial.hiv.subset)
    rv = expand.to.new.infection.proportions(jheem, sub.array)
    rv[abs(rv)<0.001] = 0

    rv
}



##----------------------##
##----------------------##
##-- HELPER FUNCTIONS --##
##----------------------##
##----------------------##

#A helper that makes an array by expanding each array in ... and multiplying it with all the other expanded marginals
create.product.array.from.marginals <- function(target.dim.names, ...)
{
    marginals = list(...)
    rv = array(1, dim=sapply(target.dim.names, length), dimnames=target.dim.names)

    for (marginal.array in marginals)
    {
        expanded.marginal.array = expand.population(marginal.array, target.dim.names)
        rv = rv * expanded.marginal.array
    }

    rv
}

