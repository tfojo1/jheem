#'@title Set whether the JHEEM should track HIV-specific and overall mortality
#'
#'@family function to set up JHEEM mortality
#'
#'@description A convenience function that bundles \code{\link{set.track.hiv.specific.mortality}} and \code{\link{set.track.overall.mortality}} into one function call
#'
#'@inheritParams set.general.mortality
#'@param track.hiv.specific.mortality,track.overall.hiv.positive.mortality,track.overall.hiv.negative.mortality Booleans indicating whether or not to track each type for mortality
#'@param dimensions.hiv.specific.mortality,dimensions.overall.hiv.positive.mortality The names of dimensions to keep when tracking hiv-positive mortality. Must be a subset of <'age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset'>
#'@param dimensions.overall.hiv.negative.mortality The names of dimensions to keep when tracking HIV-negative mortality. Must be a subset of <'age','race','subpopulation','sex','risk','non.hiv.subset'>
#'
#'@export
set.track.mortality <- function(jheem,
                                track.hiv.specific.mortality=T,
                                track.overall.hiv.positive.mortality=T,
                                track.overall.hiv.negative.mortality=T,
                                dimensions.hiv.specific.mortality = c('age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset'),
                                dimensions.overall.hiv.positive.mortality = c('age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset'),
                                dimensions.overall.hiv.negative.mortality = c('age','race','subpopulation','sex','risk','non.hiv.subset'))
{
    jheem = set.track.hiv.specific.mortality(jheem, track.hiv.specific.mortality,
                                             dimensions.hiv.specific.mortality = dimensions.hiv.specific.mortality)
    jheem = set.track.overall.mortality(jheem,
                                        track.hiv.positive = track.overall.hiv.positive.mortality,
                                        track.hiv.negative = track.overall.hiv.negative.mortality,
                                        dimensions.overall.hiv.positive.mortality = dimensions.overall.hiv.positive.mortality,
                                        dimensions.overall.hiv.negative.mortality = dimensions.overall.hiv.negative.mortality)

    jheem
}

#'@title Set whether the JHEEM should track HIV-specific mortality
#'
#'@family function to set up JHEEM mortality
#'
#'@inheritParams set.general.mortality
#'@param track Boolean indicating whether or not to track HIV-specific mortality
#'
#'@export
set.track.hiv.specific.mortality <- function(jheem, track,
                                             dimensions.hiv.specific.mortality = c('age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset'))
{
    dimensions.hiv.specific.mortality = tolower(dimensions.hiv.specific.mortality)
    allowed = names(get.dimnames.hiv(jheem))
    not.allowed = setdiff(dimensions.hiv.specific.mortality, allowed)
    if (length(not.allowed)>0)
        stop(paste0("dimensions.hiv.specific.mortality must be a subset of <",
                    paste0("'", allowed, "'", collapse=', '),
                    ">. The following are not allowed: ",
                    paste0("'", not.allowed, "'", collapse=', ')))

    jheem$parameters$hiv.specific.mortality.keep.dimensions = intersect(allowed, dimensions.hiv.specific.mortality)

    jheem$parameters$TRACK_HIV_SPECIFIC_MORTALITY = track

    jheem
}

#'@title Set whether the JHEEM should track general mortality
#'
#'@details Overall mortality, for HIV-positive, is the sum of HIV-specific mortality and non-HIV mortality
#'
#'@family function to set up JHEEM mortality
#'
#'@inheritParams set.general.mortality
#'@inheritParams set.track.mortality
#'@param track.hiv.positive,track.hiv.negative Boolean indicators for whether or not whether or not to track overall mortality for HIV-positive and HIV-negative
#'
#'@export
set.track.overall.mortality <- function(jheem, track.hiv.positive, track.hiv.negative,
                                        dimensions.overall.hiv.positive.mortality = c('age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset'),
                                        dimensions.overall.hiv.negative.mortality = c('age','race','subpopulation','sex','risk','non.hiv.subset'))
{

    dimensions.overall.hiv.positive.mortality = tolower(dimensions.overall.hiv.positive.mortality)
    allowed = names(get.dimnames.hiv(jheem))
    not.allowed = setdiff(dimensions.overall.hiv.positive.mortality, allowed)
    if (length(not.allowed)>0)
        stop(paste0("dimensions.overall.hiv.positive.mortality must be a subset of <",
                    paste0("'", allowed, "'", collapse=', '),
                    ">. The following are not allowed: ",
                    paste0("'", not.allowed, "'", collapse=', ')))

    jheem$parameters$overall.hiv.mortality.keep.dimensions = intersect(allowed, dimensions.overall.hiv.positive.mortality)

    dimensions.overall.hiv.negative.mortality = tolower(dimensions.overall.hiv.negative.mortality)
    allowed = names(get.dimnames.nonhiv(jheem))
    not.allowed = setdiff(dimensions.overall.hiv.negative.mortality, allowed)
    if (length(not.allowed)>0)
        stop(paste0("dimensions.overall.hiv.negative.mortality must be a subset of <",
                    paste0("'", allowed, "'", collapse=', '),
                    ">. The following are not allowed: ",
                    paste0("'", not.allowed, "'", collapse=', ')))

    jheem$parameters$overall.nonhiv.mortality.keep.dimensions = intersect(allowed, dimensions.overall.hiv.negative.mortality)

    jheem$parameters$TRACK_OVERALL_HIV_MORTALITY = track.hiv.positive
    jheem$parameters$TRACK_OVERALL_NONHIV_MORTALITY = track.hiv.negative

    jheem
}


#'@title Set which dimensions to keep when tracking incidence
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param dimensions The names of dimensions to keep when tracking incidence. Must be a subset of <'age','race','subpopulation','sex','risk','non.hiv.subset','continuum','cd4','hiv.subset'>
#'
#'@export
set.track.incidence.dimensions <- function(jheem,
                                           dimensions=c('age','race','subpopulation','sex','risk','non.hiv.subset','continuum','cd4','hiv.subset'))
{
    dimensions = tolower(dimensions)
    allowed = names(get.dimnames.all(jheem))
    not.allowed = setdiff(dimensions, allowed)
    if (length(not.allowed)>0)
        stop(paste0("dimensions must be a subset of <",
                    paste0("'", allowed, "'", collapse=', '),
                    ">. The following are not allowed: ",
                    paste0("'", not.allowed, "'", collapse=', ')))

    jheem$parameters$incidence.keep.dimensions = intersect(allowed, dimensions)

    jheem
}

get.collapse.indices <- function(jheem,
                                 full.dimension.names,
                                 collapsed.dimension.names)
{
    full.dimensions = get.dimnames.by.name(jheem, full.dimension.names)
    collapsed.dimension.names = get.dimnames.by.name(jheem, collapsed.dimension.names)

    collapsed.indices = array(1:prod(sapply(collapsed.dimension.names, length)),
                              dim = sapply(collapsed.dimension.names, length),
                              dimnames = collapsed.dimension.names)

    expand.population(source.population = collapsed.indices, target.dim.names = full.dimensions)
}
