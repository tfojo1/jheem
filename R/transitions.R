
##----------------------------------------------##
##-- CHECKING AND EXPANDING TRANSITION ARRAYS --##
##----------------------------------------------##

check.transition.dimension <- function(jheem, transition.dimension, hiv.negative=T,
                                       allowed.dimname.names=NULL)
{
    if (class(transition.dimension) != 'character' || length(transition.dimension) != 1)
        stop("transition.dimension must a single character value")

    if (is.null(allowed.dimname.names))
    {
        if (hiv.negative)
            allowed.dimname.names = names(get.dimnames.nonhiv(jheem))
        else
            allowed.dimname.names = names(get.dimnames.hiv(jheem))
    }

    allowed.dimname.names = setdiff(allowed.dimname.names, c('age','race'))

    if (all(allowed.dimname.names != transition.dimension))
        stop(paste0("transition.dimension must be one of: ",
                    paste0("'", allowed.dimname.names, "'", collapse=', ')))
}

get.transition.dimnames <- function(jheem, transition.dimension, hiv.negative=T,
                                    dim.names=NULL)
{
    if (is.null(dim.names))
    {
        if (hiv.negative)
            dim.names = get.dimnames.nonhiv(jheem)
        else
            dim.names = get.dimnames.hiv(jheem)
    }

    to.dim = dim.names[names(dim.names)==transition.dimension]

    names(dim.names)[names(dim.names)==transition.dimension] = paste0(transition.dimension, '.from')

    dim.names = c(dim.names, to.dim)
    names(dim.names)[length(dim.names)] = paste0(transition.dimension, '.to')

    dim.names
}

check.valid.transition.array.subset <- function(jheem, arr, transition.dimension,
                                                hiv.negative)
{
    check.transition.dimension(jheem, transition.dimension, hiv.negative)

    allowed.dim.name.names = names(get.transition.dimnames(jheem, hiv.negative = hiv.negative, transition.dimension = transition.dimension))


    if (is.null(names(dimnames(arr))))
        arr.dim.name.names = identify.dimnames.in.array(jheem, arr, single.occurence.as=T)
    else
        arr.dim.name.names = names(dimnames(arr))


    if (all(arr.dim.name.names!=paste0(transition.dimension, '.from')) ||
        all(arr.dim.name.names!=paste0(transition.dimension, '.to')))
        stop(paste0("A transition array for '", transition.dimension,
                    "' must contain dimensions '", transition.dimension,
                    ".from', and '", transition.dimension,
                    ".to'"))

    not.allowed.mask = sapply(arr.dim.name.names, function(name){
        all(allowed.dim.name.names != name)
    })

    if (any(not.allowed.mask))
        stop(paste0("Dimensions for a transition array for '",
                    transition.dimension,
                    "' must be a subset of: ",
                    paste0("'", allowed.dim.name.names, "'", collapse=', ')))
}

#'@title Expands an array into a transition array for the HIV-negative population for a given dimension
#'
#'@family Managing transition arrays
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param to.expand The array to expand. At minimum, must have .from and .to dimensions for the the transition dimension. May also have age, race, subpopulation, sex, risk, non.hiv.subset
#'@param transition.dimension The name of the dimension for which this is a transition array. Must be one of age, race, subpopulation, sex, risk, non.hiv.subset
#'
#'@export
expand.to.transition.array.hiv.negative <- function(jheem, to.expand, transition.dimension)
{
    check.valid.transition.array.subset(jheem, arr=to.expand,
                                        transition.dimension=transition.dimension,
                                        hiv.negative=T)

    expand.population(to.expand, get.transition.dimnames(jheem, transition.dimension=transition.dimension, hiv.negative=T))
}

#'@title Expands an array into a transition array for the HIV-positive population for a given dimension
#'
#'@family Managing transition arrays
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param to.expand The array to expand. At minimum, must have .from and .to dimensions for the the transition dimension. May also have age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset
#'@param transition.dimension The name of the dimension for which this is a transition array. Must be one of age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset
#'
#'@export
expand.to.transition.array.hiv.positive <- function(jheem, to.expand, transition.dimension)
{

    check.valid.transition.array.subset(jheem, arr=to.expand,
                                        transition.dimension=transition.dimension,
                                        hiv.negative=F)

    expand.population(to.expand, get.transition.dimnames(jheem, transition.dimension=transition.dimension, hiv.negative=F))
}


##-----------------------------------------##
##-- CREATING TRANSITION ARRAY SKELETONS --##
##-----------------------------------------##

#'@title Create a skeleton transition array for the HIV-negative population
#'
#'@param transition.dimension The name of the dimension for which this is a transition array. Must be one of age, race, subpopulation, sex, risk, non.hiv.subset
#'@inheritParams get.population.skeleton
#'
#'@export
get.hiv.negative.transition.array.skeleton <- function(jheem, transition.dimension, value=0)
{
    check.transition.dimension(jheem, transition.dimension=transition.dimension, hiv.negative = T)

    do.get.population.skeleton( get.transition.dimnames(jheem, transition.dimension, hiv.negative=T), value=value)
}

#'@title Create a skeleton transition array for the HIV-positive population
#'
#'@param transition.dimension The name of the dimension for which this is a transition array. Must be one of age, race, subpopulation, sex, risk, continuum, cd4, hiv.subset
#'@inheritParams get.population.skeleton
#'
#'@export
get.hiv.positive.transition.array.skeleton <- function(jheem, transition.dimension, value=0)
{
    check.transition.dimension(jheem, transition.dimension=transition.dimension, hiv.negative = F)

    do.get.population.skeleton(get.transition.dimnames(jheem, transition.dimension, hiv.negative=F), value=value)
}

#'@title Create a skeleton transition array for a general population (which applies to both HIV-negative and HIV-positive)
#'
#'@param transition.dimension The name of the dimension for which this is a transition array.
#'@inheritParams get.population.skeleton
#'
#'@export
get.general.transition.array.skeleton <- function(jheem,
                                                  transition.dimension,
                                                  value=0)
{
    dim.names = get.dimnames.general(jheem)
    check.transition.dimension(jheem, transition.dimension=transition.dimension,
                               allowed.dimname.names = names(dim.names))

    do.get.population.skeleton(get.transition.dimnames(jheem, transition.dimension, dim.names = dim.names),
                               value=value)
}

#'@title Create a skeleton transition array with arbitrary dimensions
#'
#'@param transition.dimension The name of the dimension for which this is a transition array.
#'@param age,race,subpopulation,sex,risk,non.hiv.subset,continuum,cd4,hiv.subset Whether to include these dimensions in the skeleton array
#'@inheritParams get.population.skeleton
#'
#'@export
get.transition.array.skeleton <- function(jheem,
                                          transition.dimension,
                                          value=0,
                                          age=F,
                                          race=F,
                                          subpopulation=F,
                                          sex=F,
                                          risk=F,
                                          non.hiv.subset=F,
                                          continuum=F,
                                          cd4=F,
                                          hiv.subset=F)
{
    dim.names = get.dimnames(jheem,
                             age=age, race=race,
                             subpopulation=subpopulation || transition.dimension=='subpopulation',
                             sex=sex || transition.dimension=='sex',
                             risk=risk || transition.dimension=='risk',
                             non.hiv.subset=non.hiv.subset || transition.dimension=='non.hiv.subset',
                             continuum=continuum || transition.dimension=='continuum',
                             cd4=cd4 || transition.dimension=='cd4',
                             hiv.subset=hiv.subset || transition.dimension=='hiv.subset')

    check.transition.dimension(jheem, transition.dimension=transition.dimension,
                               allowed.dimname.names = names(dim.names))

    do.get.population.skeleton(get.transition.dimnames(jheem, transition.dimension, dim.names = dim.names),
                               value=value)
}

##----------------------------------##
##-- THE TRANSITION ARRAY SETTERS --##
##----------------------------------##

#'@title Sets the array of transition rates between states (within strata of age and race) for the HIV-negative population
#'
#'@family JHEEM parameter setters
#'@family functions for building transition arrays
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param transition.array A transition array
#'@param time The time at which these transition rates apply, if time-varying
#'@inheritParams expand.to.transition.array.hiv.negative
#'
#'@inherit get.transition.array.skeleton details
#'
#'@export
set.transition.array.hiv.negative <- function(jheem,
                                              transition.array,
                                              time=-Inf)
{
    dim.name.names = names(dimnames(transition.array))
    if (!any(grepl('.from', dim.name.names)))
        stop("To be valid, transition.array must have '.to' and '.from' dimensions")
    transition.dimension = tolower(gsub('.from', '', dim.name.names[grepl('.from', dim.name.names)]))[1]
    transition.array = expand.to.transition.array.hiv.negative(jheem, to.expand=transition.array,
                                                               transition.dimension = transition.dimension)


    jheem$parameters$hiv.negative.transition.dimensions = unique(c(jheem$parameters$hiv.negative.transition.dimensions,
                                                                   transition.dimension))

    set.time.varying.parameter.value(jheem,
                                     parameter.name = paste0('HIV_NEGATIVE_TRANSITION_ARRAY_', toupper(transition.dimension)),
                                     parameter.value = transition.array,
                                     time = time)
}

#'@title Sets the array of transition rates between states (within strata of age and race) for the HIV-positive population
#'
#'@family JHEEM parameter setters
#'@family functions for building transition arrays
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param transition.array A transition array
#'@param time The time at which these transition rates apply, if time-varying
#'
#'@inherit get.transition.array.skeleton details
#'
#'@export
set.transition.array.hiv.positive <- function(jheem,
                                              transition.array,
                                              time=-Inf)
{
    dim.name.names = names(dimnames(transition.array))
    if (!any(grepl('.from', dim.name.names)))
        stop("To be valid, transition.array must have '.to' and '.from' dimensions")
    transition.dimension = tolower(gsub('.from', '', dim.name.names[grepl('.from', dim.name.names)]))[1]
    transition.array = expand.to.transition.array.hiv.positive(jheem, to.expand=transition.array,
                                                               transition.dimension = transition.dimension)


    jheem$parameters$hiv.positive.transition.dimensions = unique(c(jheem$parameters$hiv.positive.transition.dimensions,
                                                                   transition.dimension))

    set.time.varying.parameter.value(jheem,
                                     parameter.name = paste0('HIV_POSITIVE_TRANSITION_ARRAY_', toupper(transition.dimension)),
                                     parameter.value = transition.array,
                                     time = time)
}


##-----------------------------##
##-- SET TRACKED TRANSITIONS --##
##-----------------------------##

#'@title Specify an HIV-negative transition to track
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param name The name by which to refer to this tracked transition
#'@param transition.dimension The dimension to which this transition applies
#'@param from,to The names of the to and from states within the transition dimension that define the transition
#'@param age,race,subpopulation,sex,risk,non.hiv.subset,cd4,continuum,hiv.subset The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension
#'@param keep.dimensions Which dimensions to keep in the tracked array (non-keep dimensions are marginalized out)
#'
#'@export
set.tracked.transition.hiv.negative <- function(jheem, name,
                                                transition.dimension,
                                                from, to,
                                                keep.dimensions=c('age','race','subpopulation','sex','risk','non.hiv.subset'),
                                                age=NULL, race=NULL,
                                                subpopulation=NULL, sex=NULL, risk=NULL, non.hiv.subset=NULL)
{
    #Check transition.dimension
    check.transition.dimension(jheem, transition.dimension = transition.dimension, hiv.negative = T)

    #Make sure we have not already used this tracked name
    if (any(names(jheem$tracked.transitions$dimnames.for.tracked)==name))
        stop(paste0("The name '", name, "' is already in use for a tracked transition"))

    #Set up to and from dimensions to index into the arrays
    dd = list(age=age, race=race, subpopulation=subpopulation, sex=sex, risk=risk,
              non.hiv.subset=non.hiv.subset)
    dd[[paste0(transition.dimension, '.from')]] = from
    dd[[paste0(transition.dimension, '.to')]] = to
    dd[[transition.dimension]] = NULL

    #Flag which transitions we are going to track
    tracked.transition.mask = get.hiv.negative.transition.array.skeleton(jheem,
                                                                         transition.dimension = transition.dimension,
                                                                         value = F)
    access(tracked.transition.mask,
                      age=dd$age, race=dd$race,
                      subpopulation.from=dd$subpopulation.from, sex.from=dd$sex.from, risk.from=dd$risk.from, non.hiv.subset.from=dd$non.hiv.subset.from,
                      subpopulation.to=dd$subpopulation.to, sex.to=dd$sex.to, risk.to=dd$risk.to, non.hiv.subset.to=dd$non.hiv.subset.to) = T


    #Pull an array of the tracked transitions to pattern our hydrated value on
    tracked.transition.arr = access(tracked.transition.mask,
                                               age=dd$age, race=dd$race,
                                               subpopulation.from=dd$subpopulation.from, sex.from=dd$sex.from, risk.from=dd$risk.from, non.hiv.subset.from=dd$non.hiv.subset.from,
                                               subpopulation.to=dd$subpopulation.to, sex.to=dd$sex.to, risk.to=dd$risk.to, non.hiv.subset.to=dd$non.hiv.subset.to,
                                               collapse.length.one.dimensions = F)

    #Set up indices to pass to sub-function or subset
    population.indices = get.hiv.negative.population.skeleton(jheem, value='sequential')
    transition.indices = get.hiv.negative.transition.array.skeleton(jheem,
                                                                    transition.dimension=transition.dimension,
                                                                    value='sequential')

    #Get the order we'll need to put the indices into so they are hydrated properly into an array
    o = get.tracked.transition.order(tracked.transition.mask = tracked.transition.mask,
                                     tracked.transition.arr = tracked.transition.arr)

    #Pull the population indices we are tracking using the mask
    tracked.from = get.tracked.transition.from.indices(population.indices, tracked.transition.mask)[o]
    tracked.trans = as.integer(transition.indices[tracked.transition.mask])[o]

    #Save into the JHEEM
    # First, fold the indices into the parameters
    jheem$parameters$TRACKED_NONHIV_TRANSITIONS_FROM = c(jheem$parameters$TRACKED_NONHIV_TRANSITIONS_FROM, tracked.from)
    jheem$parameters$TRACKED_NONHIV_TRANSITION_INDICES = c(jheem$parameters$TRACKED_NONHIV_TRANSITION_INDICES, tracked.trans)

    # Check keep dimensions
    allowed = names(get.transition.dimnames(jheem, transition.dimension=transition.dimension, hiv.negative = T))
    not.allowed = setdiff(keep.dimensions, allowed)
    if (length(not.allowed)>0)
        stop(paste0("keep.dimensions for an HIV-negative transition array for '",
                    transition.dimension, "' must be a subset of <",
                    paste0("'", allowed, "'", collapse=', '),
                    ">. The following are not allowed: ",
                    paste0("'", not.allowed, "'", collapse=', ')))

    # Now set up the collapsed array for indices and names
    collapsed.arr = apply(tracked.transition.arr, keep.dimensions, sum)
    n.collapsed = length(collapsed.arr)
    collapsed.arr[] = 1:n.collapsed

    collapse.indices = expand.population(collapsed.arr,
                                         target.dim.names = dimnames(tracked.transition.arr))

    jheem$parameters$TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES = collapse.indices +
        jheem$parameters$N_COLLAPSED_TRACKED_NONHIV_TRANSITIONS

    jheem$parameters$N_COLLAPSED_TRACKED_NONHIV_TRANSITIONS =
        max(jheem$parameters$TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES)

    collapsed.dim.names = lapply(names(dimnames(tracked.transition.arr)), function(one.dim.names.name){
        one.dim.names = dimnames(tracked.transition.arr)[[one.dim.names.name]]
        if (length(one.dim.names)==1 || any(keep.dimensions == one.dim.names.name))
            one.dim.names
        else if (length(one.dim.names)==length(dimnames(tracked.transition.mask)[[one.dim.names.name]]))
            'all'
        else
            paste0(one.dim.names, collapse = ';')
    })
    names(collapsed.dim.names) = names(dimnames(tracked.transition.arr))

    # Then save the names for those indices into the JHEEM$tracked.transitions
    jheem$tracked.transitions$names.for.collapsed.nonhiv = c(jheem$tracked.transitions$names.for.indices.nonhiv,
                                                           rep(name, n.collapsed))
    jheem$tracked.transitions$names.for.indices.nonhiv = c(jheem$tracked.transitions$names.for.indices.nonhiv,
                                                           rep(name, length(tracked.from)))

    # Last, save the transition dimension and the dimnames for re-hydrating this tracked transition
    jheem$tracked.transitions$transition.dimensions[name] = tolower(transition.dimension)
    jheem$tracked.transitions$dimnames.for.tracked[[name]] = collapsed.dim.names


    #Return it
    jheem
}

#'@title Specify an HIV-positive transition to track
#'
#'@inheritParams set.tracked.transition.hiv.negative
#'
#'@export
set.tracked.transition.hiv.positive <- function(jheem, name,
                                                transition.dimension,
                                                from=NULL, to=NULL,
                                                keep.dimensions=c('age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset'),
                                                age=NULL, race=NULL,
                                                subpopulation=NULL, sex=NULL, risk=NULL,
                                                continuum=NULL, cd4=NULL, hiv.subset=NULL)
{
    #Check transition.dimension
    check.transition.dimension(jheem, transition.dimension = transition.dimension, hiv.negative = F)

    #Make sure we have not already used this tracked name
    if (any(names(jheem$tracked.transitions$dimnames.for.tracked)==name))
        stop(paste0("The name '", name, "' is already in use for a tracked transition"))

    #Set up to and from dimensions to index into the arrays
    dd = list(age=age, race=race, subpopulation=subpopulation, sex=sex, risk=risk,
              continuum=continuum, cd4=cd4, hiv.subset=hiv.subset)
    dd[[paste0(transition.dimension, '.from')]] = from
    dd[[paste0(transition.dimension, '.to')]] = to
    dd[[transition.dimension]] = NULL

    #Flag which transitions we are going to track
    tracked.transition.mask = get.hiv.positive.transition.array.skeleton(jheem,
                                                                         transition.dimension = transition.dimension,
                                                                         value = F)
    access(tracked.transition.mask,
           age=dd$age, race=dd$race,
           subpopulation=dd$subpopulation, sex=dd$sex, risk=dd$risk,
           continuum=dd$continuum, cd4=dd$cd4, hiv.subset=dd$hiv.subset,
           subpopulation.from=dd$subpopulation.from, sex.from=dd$sex.from, risk.from=dd$risk.from,
           continuum.from=dd$continuum.from, cd4.from=dd$cd4.from, hiv.subset.from=dd$hiv.subset.from,
           subpopulation.to=dd$subpopulation.to, sex.to=dd$sex.to, risk.to=dd$risk.to,
           continuum.to=dd$continuum.to, cd4.to=dd$cd4.to, hiv.subset.to=dd$hiv.subset.to) = T


    #Pull an array of the tracked transitions to pattern our hydrated value on
    tracked.transition.arr = access(tracked.transition.mask,
                                               age=dd$age, race=dd$race,
                                               subpopulation.from=dd$subpopulation.from, sex.from=dd$sex.from, risk.from=dd$risk.from,
                                               continuum.from=dd$continuum.from, cd4.from=dd$cd4.from, hiv.subset.from=dd$hiv.subset.from,
                                               subpopulation.to=dd$subpopulation.to, sex.to=dd$sex.to, risk.to=dd$risk.to,
                                               continuum.to=dd$continuum.to, cd4.to=dd$cd4.to, hiv.subset.to=dd$hiv.subset.to,
                                               collapse.length.one.dimensions = F)

    #Set up indices to pass to sub-function or subset
    population.indices = get.hiv.positive.population.skeleton(jheem, value='sequential')
    transition.indices = get.hiv.positive.transition.array.skeleton(jheem,
                                                                    transition.dimension = transition.dimension,
                                                                    value='sequential')

    #Get the order we'll need to put the indices into so they are hydrated properly into an array
    o = get.tracked.transition.order(tracked.transition.mask = tracked.transition.mask,
                                     tracked.transition.arr = tracked.transition.arr)

    #Pull the population indices we are tracking using the mask
    tracked.from = get.tracked.transition.from.indices(population.indices, tracked.transition.mask)[o]
    tracked.trans = as.integer(transition.indices[tracked.transition.mask])[o]

    #Save into the JHEEM
    # First, fold the indices into the parameters
    jheem$parameters$TRACKED_HIV_TRANSITIONS_FROM = c(jheem$parameters$TRACKED_HIV_TRANSITIONS_FROM, tracked.from)
    jheem$parameters$TRACKED_HIV_TRANSITION_INDICES = c(jheem$parameters$TRACKED_HIV_TRANSITION_INDICES, tracked.trans)

    # Check keep dimensions
    allowed = names(get.transition.dimnames(jheem, transition.dimension=transition.dimension, hiv.negative = F))
    not.allowed = setdiff(keep.dimensions, allowed)
    if (length(not.allowed)>0)
        stop(paste0("keep.dimensions for an HIV-positive transition array for '",
                    transition.dimension, "' must be a subset of <",
                    paste0("'", allowed, "'", collapse=', '),
                    ">. The following are not allowed: ",
                    paste0("'", not.allowed, "'", collapse=', ')))

    # Now set up the collapsed array for indices and names
    collapsed.arr = apply(tracked.transition.arr, keep.dimensions, sum)
    n.collapsed = length(collapsed.arr)
    collapsed.arr[] = 1:n.collapsed

    collapse.indices = expand.population(collapsed.arr,
                                         target.dim.names = dimnames(tracked.transition.arr))

    jheem$parameters$TRACKED_HIV_TRANSITION_COLLAPSE_INDICES = collapse.indices +
        jheem$parameters$N_COLLAPSED_TRACKED_HIV_TRANSITIONS

    jheem$parameters$N_COLLAPSED_TRACKED_HIV_TRANSITIONS =
        max(jheem$parameters$TRACKED_HIV_TRANSITION_COLLAPSE_INDICES)

    collapsed.dim.names = lapply(names(dimnames(tracked.transition.arr)), function(one.dim.names.name){
        one.dim.names = dimnames(tracked.transition.arr)[[one.dim.names.name]]
        if (length(one.dim.names)==1 || any(keep.dimensions == one.dim.names.name))
            one.dim.names
        else if (length(one.dim.names)==length(dimnames(tracked.transition.mask)[[one.dim.names.name]]))
            'all'
        else
            paste0(one.dim.names, collapse = ';')
    })
    names(collapsed.dim.names) = names(dimnames(tracked.transition.arr))

    # Then save the names for those indices into the JHEEM$tracked.transitions
    jheem$tracked.transitions$names.for.collapsed.hiv = c(jheem$tracked.transitions$names.for.indices.hiv,
                                                             rep(name, n.collapsed))
    jheem$tracked.transitions$names.for.indices.hiv = c(jheem$tracked.transitions$names.for.indices.hiv,
                                                           rep(name, length(tracked.from)))

    # Last, save the transition dimension and the dimnames for re-hydrating this tracked transition
    jheem$tracked.transitions$transition.dimensions[name] = tolower(transition.dimension)
    jheem$tracked.transitions$dimnames.for.tracked[[name]] = collapsed.dim.names


    #Return it
    jheem
}

get.tracked.transition.order <- function(tracked.transition.mask, tracked.transition.arr)
{
    target.dim.names = dimnames(tracked.transition.arr)
    counter.dim.names = lapply(names(target.dim.names), function(dim.name.name){
        if (any(names(dimnames(tracked.transition.mask))==dim.name.name))
            dimnames(tracked.transition.mask)[[dim.name.name]]
        else
            target.dim.names[[dim.name.name]]
    })
    names(counter.dim.names) = names(target.dim.names)

    counter = array(0, dim=sapply(counter.dim.names, length), dimnames=counter.dim.names)
    access(counter,
           age=target.dim.names[['age']], race=target.dim.names[['race']],
           subpopulation=target.dim.names[['subpopulation']], sex=target.dim.names[['sex']], risk=target.dim.names[['risk']],
           non.hiv.subset=target.dim.names[['non.hiv.subset']], continuum=target.dim.names[['continuum']],
           cd4=target.dim.names[['cd4']], hiv.subset=target.dim.names[['hiv.subset']],
           subpopulation.from=target.dim.names[['subpopulation.from']], sex.from=target.dim.names[['sex.from']], risk.from=target.dim.names[['risk.from']],
           non.hiv.subset.from=target.dim.names[['non.hiv.subset.from']], continuum.from=target.dim.names[['continuum.from']],
           cd4.from=target.dim.names[['cd4.from']], hiv.subset.from=target.dim.names[['hiv.subset.from']],
           subpopulation.to=target.dim.names[['subpopulation.to']], sex.to=target.dim.names[['sex.to']], risk.to=target.dim.names[['risk.to']],
           non.hiv.subset.to=target.dim.names[['non.hiv.subset.to']], continuum.to=target.dim.names[['continuum.to']],
           cd4.to=target.dim.names[['cd4.to']], hiv.subset.to=target.dim.names[['hiv.subset.to']]
    ) = 1:prod(dim(tracked.transition.arr))

    order(get.tracked.transition.from.indices(counter, tracked.transition.mask = tracked.transition.mask))
}

#A helper to set up the list elements that will be needed for tracked transitions
init.transition.tracking <- function(jheem)
{
    #The indices in parameters
    jheem$parameters$TRACKED_NONHIV_TRANSITIONS_FROM = integer()
    jheem$parameters$TRACKED_NONHIV_TRANSITION_INDICES = integer()

    jheem$parameters$TRACKED_HIV_TRANSITIONS_FROM = integer()
    jheem$parameters$TRACKED_HIV_TRANSITION_INDICES = integer()

    jheem$parameters$TRACKED_HIV_TRANSITION_COLLAPSE_INDICES = integer()
    jheem$parameters$TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES = integer()

    jheem$parameters$N_COLLAPSED_TRACKED_HIV_TRANSITIONS = 0
    jheem$parameters$N_COLLAPSED_TRACKED_NONHIV_TRANSITIONS = 0

    #The names for those indices into the JHEEM$tracked.transitions
    jheem$tracked.transitions$names.for.indices.nonhiv = character()
    jheem$tracked.transitions$names.for.indices.hiv = character()

    #The dimnames for tracked transitions
    jheem$tracked.transitions$dimnames.for.tracked = list()

    #The transition dimension for each tracked transtions
    jheem$tracked.transitions$transition.dimensions = character()

    #Return
    jheem
}

#A helper function with elements common to HIV+ and HIV- tracked transitions
get.tracked.transition.from.indices <- function(population.indices, tracked.transition.mask)
{
    #pull names
    dim.name.names = names(dimnames(population.indices))
    transition.dim.name.names = names(dimnames(tracked.transition.mask))

    #overwrite names that should be either 'from' or 'to' as 'from'
    dim.name.names.from = sapply(dim.name.names, function(name){
        if (any(transition.dim.name.names == name))
            name
        else
            paste0(name, '.from')
    })

    #expand to from+to population
    # (in which we have labeled indices corresponding to the 'from' part)
    population.indices.from = expand.population(set.dim.name.names(population.indices, dim.name.names.from),
                                                dimnames(tracked.transition.mask))

    #pull tracked indices and return
    as.integer(population.indices.from[tracked.transition.mask])
}


