

##--------------------------------------------##
##-- THE PRIMARY EXPAND POPULATION FUNCTION --##
##--    (formulated as a private helper)    --##
##--------------------------------------------##

#Helper
# Shoehorns source population array into an array with dims corresponding to target.dim.names
#If names(target.dim.names) and names(dimnames(source.population)) are both non-null, will match
# on names of the dimnames
#Otherwise, will match on vector.equals for each dimname
#'@export
expand.population <- function(source.population,
                              target.dim.names,
                              non.conforming.error=NULL)
{
    if (list.equals(dimnames(source.population), target.dim.names))
        return(source.population)

    if (is.null(dim(source.population)) && !is.null(names(source.population)))
    {
        dim.names = list(names(source.population))
        dim(source.population) = length(source.population)
        dimnames(source.population) = dim.names
    }

    #map dims from source to target
    source.dim.names = dimnames(source.population)
    source.dim.ids = sapply(source.dim.names, function(dim.names){paste0(dim.names, collapse=', ')})
    target.dim.ids = sapply(target.dim.names, function(dim.names){paste0(dim.names, collapse=', ')})
    if (!is.null(names(target.dim.names)) && !is.null(names(source.dim.names)))
    {
        source.dim.ids = paste0(names(source.dim.names), ': ', source.dim.ids)
        target.dim.ids = paste0(names(target.dim.names), ': ', target.dim.ids)
    }

    source.to.target.dims = integer()
    for (source.id in source.dim.ids)
    {
        match.indices = all.indices.of(source.id, target.dim.ids)
        unused.match.indices = setdiff(match.indices, source.to.target.dims)
        if (length(match.indices)==0)
            source.to.target.dims = c(source.to.target.dims, NA)
        else if (length(unused.match.indices)==0)
            source.to.target.dims = c(source.to.target.dims, match.indices[1])
        else
            source.to.target.dims = c(source.to.target.dims, unused.match.indices[1])
    }

    #Throw an error if any dimensions in the source are not present in the target
    if (is.null(names(source.dim.names)))
        source.names.to.print = paste0("[", source.dim.ids, "]")
    else
        source.names.to.print = names(source.dim.names)

    if (any(is.na(source.to.target.dims)))
    {
        if (is.null(non.conforming.error))
        {
            missing = is.na(source.to.target.dims)
            n.missing = sum(missing)
            non.conforming.error = paste0(n.missing, ' dimension', ifelse(n.missing==1, '', 's') ,' in the array to expand ', ifelse(n.missing==1, 'is', 'are'), ' not present in the target population dimensions: ',
                                          paste0(paste0('(', 1:n.missing, ') ', source.names.to.print[missing]), collapse = '; '))
        }
        stop(non.conforming.error)
    }

    if (is.numeric(source.population))
        rv = do_expand_population(src=source.population,
                                  target_dims=sapply(target.dim.names, length),
                                  src_to_target_dim_map=source.to.target.dims-1)
    else if (is.logical(source.population))
    {
        rv = as.logical(do_expand_population(src=source.population,
                                    target_dims=sapply(target.dim.names, length),
                                    src_to_target_dim_map=source.to.target.dims-1))
    }
    else if (is.character(source.population))
        rv = do_expand_population_character(src=source.population,
                                            target_dims=sapply(target.dim.names, length),
                                            src_to_target_dim_map=source.to.target.dims-1)
    else
        stop("source.population must contain either numeric, logical, character values")

    dim(rv) = sapply(target.dim.names, length)
    dimnames(rv) = target.dim.names

    return (rv)
}


##---------------------------------------------------------##
##--             SPECIFIC 'EXPAND' FUNCTIONS             --##
##-- (pre-set for specific populations/parameter arrays) --##
##---------------------------------------------------------##

#'@title Shoehorn a given population array into an array of the dimensions shared by both the HIV-negative and HIV-positive model populations
#'
#'@description Returns an array with dimensions and dimnames shared by both the HIV-negative and HIV-positive populations, such that every element in the returned array is drawn from the population.to.expand array
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param population.to.expand An array representing the population to expand. dimnames(population.to.expand) must be a subset of the dimnames of the HIV-negative population and of the dimnames of the HIV-positive population
#'
#'@seealso \code{\link{expand.population.to.hiv.positive}}, \code{\link{expand.population.to.hiv.negative}}
#'
#'@export
expand.population.to.general <- function(jheem, population.to.expand)
{
    expand.population(population.to.expand, get.dimnames.general(jheem))
}

#'@title Shoehorn a given population array into an array of the dimensions needed for the HIV-negative model population
#'
#'@description Returns an array with dimensions and dimnames matching the HIV negative population, such that every element in the returned array is drawn from the population.to.expand array
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param population.to.expand An array representing the population to expand. dimnames(population.to.expand) must be a subset of the dimnames of the HIV-negative population
#'
#'@seealso \code{\link{expand.population.to.hiv.positive}}
#'
#'@export
expand.population.to.hiv.negative <- function(jheem, population.to.expand)
{
    expand.population(population.to.expand, get.dimnames.nonhiv(jheem))
}

#'@title Shoehorn a given population array into an array of the dimensions needed for the HIV-positive model population
#'
#'@description Returns an array with dimensions and dimnames matching the HIV positive population, such that every element in the returned array is drawn from the population.to.expand array
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param population.to.expand An array representing the population to expand. dimnames(population.to.expand) must be a subset of the dimnames of the HIV-positive population
#'
#'@seealso \code{\link{expand.population.to.hiv.negative}}
#'
#'@export
expand.population.to.hiv.positive <- function(jheem, population.to.expand)
{
    expand.population(population.to.expand, get.dimnames.hiv(jheem))
}

#'@title Shoehorn a given population array into an array of the dimensions for new infection proportions
#'
#'@description Returns an array with dimensions and dimnames for the new infection proportions array, such that every element in the returned array is drawn from the population.to.expand array
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param to.expand An array representing the population to expand. dimnames(population.to.expand) must be a subset of the dimnames of the HIV-negative population and/or the dimnames of the HIV-positive population: [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]
#'
#'@inherit set.new.infection.proportions details
#'
#'@export
expand.to.new.infection.proportions <- function(jheem, to.expand)
{
    expand.population(to.expand, get.dimnames.all(jheem))
}

#'@title Shoehorn a given array into a transmission contact array
#'
#'@description Returns an array with dimensions and dimnames for a transmission contact array
#'
#'@param jheem An object of \code{\link{jheem-class}}
#'@param to.expand An array representing the contact array to expand. dimnames(population.to.expand) must be a subset of the dimnames for a transmission contact array: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]
#'
#'@export
expand.to.contact.array <- function(jheem, to.expand)
{
    if (is.null(dimnames(to.expand)))
        dim.name.names = character()
    else
    {
        if (is.null(names(dimnames(to.expand))))
            to.expand = set.dim.name.names(to.expand, identify.dimnames.in.array(jheem, to.expand))

        dim.name.names = names(dimnames(to.expand))
    }

    target.dim.names = get.dimnames(jheem,
                                    age.from = any(dim.name.names=='age.from'),
                                    race.from = any(dim.name.names=='race.from'),
                                    subpopulation.from = any(dim.name.names=='subpopulation.from'),
                                    sex.from = any(dim.name.names=='sex.from'),
                                    risk.from = any(dim.name.names=='risk.from'),
                                    age.to=T,
                                    race.to=T,
                                    subpopulation.to=T,
                                    sex.to=T,
                                    risk.to=T,
                                    non.hiv.subset.to=T)

    expand.population(to.expand, target.dim.names)
}

