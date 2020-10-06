

##----------------------------##
##-- FRONT-FACING FUNCTIONS --##
##----------------------------##

#'@title Extract a subset of the population from JHEEM Results
#'
#'@param results The results of a call to \code{\link{run.jheem}}
#'@param years,ages,races,subpopulations,risks,non.hiv.subsets,continuum,cd4s The elements of each of the possible result dimensions for which to extract results. If passed null, defaults to all elements of each dimension
#'@param keep.dimensions The names of which dimensions to marginalize over (all other dimensions will be summed out) - a subset of 'year', 'age', 'race', 'subpopulation', 'sex', 'risk', 'non.hiv.subset', 'continuum', 'cd4', and 'hiv.subset'. If passed null, defaults to keep the year dimension, and any other dimension for which more than one but less than all elements of that dimension are selected
#'@param denominator.dimensions If the population should be normalized, which dimensions to use in the denominator (only applies if per.population is not NA). This must be a subset of the keep.dimensions. If NULL, the results are reported as fractions of the total population
#'@param include.hiv.positive,include.hiv.negative Whether to include the hiv.positive and hiv.negative subpopulations (at least one of these two must be set to true)
#'@param per.population The unit of the denominator if the return value should be a proportion of the population. If NA, the absolute number is returned without dividing by the population size
#'
#'@export
extract.population.subset <- function(results,
                                      years=NULL,
                                      ages=NULL,
                                      races=NULL,
                                      subpopulations=NULL,
                                      sexes=NULL,
                                      risks=NULL,
                                      non.hiv.subsets=NULL,
                                      continuum=NULL,
                                      cd4s=NULL,
                                      hiv.subsets=NULL,
                                      include.hiv.positive=T,
                                      include.hiv.negative=T,
                                      keep.dimensions=NULL,
                                      denominator.dimensions='year',
                                      per.population=NA,
                                      transformation.fn=NULL,
                                      use.cdc.categorizations=F
)
{
    #Check arguments
    if (!include.hiv.positive && !include.hiv.negative)
        stop("Extracted data must include at least one of HIV-positive or HIV-negative")

    ALLOWED.DIMENSIONS = c('year','age','race','subpopulation','sex','risk','non.hiv.subset','continuum','cd4','hiv.subset')
    invalid.keep.dimensions = setdiff(keep.dimensions, ALLOWED.DIMENSIONS)
    if (!is.null(keep.dimensions) && length(invalid.keep.dimensions)>0)
        stop(paste0("keep.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.keep.dimensions, "'"), collapse=',')))
    invalid.denominator.dimensions = setdiff(keep.dimensions, ALLOWED.DIMENSIONS)
    if (!is.na(per.population) && length(invalid.denominator.dimensions)>0)
        stop(paste0("denominator.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.denominator.dimensions, "'"), collapse=',')))

    in.denominator.but.not.keep = setdiff(denominator.dimensions, keep.dimensions)
    if (!is.na(per.population) && !is.null(keep.dimensions) && length(in.denominator.but.not.keep)>0)
        stop("denominator.dimensions must be a subset of keep.dimensions")

    if (include.hiv.positive && (any(keep.dimensions=='non.hiv.subset') ||
                                 any(denominator.dimensions=='non.hiv.subset')))
        stop("If including hiv.positive, 'non.hiv.subset' cannot be one of the keep.dimensions or denominator.dimensions")
    if (include.hiv.negative && (any(keep.dimensions=='continuum') ||
                                 any(keep.dimensions=='cd4') ||
                                 any(keep.dimensions=='hiv.subset') ||
                                 any(denominator.dimensions=='continuum') ||
                                 any(denominator.dimensions=='cd4') ||
                                 any(denominator.dimensions=='hiv.subset')))
        stop("If including hiv.negative, neither 'continuum', 'cd4', nor 'hiv.subset' can be one of the keep.dimensions or denominator.dimensions")

    #Fill in missing dimensions
    if (is.null(years))
        years = results$years
    if (is.null(ages))
        ages = results$ages
    if (is.null(races))
        races = results$races
    if (is.null(subpopulations))
        subpopulations = results$subpopulations
    if (is.null(sexes))
    {
        if (use.cdc.categorizations)
            sexes = results$cdc.sexes
        else
            sexes = results$sexes
    }
    if (is.null(risks))
    {
        if (use.cdc.categorizations)
            risks = results$cdc.risks
        else
            risks = results$risks
    }
    if (is.null(non.hiv.subsets))
        non.hiv.subsets = results$non.hiv.subsets
    if (is.null(continuum))
        continuum = results$continuum
    if (is.null(cd4s))
        cd4s = results$cd4
    if (is.null(hiv.subsets))
        hiv.subsets = results$hiv.subsets

    years = as.character(years)

    #Set up keep dimensions
    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=non.hiv.subsets,
                                                      continuum=continuum,
                                                      cd4s=cd4s,
                                                      hiv.subsets=hiv.subsets,
                                                      use.cdc.categorizations=use.cdc.categorizations,
                                                      at.minimum=denominator.dimensions)

    #set CDC or non-CDC
    if (use.cdc.categorizations)
    {
        hiv.negative = results$hiv.negative.cdc
        hiv.positive = results$hiv.positive.cdc
    }
    else
    {
        hiv.negative = results$hiv.negative
        hiv.positive = results$hiv.positive
    }

    #Pull the numerators
    numerators.hiv.negative = numerators.hiv.positive = 0
    if (include.hiv.negative)
    {
        if (length(keep.dimensions)==0)
            numerators.hiv.negative = sum(hiv.negative[years,ages,races,subpopulations,sexes,risks,non.hiv.subsets])
        else
        {
            dim.names = list(year=years, age=ages, race=races, subpopulation=subpopulations, sex=sexes, risk=risks, non.hiv.subset=non.hiv.subsets)
            numerators.hiv.negative = hiv.negative[years,ages,races,subpopulations,sexes,risks,non.hiv.subsets]
            dim(numerators.hiv.negative) = sapply(dim.names, length)
            dimnames(numerators.hiv.negative) = dim.names
            numerators.hiv.negative = apply(numerators.hiv.negative, keep.dimensions, sum)
        }
    }
    if (include.hiv.positive)
    {
        if (length(keep.dimensions)==0)
            numerators.hiv.positive = sum(hiv.positive[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets])
        else
        {
            dim.names = list(year=years, age=ages, race=races, subpopulation=subpopulations, sex=sexes, risk=risks, continuum=continuum, cd4=cd4s, hiv.subset=hiv.subsets)
            numerators.hiv.positive = hiv.positive[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets]
            dim(numerators.hiv.positive) = sapply(dim.names, length)
            dimnames(numerators.hiv.positive) = dim.names
            numerators.hiv.positive = apply(numerators.hiv.positive, keep.dimensions, sum)
        }
    }

    rv = numerators.hiv.negative + numerators.hiv.positive


    #If one-dimensional, hydrate up to array
    if (length(keep.dimensions)>1 && is.null(dim(rv)))
    {
        dim.names = list(names(rv))
        names(dim.names) = keep.dimensions

        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }


    #If requested, divide by denominator and normalize
    if (!is.na(per.population))
    {
        if (is.null(denominator.dimensions))
            denominator = sum(rv)
        else
        {
            if (length(denominator.dimensions)==0)
                denominator = sum(rv)
            else
                denominator = apply(rv, denominator.dimensions, sum)

            if (is.null(dim(rv)))
                denominator = rep(denominator, length(rv))
            else
            {
                if (is.null(dim(denominator)))
                {
                    dim.names = list(names(denominator))
                    names(dim.names) = denominator.dimensions

                    dim(denominator) = sapply(dim.names, length)
                    dimnames(denominator) = dim.names
                }
                denominator = expand.population(denominator,
                                                target.dim.names = dimnames(rv))
            }
        }

        rv = rv / denominator * per.population
    }

    rv
}

#'@title Extract the subset of the HIV-positive population which has been diagnosed from JHEEM Results
#'
#'@param results The results of a call to \code{\link{run.jheem}}
#'@param years,ages,races,subpopulations,risks,cd4s The elements of each of the possible result dimensions for which to extract results. If passed null, defaults to all elements of each dimension
#'@param keep.dimensions The names of which dimensions to marginalize over (all other dimensions will be summed out) - a subset of 'year', 'age', 'race', 'subpopulation', 'sex', 'risk', 'non.hiv.subset', 'continuum', 'cd4', and 'hiv.subset'. If passed null, defaults to keep the year dimension, and any other dimension for which more than one but less than all elements of that dimension are selected
#'@param per.population The unit of the denominator if the return value should be a proportion of all of those with HIV If NA, the absolute number is returned without dividing by the population size
#'
#'@export
extract.diagnosed.hiv <- function(results,
                                      years=NULL,
                                      ages=NULL,
                                      races=NULL,
                                      subpopulations=NULL,
                                      sexes=NULL,
                                      risks=NULL,
                                      cd4s=NULL,
                                      hiv.subsets=NULL,
                                      keep.dimensions='year',
                                      per.population=1,
                                      use.cdc.categorizations=F
)
{
    if (is.null(results$diagnosed.continuum.states))
        stop("Which continuum states represent diagnosed states has not been indicated for this JHEEM. Unable to report diagnosed HIV")
    keep.dimensions = setdiff(keep.dimensions, 'continuum')
    rv = extract.population.subset(results=results,
                              continuum=results$diagnosed.continuum.states,
                              include.hiv.negative=F,
                              keep.dimensions=keep.dimensions,
                              years=years,
                              ages=ages,
                              races=races,
                              subpopulations=subpopulations,
                              sexes=sexes,
                              risks=risks,
                              cd4s=cd4s,
                              hiv.subsets=hiv.subsets,
                              per.population=NA,
                              use.cdc.categorizations=use.cdc.categorizations)

    if (!is.na(per.population))
    {
        denominators = extract.population.subset(results=results,
                                                 continuum=NULL,
                                                 include.hiv.negative=F,
                                                 keep.dimensions=keep.dimensions,
                                                 years=years,
                                                 ages=ages,
                                                 races=races,
                                                 subpopulations=subpopulations,
                                                 sexes=sexes,
                                                 risks=risks,
                                                 cd4s=cd4s,
                                                 hiv.subsets=hiv.subsets,
                                                 per.population=NA,
                                                 use.cdc.categorizations=use.cdc.categorizations)

        rv = rv / denominators * per.population
    }

    rv
}
#'@title Extract HIV incidence from JHEEM Results
#'
#'@inheritParams extract.population.subset
#'@param include.hiv.positive.in.denominator Whether to include HIV positive individuals in the denominator (only applies if per.population is set)
#'
#'@export
extract.incidence <- function(results,
                              years=NULL,
                              ages=NULL,
                              races=NULL,
                              subpopulations=NULL,
                              sexes=NULL,
                              risks=NULL,
                              non.hiv.subsets=NULL,
                              continuum=NULL,
                              cd4s=NULL,
                              hiv.subsets=NULL,
                              keep.dimensions=NULL,
                              include.hiv.positive.in.denominator=T,
                              per.population=100000,
                              use.cdc.categorizations=F
)
{
    #Check arguments
    ALLOWED.KEEP.DIMENSIONS = c('year','age','race','subpopulation','sex','risk')#c('year','age','race','subpopulation','sex','risk','continuum.from','cd4','hiv.subset','continuum.to')
    if (!include.hiv.positive.in.denominator)
        ALLOWED.KEEP.DIMENSIONS = c(ALLOWED.KEEP.DIMENSIONS, 'non.hiv.subset')
    invalid.keep.dimensions = setdiff(keep.dimensions, ALLOWED.KEEP.DIMENSIONS)
    if (!is.null(keep.dimensions) && length(invalid.keep.dimensions)>0)
        stop(paste0("When ", ifelse(include.hiv.positive.in.denominator, '', 'NOT'),
                    "including HIV-positive in the denominator, keep.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.KEEP.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.keep.dimensions, "'"), collapse=',')))

    #Fill in missing dimensions
    if (is.null(years))
        years = dimnames(results$incidence)$year
    if (is.null(ages))
        ages = dimnames(results$incidence)$age
    if (is.null(races))
        races = dimnames(results$incidence)$race
    if (is.null(subpopulations))
        subpopulations = dimnames(results$incidence)$subpopulation
    if (is.null(sexes))
    {
        if (use.cdc.categorizations)
            sexes = dimnames(results$incidence.cdc)$sex
        else
            sexes = dimnames(results$incidence)$sex
    }
    if (is.null(risks))
    {
        if (use.cdc.categorizations)
            risks = dimnames(results$incidence.cdc)$risk
        else
            risks = dimnames(results$incidence)$risk
    }
    if (is.null(non.hiv.subsets))
        non.hiv.subsets = dimnames(results$incidence)$non.hiv.subset
    if (is.null(continuum))
        continuum = dimnames(results$incidence)$continuum
    if (is.null(cd4s))
        cd4s = dimnames(results$incidence)$cd4
    if (is.null(hiv.subsets))
        hiv.subsets = dimnames(results$incidence)$hiv.subset

    years = as.character(years)

    #Set up keep dimensions
    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=if (!include.hiv.positive.in.denominator) NULL else hiv.subsets,
                                                      continuum=NULL,
                                                      cd4s=NULL,
                                                      hiv.subsets=NULL,
                                                      use.cdc.categorizations=use.cdc.categorizations)
    #set CDC or non-CDC
    if (use.cdc.categorizations)
        incidence = results$incidence.cdc
    else
        incidence = results$incidence

    #Pull the numerators
    if (length(keep.dimensions)==0)
        rv = sum(incidence[years,ages,races,subpopulations,sexes,risks,non.hiv.subsets,continuum,cd4s,hiv.subsets])
    else
    {
        dim.names = list(year=years, age=ages, race=races, subpopulation=subpopulations, sex=sexes, risk=risks,
                         non.hiv.subset=non.hiv.subsets, continuum=continuum, cd4=cd4s, hiv.subset=hiv.subsets)
        rv = incidence[years,ages,races,subpopulations,sexes,risks,non.hiv.subsets,continuum,cd4s,hiv.subsets]
        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
        rv = apply(rv, keep.dimensions, sum)
    }

    #If one-dimensional, hydrate up to array
    if (length(keep.dimensions)>1 && is.null(dim(rv)))
    {
        dim.names = list(names(rv))
        names(dim.names) = keep.dimensions

        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }


    #If requested, divide by denominator and normalize
    if (!is.na(per.population))
    {
        denominator = extract.population.subset(results,
                                                years=years,
                                                ages=ages,
                                                races=races,
                                                subpopulations=subpopulations,
                                                sexes=sexes,
                                                risks=risks,
                                                non.hiv.subsets=NULL,
                                                continuum=NULL,
                                                cd4s=NULL,
                                                hiv.subsets=NULL,
                                                keep.dimensions=keep.dimensions,
                                                per.population=NA,
                                                include.hiv.negative = T,
                                                include.hiv.positive = include.hiv.positive.in.denominator,
                                                use.cdc.categorizations = use.cdc.categorizations)

        rv = rv / denominator * per.population
    }

    rv
}


#'@title Extract HIV prevalence from JHEEM Results
#'
#'@inheritParams extract.incidence
#'
#'@export
extract.prevalence <- function(results,
                               years=NULL,
                               ages=NULL,
                               races=NULL,
                               subpopulations=NULL,
                               sexes=NULL,
                               risks=NULL,
                               continuum=NULL,
                               cd4s=NULL,
                               hiv.subsets=NULL,
                               keep.dimensions=NULL,
                               per.population=100,
                               transformation.fn=NULL,
                               use.cdc.categorizations=F
)
{
    #Check arguments
    ALLOWED.KEEP.DIMENSIONS = c('year','age','race','subpopulation','sex','risk')#c('year','age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset')
    invalid.keep.dimensions = setdiff(keep.dimensions, ALLOWED.KEEP.DIMENSIONS)
    if (!is.null(keep.dimensions) && length(invalid.keep.dimensions)>0)
        stop(paste0("keep.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.KEEP.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.keep.dimensions, "'"), collapse=',')))

#    ALLOWED.DENOMINATOR.DIMENSIONS = c('year','age','race','subpopulation','sex','risk')
#    invalid.denominator.dimensions = setdiff(keep.dimensions, ALLOWED.DENOMINATOR.DIMENSIONS)
#    if (!is.na(per.population) && length(invalid.denominator.dimensions)>0)
#        stop(paste0("denominator.dimensions must be a subset of ",
#                    paste0(paste0("'", ALLOWED.DENOMINATOR.DIMENSIONS, "'"), collapse=','),
#                    ". The following dimensions were given but not allowed: ",
#                    paste0(paste0("'", invalid.denominator.dimensions, "'"), collapse=',')))

#    in.denominator.but.not.keep = setdiff(denominator.dimensions, keep.dimensions)
#    if (!is.null(keep.dimensions) && length(in.denominator.but.not.keep)>0)
#        stop("denominator.dimensions must be a subset of keep.dimensions")

    #Fill in missing dimensions
    if (is.null(years))
        years = results$years
    if (is.null(ages))
        ages = results$ages
    if (is.null(races))
        races = results$races
    if (is.null(subpopulations))
        subpopulations = results$subpopulations
    if (is.null(sexes))
    {
        if (use.cdc.categorizations)
            sexes = results$cdc.sexes
        else
            sexes = results$sexes
    }
    if (is.null(risks))
    {
        if (use.cdc.categorizations)
            risks = results$cdc.risks
        else
            risks = results$risks
    }
    if (is.null(continuum))
        continuum = results$continuum
    if (is.null(cd4s))
        cd4s = results$cd4
    if (is.null(hiv.subsets))
        hiv.subsets = results$hiv.subsets

    years = as.character(years)

    #Set up keep dimensions
    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=NULL,
                                                      continuum=NULL,
                                                      cd4s=NULL,
                                                      hiv.subsets=NULL,
                                                      use.cdc.categorizations=use.cdc.categorizations)


    #set CDC or non-CDC
    if (use.cdc.categorizations)
        hiv.positive = results$hiv.positive.cdc
    else
        hiv.positive = results$hiv.positive

    #Pull the numerators
    if (length(keep.dimensions)==0)
        rv = sum(hiv.positive[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets])
    else
    {
        dim.names = list(year=years, age=ages, race=races, subpopulation=subpopulations, sex=sexes, risk=risks,
                         continuum=continuum, cd4=cd4s, hiv.subset=hiv.subsets)
        rv = hiv.positive[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets]
        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
        rv = apply(rv, keep.dimensions, sum)
    }

    #If one-dimensional, hydrate up to array
    if (length(keep.dimensions)>1 && is.null(dim(rv)))
    {
        dim.names = list(names(rv))
        names(dim.names) = keep.dimensions

        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }


    #If requested, divide by denominator and normalize
    if (!is.na(per.population))
    {
        denominator = extract.population.subset(results,
                                                years=years,
                                                ages=ages,
                                                races=races,
                                                subpopulations=subpopulations,
                                                sexes=sexes,
                                                risks=risks,
                                                non.hiv.subsets=NULL,
                                                continuum=NULL,
                                                cd4s=NULL,
                                                hiv.subsets=NULL,
                                                keep.dimensions=keep.dimensions,
                                                per.population=NA,
                                                include.hiv.negative = T,
                                                include.hiv.positive = T,
                                                use.cdc.categorizations = use.cdc.categorizations)

        rv = rv / denominator * per.population
    }

    rv
}

#'@title Extract new HIV diagnoses from JHEEM Results
#'
#'@inheritParams extract.incidence
#'
#'@export
extract.new.diagnoses <- function(results,
                                  years=NULL,
                                  ages=NULL,
                                  races=NULL,
                                  subpopulations=NULL,
                                  sexes=NULL,
                                  risks=NULL,
                                  continuum.from=NULL,
                                  cd4s=NULL,
                                  hiv.subsets=NULL,
                                  keep.dimensions=NULL,
                                  continuum.to=NULL,
                                  include.hiv.positive.in.denominator=T,
                                  per.population=100000,
                                  transformation.fn=NULL,
                                  use.cdc.categorizations=F
)
{
    #Check arguments
    ALLOWED.KEEP.DIMENSIONS = c('year','age','race','subpopulation','sex','risk')#c('year','age','race','subpopulation','sex','risk','continuum.from','cd4','hiv.subset','continuum.to')
    invalid.keep.dimensions = setdiff(keep.dimensions, ALLOWED.KEEP.DIMENSIONS)
    if (!is.null(keep.dimensions) && length(invalid.keep.dimensions)>0)
        stop(paste0("keep.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.KEEP.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.keep.dimensions, "'"), collapse=',')))

#    ALLOWED.DENOMINATOR.DIMENSIONS = c('year','age','race','subpopulation','sex','risk')
#    invalid.denominator.dimensions = setdiff(keep.dimensions, ALLOWED.DENOMINATOR.DIMENSIONS)
#    if (!is.na(per.population) && length(invalid.denominator.dimensions)>0)
#        stop(paste0("denominator.dimensions must be a subset of ",
#                    paste0(paste0("'", ALLOWED.DENOMINATOR.DIMENSIONS, "'"), collapse=','),
#                    ". The following dimensions were given but not allowed: ",
#                    paste0(paste0("'", invalid.denominator.dimensions, "'"), collapse=',')))

#    in.denominator.but.not.keep = setdiff(denominator.dimensions, keep.dimensions)
#    if (!is.null(keep.dimensions) && length(in.denominator.but.not.keep)>0)
#        stop("denominator.dimensions must be a subset of keep.dimensions")

    #Fill in missing dimensions
    if (is.null(years))
        years = dimnames(results$new.diagnoses)$year
    if (is.null(ages))
        ages = dimnames(results$new.diagnoses)$age
    if (is.null(races))
        races = dimnames(results$new.diagnoses)$race
    if (is.null(subpopulations))
        subpopulations = dimnames(results$new.diagnoses)$subpopulation
    if (is.null(sexes))
    {
        if (use.cdc.categorizations)
            sexes = dimnames(results$new.diagnoses.cdc)$sex
        else
            sexes = dimnames(results$new.diagnoses)$sex
    }
    if (is.null(risks))
    {
        if (use.cdc.categorizations)
            risks = dimnames(results$new.diagnoses.cdc)$risk
        else
            risks = dimnames(results$new.diagnoses)$risk
    }
    if (is.null(continuum.from))
        continuum.from = dimnames(results$new.diagnoses)$continuum.from
    if (is.null(cd4s))
        cd4s = dimnames(results$new.diagnoses)$cd4
    if (is.null(hiv.subsets))
        hiv.subsets = dimnames(results$new.diagnoses)$hiv.subset
    if (is.null(continuum.to))
        continuum.to = dimnames(results$new.diagnoses)$continuum.to

    years = as.character(years)

    #Set up keep dimensions
    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=NULL,
                                                      continuum=NULL,
                                                      cd4s=NULL,
                                                      hiv.subsets=NULL,
                                                      use.cdc.categorizations=use.cdc.categorizations)
    #set CDC or non-CDC
    if (use.cdc.categorizations)
        new.diagnoses = results$new.diagnoses.cdc
    else
        new.diagnoses = results$new.diagnoses

    #Pull the numerators
    if (length(keep.dimensions)==0)
        rv = sum(new.diagnoses[years,ages,races,subpopulations,sexes,risks,continuum.from,cd4s,hiv.subsets,continuum.to])
    else
    {
        dim.names = list(year=years, age=ages, race=races, subpopulation=subpopulations, sex=sexes, risk=risks,
                         continuum.from=continuum.from, cd4=cd4s, hiv.subsets=hiv.subsets, continuum.to=continuum.to)
        rv = new.diagnoses[years,ages,races,subpopulations,sexes,risks,continuum.from,cd4s,hiv.subsets,continuum.to]
        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
        rv = apply(rv, keep.dimensions, sum)
    }


    #If one-dimensional, hydrate up to array
    if (length(keep.dimensions)>1 && is.null(dim(rv)))
    {
        dim.names = list(names(rv))
        names(dim.names) = keep.dimensions

        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }


    #If requested, divide by denominator and normalize
    if (!is.na(per.population))
    {
        denominator = extract.population.subset(results,
                                                years=years,
                                                ages=ages,
                                                races=races,
                                                subpopulations=subpopulations,
                                                sexes=sexes,
                                                risks=risks,
                                                non.hiv.subsets=NULL,
                                                continuum=NULL,
                                                cd4s=NULL,
                                                hiv.subsets=NULL,
                                                keep.dimensions=keep.dimensions,
                                                per.population=NA,
                                                include.hiv.negative = T,
                                                include.hiv.positive = include.hiv.positive.in.denominator,
                                                use.cdc.categorizations = use.cdc.categorizations)

        rv = rv / denominator * per.population
    }

    rv
}

#'@title Extract a tracked transition from JHEEM Results
#'
#'@inheritParams extract.incidence
#'
#'@export
extract.tracked.transition <- function(results,
                                       transition.name,
                                       years=NULL,
                                       ages=NULL,
                                       races=NULL,
                                       subpopulations=NULL,
                                       sexes=NULL,
                                       risks=NULL,
                                       non.hiv.subsets=NULL,
                                       continuum=NULL,
                                       cd4s=NULL,
                                       hiv.subsets=NULL,
                                       keep.dimensions=NULL,
                                       include.hiv.negative.in.denominator=T,
                                       include.hiv.positive.in.denominator=T,
                                       per.population=100000,
                                       transformation.fn=NULL
)
{
    stop("This function is not implented correctly now - needs a rewrite")

    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=non.hiv.subsets,
                                                      continuum=continuum,
                                                      cd4s=cd4s,
                                                      hiv.subsets=hiv.subsets,
                                                      use.cdc.categorizations=use.cdc.categorizations)


    if (!is.na(per.population) && !include.hiv.positive.in.denominator && !include.hiv.negative.in.denominator)
        stop("Extracted data must include at least one of HIV-positive or HIV-negative in the denominator")

    denominators = list(results$hiv.negative, results$hiv.positive, NULL)[c(include.hiv.negative.in.denominator, include.hiv.positive.in.denominator, T)]


    do.extract.results(numerators=results$tracked.transitions[[transition.name]],
                       denominators.1=denominators[[1]],
                       denominators.2=denominators[[2]],
                       years=years,
                       ages=ages,
                       races=races,
                       subpopulations=subpopulations,
                       sexes=sexes,
                       risks=risks,
                       non.hiv.subsets=non.hiv.subsets,
                       continuum=continuum,
                       cd4s=cd4s,
                       hiv.subsets=hiv.subsets,
                       keep.dimensions=keep.dimensions,
                       per.population=per.population,
                       transformation.fn=transformation.fn)
}

#'@title Extract HIV-specific mortality from JHEEM Results
#'
#'@inheritParams extract.incidence
#'@param include.hiv.negative.in.denominator Whether to include HIV negative individuals in the denominator (only applies if per.population is set)
#'
#'@export
extract.hiv.specific.mortality <- function(results,
                                           years=NULL,
                                           ages=NULL,
                                           races=NULL,
                                           subpopulations=NULL,
                                           sexes=NULL,
                                           risks=NULL,
                                           continuum=NULL,
                                           cd4s=NULL,
                                           hiv.subsets=NULL,
                                           keep.dimensions=NULL,
                                           include.hiv.negative.in.denominator=T,
                                           per.population=100000,
                                           transformation.fn=NULL,
                                           use.cdc.categorizations=F
)
{
stop("This needs to be cleaned up and merged with overall hiv mortality")
    #Check arguments
    ALLOWED.DIMENSIONS = c('year','age','race','subpopulation','sex','risk','non.hiv.subset','continuum','cd4','hiv.subset')
    invalid.keep.dimensions = setdiff(keep.dimensions, ALLOWED.DIMENSIONS)
    if (!is.null(keep.dimensions) && length(invalid.keep.dimensions)>0)
        stop(paste0("keep.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.keep.dimensions, "'"), collapse=',')))
    invalid.denominator.dimensions = setdiff(keep.dimensions, ALLOWED.DIMENSIONS)
    if (!is.na(per.population) && length(invalid.denominator.dimensions)>0)
        stop(paste0("denominator.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.denominator.dimensions, "'"), collapse=',')))

    in.denominator.but.not.keep = setdiff(denominator.dimensions, keep.dimensions)
    if (!is.null(keep.dimensions) && length(in.denominator.but.not.keep)>0)
        stop("denominator.dimensions must be a subset of keep.dimensions")

    if (include.hiv.negative.in.denominator && (any(denominator.dimensions=='continuum') ||
                                                any(denominator.dimensions=='cd4') ||
                                                any(denominator.dimensions=='hiv.subset')))
        stop("If including hiv.negative, neither 'continuum', 'cd4', nor 'hiv.subset' can be one of the denominator.dimensions")


    #Fill in missing dimensions
    if (is.null(years))
        years = results$years
    if (is.null(ages))
        ages = results$ages
    if (is.null(races))
        races = results$races
    if (is.null(subpopulations))
        subpopulations = results$subpopulations
    if (is.null(sexes))
    {
        if (use.cdc.categorizations)
            sexes = results$cdc.sexes
        else
            sexes = results$sexes
    }
    if (is.null(risks))
    {
        if (use.cdc.categorizations)
            risks = results$cdc.risks
        else
            risks = results$risks
    }
    if (is.null(continuum))
        continuum = results$continuum
    if (is.null(cd4s))
        cd4s = results$cd4
    if (is.null(hiv.subsets))
        hiv.subsets = results$hiv.subsets

    years = as.character(years)

    #Set up keep dimensions
    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=NULL,
                                                      continuum=continuum,
                                                      cd4s=cd4s,
                                                      hiv.subsets=hiv.subsets,
                                                      at.minimum=denominator.dimensions,
                                                      use.cdc.categorizations=use.cdc.categorizations)
    #set CDC or non-CDC
    if (use.cdc.categorizations)
        hiv.specific.mortality = results$hiv.specific.mortality.cdc
    else
        hiv.specific.mortality = results$hiv.specific.mortality

    #Pull the numerators
    rv = apply(hiv.specific.mortality[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets], keep.dimensions, sum)


    #If one-dimensional, hydrate up to array
    if (is.null(dim(rv)))
    {
        dim.names = list(names(rv))
        names(dim.names) = keep.dimensions

        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }


    #If requested, divide by denominator and normalize
    if (!is.na(per.population))
    {
        denominator = extract.population.subset(results,
                                                years=years,
                                                ages=ages,
                                                races=races,
                                                subpopulations=subpopulations,
                                                sexes=sexes,
                                                risks=risks,
                                                non.hiv.subsets=NULL,
                                                continuum=NULL,
                                                cd4s=NULL,
                                                hiv.subsets=NULL,
                                                keep.dimensions=c('year','age','race','subpopulation','sex','risk'),
                                                include.hiv.negative = include.hiv.negative.in.denominator,
                                                include.hiv.positive = T,
                                                use.cdc.categorizations = use.cdc.categorizations)

        if (is.null(denominator.dimensions))
            denominator = sum(rv)
        else
        {
            denominator = apply(rv, denominator.dimensions, sum)
            if (is.null(dim(denominator)))
            {
                dim.names = list(names(denominator))
                names(dim.names) = denominator.dimensions

                dim(denominator) = sapply(dim.names, length)
                dimnames(denominator) = dim.names
            }
            denominator = expand.population(denominator,
                                            target.dim.names = dimnames(rv))
        }

        rv = rv / denominator * per.population
    }

    rv
}

#'@title Extract Overall Mortality in HIV-positive from JHEEM Results
#'
#'@inheritParams extract.incidence
#'@param include.hiv.negative.in.denominator Whether to include HIV negative individuals in the denominator (only applies if per.population is set)
#'
#'@details Overall mortality for HIV-positive is HIV-specific mortality plus non-HIV mortality
#'
#'@export
extract.overall.hiv.mortality <- function(results,
                                  years=NULL,
                                  ages=NULL,
                                  races=NULL,
                                  subpopulations=NULL,
                                  sexes=NULL,
                                  risks=NULL,
                                  continuum=NULL,
                                  cd4s=NULL,
                                  hiv.subsets=NULL,
                                  keep.dimensions=NULL,
                                  include.hiv.negative.in.denominator=T,
                                  per.population=100000,
                                  transformation.fn=NULL,
                                  use.cdc.categorizations=F
)
{
    #Check arguments
    if (include.hiv.negative.in.denominator)
        ALLOWED.KEEP.DIMENSIONS = c('year','age','race','subpopulation','sex','risk')#c('year','age','race','subpopulation','sex','risk','continuum.from','cd4','hiv.subset','continuum.to')
    else
        ALLOWED.KEEP.DIMENSIONS = c('year','age','race','subpopulation','sex','risk','continuum','cd4','hiv.subset')
    invalid.keep.dimensions = setdiff(keep.dimensions, ALLOWED.KEEP.DIMENSIONS)
    if (!is.null(keep.dimensions) && length(invalid.keep.dimensions)>0)
        stop(paste0(ifelse(include.hiv.negative.in.denominator,
                           "When including HIV-negative in the denominator, ",
                           "When excluding HIV-negative from the denominator, "),
                    "keep.dimensions must be a subset of ",
                    paste0(paste0("'", ALLOWED.KEEP.DIMENSIONS, "'"), collapse=','),
                    ". The following dimensions were given but not allowed: ",
                    paste0(paste0("'", invalid.keep.dimensions, "'"), collapse=',')))


    if (include.hiv.negative.in.denominator && !is.na(per.population))
    {
        if (!is.null(continuum) || !is.null(cd4s) || !is.null(hiv.subsets))
            stop("If including HIV-positive in the denominator, continuum, cd4s, and hiv.subsets must be set to NULL")
    }

    may.be.null.continuum = continuum
    may.be.null.cd4s = cd4s
    may.be.null.hiv.subsets = hiv.subsets
    #Fill in missing dimensions
    if (is.null(years))
        years = dimnames(results$hiv.positive.mortality)$year
    if (is.null(ages))
        ages = dimnames(results$hiv.positive.mortality)$age
    if (is.null(races))
        races = dimnames(results$hiv.positive.mortality)$race
    if (is.null(subpopulations))
        subpopulations = dimnames(results$hiv.positive.mortality)$subpopulation
    if (is.null(sexes))
    {
        if (use.cdc.categorizations)
            sexes = dimnames(results$hiv.positive.mortality.cdc)$sex
        else
            sexes = dimnames(results$hiv.positive.mortality)$sex
    }
    if (is.null(risks))
    {
        if (use.cdc.categorizations)
            risks = dimnames(results$hiv.positive.mortality.cdc)$risk
        else
            risks = dimnames(results$hiv.positive.mortality)$risk
    }
    if (is.null(continuum))
        continuum = dimnames(results$hiv.positive.mortality)$continuum
    if (is.null(cd4s))
        cd4s = dimnames(results$hiv.positive.mortality)$cd4
    if (is.null(hiv.subsets))
        hiv.subsets = dimnames(results$hiv.positive.mortality)$hiv.subset

    years = as.character(years)

    #Set up keep dimensions
    if (is.null(keep.dimensions))
        keep.dimensions = get.default.keep.dimensions(results,
                                                      years=years,
                                                      ages=ages,
                                                      races=races,
                                                      subpopulations=subpopulations,
                                                      sexes=sexes,
                                                      risks=risks,
                                                      non.hiv.subsets=NULL,
                                                      continuum=may.be.null.continuum,
                                                      cd4s=may.be.null.cd4s,
                                                      hiv.subsets=may.be.null.hiv.subsets,
                                                      use.cdc.categorizations=use.cdc.categorizations)
    #set CDC or non-CDC
    if (use.cdc.categorizations)
        mortality = results$hiv.positive.mortality.cdc
    else
        mortality = results$hiv.positive.mortality

    #Pull the numerators
    if (length(keep.dimensions)==0)
        rv = sum(mortality[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets])
    else
    {
        dim.names = list(year=years, age=ages, race=races, subpopulation=subpopulations, sex=sexes, risk=risks,
                         continuum=continuum, cd4=cd4s, hiv.subset=hiv.subsets)
        rv = mortality[years,ages,races,subpopulations,sexes,risks,continuum,cd4s,hiv.subsets]
        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
        rv = apply(rv, keep.dimensions, sum)
    }

    #If one-dimensional, hydrate up to array
    if (length(keep.dimensions)>1 && is.null(dim(rv)))
    {
        dim.names = list(names(rv))
        names(dim.names) = keep.dimensions

        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }


    #If requested, divide by denominator and normalize
    if (!is.na(per.population))
    {
        denominator = extract.population.subset(results,
                                                years=years,
                                                ages=ages,
                                                races=races,
                                                subpopulations=subpopulations,
                                                sexes=sexes,
                                                risks=risks,
                                                non.hiv.subsets=NULL,
                                                continuum=may.be.null.hiv.subsets,
                                                cd4s=may.be.null.cd4s,
                                                hiv.subsets=may.be.null.hiv.subsets,
                                                keep.dimensions=keep.dimensions,
                                                per.population=NA,
                                                include.hiv.negative = include.hiv.negative.in.denominator,
                                                include.hiv.positive = T,
                                                use.cdc.categorizations = use.cdc.categorizations)

        rv = rv / denominator * per.population
    }

    rv
}


##-------------------------------------------##
##-- INTERNAL FUNCTIONS USED BY OTHER CODE --##
##-------------------------------------------##

#export for now only, for testing
#'@export
do.extract.population.subset <- function(population.1,
                                         population.2=NULL,
                                         years=NULL,
                                         ages=NULL,
                                         races=NULL,
                                         subpopulations=NULL,
                                         sexes=NULL,
                                         risks=NULL,
                                         non.hiv.subsets=NULL,
                                         continuum=NULL,
                                         cd4s=NULL,
                                         hiv.subsets=NULL,
                                         keep.dimensions,
                                         denominator.dimensions=NULL,
                                         per.population=1,
                                         transformation.fn=NULL
)
{
    if (!is.null(transformation.fn))
    {
        if (!is.null(population.1))
            population.1 = transformation.fn(population.1, years=years)
        if (!is.null(denominators.2))
            population.2 = transformation.fn(population.2, years=years)
    }

    if (!is.null(years))
        years = as.character(years)

    #-- Check dimensions to make sure they're in the given arrays --#
    if (is.null(names(dim(population.1))) ||
        (!is.null(population.2) && is.null(names(dim(population.2)))))
        stop("The population arrays must have named dimensions to extract results")

    missing.dims.from.subset.1 = setdiff(keep.dimensions, names(dim(population.1)))
    if (length(missing.dims.from.subset.1)>0)
        stop(paste0("The following dimensions listed as keep dimensions are not present in the given population.1 array: "),
             paste0(missing.dims.from.subset.1, collapse=', '))

    if (!is.null(population.2))
    {
        missing.dims.from.subset.2 = setdiff(keep.dimensions, names(dim(population.2)))
        if (length(missing.dims.from.subset.2)>0)
            stop(paste0("The following dimensions listed as keep dimensions are not present in the given population.2 array: "),
                 paste0(missing.dims.from.subset.2, collapse=', '))
    }

    if (!is.na(per.population) && !is.null(denominator.dimensions))
    {
        missing.denominator.dimensions = setdiff(denominator.dimensions, keep.dimensions)
        if (length(missing.denominator.dimensions)>0)
            stop(paste0("denominator.dimensions must be a subset of keep.dimensions. The following are present in denominator.dimensions but not in keep.dimensions: ",
                        paste0(missing.denominator.dimensions, collapse=', ')))
    }

    #-- PULL FROM FIRST SUBSET --#
    rv = access(population.1,
                year=years,
                age=ages,
                race=races,
                subpopulation=subpopulations,
                sex=sexes,
                risk=risks,
                continuum=continuum,
                cd4=cd4s,
                hiv.subset=hiv.subsets,
                collapse.length.one.dimensions = F)
    accessed.dimnames = dimnames(rv)
    rv = apply(rv, keep.dimensions, sum)

    #-- PULL FROM SECOND SUBSET (if non-null) --#
    if (!is.null(population.2))
        rv = rv + apply(access(population.2,
                               year=years,
                               age=ages,
                               race=races,
                               subpopulation=subpopulations,
                               sex=sexes,
                               risk=risks,
                               non.hiv.subset=non.hiv.subsets,
                               collapse.length.one.dimensions = F),
                        keep.dimensions, sum)

    #-- IF ONE-DIMENSIONAL, HYDRATE UP TO ARRAY --#
    if (is.null(dim(rv)))
    {
        dim(rv) = length(accessed.dimnames[[keep.dimensions]])
        names(dim(rv)) = names(accessed.dimnames[keep.dimensions])
        dimnames(rv) = accessed.dimnames[keep.dimensions]
    }

    #-- IF REQUESTED, NORMALIZE TO POPULATION --#

    if (!is.na(per.population))
    {
        if (is.null(denominator.dimensions))
            denominator = sum(rv)
        else
        {
            denominator = apply(rv, denominator.dimensions, sum)
            if (is.null(dim(denominator)))
            {
                dim(denominator) = length(accessed.dimnames[[denominator.dimensions]])
                names(dim(denominator)) = names(accessed.dimnames[denominator.dimensions])
                dimnames(denominator) = accessed.dimnames[denominator.dimensions]
            }
            denominator = expand.population(denominator,
                                            target.dim.names = dimnames(rv))
        }

        rv = rv / denominator * per.population
    }

    #-- RETURN --#
    rv
}

do.extract.results <- function(numerators,
                               denominators.1,
                               denominators.2,
                               years=NULL,
                               ages=NULL,
                               races=NULL,
                               subpopulations=NULL,
                               sexes=NULL,
                               risks=NULL,
                               non.hiv.subsets=NULL,
                               continuum=NULL,
                               cd4s=NULL,
                               hiv.subsets=NULL,
                               keep.dimensions,
                               include.hiv.positive.in.denominator=T,
                               per.population=100000,
                               transformation.fn=NULL
                               )
{
    if (!is.null(transformation.fn))
    {
        numerators = transformation.fn(numerators, years=years)
        if (!is.null(denominators.1))
            denominators.1 = transformation.fn(denominators.1, years=years)
        if (!is.null(denominators.2))
            denominators.2 = transformation.fn(denominators.2, years=years)
    }

    if (!is.null(years))
        years = as.character(years)

    #-- Check dimensions to make sure they're in the given array --#
    if (is.null(names(dim(numerators))))
        stop("The arrays must have named dimensions to extract results")

    missing.dims.from.numerators = setdiff(keep.dimensions, names(dim(numerators)))
    if (length(missing.dims.from.numerators)>0)
        stop(paste0("The following dimensions listed as keep dimensions are not present in the given numerators array: "),
             paste0(missing.dims.from.numerators, collapse=', '))


    #-- PULL FROM NUMERATORS --#
    rv = access(numerators,
                year=years,
                age=ages,
                race=races,
                subpopulation=subpopulations,
                sex=sexes,
                risk=risks,
                continuum=continuum,
                cd4=cd4s,
                hiv.subset=hiv.subsets,
                collapse.length.one.dimensions = F)
    accessed.dimnames = dimnames(rv)
    rv = apply(rv, keep.dimensions, sum)

    #-- IF ONE-DIMENSIONAL, HYDRATE UP TO ARRAY --#
    if (is.null(dim(rv)))
    {
        dim(rv) = length(accessed.dimnames[[keep.dimensions]])
        dimnames(rv) = accessed.dimnames[keep.dimensions]
    }

    #-- NORMALIZE TO DENOMINATORS, IF REQUESTED --#

    if (!is.na(per.population))
    {
        if (is.null(denominators.1) && is.null(denominators.2))
            denominator = 1
        else
            denominator = do.extract.population.subset(population.1 = denominators.1,
                                                       population.2 = denominators.2,
                                                       years=years,
                                                       ages=ages,
                                                       races=races,
                                                       subpopulations=subpopulations,
                                                       sexes=sexes,
                                                       risks=risks,
                                                       non.hiv.subsets=non.hiv.subsets,
                                                       continuum=continuum,
                                                       cd4s=cd4s,
                                                       hiv.subsets=hiv.subsets,
                                                       keep.dimensions=keep.dimensions,
                                                       per.population = NA)

        rv = rv / denominator * per.population
    }

    #-- RETURN --#
    rv
}


##-----------------------##
##-- LOW-LEVEL HELPERS --##
##-----------------------##

get.default.keep.dimensions <- function(results,
                                        years,
                                        ages,
                                        races,
                                        subpopulations,
                                        sexes,
                                        risks,
                                        non.hiv.subsets,
                                        continuum,
                                        cd4s,
                                        hiv.subsets,
                                        use.cdc.categorizations=F,
                                        at.minimum=NULL
)
{
    aggregate.years=F
    aggregate.ages=is.null(ages) || length(ages)==1 || length(ages)==length(results$ages)
    aggregate.races=is.null(races) || length(races)==1 || length(races)==length(results$races)
    aggregate.subpopulations=is.null(subpopulations) || length(subpopulations)==1 || length(subpopulations)==length(results$subpopulations)
    aggregate.sexes=is.null(sexes) || length(sexes)==1 || (use.cdc.categorizations && length(sexes)==length(results$cdc.sexes)) ||
        (!use.cdc.categorizations && length(sexes)==length(results$sexes))
    aggregate.risks=is.null(risks) || length(risks)==1 || (use.cdc.categorizations && length(risks)==length(results$cdc.risks)) ||
        (!use.cdc.categorizations && length(risks)==length(results$risks))
    aggregate.non.hiv.subsets=is.null(non.hiv.subsets) || length(non.hiv.subsets)==1 || length(non.hiv.subsets)==length(results$non.hiv.subsets)
    aggregate.continuum=is.null(continuum) || length(continuum)==1 || length(continuum)==length(results$continuum)
    aggregate.cd4=is.null(cd4s) || length(cd4s)==1 || length(cd4s)==length(results$cd4s)
    aggregate.hiv.subsets=is.null(hiv.subsets) || length(hiv.subsets==1) || length(hiv.subsets)==length(results$hiv.subsets)

    if (!is.null(at.minimum))
    {
        aggregate.years = aggregate.years & !any(at.minimum=='year')
        aggregate.ages = aggregate.years & !any(at.minimum=='age')
        aggregate.races = aggregate.years & !any(at.minimum=='race')
        aggregate.subpopulations = aggregate.years & !any(at.minimum=='subpopulation')
        aggregate.sexes = aggregate.years & !any(at.minimum=='sex')
        aggregate.risks = aggregate.years & !any(at.minimum=='risk')
        aggregate.non.hiv.subsets = aggregate.years & !any(at.minimum=='non.hiv.subset')
        aggregate.continuum = aggregate.years & !any(at.minimum=='continuum')
        aggregate.cd4 = aggregate.years & !any(at.minimum=='cd4')
        aggregate.hiv.subsets = aggregate.years & !any(at.minimum=='hiv.subset')
    }

    keep.dimensions = c('year','age','race','subpopulation','sex','risk','non.hiv.subset','continuum','cd4','hiv.subset')
    keep.dimensions = keep.dimensions[c(!aggregate.years,
                                        !aggregate.ages,
                                        !aggregate.races,
                                        !aggregate.subpopulations,
                                        !aggregate.sexes,
                                        !aggregate.risks,
                                        !aggregate.non.hiv.subsets,
                                        !aggregate.continuum,
                                        !aggregate.cd4,
                                        !aggregate.hiv.subsets)]

    keep.dimensions
}
