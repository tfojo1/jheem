
#result.year.minus.ode.year - eg, 'year end 2008' is actually represented by the ODE as 2009 (year start) -> result.year.minus.ode.year = -1
hydrate.ode.results <- function(ode.results,
                                jheem,
                                parameters,
                                keep.years,
                                result.year.minus.ode.year=-1)
{
    ##-- SET UP THE RV and INDEX --##
    rv = list()
    index = 1

    ode.results = as.matrix(ode.results)

    ##-- SET UP THE YEARS --##
    years = as.character(ode.results[,1]+result.year.minus.ode.year)
    num.years = length(years)
    keep.years = intersect(as.character(keep.years), years)
    keep.mask = sapply(years, function(year){any(year==keep.years)})
    num.keep = length(keep.years)

    ##-- SET UP DIMENSION NAMES in RV --#
    rv$years = as.numeric(keep.years)
    rv$ages = jheem$age$labels
    rv$races = jheem$race
    rv$subpopulations = jheem$subpopulations
    rv$sexes = jheem$sex
    rv$risks = jheem$risk.strata
    rv$non.hiv.subsets = jheem$nonhiv.subsets
    rv$continuum = jheem$continuum.of.care
    rv$cd4 = jheem$cd4.strata
    rv$hiv.subsets = jheem$hiv.subsets
    rv$diagnosed.continuum.states = jheem$DIAGNOSED_CONTINUUM_STATES

    ##-- SET UP DIMNAMES --##
    dimnames.hiv.positive = c(list(year=keep.years), get.dimnames.hiv(jheem))
    dimnames.hiv.negative = c(list(year=keep.years), get.dimnames.nonhiv(jheem))
    shared.dims = length(get.dimnames.general(jheem))
    dimnames.all = c(list(year=keep.years), get.dimnames.all(jheem))

    ##-- HIV POSITIVE --##
    rv$hiv.positive = array(ode.results[keep.mask,index + 1:parameters$NUM_HIV_STATES],
                            dim = sapply(dimnames.hiv.positive, length),
                            dimnames = dimnames.hiv.positive)
    index = index + parameters$NUM_HIV_STATES

    ##-- HIV NEGATIVE --##
    rv$hiv.negative = array(ode.results[keep.mask,index + 1:parameters$NUM_NONHIV_STATES],
                            dim = sapply(dimnames.hiv.negative, length),
                            dimnames = dimnames.hiv.negative)
    index = index + parameters$NUM_NONHIV_STATES

    ##-- ADD THEM UP TO GET TOTAL POPULATION --##
#    rv$total.population = rowSums(rv$hiv.positive, dims=shared.dims) + rowSums(rv$hiv.negative, dims=shared.dims)

    ##-- NEW CASES --##
    new.cases.dim.names = get.collapsed.dimnames(full.dimnames=dimnames.all,
                                                 collapsed.dimname.names = c('year', jheem$parameters$incidence.keep.dimensions))

    n.new.cases = prod(sapply(new.cases.dim.names[-1], length))

    new.cases = ode.results[,index + 1:n.new.cases]
    new.cases[-1,] = new.cases[-1,] - new.cases[-num.years,]
    new.cases[1,] = 0
    rv$incidence = array(new.cases[keep.mask,],
                         dim = sapply(new.cases.dim.names, length),
                         dimnames = new.cases.dim.names)
    index = index + n.new.cases

    ##-- HIV SPECIFIC MORTALITY --##
    if (parameters$TRACK_HIV_SPECIFIC_MORTALITY)
    {
        hiv.specific.mortality.dimnames = get.collapsed.dimnames(full.dimnames=dimnames.hiv.positive,
                                                                 collapsed.dimname.names = c('year', jheem$parameters$hiv.specific.mortality.keep.dimensions))
        n.hiv.specific.mortality = prod(sapply(hiv.specific.mortality.dimnames[-1], length))

        mortality = ode.results[,index + 1:n.hiv.specific.mortality]
        mortality[-1,] = mortality[-1,] - mortality[-num.years,]
        mortality[1,] = 0
        rv$hiv.specific.mortality = array(mortality[keep.mask,],
                                          dim = sapply(hiv.specific.mortality.dimnames, length),
                                          dimnames = hiv.specific.mortality.dimnames)
        index = index + n.hiv.specific.mortality
    }

    ##-- HIV POSITIVE MORTALITY --##
    if (parameters$TRACK_OVERALL_HIV_MORTALITY)
    {
        overall.hiv.mortality.dimnames = get.collapsed.dimnames(full.dimnames=dimnames.hiv.positive,
                                                                 collapsed.dimname.names = c('year', jheem$parameters$overall.hiv.mortality.keep.dimensions))
        n.overall.hiv.mortality = prod(sapply(overall.hiv.mortality.dimnames[-1], length))

        mortality = ode.results[,index + 1:n.overall.hiv.mortality]
        mortality[-1,] = mortality[-1,] - mortality[-num.years,]
        mortality[1,] = 0
        rv$hiv.positive.mortality = array(mortality[keep.mask,],
                                          dim = sapply(overall.hiv.mortality.dimnames, length),
                                          dimnames = overall.hiv.mortality.dimnames)
        index = index + n.overall.hiv.mortality
    }

    ##-- HIV NEGATIVE MORTALITY --##
    if (parameters$TRACK_OVERALL_NONHIV_MORTALITY)
    {
        nonhiv.mortality.dimnames = get.collapsed.dimnames(full.dimnames=dimnames.hiv.negative,
                                                                collapsed.dimname.names = c('year', jheem$parameters$overall.nonhiv.mortality.keep.dimensions))
        n.nonhiv.mortality = prod(sapply(nonhiv.mortality.dimnames[-1], length))

        mortality = ode.results[,index + 1:n.nonhiv.mortality]
        mortality[-1,] = mortality[-1,] - mortality[-num.years,]
        mortality[1,] = 0
        rv$hiv.negative.mortality = array(mortality[keep.mask,],
                                          dim = sapply(nonhiv.mortality.dimnames, length),
                                          dimnames = nonhiv.mortality.dimnames)
        index = index + parameters$NUM_NONHIV_STATES
    }

    ##-- TRACKED TRANSITIONS --##

    transition.names = c(jheem$tracked.transitions$names.for.collapsed.hiv, jheem$tracked.transitions$names.for.collapsed.nonhiv)
    if (length(transition.names)>0)
    {
        raw.transitions = ode.results[, index + (1:length(transition.names))]
        raw.transitions[-1,] = raw.transitions[-1,] - raw.transitions[-num.years,]
        raw.transitions[1,] = 0

        rv$tracked.transitions = lapply(names(jheem$tracked.transitions$dimnames.for.tracked), function(name){
            dimnames.for.transition = c(list(year=keep.years), jheem$tracked.transitions$dimnames.for.tracked[[name]])
            name.mask = transition.names == name
            array(raw.transitions[keep.mask, name.mask],
                  dim = sapply(dimnames.for.transition, length),
                  dimnames = dimnames.for.transition)
        })
        names(rv$tracked.transitions) = names(jheem$tracked.transitions$dimnames.for.tracked)

        index = index + length(transition.names)

        ##-- Pull out new diagnoses as a special tracked transition --##
        if (any(names(rv$tracked.transitions)=='new_diagnoses'))
        {
            rv$new.diagnoses = rv$tracked.transitions[['new_diagnoses']]
            rv$tracked.transitions[['new_diagnoses']] = NULL
        }
    }

    ##-- ADD CDC CATEGORIZATIONS --33
    #rv = recategorize.for.cdc(rv)

    ##-- RETURN IT --##

    class(rv) = 'jheem.results'
    rv
}

CDC.SEX = c('male','female')
CDC.RISK = c('msm','idu', 'msm_idu', 'heterosexual')
#'@export
recategorize.for.cdc <- function(results)
{
    results$cdc.sexes = CDC.SEX
    results$cdc.risks = CDC.RISK

    for (elem.name in names(results))
    {
     #   print(elem.name)
        if (class(results[[elem.name]])=='array')
            results[[paste0(elem.name, '.cdc')]] = convert.to.cdc.array(results[[elem.name]])
    }

    for (elem.name in names(results$tracked.transitions))
    {
    #    print(elem.name)
        if (class(results$tracked.transitions[[elem.name]])=='array')
            results$tracked.transitions[[paste0(elem.name, '.cdc')]] = convert.to.cdc.array(results$tracked.transitions[[elem.name]])
    }

    results
}

convert.to.cdc.array <- function(arr)
{
    dim.names = dimnames(arr)
    if (any(dim.names[['sex']]!='all') || any(dim.names[['risk']]!='all'))
    {
        dim.names[['sex']] = CDC.SEX
        dim.names[['risk']] = CDC.RISK

        pre.risk.dimnames.old = dimnames(arr)[1:6]
        pre.risk.dimnames.new = dim.names[1:6]
        post.risk.dimnames = dim.names[7:length(dim.names)]
        post.risk.dim = prod(sapply(post.risk.dimnames, length))

        #set up collapsed rv
        rv = array(0, dim=c(sapply(pre.risk.dimnames.new, length), other=post.risk.dim),
                   dimnames=c(pre.risk.dimnames.new, list(other=NULL)))

        #collapse arr
        dim(arr) = c(sapply(pre.risk.dimnames.old, length), other=post.risk.dim)
        dimnames(arr) = c(pre.risk.dimnames.old, list(other=NULL))

        #pull msm
        rv[,,,,'male','msm',] = arr[,,,,'msm','never_IDU',]

        #pull idu
        rv[,,,,'female','idu',] = arr[,,,,'female','active_IDU',] + arr[,,,,'female','IDU_in_remission',]
        rv[,,,,'male','idu',] = arr[,,,,'heterosexual_male','active_IDU',] + arr[,,,,'heterosexual_male','IDU_in_remission',]

        #pull msm+idu
        rv[,,,,'male','msm_idu',] = arr[,,,,'msm','active_IDU',] + arr[,,,,'msm','IDU_in_remission',]

        #pull heterosexual
        rv[,,,,'female','heterosexual',] = arr[,,,,'female','never_IDU',]
        rv[,,,,'male','heterosexual',] = arr[,,,,'heterosexual_male','never_IDU',]

        #hydrate and return
        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
        rv
    }
    else
        arr
}


#'@title Combine Two JHEEM Results Objects into One by year
#'
#'@param results1,results2 Two objects of class jheem.results. If years from the two overlap, data from results1 will be used
#'
#'@export
combine.jheem.results <- function(results1, results2)
{
    years = unique(sort(c(results1$years, results2$years)))

    conversion.fn <- function(r1, r2)
    {
        rv=lapply(names(r1), function(elem.name){
            elem = r1[[elem.name]]

            if (elem.name == 'years')
                years
            else if (any(class(elem)=='array'))
            {
                #grab elem1 and elem2
                elem1 = elem
                years1 = dimnames(elem1)[['year']]
                elem2 = r2[[elem.name]]
                years2 = dimnames(elem2)[['year']]

                #smush elem1 and elem2 to 2d
                dim(elem1) = c(dim(elem1)[1], prod(dim(elem1)[-1]))
                dim(elem2) = c(dim(elem2)[1], prod(dim(elem2)[-1]))

                #set up a (smushed) rv
                dim.names = dimnames(elem)
                dim.names[['year']] = as.character(years)
                elem = array(0, dim=c(length(dim.names[[1]]), prod(sapply(dim.names[-1], length))),
                             dimnames=list(dim.names[[1]], NULL))

#                access(elem, year=dimnames(elem2)[['year']]) = elem2
#                access(elem, year=dimnames(elem1)[['year']]) = elem1

                #pull into the rv
                elem[years1,] = elem1
                elem[years2,] = elem2

                #hydrate dimensions and return
                dim(elem) = sapply(dim.names, length)
                dimnames(elem) = dim.names

                elem
            }
            else if (any(class(elem)=='list'))
                conversion.fn(elem, r2[[elem.name]])
            else
                elem
        })
        names(rv) = names(r1)
        rv
    }

    rv = conversion.fn(results1, results2)
    class(rv)='jheem.results'
    rv
}

#'@title Access a Subset of JHEEM Results by Year
#'
#'@param results An object of class jheem.results to subset
#'@param years The years to subset
#'
#'@export
subset.jheem.results <- function(results, years)
{
    conversion.fn <- function(rr)
    {
        rv = lapply(names(rr), function(elem.name){
            elem = rr[[elem.name]]

            if (elem.name == 'years')
                years
            else if (any(class(elem)=='array'))
            {
#                access(elem, year=as.character(years), collapse.length.one.dimensions = F)
                orig.dim.names = dimnames(elem)
                new.dim.names = orig.dim.names
                new.dim.names[['year']] = as.character(years)

                dim(elem) = c(dim(elem)[1], prod(dim(elem)[-1]))
                dimnames(elem)[[1]] = orig.dim.names[[1]]

                elem = elem[as.character(years),]
                dim(elem) = sapply(new.dim.names, length)
                dimnames(elem) = new.dim.names

                elem
            }
            else if (any(class(elem)=='list'))
                conversion.fn(elem)
            else
                elem
        })
        names(rv) = names(rr)
        rv
    }

    rv = conversion.fn(results)
    class(rv)='jheem.results'
    rv
}

#'@title Expand JHEEM Results to Contain New Dimensions
#'
#'@description Expands the components of a JHEEM Results object to have new dimensions (the values in these new dimensions are set to zero)
#'
#'@param results An object of class jheem.results
#'@param new.jheem An object of class jheem that describes the new dimensions. Each dimension (age, race, etc.) must be a superset of the dimensions in the results object
#'
#'@export
expand.jheem.results <- function(results, new.jheem)
{
    all.dimnames = all.dimnames.from = all.dimnames.to =  get.dimnames.all(new.jheem)
    names(all.dimnames.from) = paste0(names(all.dimnames), '.from')
    names(all.dimnames.to) = paste0(names(all.dimnames), '.to')
    all.dimnames = c(all.dimnames, all.dimnames.from, all.dimnames.to)
    all.dimnames = c(list(year=as.character(results$year)), all.dimnames)


    conversion.fn <- function(rr, is.transition=F)
    {
        rv = lapply(names(rr), function(elem.name){
            elem = rr[[elem.name]]

            if (any(elem.name == names(all.dimnames)))
                all.dimnames[[elem.name]]
            else if (elem.name=='new.diagnoses')
                NULL #we will fill this in later, from the tracked transition
            else if (any(elem.name == paste0(names(all.dimnames), 's')))
                all.dimnames[[substr(elem.name, 1, nchar(elem.name)-1)]]
            else if (any(class(elem)=='array'))
            {
                if (is.transition)
                    dimnames.superset = c(list(year=as.character(results$year)), new.jheem$tracked.transitions$dimnames.for.tracked[[elem.name]])
                else
                    dimnames.superset = all.dimnames

                old.dimnames = dimnames(elem)
                new.dimnames.for.elem = dimnames.superset[names(old.dimnames)]
                new.elem = array(0, dim=sapply(new.dimnames.for.elem, length), dimnames=new.dimnames.for.elem)

                access(new.elem,
                       year=old.dimnames[['year']],
                       age=old.dimnames[['age']],
                       race=old.dimnames[['race']],
                       subpopulation=old.dimnames[['subpopulation']],
                       sex=old.dimnames[['sex']],
                       risk=old.dimnames[['risk']],
                       non.hiv.subset=old.dimnames[['non.hiv.subset']],
                       continuum=old.dimnames[['continuum']],
                       cd4=old.dimnames[['cd4']],
                       hiv.subset=old.dimnames[['hiv.subset']],
                       age.from=old.dimnames[['age.from']],
                       race.from=old.dimnames[['race.from']],
                       subpopulation.from=old.dimnames[['subpopulation.from']],
                       sex.from=old.dimnames[['sex.from']],
                       risk.from=old.dimnames[['risk.from']],
                       non.hiv.subset.from=old.dimnames[['non.hiv.subset.from']],
                       continuum.from=old.dimnames[['continuum.from']],
                       cd4.from=old.dimnames[['cd4.from']],
                       hiv.subset.from=old.dimnames[['hiv.subset.from']],
                       age.to=old.dimnames[['age.to']],
                       race.to=old.dimnames[['race.to']],
                       subpopulation.to=old.dimnames[['subpopulation.to']],
                       sex.to=old.dimnames[['sex.to']],
                       risk.to=old.dimnames[['risk.to']],
                       non.hiv.subset.to=old.dimnames[['non.hiv.subset.to']],
                       continuum.to=old.dimnames[['continuum.to']],
                       cd4.to=old.dimnames[['cd4.to']],
                       hiv.subset.to=old.dimnames[['hiv.subset.to']]) = elem

                new.elem
            }
            else if (any(class(elem)=='list'))
                conversion.fn(elem, elem.name=='tracked.transitions')
            else
                elem
        })
        names(rv) = names(rr)
        rv
    }

    rv = conversion.fn(results)
    rv$new.diagnoses = rv$tracked.transitions$new_diagnoses
    class(rv)='jheem.results'
    rv
}


#A helper to make a set of dimnames where some dimensions have been summed out
get.collapsed.dimnames <- function(full.dimnames,
                                    collapsed.dimname.names)
{
    collapsed.dimnames = full.dimnames
    summed.out.mask = sapply(names(full.dimnames), function(name){
        all(name != collapsed.dimname.names)
    })
    collapsed.dimnames[summed.out.mask] = 'all'
    collapsed.dimnames
}
