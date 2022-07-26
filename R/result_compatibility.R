

get.jheem.result.dim.names <- function(results)
{
    if (!is(results, 'jheem.results'))
        stop("results must be an object of class 'jheem.results'")

    list(
        age = result$ages,
        race = result$races,
        #location = result$locations,
        subpopulations = jheem$subpopulations,
        sexes = result$sexes,
        risks = result$risks,
        non.hiv.subsets = result$nonhiv.subsets,
        continuum = result$continuum,
        cd4 = result$cd4,
        hiv.subsets = result$hiv.subsets
    )
}

jheem.dimensions.match.results <- function(results, jheem)
{
    results.dim.names = get.jheem.result.dim.names(results)
    jheem.dim.names = get.dimnames(jheem,
                                   age=T, race=T,
                                   #location=T,
                                   subpopulation=T,
                                   sex=T, risk=T,
                                   non.hiv.subset=T,
                                   continuum.of.care=T, cd4=T, hiv.subset=T)

    all(sapply(names(results.dim.names), function(name){
        length(results.dim.names[[name]]==length(jheem.dim.names[[name]]) &&
                   all(results.dim.names[[name]]==jheem.dim.names[[name]]))
    }))
}

jheem.result.dimensions.are.subset <- function(results, jheem)
{
    results.dim.names = get.jheem.result.dim.names(results)
    jheem.dim.names = get.dimnames(jheem,
                                   age=T, race=T,
                                   #location=T,
                                   subpopulation=T,
                                   sex=T, risk=T,
                                   non.hiv.subset=T,
                                   continuum.of.care=T, cd4=T, hiv.subset=T)

    all(sapply(names(results.dim.names), function(name){
        length(setdiff(jheem.dim.names[[name]], results.dim.names[[name]]))==0
    }))
}

#'@description Expand the results arrays in a jheem.results object to accomodate a more expansive JHEEM instantiation
#'
#'@export
expand.jheem.results <- function(results, jheem)
{
    to.expand.names = names(results)[sapply(results, is.array)]

    rv = initialize.jheem.results(jheem, years=results$years)
    for (one.name.to.expand in to.expand.names)
        rv[[one.name.to.expand]] = expand.jheem.component(results[[one.name.to.expand]], jheem)

    rv
}

expand.jheem.component <- function(comp, jheem, default.value=0)
{
    old.dim.names = dimnames(comp)
    new.dim.names = get.dimnames(jheem,
                                   age=!is.null(old.dim.names$age),
                                   race=!is.null(old.dim.names$race),
                                   #location=!is.null(old.dim.names$location),
                                   subpopulation=!is.null(old.dim.names$subpopulation),
                                   sex=!is.null(old.dim.names$sex),
                                   risk=!is.null(old.dim.names$risk),
                                   non.hiv.subset=!is.null(old.dim.names$non.hiv.subset),
                                   continuum.of.care=!is.null(old.dim.names$continuum),
                                   cd4=!is.null(old.dim.names$cd4),
                                   hiv.subst=!is.null(old.dim.names$hiv.subset))
    if (!is.null(old.dim.names$year))
        new.dim.names = c(old.dim.names['year'],
                          new.dim.names)

    rv = array(default.value,
               dim=sapply(new.dim.names, length),
               dimnames=new.dim.names)

    access(rv,
           year=old.dim.names$year,
           age=old.dim.names$age,
           race=old.dim.names$race,
           #location=old.dim.names$location,
           subpopulation=old.dim.names$subpopulation,
           sex=old.dim.names$sex,
           risk=old.dim.names$risk,
           non.hiv.subset=old.dim.names$non.hiv.subset,
           continuum=old.dim.names$continuum,
           cd4=old.dim.names$cd4,
           hiv.subset=old.dim.names$hiv.subset) = comp

    rv
}
