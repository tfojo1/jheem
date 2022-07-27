

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


