

##--------------------------------------------##
##-- THE MAIN/GENERAL GET DIMNAMES FUNCTION --##
##--------------------------------------------##

#The main/general helper
get.dimnames <- function(jheem,
                         age=F,
                         race=F,
                         subpopulation=F,
                         sex=F,
                         risk=F,
                         non.hiv.subset=F,
                         continuum.of.care=F,
                         cd4=F,
                         hiv.subset=F,
                         age.from=F,
                         race.from=F,
                         subpopulation.from=F,
                         sex.from=F,
                         risk.from=F,
                         non.hiv.subset.from=F,
                         continuum.of.care.from=F,
                         cd4.from=F,
                         hiv.subset.from=F,
                         age.to=F,
                         race.to=F,
                         subpopulation.to=F,
                         sex.to=F,
                         risk.to=F,
                         non.hiv.subset.to=F,
                         continuum.of.care.to=F,
                         cd4.to=F,
                         hiv.subset.to=F)
{
    rv = list(age=jheem$age$labels, race=jheem$race, subpopulation=jheem$subpopulations,
              sex=jheem$sex, risk=jheem$risk.strata, non.hiv.subset=jheem$nonhiv.subsets,
              continuum=jheem$continuum.of.care, cd4=jheem$cd4.strata, hiv.subset=jheem$hiv.subsets,
              age.from=jheem$age$labels, race.from=jheem$race, subpopulation.from=jheem$subpopulations,
              sex.from=jheem$sex, risk.from=jheem$risk.strata, non.hiv.subset.from=jheem$nonhiv.subsets,
              continuum.from=jheem$continuum.of.care, cd4.from=jheem$cd4.strata, hiv.subset.from=jheem$hiv.subsets,
              age.to=jheem$age$labels, race.to=jheem$race, subpopulation.to=jheem$subpopulations,
              sex.to=jheem$sex, risk.to=jheem$risk.strata, non.hiv.subset.to=jheem$nonhiv.subsets,
              continuum.to=jheem$continuum.of.care, cd4.to=jheem$cd4.strata, hiv.subset.to=jheem$hiv.subsets)
    mask = c(age, race, subpopulation, sex, risk, non.hiv.subset, continuum.of.care, cd4, hiv.subset,
             age.from, race.from, subpopulation.from, sex.from, risk.from, non.hiv.subset.from, continuum.of.care.from, cd4.from, hiv.subset.from,
             age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to, continuum.of.care.to, cd4.to, hiv.subset.to)

    rv[mask]
}



get.dimnames.by.name <- function(jheem,
                                 dim.name.names)
{
    get.dimnames(jheem,
                 age=any(dim.name.names=='age'),
                 race=any(dim.name.names=='race'),
                 subpopulation=any(dim.name.names=='subpopulation'),
                 sex=any(dim.name.names=='sex'),
                 risk=any(dim.name.names=='risk'),
                 non.hiv.subset=any(dim.name.names=='non.hiv.subset'),
                 continuum.of.care=any(dim.name.names=='continuum.of.care') || any(dim.name.names=='continuum'),
                 cd4=any(dim.name.names=='cd4'),
                 hiv.subset=any(dim.name.names=='hiv.subset'),
                 age.from=any(dim.name.names=='age.from'),
                 race.from=any(dim.name.names=='race.from'),
                 subpopulation.from=any(dim.name.names=='subpopulation.from'),
                 sex.from=any(dim.name.names=='sex.from'),
                 risk.from=any(dim.name.names=='risk.from'),
                 non.hiv.subset.from=any(dim.name.names=='non.hiv.subset.from'),
                 continuum.of.care.from=any(dim.name.names=='continuum.of.care.from') || any(dim.name.names=='continuum.from'),
                 cd4.from=any(dim.name.names=='cd4.from'),
                 hiv.subset.from=any(dim.name.names=='hiv.subset.from'),
                 age.to=any(dim.name.names=='age.to'),
                 race.to=any(dim.name.names=='race.to'),
                 subpopulation.to=any(dim.name.names=='subpopulation.to'),
                 sex.to=any(dim.name.names=='sex.to'),
                 risk.to=any(dim.name.names=='risk.to'),
                 non.hiv.subset.to=any(dim.name.names=='non.hiv.subset.to'),
                 continuum.of.care.to=any(dim.name.names=='continuum.of.care.to') || any(dim.name.names=='continuum.to'),
                 cd4.to=any(dim.name.names=='cd4.to'),
                 hiv.subset.to=any(dim.name.names=='hiv.subset.to'))
}

##-----------------------------------------------------##
##-- GET DIMNAMES FUNCTIONS FOR SPECIFIC POPULATIONS --##
##-----------------------------------------------------##


get.dimnames.nonhiv <- function(jheem)
{
    get.dimnames(jheem, age = T, race = T, subpopulation = T, sex = T, risk = T,
                 non.hiv.subset = T, continuum.of.care = F, cd4 = F, hiv.subset = F)
}

get.dimnames.hiv <- function(jheem)
{
    get.dimnames(jheem, age = T, race = T, subpopulation = T, sex = T, risk = T,
                 non.hiv.subset = F, continuum.of.care = T, cd4 = T, hiv.subset = T)
}

get.dimnames.general <- function(jheem)
{
    get.dimnames(jheem, age = T, race = T, subpopulation = T, sex = T, risk = T,
                 non.hiv.subset = F, continuum.of.care = F, cd4 = F, hiv.subset = F)
}

get.dimnames.all <- function(jheem)
{
    get.dimnames(jheem, age = T, race = T, subpopulation = T, sex = T, risk = T,
                 non.hiv.subset = T, continuum.of.care = T, cd4 = T, hiv.subset = T)
}

get.dimnames.contacts <- function(jheem)
{
    get.dimnames(jheem,
                 age.from = T, race.from = T, subpopulation.from = T, sex.from = T, risk.from = T,
                 age.to = T, race.to = T, subpopulation.to = T, sex.to = T, risk.to = T,
                 non.hiv.subset.to=T)
}

get.dimnames.transitions <- function(jheem, hiv=F)
{
    get.dimnames(jheem,
                 age = T, race = T,
                 subpopulation.from = T, sex.from = T, risk.from = T,
                 non.hiv.subset.from = !hiv,
                 continuum.of.care.from = hiv, cd4.from = hiv, hiv.subset.from = hiv,
                 subpopulation.to = T, sex.to = T, risk.to = T,
                 non.hiv.subset.to = !hiv,
                 continuum.of.care.to = hiv, cd4.to = hiv, hiv.subset.to = hiv)
}

get.dimnames.birth.proportions <- function(jheem, from.hiv=F, to.hiv=from.hiv)
{
    get.dimnames(jheem,
                 age.from = T, race.from = T,
                 subpopulation.from = T, sex.from = T, risk.from = T,
                 non.hiv.subset.from = !from.hiv,
                 continuum.of.care.from = from.hiv, cd4.from = from.hiv, hiv.subset.from = from.hiv,
                 subpopulation.to = T, sex.to = T, risk.to = T,
                 non.hiv.subset.to = !to.hiv,
                 continuum.of.care.to = to.hiv, cd4.to = to.hiv, hiv.subset.to = to.hiv)
}


##-------------##
##-- HELPERS --##
##-------------##


identify.dimnames.in.array <- function(jheem, arr,
                                       single.occurence.as='.to',
                                       error.if.no.match=T,
                                       age.single.as=single.occurence.as,
                                       race.single.as=single.occurence.as,
                                       subpopulation.single.as=single.occurence.as,
                                       sex.single.as=single.occurence.as,
                                       risk.single.as=single.occurence.as,
                                       non.hiv.subset.single.as=single.occurence.as,
                                       continuum.single.as=single.occurence.as,
                                       cd4.single.as=single.occurence.as,
                                       hiv.subset.single.as=single.occurence.as)
{
    if (is.null(dimnames(arr)))
        return (character())

    single.occurence.suffixes = c(age=age.single.as, race=race.single.as, subpopulation=subpopulation.single.as,
                                  sex=sex.single.as, risk=risk.single.as, non.hiv.subset=non.hiv.subset.single.as,
                                  continuum=continuum.single.as, cd4=cd4.single.as, hiv.subset=hiv.subset.single.as)

    if (is.null(names(dimnames(arr))))
    {
        #Figure out which dimension (irrespective of to/from) they match to
        all.possible.dim.names = get.dimnames.all(jheem)
        matched.names = sapply(dimnames(arr), function(one.dim.name.arr){
            match.mask = sapply(all.possible.dim.names, vector.equals, one.dim.name.arr)

            if (!any(match.mask))
            {
                if (error.if.no.match)
                    stop("The dimension names: [", paste0("'", one.dim.name.arr, "'", collapse=", "),
                         "] have no match in the given JHEEM")
                else
                    NA
            }
            else
                names(all.possible.dim.names)[match.mask][1]
        })
    }
    else #Just go by the names(dimnames(arr))
        matched.names = names(dimnames(arr))


    #Classify as to/from
    tabled.matched.names = table(matched.names)
    to.from.mask = grepl('.to', names(tabled.matched.names)) | grepl('.from', names(tabled.matched.names))

    #If there are more than two occurrences of any name, throw an error
    if (any(tabled.matched.names>2 & !to.from.mask))
        stop("A dimension can appear at most two times in the array. The following dimensions appear more than twice: ",
             paste0("'", names(tabled.matched.names)[tabled.matched.names>2], "'", collapse=", "))

    #The names with one match get single.occurrence.as appended
    for (name in names(tabled.matched.names)[tabled.matched.names==1 & !to.from.mask])
        matched.names[matched.names==name] = paste0(name, single.occurence.suffixes[name])

    #The names with two matches get from and to appended
    for (name in names(tabled.matched.names)[tabled.matched.names==2 & !to.from.mask])
        matched.names[matched.names==name] = paste0(name, c(".from", ".to"))

    matched.names
}

set.dim.name.names <- function(arr, dim.name.names)
{
    dim.names = dimnames(arr)
    names(dim.names) = dim.name.names

    dim(arr) = sapply(dim.names, length)
    dimnames(arr) = dim.names

    arr

}


#if the given array has dimensions matching the dimnames
array.matches.dimnames <- function(arr, dim.names, allow.na=F, allow.missing.length.one.dims=F)
{
    arr.dim.names = dimnames(arr)
    arr.dim = dim(arr)

    target.dim.names = dim.names
    target.dim = sapply(target.dim.names, length)

    if (is.null(dim(arr)))
        return(F)

    if (allow.missing.length.one.dims)
    {
        arr.length.one = arr.dim == 1
        arr.dim.names = arr.dim.names[!arr.length.one]
        arr.dim = arr.dim[!arr.length.one]

        target.length.one = target.dim == 1
        target.dim.names = target.dim.names[!target.length.one]
        target.dim = target.dim[!target.length.one]
    }

    length(arr.dim) == length(target.dim) &&
        all(arr.dim == target.dim) &&
        all(sapply(1:length(target.dim.names), function(i){
            all(target.dim.names[[i]]==arr.dim.names[[i]])
        })) &&
        (allow.na || all(!is.na(arr)))
}


#returns the numeric indices (with names) of the 'x.from' dimensions present in the given contact array
get.contact.from.dim.indices <- function(jheem, contact.array)
{
    possible.from.names = paste0(names(get.dimnames.general(jheem)), '.from')
    from.dims = 1:length(possible.from.names)
    names(from.dims) = possible.from.names
    names.of.dims = names(dimnames(contact.array))

    from.dims = from.dims[sapply(possible.from.names, function(name){
        any(names.of.dims == name)
    })]

    from.dims
}

#returns a boolean vector indicating which of the standard dimensions
# for the JHEEM are included in the given array
#Matches by dimnames *NOT by names(dimnames)
dims.present <- function(jheem, arr)
{
    possible.dims = get.dimnames.all(jheem)

    rv = sapply(possible.dims, function(x.elem){
        any(sapply(dimnames(arr), function(y.elem){
            length(x.elem) == length(y.elem) && all(x.elem == y.elem)
        }))
    })

    names(rv) = names(possible.dims)

    rv
}



##----------------------##
##-- RENAMING HELPERS --##
##----------------------##

# Returns a list of names such that repeats in the list are appended with the suffix
# ie
# if (all(dim.names[[i]]==dim.names[[j]])) for any j<i
# dim.names[[i]] <- paste0(dim.names[[i]], suffix)
append.repeat.names <- function(dim.names, suffix)
{
    if (length(dim.names)==1)
        return(dim.names)

    for (i in 2:length(dim.names))
    {
        is.repeat = any(sapply(1:(i-1), function(j){
            length(dim.names[[i]]) == length(dim.names[[j]]) && all(dim.names[[i]] == dim.names[[j]])
        }))

        if (is.repeat)
            dim.names[[i]] = paste0(dim.names[[i]], suffix)
    }

    dim.names
}
