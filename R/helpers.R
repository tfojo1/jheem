

##------------------------##
##------------------------##
##-- HIGH-LEVEL HELPERS --##
##------------------------##
##------------------------##

##-----------------------##
##-----------------------##
##-- MID-LEVEL HELPERS --##
##-----------------------##
##-----------------------##

#makes a string of the format init[repeater]*
#such that the rv is not a subset of any of the strings in the strings argument
get.non.subset.string <- function(strings, init='_2', repeater='z')
{
    rv = init
    while (any(grepl(rv, strings)))
        rv = paste0(init, repeater)

    rv
}


##-----------------------##
##-----------------------##
##-- LOW-LEVEL HELPERS --##
##-----------------------##
##-----------------------##

#returns true if to.check is a scalar
# if specific.value is not NA, also checks whether to.check equals specific.value
is.scalar <- function(to.check, specific.value=NA)
{
    if ((class(to.check) != 'numeric' && class(to.check) != 'integer') ||
        length(to.check) != 1)
        F
    else
    {
        if (!is.na(specific.value) && !is.scalar(to.check=specific.value, specific.value=NA))
            stop("specific.value must be a scalar numeric or integer")
        else if (is.na(specific.value))
            T
        else
            to.check == specific.value
    }
}

index.of <- function(needle, haystack)
{
    sapply(needle, function(x){
        match = haystack==x
        if (any(match))
            (1:length(haystack))[match][1]
        else
            NA
    })
}

all.indices.of <- function(needle, haystack)
{
    rv = lapply(needle, function(x){
        match = haystack==x
        if (any(match))
            (1:length(haystack))[match]
        else
            NULL
    })

    if (length(needle)==1)
        rv[[1]]
    else
        rv
}

setdiff.list <- function(x, y)
{
    x[sapply(x, function(x.elem){
        all(sapply(y, function(y.elem){
            length(x.elem) != length(y.elem) || any(x.elem != y.elem)
        }))
    })]
}

intersect.list <- function(x, y)
{
    x[sapply(x, function(x.elem){
        any(sapply(y, function(y.elem){
            length(x.elem) == length(y.elem) && all(x.elem == y.elem)
        }))
    })]
}

list.equals <- function(x, y)
{
    if (is.null(x) && is.null(y))
        T
    else
        !is.null(x) && !is.null(y) &&
        length(x)==length(y) &&
        all(sapply(1:length(x), function(i){
            length(x[[i]]) == length(y[[i]]) &&
                all(x[[i]] == y[[i]])
        }))
}

vector.equals <- function(x, y)
{
    (is.null(x) && is.null(y)) ||
        (!is.null(x) && !is.null(y) && length(x)==length(y) && all(x==y))
}

get.ordinal <- function(nums)
{
    sapply(nums, function(num){

        two.digits = num %% 100
        one.digit = num %% 10

        if (two.digits > 10 && two.digits <20)
            paste0(num, 'th')
        else if (one.digit==1)
            paste0(num, 'st')
        else if (one.digit==2)
            paste0(num, 'nd')
        else if (one.digit==3)
            paste0(num, 'rd')
        else
            paste0(num, 'th')

    })
}

get.character.list <- function(values,
                               pre.last.sep=', ',
                               last.sep.if.two = ' and ',
                               last.sep.if.three.or.more = ', and ')
{
    n = length(values)
    if (n==0)
        ''
    else if (n==1)
        values
    else if (n==2)
        paste0(values[1], last.sep.if.two, values[2])
    else
    {
        paste0(paste0(values[-n], collapse=pre.last.sep),
               last.sep.if.three.or.more,
               values[n])
    }
}
