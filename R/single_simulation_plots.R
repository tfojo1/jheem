
RED = '#C75D4D'
GREEN = '#458B00'
#GREEN = 'darkgreen'
BLUE = '#3278FA'
ORANGE = 'darkorange3'

DEFAULT.SHAPES = c(21,24,23,22,25)

#'@param transformation.fn Applied only to the JHEEM
#'
#'@export
jplot.subpopulation <- function(jheem.results,
                                reference=NULL,
                                years,
                                split.by.dimensions=NULL,
                                facet.by.dimensions=NULL,
                                jheem.color=BLUE,
                                reference.color=GREEN,
                                show.labels=F,
                                label.size=5,
                                line.size=2,
                                point.size=5,
                                jheem.name='Simulated',
                                reference.name='Actual',
                                per.population=1,
                                transformation.fn=NULL,
                                as.line=F,
                                shapes=DEFAULT.SHAPES
                                )
{
    if (as.line)
    {
        split.by.dimensions = setdiff(split.by.dimensions, 'year')
        facet.by.dimensions = setdiff(facet.by.dimensions, 'year')

        non.split.keep.dimensions = c(facet.by.dimensions, 'year')
    }
    else
        non.split.keep.dimensions = facet.by.dimensions

    df = NULL
    if (!is.null(jheem.results))
    {
        jheem.years = intersect(as.numeric(years), as.numeric(jheem.results$years))
        df = prepare.population.subset.df(arr=jheem.results,
                                              years=jheem.years,
                                              split.by.dimensions = split.by.dimensions,
                                              facet.by.dimensions = non.split.keep.dimensions,
                                              per.population = per.population,
                                          transformation.fn=transformation.fn)
        df$Population=jheem.name
    }

    if (!is.null(reference))
    {
        ref.years = intersect(as.numeric(years), as.numeric(dimnames(reference)[['year']]))
        ref.df = prepare.population.subset.df(arr=reference,
                                              years=ref.years,
                                              split.by.dimensions = split.by.dimensions,
                                              facet.by.dimensions = non.split.keep.dimensions,
                                              per.population = per.population)
        ref.df$Population=reference.name

        df = rbind(df, ref.df)
    }

    if (!is.na(per.population) & per.population==1)
        df$label = paste0(round(100*df$value), '%')
    else if (!is.na(per.population) & per.population==100)
        df$label = paste0(round(df$value), '%')
    else
        df$label = df$value

    colors = c(jheem.color, reference.color)
    names(colors) = c(jheem.name, reference.name)


    position = position_dodge(width = 1)

    if (as.line)
        rv = ggplot(df, aes(x=year, y=value, fill=Population, shape=subset, linetype=subset)) +
            geom_line(aes(color=Population), size=line.size) +
            geom_point(size=point.size)
    else
        rv = ggplot(df, aes(x=subset, y=value, fill=Population)) +
            geom_col(position=position) +
            theme(axis.text.x = element_text(angle = 30))

    if (show.labels & !as.line)
        rv = rv + geom_text(aes(label=label),
                            position=position,
                            size=label.size, vjust='bottom')

    if (!is.null(facet.by.dimensions))
        rv = rv + facet_wrap(as.formula(paste0('~',paste0(facet.by.dimensions, collapse='+'))))

    rv
}

#'@export
jplot.incidence <- function(jheem.results,
                            numerators.ref=NULL,
                            denominators.ref=NULL,
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
                            split.by.dimensions=NULL,
                            facet.by.dimensions=NULL,
                            include.hiv.positive.in.denominator=T,
                            per.population=100000,
                            transformation.fn=NULL,
                            split.string=', ',
                            jheem.color=BLUE,
                            reference.color=GREEN,
                            show.labels=F,
                            label.size=5,
                            line.size=2,
                            point.size=5,
                            jheem.name='Simulated',
                            reference.name='Actual',
                            shapes=DEFAULT.SHAPES)
{
    dimensions = clean.keep.dimensions(split.by.dimensions = split.by.dimensions,
                                       facet.by.dimensions = facet.by.dimensions)

    if (is.null(years))
        jheem.years = jheem.results$years
    else
        jheem.years = intersect(as.character(years), as.character(jheem.results$years))

    extracted.jheem = extract.incidence(jheem.results,
                                        years=jheem.years,
                                        ages=ages,
                                        races=races,
                                        subpopulations=subpopulations,
                                        sexes=sexes,
                                        risks=risks,
                                        non.hiv.subsets=non.hiv.subsets,
                                        continuum=continuum,
                                        cd4s=cd4s,
                                        hiv.subsets=hiv.subsets,
                                        keep.dimensions=dimensions$keep.dimensions,
                                        per.population=per.population
                                        )
    jplot.ratios(extracted.jheem,
                 numerators.ref=numerators.ref,
                 denominators.ref=denominators.ref,
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
                 per.population=per.population,
                 split.by.dimensions = split.by.dimensions,
                 facet.by.dimensions = facet.by.dimensions,
                 subset.split.str = split.string,
                 jheem.color=jheem.color,
                 reference.color=reference.color,
                 show.labels=show.labels,
                 label.size=label.size,
                 line.size=line.size,
                 point.size=point.size,
                 jheem.name=jheem.name,
                 reference.name=reference.name,
                 shapes=shapes
    )
}

#'@export
jplot.new.diagnoses <- function(jheem.results,
                                numerators.ref=NULL,
                                denominators.ref=NULL,
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
                                split.by.dimensions=NULL,
                                facet.by.dimensions=NULL,
                                include.hiv.positive.in.denominator=T,
                                per.population=100000,
                                transformation.fn=NULL,
                                split.string=', ',
                                jheem.color=BLUE,
                                reference.color=GREEN,
                                show.labels=F,
                                label.size=5,
                                line.size=2,
                                point.size=5,
                                jheem.name='Simulated',
                                reference.name='Actual',
                                shapes=DEFAULT.SHAPES)
{
    dimensions = clean.keep.dimensions(split.by.dimensions = split.by.dimensions,
                                       facet.by.dimensions = facet.by.dimensions)

    if (is.null(years))
        jheem.years = jheem.results$years
    else
        jheem.years = intersect(as.character(years), as.character(jheem.results$years))

    extracted.jheem = extract.new.diagnoses(jheem.results,
                                            years=jheem.years,
                                            ages=ages,
                                            races=races,
                                            subpopulations=subpopulations,
                                            sexes=sexes,
                                            risks=risks,
                                            non.hiv.subsets=non.hiv.subsets,
                                            continuum=continuum,
                                            cd4s=cd4s,
                                            hiv.subsets=hiv.subsets,
                                            keep.dimensions=dimensions$keep.dimensions,
                                            per.population=per.population,
                                            transformation.fn=transformation.fn
    )

    jplot.ratios(extracted.jheem,
                 numerators.ref=numerators.ref,
                 denominators.ref=denominators.ref,
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
                 per.population=per.population,
                 split.by.dimensions = split.by.dimensions,
                 facet.by.dimensions = facet.by.dimensions,
                 subset.split.str = split.string,
                 jheem.color=jheem.color,
                 reference.color=reference.color,
                 show.labels=show.labels,
                 label.size=label.size,
                 line.size=line.size,
                 point.size=point.size,
                 jheem.name=jheem.name,
                 reference.name=reference.name,
                 shapes=shapes
    )
}

#'@export
jplot.prevalence <- function(jheem.results,
                             numerators.ref=NULL,
                             denominators.ref=NULL,
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
                             split.by.dimensions=NULL,
                             facet.by.dimensions=NULL,
                             include.hiv.positive.in.denominator=T,
                             per.population=100000,
                             transformation.fn=NULL,
                             split.string=', ',
                             jheem.color=BLUE,
                             reference.color=GREEN,
                             show.labels=F,
                             label.size=5,
                             line.size=2,
                             point.size=5,
                             jheem.name='Simulated',
                             reference.name='Actual',
                             shapes=DEFAULT.SHAPES)
{
    dimensions = clean.keep.dimensions(split.by.dimensions = split.by.dimensions,
                                       facet.by.dimensions = facet.by.dimensions)

    if (is.null(years))
        jheem.years = jheem.results$years
    else
        jheem.years = intersect(as.character(years), as.character(jheem.results$years))

    extracted.jheem = extract.prevalence(jheem.results,
                                            years=jheem.years,
                                            ages=ages,
                                            races=races,
                                            subpopulations=subpopulations,
                                            sexes=sexes,
                                            risks=risks,
                                            continuum=continuum,
                                            cd4s=cd4s,
                                            hiv.subsets=hiv.subsets,
                                            keep.dimensions=dimensions$keep.dimensions,
                                            per.population=per.population,
                                            transformation.fn=transformation.fn
    )

    jplot.ratios(extracted.jheem,
                 numerators.ref=numerators.ref,
                 denominators.ref=denominators.ref,
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
                 per.population=per.population,
                 split.by.dimensions = split.by.dimensions,
                 facet.by.dimensions = facet.by.dimensions,
                 subset.split.str = split.string,
                 jheem.color=jheem.color,
                 reference.color=reference.color,
                 show.labels=show.labels,
                 label.size=label.size,
                 line.size=line.size,
                 point.size=point.size,
                 jheem.name=jheem.name,
                 reference.name=reference.name,
                 shapes=shapes
    )
}

#'@export
jplot.hiv.specific.mortality <- function(jheem.results,
                             numerators.ref=NULL,
                             denominators.ref=NULL,
                             years=NULL,
                             ages=NULL,
                             races=NULL,
                             subpopulations=NULL,
                             sexes=NULL,
                             risks=NULL,
                             continuum=NULL,
                             cd4s=NULL,
                             hiv.subsets=NULL,
                             split.by.dimensions=NULL,
                             facet.by.dimensions=NULL,
                             include.hiv.negative.in.denominator=F,
                             per.population=100000,
                             transformation.fn=NULL,
                             split.string=', ',
                             jheem.color=BLUE,
                             reference.color=GREEN,
                             show.labels=F,
                             label.size=5,
                             line.size=2,
                             point.size=5,
                             jheem.name='Simulated',
                             reference.name='Actual',
                             shapes=DEFAULT.SHAPES)
{
    dimensions = clean.keep.dimensions(split.by.dimensions = split.by.dimensions,
                                       facet.by.dimensions = facet.by.dimensions)

    if (is.null(years))
        jheem.years = jheem.results$years
    else
        jheem.years = intersect(as.character(years), as.character(jheem.results$years))

    extracted.jheem = extract.hiv.specific.mortality(jheem.results,
                                                     years=jheem.years,
                                                     ages=ages,
                                                     races=races,
                                                     subpopulations=subpopulations,
                                                     sexes=sexes,
                                                     risks=risks,
                                                     continuum=continuum,
                                                     cd4s=cd4s,
                                                     hiv.subsets=hiv.subsets,
                                                     include.hiv.negative.in.denominator=include.hiv.negative.in.denominator,
                                                     keep.dimensions=dimensions$keep.dimensions,
                                                     per.population=per.population,
                                                     transformation.fn=transformation.fn
    )

    jplot.ratios(extracted.jheem,
                 numerators.ref=numerators.ref,
                 denominators.ref=denominators.ref,
                 years=years,
                 ages=ages,
                 races=races,
                 subpopulations=subpopulations,
                 sexes=sexes,
                 risks=risks,
                 continuum=continuum,
                 cd4s=cd4s,
                 hiv.subsets=hiv.subsets,
                 per.population=per.population,
                 split.by.dimensions = split.by.dimensions,
                 facet.by.dimensions = facet.by.dimensions,
                 subset.split.str = split.string,
                 jheem.color=jheem.color,
                 reference.color=reference.color,
                 show.labels=show.labels,
                 label.size=label.size,
                 line.size=line.size,
                 point.size=point.size,
                 jheem.name=jheem.name,
                 reference.name=reference.name,
                 shapes=shapes
    )
}

#'@export
jplot.overall.hiv.mortality <- function(jheem.results,
                                         numerators.ref=NULL,
                                         denominators.ref=NULL,
                                         years=NULL,
                                         ages=NULL,
                                         races=NULL,
                                         subpopulations=NULL,
                                         sexes=NULL,
                                         risks=NULL,
                                         continuum=NULL,
                                         cd4s=NULL,
                                         hiv.subsets=NULL,
                                         split.by.dimensions=NULL,
                                         facet.by.dimensions=NULL,
                                         include.hiv.negative.in.denominator=F,
                                         per.population=100000,
                                         transformation.fn=NULL,
                                         split.string=', ',
                                         jheem.color=BLUE,
                                         reference.color=GREEN,
                                         show.labels=F,
                                         label.size=5,
                                         line.size=2,
                                         point.size=5,
                                         jheem.name='Simulated',
                                         reference.name='Actual',
                                         shapes=DEFAULT.SHAPES)
{
    dimensions = clean.keep.dimensions(split.by.dimensions = split.by.dimensions,
                                       facet.by.dimensions = facet.by.dimensions)

    if (is.null(years))
        jheem.years = jheem.results$years
    else
        jheem.years = intersect(as.character(years), as.character(jheem.results$years))

    extracted.jheem = extract.overall.hiv.mortality(jheem.results,
                                                     years=jheem.years,
                                                     ages=ages,
                                                     races=races,
                                                     subpopulations=subpopulations,
                                                     sexes=sexes,
                                                     risks=risks,
                                                     continuum=continuum,
                                                     cd4s=cd4s,
                                                     hiv.subsets=hiv.subsets,
                                                     include.hiv.negative.in.denominator=include.hiv.negative.in.denominator,
                                                     keep.dimensions=dimensions$keep.dimensions,
                                                     per.population=per.population,
                                                     transformation.fn=transformation.fn
    )

    jplot.ratios(extracted.jheem,
                 numerators.ref=numerators.ref,
                 denominators.ref=denominators.ref,
                 years=years,
                 ages=ages,
                 races=races,
                 subpopulations=subpopulations,
                 sexes=sexes,
                 risks=risks,
                 continuum=continuum,
                 cd4s=cd4s,
                 hiv.subsets=hiv.subsets,
                 per.population=per.population,
                 split.by.dimensions = split.by.dimensions,
                 facet.by.dimensions = facet.by.dimensions,
                 subset.split.str = split.string,
                 jheem.color=jheem.color,
                 reference.color=reference.color,
                 show.labels=show.labels,
                 label.size=label.size,
                 line.size=line.size,
                 point.size=point.size,
                 jheem.name=jheem.name,
                 reference.name=reference.name,
                 shapes=shapes
    )
}


#The big helper for jplots (besides proportions)
jplot.ratios <- function(extracted.jheem,
                         numerators.ref,
                         denominators.ref,
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
                         split.by.dimensions=NULL,
                         facet.by.dimensions=NULL,
                         include.hiv.positive.in.denominator=T,
                         per.population=100000,
                         subset.split.str=', ',
                         jheem.color=BLUE,
                         reference.color=GREEN,
                         show.labels=F,
                         label.size=5,
                         line.size=2,
                         point.size=5,
                         jheem.name='Simulated',
                         reference.name='Actual',
                         shapes=DEFAULT.SHAPES)
{
    dimensions = clean.keep.dimensions(split.by.dimensions = split.by.dimensions,
                                       facet.by.dimensions = facet.by.dimensions)

    df = NULL
    if (!is.null(extracted.jheem))
    {
        df = convert.ratio.arr.to.df(extracted.jheem,
                                     split.by.dimensions = dimensions$split.by.dimensions,
                                     subset.split.str = subset.split.str)
        df$Population = jheem.name
    }

    if (!is.null(numerators.ref))
    {
        if (is.null(years))
            years.ref = dimnames(numerators.ref)[['year']]
        else
            years.ref = intersect(as.character(years), dimnames(numerators.ref)[['year']])
        if (!is.null(denominators.ref))
            years.ref = intersect(years.ref, dimnames(denominators.ref)[['year']])
        extracted.ref = do.extract.results(numerators=numerators.ref,
                                           denominators.1=denominators.ref,
                                           denominators.2=NULL,
                                           years=years.ref,
                                           ages=ages,
                                           races=races,
                                           subpopulations=subpopulations,
                                           sexes=sexes,
                                           risks=risks,
                                           non.hiv.subsets=non.hiv.subsets,
                                           continuum=continuum,
                                           cd4s=cd4s,
                                           hiv.subsets=hiv.subsets,
                                           keep.dimensions=dimensions$keep.dimensions,
                                           per.population=per.population)

        ref.df = convert.ratio.arr.to.df(extracted.ref,
                                         split.by.dimensions = dimensions$split.by.dimensions,
                                         subset.split.str = subset.split.str)
        ref.df$Population = reference.name

        df = rbind(df,
                   ref.df)
    }

    colors = c(jheem.color, reference.color)
    names(colors) = c(jheem.name, reference.name)
    subsets = unique(df$subset)
    shapes = shapes[1:length(subsets)]
    names(shapes) = subsets

    rv = ggplot(df, aes(x=year, y=value, fill=Population, shape=subset, linetype=subset)) +
        geom_line(aes(color=Population), size=line.size) +
        geom_point(size=point.size)

    if (!is.null(facet.by.dimensions))
        rv = rv + facet_wrap(as.formula(paste0('~',paste0(facet.by.dimensions, collapse='+'))))

    rv = rv + scale_fill_manual(values=colors) +
        scale_color_manual(values=colors) +
        scale_shape_manual(values=shapes)

    rv

}

##-------------##
##-- HELPERS --##
##-------------##

prepare.population.subset.df <- function(arr,
                                         years,
                                         split.by.dimensions,
                                         facet.by.dimensions,
                                         per.population=1,
                                         subset.split.str=', ',
                                         transformation.fn=NULL)
{
    if (class(arr)=='jheem.results')
        extracted = extract.population.subset(arr,
                                              years=years,
                                              keep.dimensions = c(split.by.dimensions, facet.by.dimensions),
                                              denominator.dimensions = facet.by.dimensions,
                                              per.population = per.population,
                                              transformation.fn=transformation.fn)
    else
        extracted = do.extract.population.subset(population.1 = arr,
                                                 years=years,
                                                 keep.dimensions = c(split.by.dimensions, facet.by.dimensions),
                                                 denominator.dimensions = facet.by.dimensions,
                                                 per.population = per.population,
                                                 transformation.fn=transformation.fn)

    df = melt(extracted)

    if (!is.null(split.by.dimensions))
    {
        df$subset = df[,split.by.dimensions[1]]
        if (length(split.by.dimensions)>1)
        {
            for (i in 2:length(split.by.dimensions))
                df$subset = paste0(df$subset, subset.split.str, df[,split.by.dimensions[i]])
        }
    }
    else
        df$subset = 'All'

    df
}

convert.ratio.arr.to.df <- function(extracted,
                                    split.by.dimensions,
                                    subset.split.str=', ')
{
    df = melt(extracted)

    if (!is.null(split.by.dimensions))
    {
        df$subset = df[,split.by.dimensions[1]]
        if (length(split.by.dimensions)>1)
        {
            for (i in 2:length(split.by.dimensions))
                df$subset = paste0(df$subset, subset.split.str, df[,split.by.dimensions[i]])
        }
    }
    else
        df$subset = 'All'

    df
}

clean.keep.dimensions <- function(split.by.dimensions,
                                  facet.by.dimensions,
                                  exclude.year=T)
{
    if (exclude.year)
    {
        split.by.dimensions = setdiff(split.by.dimensions, 'year')
        facet.by.dimensions = setdiff(facet.by.dimensions, 'year')
    }

    facet.by.dimensions = setdiff(facet.by.dimensions, split.by.dimensions)
    keep.dimensions = c(split.by.dimensions, facet.by.dimensions)

    if (exclude.year)
        keep.dimensions = c('year', keep.dimensions)

    list(split.by.dimensions=split.by.dimensions,
         facet.by.dimensions=facet.by.dimensions,
         keep.dimensions=keep.dimensions)
}
