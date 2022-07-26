#* Potential optimizations
#
# - block or sparse transition matrices
# - sum out dimensions for fertility or birth proportions before multiplying


##-- OVERVIEW OF INDEXING
##   We're going to split the population by HIV+ and HIV-
##   Both subpopulations will be indexed by age, race, subpopulation, sex, and risk, first
##   HIV+ are further indexed by continuum of care status, CD4 count, and hiv subset
##   HIV- are further indexed by non-hiv subset

#TOLERANCE - below this value, we treat as zero (to avoid extreme small number errors)
TOLERANCE = 1e-11


model.dx.Rcpp <- function(time, y, parameters)
{
    ##--------------------------------------##
    ##-- IF WE ARE TAKING TOO LONG, ABORT --##
    ##--------------------------------------##

    if ((as.numeric(Sys.time()) - as.numeric(parameters$START_TIME)) > parameters$MAX_RUN_TIME)
    {
        return(rep(1,length(y)))
    }

    ##-------------------------------------------------##
    ##-- CUT ANY NUMBERS BELOW THE TOLERANCE TO ZERO --##
    ##-------------------------------------------------##

    y[y<TOLERANCE] = 0


    ##---------------------------------------##
    ##-- CALCULATE TIME-VARYING PARAMETERS --##
    ##---------------------------------------##

    ## In general, for each parameter, we are going to find the closest anchor time before
    ## and the closest anchor time after, and linearly interpolate between the two
    ## If there is no anchor time before, then the parameter is just the first value
    ## If there is no anchor time after, then the parameter is just the last value

    #The pp list holds the values of all the time-varying parameters

    just.one.param.value = sapply(parameters$time.varying.parameters, function(param.set){
        length(param.set$values)==1
    })

    pp.one.value = lapply(parameters$time.varying.parameters[just.one.param.value], function(param.set){
        param.set$values[[1]]
    })
    names(pp.one.value) = names(parameters$time.varying.parameters)[just.one.param.value]

    pp.interpolated = pull_time_varying_parameters(parameters$time.varying.parameters[!just.one.param.value], time)
    names(pp.interpolated) = names(parameters$time.varying.parameters)[!just.one.param.value]


#    pp.interpolated = lapply(parameters$time.varying.parameters[!just.one.param.value], function(param.set){
#        n.times = length(param.set$times)
#        if (n.times==1)
#            return (param.set$values[[1]])

#        index.before = (1:n.times)[param.set$times <= time]
#        index.before = index.before[length(index.before)]
#        index.after = (1:n.times)[param.set$times > time][1]

#        if (length(index.before)==0)
#           param.set$values[[index.after]]
#        else if (param.set$stepwise || is.na(index.after) || is.infinite(param.set$times[index.before]))
#            param.set$values[[index.before]]
#        else if (is.infinite(param.set$times[index.after]))
#            param.set$values[[index.after]]
#        else
#            (param.set$values[[index.before]] * (param.set$times[index.after] - time) +
#                 param.set$values[[index.after]] * (time - param.set$times[index.before])) /
#            (param.set$times[index.after] - param.set$times[index.before])
#    })

    pp = c(pp.one.value, pp.interpolated)

    ##-----------------------------------##
    ##-- PASS THE CALL THROUGH TO RCPP --##
    ##-----------------------------------##

    #The constant for fixed strata size
    FIXED_SIZE_BY_STRATA = 1

    if (is.null(pp$BIRTH_PROPORTIONS_HIV_TO_HIV))
        pp$BIRTH_PROPORTIONS_HIV_TO_HIV = 0

    delta = dx_function(hiv_positive = y[1:parameters$NUM_HIV_STATES],
                        hiv_negative = y[parameters$NUM_HIV_STATES + 1:parameters$NUM_NONHIV_STATES],
                        DIMS_HIV = parameters$DIMS_HIV,
                        DIMS_NONHIV = parameters$DIMS_NONHIV,
                        AGING_RATES_HIV = pp$AGING_RATES_HIV,
                        AGING_RATES_NONHIV = pp$AGING_RATES_NONHIV,
                        FERTILITY_RATES_HIV = pp$FERTILITY_RATES_HIV,
                        FERTILITY_RATES_NONHIV = pp$FERTILITY_RATES_NONHIV,
                        FRACTION_BIRTHS_INFECTED = pp$FRACTION_BIRTHS_INFECTED,
                        BIRTHS_DIMMASK_NONHIV=parameters$BIRTHS.DIMMASK.NONHIV,
                        BIRTH_PROPORTIONS_NONHIV = pp$BIRTH_PROPORTIONS_NONHIV,
                        BIRTHS_DIMMASK_HIV_TO_NONHIV = parameters$BIRTHS.DIMMASK.HIV.TO.NONHIV,
                        BIRTH_PROPORTIONS_HIV_TO_NONHIV = pp$BIRTH_PROPORTIONS_HIV_TO_NONHIV,
                        BIRTHS_DIMMASK_HIV_TO_HIV = parameters$BIRTHS.DIMMASK.HIV.TO.HIV,
                        BIRTH_PROPORTIONS_HIV_TO_HIV = pp$BIRTH_PROPORTIONS_HIV_TO_HIV,
                        NUM_TRANSMISSION_ROUTES = parameters$NUM_TRANSMISSION_ROUTES,
                        CONTACT_MATRICES = pp[sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){paste0('CONTACT_', i)})],
                        TRANSMISSIBILITIES = pp[sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){paste0('TRANSMISSIBILITY_', i)})],
                        SUSCEPTIBILITIES = pp[sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){paste0('SUSCEPTIBILITY_', i)})],
                        GLOBAL_TRANSMISSION_RATES = sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){pp[[paste0('GLOBAL_TRANSMISSION_RATE_', i)]]}),
                        CONTACT_DIMMASKS = parameters$TRANSMISSION_ROUTE_DIMMASKS,
                        NEW_INFECTION_PROPORTIONS = pp$NEW_INFECTION_PROPORTIONS,
                        COLLAPSE_INCIDENCE_INDICES = parameters$collapse.incidence.indices-1,
                        HIV_TRANSITION_ARRAYS = pp[parameters$transition.mappings$hiv.array.names],
                        HIV_TRANSITION_DIMENSIONS = parameters$transition.mappings$dimensions.for.hiv.arrays,
                        NONHIV_TRANSITION_ARRAYS = pp[parameters$transition.mappings$nonhiv.array.names],
                        NONHIV_TRANSITION_DIMENSIONS = parameters$transition.mappings$dimensions.for.nonhiv.arrays,
                        TRACKED_HIV_TRANSITION_INDICES = parameters$TRACKED_HIV_TRANSITION_INDICES-1,
                        TRACKED_HIV_TRANSITION_FROM = parameters$TRACKED_HIV_TRANSITIONS_FROM-1,
                        TRACKED_HIV_TRANSITION_ARRAY_INDICES = parameters$transition.mappings$tracked.hiv.array.indices,
                        TRACKED_HIV_TRANSITION_COLLAPSE_INDICES = parameters$TRACKED_HIV_TRANSITION_COLLAPSE_INDICES-1,
                        TRACKED_NONHIV_TRANSITION_INDICES = parameters$TRACKED_NONHIV_TRANSITION_INDICES-1,
                        TRACKED_NONHIV_TRANSITION_FROM = parameters$TRACKED_NONHIV_TRANSITIONS_FROM-1,
                        TRACKED_NONHIV_TRANSITION_ARRAY_INDICES = parameters$transition.mappings$tracked.nonhiv.array.indices,
                        TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES = parameters$TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES-1,
                        HIV_SPECIFIC_MORTALITY_RATES = pp$HIV_SPECIFIC_MORTALITY,
                        GENERAL_MORTALITY_RATES_HIV = pp$GENERAL_MORTALITY_FOR_HIV_POSITIVE,
                        GENERAL_MORTALITY_RATES_NONHIV = pp$GENERAL_MORTALITY_FOR_HIV_NEGATIVE,
                        FIXED_SIZE=if (pp$KEEP_STRATA_SIZE_CONSTANT) FIXED_SIZE_BY_STRATA else 0,
                        FIXED_SIZE_STRATA_DIMMASK = parameters$FIXED_SIZE_STRATA_DIMMASK,
                        MODEL_BIRTHS = pp$MODEL_BIRTHS,
                        MODEL_MATERNAL_TRANSMISSION = parameters$MODEL_MATERNAL_TRANSMISSION,
                        TRACK_HIV_SPECIFIC_MORTALITY = parameters$TRACK_HIV_SPECIFIC_MORTALITY,
                        COLLAPSE_DIMMASK_HIV_SPECIFIC_MORTALITY = parameters$COLLAPSE_DIMMASK_HIV_SPECIFIC_MORTALITY,
                        TRACK_HIV_OVERALL_MORTALITY = parameters$TRACK_OVERALL_HIV_MORTALITY,
                        COLLAPSE_DIMMASK_HIV_OVERALL_MORTALITY = parameters$COLLAPSE_DIMMASK_HIV_OVERALL_MORTALITY,
                        TRACK_NONHIV_MORTALITY = parameters$TRACK_OVERALL_NONHIV_MORTALITY,
                        COLLAPSE_DIMMASK_NONHIV_MORTALITY = parameters$COLLAPSE_DIMMASK_NONHIV_MORTALITY)


    # Return
    # list(delta)
    delta
}

#'@title Run a JHEEM Model
#'
#'@param jheem A fully initialized and set up jheem object
#'@param start.year,end.year The years between which to run the model (inclusive)
#'@param keep.years Which years to save results from
#'@param prior.run.results An object of class jheem.results, if starting values for this run should be taken from a prior run. If NULL, starting values will be taken from the initial state specified by the JHEEM object
#'@param verbose Whether to print detailed status updates
#'
#'@return An object of class jheem.results
#'
#'@export
run.jheem <- function(jheem,
                      start.year,
                      end.year,
                      keep.years = start.year:end.year,
                      prior.run.results=NULL,
                      verbose=T,
                      print.warnings=T,
                      max.run.time.seconds=Inf,
                      atol=1e-06,
                      rtol=1e-06
                      )
{
    #-- Check years --#
    if (!(end.year > start.year))
        stop(paste0("The end year (", end.year, ") must be AFTER the start year (", start.year, ")"))
    if (is.null(prior.run.results) && length(setdiff(keep.years, start.year:end.year))>0)
        stop(paste0("The keep.years must all fall between ", start.year, " and ", end.year))
    else if (!is.null(prior.run.results) && length(setdiff(setdiff(keep.years, start.year:end.year), prior.run.results$years))>0)
        stop(paste0("The keep.years must all either be contained in prior.run.results or fall between ", start.year, " and ", end.year))

    #-- Pull initial state components --#
    # either from a prior run or from the JHEEM itself (if no prior run results)
    if (is.null(prior.run.results))
    {
        #Check initial states in JHEEM
        if (is.null(jheem$initial.state.hiv.positive) && is.null(jheem$initial.state.hiv.negative))
            stop("The initial states for HIV-negative and HIV-positive populations have not been set for the JHEEM")
        else if (is.null(jheem$initial.state.hiv.negative))
            stop("The initial state for the HIV-negative populations has not been set for the JHEEM")
        else if (is.null(jheem$initial.state.hiv.positive))
            stop("The initial state for the HIV-positive populations has not been set for the JHEEM")

        #Load initial states
        hiv.negative.population = jheem$initial.state.hiv.negative
        hiv.positive.population = jheem$initial.state.hiv.positive
    }
    else
    {
        # make sure the prior run includes the new start year
        if (!any(prior.run.results$years==(start.year-1)))
            stop(paste0("In order to start the new model at ", start.year,
                        " based off of previous run results, the previous run must include the year ",
                        (start.year-1), " (ie, the year leading into ", start.year, "), but it does not"))

        # check dimension compatibility

        #Load initial states
        hiv.negative.population = prior.run.results$hiv.negative[as.character(start.year-1),,,,,,]
        hiv.positive.population = prior.run.results$hiv.positive[as.character(start.year-1),,,,,,,,]
    }

    #-- Get the parameters --#
    parameters = get.jheem.parameters(jheem,
                                      years=start.year:(end.year+1),
                                      initial.hiv.negative = hiv.negative.population,
                                      initial.hiv.positive = hiv.positive.population,
                                      verbose=verbose,
                                      print.warnings=print.warnings,
                                      use.rcpp=use.rcpp)


    #-- Set up the cumulative states --##
    #   Since these are counts that are subtracted one year to the next
    #   (and since we always throw away the first year of our model run)
    #   these can all be set to zero

    #Load (empty) incident cases
#    new.cases = do.get.population.skeleton(get.dimnames.all(jheem), 0)
    new.cases = numeric(parameters$N.INCIDENT.CASES)

    #Load (empty) mortality
    if (parameters$TRACK_HIV_SPECIFIC_MORTALITY)
    {
        dim.names = get.dimnames.by.name(jheem, jheem$parameters$hiv.specific.mortality.keep.dimensions)
        hiv.specific.mortality = numeric(prod(sapply(dim.names, length)))
    }

    if (parameters$TRACK_OVERALL_HIV_MORTALITY)
    {
        dim.names = get.dimnames.by.name(jheem, jheem$parameters$overall.hiv.mortality.keep.dimensions)
        hiv.positive.mortality = numeric(prod(sapply(dim.names, length)))
    }

    if (parameters$TRACK_OVERALL_NONHIV_MORTALITY)
    {
        dim.names = get.dimnames.by.name(jheem, jheem$parameters$overall.nonhiv.mortality.keep.dimensions)
        hiv.negative.mortality = numeric(prod(sapply(dim.names, length)))
    }

    #Load (empty) tracked transitions
    tracked.transitions.hiv.positive = rep(0, jheem$parameters$N_COLLAPSED_TRACKED_HIV_TRANSITIONS)
    tracked.transitions.hiv.negative = rep(0, jheem$parameters$N_COLLAPSED_TRACKED_NONHIV_TRANSITIONS)


    #-- Flatten the initial state --#
    mortality = numeric()
    if (parameters$TRACK_HIV_SPECIFIC_MORTALITY)
        mortality = c(mortality, hiv.specific.mortality)
    if (parameters$TRACK_OVERALL_HIV_MORTALITY)
        mortality = c(mortality, hiv.positive.mortality)
    if (parameters$TRACK_OVERALL_NONHIV_MORTALITY)
        mortality = c(mortality, hiv.negative.mortality)

    init = c(as.numeric(hiv.positive.population),
             as.numeric(hiv.negative.population),
             as.numeric(new.cases),
             mortality,
             tracked.transitions.hiv.positive,
             tracked.transitions.hiv.negative)

    #-- Run ODE solver --#
    if (verbose)
        print(paste0("Running model from ", start.year, " to ", end.year, "..."))

    parameters$MAX_RUN_TIME = max.run.time.seconds
    start.time = parameters$START_TIME = Sys.time()
#    ode.results = deSolve::ode(y=init,
#                               times=start.year:(end.year+1),
#                               func=fn,
#                               parms=parameters)

#    ode.results = diffeqr::ode.solve(f=fn,
#                                     u0=init,
#                                     p=parameters,
#                                     tspan = list(start.year, end.year+1),
#                                     saveat = start.year:(end.year+1))

    ode.results = odeintr::integrate_sys(sys=function(x,t){model.dx.Rcpp(time=t,y=x,parameters=parameters)},
                                init=init,
                                start=start.year,
                                duration=end.year+1-start.year,
                                step_size=1, atol = atol, rtol = rtol)
    end.time = Sys.time()

    #see if the check bit was set
#    terminated = any(ode.results[,dim(ode.results)[2]]) != 0
    terminated = (as.numeric(end.time) - as.numeric(start.time)) > max.run.time.seconds

    if (verbose)
    {
        if (terminated)
            print("Model was terminated for long run time. Runtime was:")
        else
            print("Done running model. Runtime was:")

        print(end.time-start.time)
    }

    #-- Hydrate the results and return --#

    if (terminated)
    {
        rv = list(terminated=T)
    }
    else
    {
        if (is.null(prior.run.results))
            years.to.keep.from.this.run = keep.years
        else
            years.to.keep.from.this.run = intersect(keep.years, start.year:end.year)

        rv = hydrate.ode.results(ode.results,
                                 jheem=jheem,
                                 parameters=parameters,
                                 keep.years=years.to.keep.from.this.run,
                                 result.year.minus.ode.year = -1)

        if (!is.null(prior.run.results) && length(setdiff(keep.years, start.year:end.year))>0)
        {
            if (verbose)
                print('Combining results with prior run...')
            rv = combine.jheem.results(prior.run.results, rv)
            rv = subset.jheem.results(rv, keep.years)
        }

        rv$terminated = F
    }
    rv$times = list(start=start.time,
                    end=end.time,
                    run=end.time-start.time)

    rv
}

#A helper to check to make sure that all needed parameters are set for the JHEEM run
# Will set some default parameters if parameters with defaults are missing

#exporting for now to help testing new transitions
#'@export
get.jheem.parameters <- function(jheem,
                                 years,
                                 initial.hiv.negative,
                                 initial.hiv.positive,
                                 verbose=T,
                                 print.warnings=T,
                                 use.rcpp)
{
    # Check the parameters one by one...

    #-- Births vs. Fixed Strata Sizes --#
    # At least one of fix strata sizes or model births must be true (or both) at all times
    # If any fix strata sizes is true, we must have strata sizes set
    # If any model births is true, we must have fertility rates and birth proportions set

    # Set MODEL_BIRTHS and FIX_STRATA_SIZES if never set
    if (is.null(jheem$parameters$time.varying.parameters$MODEL_BIRTHS) &&
        is.null(jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES))
    {
        # if both fixed size and fertility rates have been set, throw an error
        # if fixed population size but not fertility rates have been set, fix population
        # if fertility rates but not fixed population size have been set, model births
        # if neither has been set, throw error
        if (!is.null(jheem$parameters$NUM_FIXED_STRATA) &&
            !is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_NONHIV))
        {
            stop("You have not set whether to fix strata sizes of model births explicitly. Information has been provided to the JHEEM about both fixed population sizes and fertility rates. Please specify when to use fixed sizes vs. model births using 'set.use.fixed.population")
        }
        else if (!is.null(jheem$parameters$NUM_FIXED_STRATA))
        {
#            print("You have not set whether to fix strata sizes or model births explicitly. Since you have set what sizes population strata should be fixed to (and have not set fertility rates), the JHEEM will be set to always fix population strata sizes by default. To set this manually, use 'set.use.fixed.population")
            jheem = set.use.fixed.population(jheem, T)
        }
        else if (!is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_NONHIV))
        {
#            print("You have not set whether to fix strata sizes or model births explicitly. Since you have set fertility rates (and have not set what sizes population strata should be fixed to), the JHEEM will be set to always fix population strata sizes by default. To set this manually, use 'set.use.fixed.population")
            jheem = set.use.fixed.population(jheem, F)
        }
        else
        {
            stop("No information has been provided to the JHEEM about either fixing population strata or modeling births. One or the other must be set, using 'set.fertility' or the set of functions around 'set.fixed.population.strata")
        }
    }
    else if (is.null(jheem$parameters$time.varying.parameters$MODEL_BIRTHS))
    {
        for (i in 1:length(jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES$values))
            jheem = set.time.varying.parameter.value(jheem,
                                                     parameter.name = 'MODEL_BIRTHS',
                                                     parameter.value = !jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES$values[[i]],
                                                     time = jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES$time[[i]],
                                                     stepwise=T)
    }
    else if (is.null(jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES))
    {
        for (i in 1:length(jheem$parameters$time.varying.parameters$MODEL_BIRTHS$values))
            jheem = set.time.varying.parameter.value(jheem,
                                                     parameter.name = 'FIX_STRATA_SIZES',
                                                     parameter.value = !jheem$parameters$time.varying.parameters$MODEL_BIRTHS$values[[i]],
                                                     time = jheem$parameters$time.varying.parameters$MODEL_BIRTHS$time[i],
                                                     stepwise=T)
    }

    #Check aging rates
    if (is.null(jheem$parameters$time.varying.parameters$AGING_RATES_HIV))
    {
        if (verbose)
            print("No aging rates set for HIV-positive. Assuming aging rates reflect a constant distribution of age within age strata")
        jheem = set.aging.hiv.positive(jheem, get.aging.skeleton.hiv.positive(jheem))
    }

    if (is.null(jheem$parameters$time.varying.parameters$AGING_RATES_NONHIV))
    {
        if (verbose)
            print("No aging rates set for HIV-negative. Assuming aging rates reflect a constant distribution of age within age strata")
        jheem = set.aging.hiv.negative(jheem, get.aging.skeleton.hiv.negative(jheem))
    }

    #Check constant strata sizes
    if (is.null(jheem$parameters$time.varying.parameters$KEEP_STRATA_SIZE_CONSTANT) ||
        !any(unlist(jheem$parameters$time.varying.parameters$KEEP_STRATA_SIZE_CONSTANT$values)))
    {
        if (!is.null(jheem$parameters$CONSTANT_SIZE_DIMENSIONS) &&
            (verbose || print.warnings))
            print("WARNING: You have set dimensions by which strata sizes are kept constant, but the JHEEM is not actually set to keep strata sizes constant at any time. Use 'set.keep.strata.sizes.constant' to indicate when the JHEEM should keep strata sizes constant.")
    }
    else
    {
        if (is.null(jheem$parameters$CONSTANT_SIZE_DIMENSIONS))
            stop("You have set the JHEEM to keep strata sizes constant, but you have not specified which dimensions define the strata to keep constant. Use 'set.fixed.size.strata' to do so.")

        if (any(!unlist(jheem$parameters$time.varying.parameters$MODEL_BIRTHS$values)) &&
            (verbose || print.warnings))
            print("WARNING: You have set the JHEEM to keep strata sizes constant for at least some time, and you have set the JHEEM to NOT model births for at least some time. Keeping strata sizes constant WITHOUT modeling births at the time will cause erroneous behavior.")

        if (any(unlist(jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES$values)) &&
            (verbose || print.warnings))
            print("WARNING: You have set the JHEEM to keep strata sizes constant for at least some time, and you have set the JHEEM to fix strata sizes for at least some time. Fixing strata sizes and keeping them constant will cause unpredictable behavior.")
    }

    if (is.null(jheem$parameters$time.varying.parameters$KEEP_STRATA_SIZE_CONSTANT))
        jheem = set.time.varying.parameter.value(jheem,
                                                 parameter.name = 'KEEP_STRATA_SIZE_CONSTANT',
                                                 parameter.value = F,
                                                 time = -Inf,
                                                 stepwise = T)


    #Check to make sure we have not plugged in model births info but are not set to model births
    if (!any(unlist(jheem$parameters$time.varying.parameters$MODEL_BIRTHS$values)))
    {
        if (((!is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_NONHIV) || !is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_HIV)) &&
            (!is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_NONHIV) || !is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_NONHIV) || !is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_HIV))) &&
            (verbose || print.warnings))
            print("WARNING: You have set fertility rates and birth proportions, but the JHEEM is not set to actually model births at any time point (ie, the fertility rates and birth proportions will never be used)")
        else if ((!is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_NONHIV) || !is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_HIV))  &&
            (verbose || print.warnings))
            print("WARNING: You have set fertility rates, but the JHEEM is not set to actually model births at any time point (ie, the fertility rates will never be used)")
        else if ((!is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_NONHIV) ||
                 !is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_NONHIV) ||
                 !is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_HIV)) &&
                 (verbose || print.warnings))
            print("WARNING: You have set birth proportions, but the JHEEM is not set to actually model births at any time point (ie, the birth proportions will never be used)")
    }


    #Check to make sure we have not plugged in fixed strata size info but are not set to fix strata sizes
    if (!any(unlist(jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES$values)))
    {
        if (!is.null(jheem$parameters$time.varying.parameters$NUM_FIXED_STRATA) &&
            (verbose || print.warnings))
            print("WARNING: You have set information on how to fix population strata sizes, but the JHEEM is not set to actually use fixed strata at any time point (ie, the fixed population strata settings will never be used)")
    }

    #Check required parameters to model births (FERTILITY_RATES and BIRTH_PROPORTIONS)
    if (any(unlist(jheem$parameters$time.varying.parameters$MODEL_BIRTHS$values)))
    {
        if (is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_HIV) && is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_NONHIV))
            stop("Fertility rates for HIV-positive and HIV-negative populations have not been set, but the JHEEM is set to model births. The JHEEM cannot model births without fertility rates. Use 'set.fertility' to set these rates")
        else if (is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_HIV))
            stop("Fertility rates for the HIV-positive population have not been set, but the JHEEM is set to model births. The JHEEM cannot model births without fertility rates. Use 'set.fertility.hiv.positive' to set these rates")
        else if (is.null(jheem$parameters$time.varying.parameters$FERTILITY_RATES_NONHIV))
            stop("Fertility rates for HIV-positive and HIV-negative populations have not been set, but the JHEEM is set to model births. The JHEEM cannot model births without fertility rates. Use 'set.fertility.hiv.negative' to set these rates")


        if (is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_NONHIV))
        {
            if (is.null(initial.hiv.negative))
                stop("In order to model births, you must either indicate the proportions according to which new births from HIV-negative mothers are distributed across the population using 'set.birth.proportions.hiv.negative', or set the initial HIV-negative population using 'set.initial.population.hiv.negative'")
            else
            {
                if (verbose)
                    print("The proportions according to which new births from HIV-negative mothers are distributed across the population were not set. By default, these will proportions will be set to the proportions across the initial HIV-negative population. Use 'set.birth.proportions.hiv.negative', to set this manually")
                birth.proportions = create.birth.proportions.from.population(jheem, initial.hiv.negative)
                jheem = set.birth.proportions.hiv.negative(jheem, birth.proportions)
            }
        }

        if (is.null(jheem$parameters$time.varying.parameters$FRACTION_BIRTHS_INFECTED))
        {
            if (verbose)
                print("The fraction of HIV-infected births to HIV-positive mothers has not been set. By default, this will be set to zero (ie, no maternal-fetal transmission). To set this manually, use 'set.birth.proportions.hiv.positive'")
            jheem = set.birth.proportions.hiv.positive(jheem, fraction.births.infected=0, proportions.for.uninfected.births = NULL)
        }

        any.maternal.transmission = any(sapply(jheem$parameters$time.varying.parameters$FRACTION_BIRTHS_INFECTED$values, function(fraction){
            any(fraction>0)
        }))
        jheem$parameters$MODEL_MATERNAL_TRANSMISSION = any.maternal.transmission
        if (is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_HIV) && any.maternal.transmission)
                stop("The proportions according to which HIV-positive new births from HIV-positive mothers (ie, maternal-fetal transmission) are distributed across the population were not set. These must be set using 'set.birth.proportions.hiv.positive' when maternal-fetal transmission is modeled.")

        if (!any.maternal.transmission)
        {
            jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_HIV = NULL
            jheem$parameters$time.varying.parameters$FRACTION_BIRTHS_INFECTED = NULL
            jheem = set.time.varying.parameter.value(jheem,
                                                     parameter.name='FRACTION_BIRTHS_INFECTED',
                                                     parameter.value=0,
                                                     times=-Inf)
        }

        if (is.null(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_HIV_TO_NONHIV))
        {
            if (verbose)
                print("The proportions according to which new HIV-negative births from HIV-positive mothers are distributed across the population were not set. By default, these will proportions will be set to match birth proportions from HIV-negative mothers. Use 'set.birth.proportions.hiv.positive', to set this manually")
            for (i in 1:length(jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_NONHIV$values))
            {
                birth.proportions = jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_NONHIV$values[[i]]
                birth.proportions = apply(birth.proportions, c(1:5,7:10), sum)
                jheem = set.birth.proportions.hiv.positive(jheem,
                                                           fraction.births.infected = NULL,
                                                           proportions.for.uninfected.births = birth.proportions,
                                                           time=jheem$parameters$time.varying.parameters$BIRTH_PROPORTIONS_NONHIV$times[i])
            }
        }

    }

    # Check required parameters to fix strata
    if (any(unlist(jheem$parameters$time.varying.parameters$FIX_STRATA_SIZES$values)))
    {
        if (is.null(jheem$parameters$NUM_FIXED_STRATA))
            stop("In order to use fixed population sizes, you must indicate which strata are fixed using 'set.fixed.population.strata'")

        if (is.null(jheem$parameters$time.varying.parameters$TARGET_STRATUM_SIZE))
        {
            if (is.null(initial.hiv.negative) || is.null(initial.hiv.positive))
                stop("In order to use fixed population sizes, you must either indicate the sizes to which to fix the population using 'set.fixed.population.size', or set initial HIV-negative and HIV-positive populations using 'set.initial.populations'")
            else
            {
                if (verbose)
                    print("The sizes to which to fix the population when using fixed population strata was not set. By default, these will be pegged to the initial HIV-negative + HIV-positive poulations. Use 'set.fixed.population.size', to set this manually")
                jheem = set.fixed.population.size(jheem,
                                                  hiv.negative.population = initial.hiv.negative,
                                                  hiv.positive.population = initial.hiv.positive)
            }
        }
    }


    #-- Mortality --#

    if (is.null(jheem$parameters$time.varying.parameters$GENERAL_MORTALITY_FOR_HIV_NEGATIVE))
        stop("General mortality rates for the HIV-negative population have not been set. The JHEEM cannot run without these. Use 'set.general.mortality' or 'set.general.mortality.hiv.negative' to set these rates")

    if (is.null(jheem$parameters$time.varying.parameters$GENERAL_MORTALITY_FOR_HIV_POSITIVE))
        stop("General mortality rates for the HIV-positive population have not been set. The JHEEM cannot run without these. Use 'set.general.mortality' or 'set.general.mortality.hiv.positive' to set these rates")

    if (is.null(jheem$parameters$time.varying.parameters$HIV_SPECIFIC_MORTALITY))
        stop("HIV-specific mortality rates have not been set. The JHEEM cannot run without these. Use 'set.hiv.specific.mortality' to set these rates")


    #-- Immigration and Emigration --#

    if (is.null(jheem$parameters$time.varying.parameters$EMIGRATION_FOR_HIV_NEGATIVE))
    {
        if ((verbose || print.warnings))
            print("WARNING: Emigration rates for the HIV-negative population have not been set. Assuming zero emigration from HIV-negative")
        jheem = set.emigration.hiv.negative(jheem, 0)
    }

    if (is.null(jheem$parameters$time.varying.parameters$EMIGRATION_FOR_HIV_POSITIVE))
    {
        if ((verbose || print.warnings))
            print("WARNING: Emigration rates for the HIV-positive population have not been set. Assuming zero emigration from HIV-positive")
        jheem = set.emigration.hiv.positive(jheem, 0)
    }

    if (is.null(jheem$parameters$time.varying.parameters$IMMIGRATION_RATES))
    {
        if (is.null(jheem$parameters$time.varying.parameters$USE_SET_IMMGRATION_PROPORTIONS) &&
            is.null(jheem$parameters$IMMIGRATION_FROM_DIMENSIONS))
        {
            if ((verbose || print.warnings))
                print("WARNING: Immigration rates have not been set. Assuming zero immigration")
            jheem = set.immigration.rates(jheem, 0)
            jheem = set.immigration.proportions(jheem, use.population.proportions=T)
        }
        else if (!is.null(jheem$parameters$time.varying.parameters$USE_SET_IMMGRATION_PROPORTIONS))
            stop("Immgration proportions have been set, but immigration rates have not been set. Use set.immigration.rates")
        else
            stop("Dimensions for immigration have been indicated, but, but immigration rates and proportions have not been set. Use set.immigration.rates and set.immigration.proportions")
    }
    else if (is.null(jheem$parameters$time.varying.parameters$USE_SET_IMMGRATION_PROPORTIONS))
        stop("Immigration rates have been set, but the proportions according to which immigrants are distributed have not been set. Use set.immigration.proportions")


    #-- Transmission Parameters --#

    for (i in 1:length(jheem$transmission.routes))
    {
        route = jheem$transmission.routes[i]

        #-- Transmissibility, Susceptibility, and Global Transmission Rate --#

        if (is.null(jheem$parameters$time.varying.parameters[[paste0('TRANSMISSIBILITY_', i)]]))
        {
            if (verbose)
                print(paste0("A transmissibility array has not been set for the '", route, "' transmission route. ",
                             "By default, all HIV-positive strata will have (equal) transmissibility set to 1. ",
                             "Use 'set.transmissibility' to set this manually"))

            jheem = set.transmissibility(jheem, get.hiv.positive.population.skeleton(jheem, value=1), transmission.route.names = route)
        }

        if (is.null(jheem$parameters$time.varying.parameters[[paste0('SUSCEPTIBILITY_', i)]]))
        {
            if (verbose)
                print(paste0("A susceptibility array has not been set for the '", route, "' transmission route. ",
                             "By default, all HIV-negative strata will have (equal) susceptibility set to 1. ",
                             "Use 'set.susceptibility' to set this manually"))

            jheem = set.susceptibility(jheem, get.hiv.negative.population.skeleton(jheem, value=1), transmission.route.names = route)
        }

        if (is.null(jheem$parameters$time.varying.parameters[[paste0('GLOBAL_TRANSMISSION_RATE_', i)]]))
        {
            if (verbose)
                print(paste0("A global transmission rate has not been set for the '", route, "' transmission route. ",
                             "By default, this will be set to (a constant) 1. ",
                             "Use 'set.global.transmission.rate' to set this manually"))

            jheem = set.global.transmission.rate(jheem, 1, transmission.route.names = route)
        }


        #-- Contact Matrices --#

        if (is.null(jheem$parameters$time.varying.parameters[[paste0('CONTACT_', i)]]))
        {
            stop(paste0("A transmission contact array has not been set for the '", route, "' transmission route. ",
                        "Every route of transmission must have a contact array set explicity using 'set.transmission.contact.array'"))
        }

    }

    #-- New Infection Proportions --#

    if (is.null(jheem$parameters$time.varying.parameters$NEW_INFECTION_PROPORTIONS))
    {
        if (verbose)
            print(paste0("The proportions by which new infections are distributed across the population has not been set. ",
                         "By default, all new HIV infections will fall into the '",
                         jheem$continuum.of.care[1], "' continuum state, the '",
                         jheem$cd4.strata[1], "' cd4 stratum, and the '",
                         jheem$hiv.subsets[1], "' hiv subset. ",
                         "You can use 'set.new.infection.proportions' to set this manually"))

        default.proportions = get.uniform.new.infection.proportions(jheem,
                                                                   initial.continuum=jheem$continuum.of.care[1],
                                                                   initial.cd4=jheem$cd4.strata[1],
                                                                   initial.hiv.subset=jheem$hiv.subsets[1])
        jheem = set.new.infection.proportions(jheem, default.proportions)
    }


    #-- Whether to Track Mortality --#

    if (is.null(jheem$parameters$TRACK_HIV_SPECIFIC_MORTALITY))
    {
        if (verbose)
            print("Whether to track HIV-specific mortality has not been set. Setting to FALSE by default")
        jheem = set.track.hiv.specific.mortality(jheem, F)
    }

    if (is.null(jheem$parameters$TRACK_OVERALL_HIV_MORTALITY) || is.null(jheem$parameters$TRACK_OVERALL_NONHIV_MORTALITY))
    {
        if (verbose)
            print("Whether to track overall mortality has not been set. Setting to TRUE for HIV-positive and FALSE for hiv-negative by default")
        jheem = set.track.overall.mortality(jheem, track.hiv.positive=T, track.hiv.negative=F)
    }


    #-- Set up parameter versions for Rcpp --#

    # Pull dimensions
    jheem$parameters$DIMS_HIV = sapply(get.dimnames.hiv(jheem), length)
    jheem$parameters$DIMS_NONHIV = sapply(get.dimnames.nonhiv(jheem), length)

    # Pull dimension masks for birth proportions
    jheem$parameters$BIRTHS.DIMMASK.HIV.TO.HIV = rep(F, NUM_HIV_STRATA)
    if (!is.null(jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS$HIV_TO_HIV))
        jheem$parameters$BIRTHS.DIMMASK.HIV.TO.HIV[jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS$HIV_TO_HIV] = T

    jheem$parameters$BIRTHS.DIMMASK.NONHIV = rep(F, NUM_NONHIV_STRATA)
    jheem$parameters$BIRTHS.DIMMASK.NONHIV[jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS$NONHIV] = T

    jheem$parameters$BIRTHS.DIMMASK.HIV.TO.NONHIV = rep(F, NUM_HIV_STRATA)
    jheem$parameters$BIRTHS.DIMMASK.HIV.TO.NONHIV[jheem$parameters$BIRTH_PROPORTIONS_FROM_DIMENSIONS$HIV_TO_NONHIV] = T

    # Set up dimension masks for transmission routes
    jheem$parameters$TRANSMISSION_ROUTE_DIMMASKS = lapply(jheem$parameters$CONTACT_TYPE_FROM_DIMENSIONS, function(from){
        rv = rep(F, NUM_GENERAL_STRATA)
        rv[from] = T
        rv
    })

    jheem$parameters$FIXED_SIZE_STRATA_DIMMASK = rep(F,5)
    if (!is.null(jheem$parameters$CONSTANT_SIZE_DIMENSIONS))
        jheem$parameters$FIXED_SIZE_STRATA_DIMMASK[jheem$parameters$CONSTANT_SIZE_DIMENSIONS] = T


    # Set up dimension masks for collapsing mortality
    jheem$parameters$COLLAPSE_DIMMASK_HIV_SPECIFIC_MORTALITY = sapply(names(get.dimnames.hiv(jheem)), function(name){
        any(name == jheem$parameters$hiv.specific.mortality.keep.dimensions)
    })
    jheem$parameters$COLLAPSE_DIMMASK_HIV_OVERALL_MORTALITY = sapply(names(get.dimnames.hiv(jheem)), function(name){
        any(name == jheem$parameters$overall.hiv.mortality.keep.dimensions)
    })
    jheem$parameters$COLLAPSE_DIMMASK_NONHIV_MORTALITY = sapply(names(get.dimnames.nonhiv(jheem)), function(name){
        any(name == jheem$parameters$overall.nonhiv.mortality.keep.dimensions)
    })

    #-- Set up Dimension Mappings --#

    #Set up mappings for transition dimensions
    jheem = prepare.transition.mappings(jheem)

    #Set up collapse dimensions for incidence tracking
    jheem$parameters$collapse.incidence.indices = get.collapse.indices(jheem,
                                                                       full.dimension.names=names(get.dimnames.all(jheem)),
                                                                       collapsed.dimension.names=jheem$parameters$incidence.keep.dimensions)

    jheem$parameters$N.INCIDENT.CASES = max(jheem$parameters$collapse.incidence.indices)

    #-- Check for NA's --#
    check.na.parameters(jheem$parameters)

    #-- Return it --#
    jheem$parameters
}

check.na.parameters <- function(parameters)
{
    #Check time-varying parameters
    sapply(names(parameters$time.varying.parameters), function(name){

        na.times = is.na(parameters$time.varying.parameters)
        if (any(na.times))
            stop(paste0("The ", get.character.list(get.ordinal((1:length(na.times))[na.times])),
                        " time", ifelse(sum(na.times)==1, '', 's'),
                        " for time-varying parameter '", name, "' are NA"
            ))

        n.val = length(parameters$time.varying.parameters[[name]]$values)
        for (i in 1:n.val)
        {
            val = parameters$time.varying.parameters[[name]]$values[[i]]
            is.multiple = class(val)=='list' || is.null(length(val)) || length(val) > 1

            if (any(is.na(val)))
                stop("The ",
                     ifelse(n.val==1, 'value',
                            paste0(get.ordinal(i), " value (at time=",
                                    parameters$time.varying.parameters[[name]]$times[i],
                                    ")")),
                     " for time-varying parameter '", name,
                     "' ",
                     ifelse(is.multiple, "has components that are", "is"),
                     " NA")
        }
    })

    #Check all other parameters
    non.time.varying.param.names = names(parameters)[names(parameters)!='time.varying.parameters']
    sapply(non.time.varying.param.names, function(name){
        val = parameters[[name]]
        is.multiple = class(val)=='list' || is.null(length(val)) || length(val) > 1

        if (any(is.na(unlist(val))))
            stop("The parameter '", name, "' ",
                 ifelse(is.multiple, "has components that are", "is"),
                 " NA")
    })


    #Return success
    T
}

prepare.transition.mappings <- function(jheem)
{
    jheem$parameters$transition.mappings = list()

    #-- Set up name to index mappings --#
    N.HIV = length(jheem$parameters$hiv.positive.transition.dimensions)
    hiv.name.to.index = 0:(N.HIV-1)
    names(hiv.name.to.index) = jheem$parameters$hiv.positive.transition.dimensions

    N.NONHIV = length(jheem$parameters$hiv.negative.transition.dimensions)
    nonhiv.name.to.index = 0:(N.NONHIV-1)
    names(nonhiv.name.to.index) = jheem$parameters$hiv.negative.transition.dimensions

    #-- Set up name to dimension mappings --#
    hiv.dimname.names = names(get.dimnames.hiv(jheem))
    hiv.name.to.dimension = 0:(length(hiv.dimname.names)-1)
    names(hiv.name.to.dimension) = hiv.dimname.names

    nonhiv.dimname.names = names(get.dimnames.nonhiv(jheem))
    nonhiv.name.to.dimension = 0:(length(nonhiv.dimname.names)-1)
    names(nonhiv.name.to.dimension) = nonhiv.dimname.names

    #-- Prepare the transition array names --#
    jheem$parameters$transition.mappings$hiv.array.names = paste0('HIV_POSITIVE_TRANSITION_ARRAY_', toupper(jheem$parameters$hiv.positive.transition.dimensions))
    jheem$parameters$transition.mappings$nonhiv.array.names = paste0('HIV_NEGATIVE_TRANSITION_ARRAY_', toupper(jheem$parameters$hiv.negative.transition.dimensions))

    #-- Prepare the (numeric) dimensions for transition arrays --#
    jheem$parameters$transition.mappings$dimensions.for.hiv.arrays = hiv.name.to.dimension[jheem$parameters$hiv.positive.transition.dimensions]
    jheem$parameters$transition.mappings$dimensions.for.nonhiv.arrays = hiv.name.to.dimension[jheem$parameters$hiv.negative.transition.dimensions]

    #-- Prepare the (numeric) dimensions for tracked transition arrays --#
    jheem$parameters$transition.mappings$tracked.hiv.array.indices = hiv.name.to.index[jheem$tracked.transitions$transition.dimensions[jheem$tracked.transitions$names.for.indices.hiv]]
    jheem$parameters$transition.mappings$tracked.nonhiv.array.indices = nonhiv.name.to.index[jheem$tracked.transitions$transition.dimensions[jheem$tracked.transitions$names.for.indices.nonhiv]]

    #-- Return --#
    jheem
}
